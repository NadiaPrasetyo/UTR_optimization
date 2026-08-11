#!/bin/bash
#SBATCH --job-name=af3_gen0081_b000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=aoraki_gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-9
#SBATCH --output=optimization_a7/absolute_optimization/gen_0081/af3_jobs/batch_0000/logs/af3_%A_%a.log

module purge
module load alphafold3/3.0.3

export AF3_DB=/opt/alphafold3/databases
export AF3_WD=optimization_a7/absolute_optimization/gen_0081
export AF3_MODELS=/projects/health_sciences/bms/biochemistry/gardner_group/alphafold3-weights

# --- Sanitize LD_PRELOAD -------------------------------------------------
# Some module environments inject LD_PRELOAD entries (e.g. a distro
# libstdc++.so.6 path) that don't exist on every node type in a
# heterogeneous cluster. When that happens, ld.so silently ignores the
# preload ("cannot be preloaded ... ignored"), but the broken/partial
# preload can still corrupt symbol resolution for the CUDA driver
# libraries and manifest downstream as jax/cuInit failures
# (CUDA_ERROR_NO_DEVICE) that have nothing to do with the GPU itself.
# Strip any entries that don't actually exist on this node before running
# anything.
if [ -n "$LD_PRELOAD" ]; then
    _valid_preload=""
    for _lib in $(echo "$LD_PRELOAD" | tr ':' ' '); do
        if [ -f "$_lib" ]; then
            _valid_preload="${_valid_preload:+$_valid_preload:}$_lib"
        else
            echo "WARNING: dropping missing LD_PRELOAD entry on $(hostname): $_lib"
        fi
    done
    export LD_PRELOAD="$_valid_preload"
fi
# --------------------------------------------------------------------------

PER_GPU=10
SLOT_START=$((SLURM_ARRAY_TASK_ID * PER_GPU))
SLOT_END=$((SLOT_START + PER_GPU - 1))
LINES=$(awk -F'\t' -v s="$SLOT_START" -v e="$SLOT_END" '$1 >= s && $1 <= e' "optimization_a7/absolute_optimization/gen_0081/af3_jobs/batch_0000/jobs.tsv")

echo "Slot $SLURM_ARRAY_TASK_ID: $(echo "$LINES" | wc -l) job(s) running sequentially on this GPU."
echo "Slot $SLURM_ARRAY_TASK_ID: running on node $(hostname), CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

# --- GPU health check -------------------------------------------------
# nvidia-smi confirms the driver/device is visible to this allocation.
# It does NOT guarantee jax/CUDA can actually cuInit() (that can still
# fail even when nvidia-smi succeeds), so we also do a lightweight
# python probe using the same jax the af3 module provides.
gpu_healthy() {
    if ! nvidia-smi --query-gpu=index,name,pstate,memory.used,memory.total \
            --format=csv,noheader > "optimization_a7/absolute_optimization/gen_0081/af3_jobs/batch_0000/logs/nvidia_smi_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" 2>&1; then
        echo "GPU health check: nvidia-smi failed — GPU not visible to this allocation."
        return 1
    fi
    probe_out=$(python3 -c "
import sys
try:
    import jax
except ImportError as e:
    print(f'SKIP: jax not importable here ({e}) — cannot confirm via jax, relying on nvidia-smi only', file=sys.stderr)
    sys.exit(2)
try:
    d = jax.local_devices(backend='gpu')
    sys.exit(0 if d else 1)
except Exception as e:
    print(f'jax GPU probe failed: {e}', file=sys.stderr)
    sys.exit(1)
" 2>&1)
    rc=$?
    echo "$probe_out" >> "optimization_a7/absolute_optimization/gen_0081/af3_jobs/batch_0000/logs/nvidia_smi_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
    if [ $rc -eq 1 ]; then
        echo "GPU health check: jax cannot see a GPU device (cuInit likely failing)."
        return 1
    fi
    # rc==0 (jax confirms GPU) or rc==2 (jax unavailable, trust nvidia-smi) both pass.
    return 0
}

# Recognize the specific GPU-init failure signature in a completed job's
# af3 log. This exists because we've observed gpu_healthy() pass
# immediately before a job, and af3 still fail with CUDA_ERROR_NO_DEVICE —
# i.e. the lightweight pre-check doesn't always catch it. Once a GPU is in
# this state it tends to stay broken for the rest of the slot, so detecting
# it from the actual failure lets us stop wasting time retrying every
# remaining individual on the same dead GPU.
is_gpu_init_failure() {
    local logfile="$1"
    grep -qE 'CUDA_ERROR_NO_DEVICE|cuInit\(0\) failed|no platforms that are instances of gpu are present' "$logfile" 2>/dev/null
}

mark_slot_gpu_unavailable() {
    # Distinguish infra failures from real per-individual failures so
    # downstream GA logic doesn't penalize an individual's fitness for
    # a problem that wasn't caused by its input.
    while IFS=$'\t' read -r idx job_name json_path out_dir; do
        [ -z "$job_name" ] && continue
        mkdir -p "$out_dir"
        echo "gpu_unavailable" > "$out_dir/.failed"
        echo "[$job_name] SKIPPED — GPU unavailable on node $(hostname) for this slot."
    done <<< "$LINES"
}

if ! gpu_healthy; then
    echo "Slot $SLURM_ARRAY_TASK_ID: GPU unhealthy at slot start on node $(hostname) — failing fast instead of burning $(echo "$LINES" | wc -l) doomed job(s)."
    mark_slot_gpu_unavailable
    exit 1
fi
# ------------------------------------------------------------------------

run_one() {
    local job_name="$1" json_path="$2" out_dir="$3"
    mkdir -p "$out_dir"
    rm -f "$out_dir/.failed"

    # Retry budget for THIS individual only. The GPU-init failure
    # (CUDA_ERROR_NO_DEVICE / cuInit failing) has been observed to be
    # transient/intermittent rather than a permanently dead GPU — other
    # individuals in the same slot, before and after, can succeed just
    # fine. So on this specific signature we retry the same job with
    # backoff rather than giving up on the rest of the slot.
    local max_attempts=3
    local attempt=1
    local task_log="optimization_a7/absolute_optimization/gen_0081/af3_jobs/batch_0000/logs/af3_task_${job_name}.log"

    while [ $attempt -le $max_attempts ]; do
        if ! gpu_healthy; then
            echo "[$job_name] GPU unhealthy before attempt $attempt/$max_attempts at $(date) on node $(hostname)."
            if [ $attempt -lt $max_attempts ]; then
                sleep $((attempt * 10))
                attempt=$((attempt + 1))
                continue
            fi
            echo "gpu_unavailable" > "$out_dir/.failed"
            echo "[$job_name] SKIPPED — GPU still unhealthy after $max_attempts attempts on node $(hostname)."
            return 1
        fi

        echo "[$job_name] starting (attempt $attempt/$max_attempts) at $(date) on node $(hostname)"
        af3 "$json_path" "$out_dir" --run_data_pipeline=true --run_inference=true \
            > "$task_log" 2>&1
        local rc=$?

        if [ $rc -eq 0 ]; then
            local n_cifs
            n_cifs=$(find "$out_dir" -path "$out_dir/**/*_model.cif" 2>/dev/null | wc -l)
            if [ "$n_cifs" -gt 0 ]; then
                echo "[$job_name] finished at $(date), $n_cifs structure file(s) confirmed."
                return 0
            fi
            echo "[$job_name] exited 0 but no output matching '**/*_model.cif' found — treating as FAILED (likely OOM)."
            echo "no_output_despite_exit_0" > "$out_dir/.failed"
            return 1
        fi

        echo "[$job_name] FAILED at $(date) (exit $rc, attempt $attempt/$max_attempts) — see $task_log"
        # Capture GPU state at time of failure to help distinguish OOM vs
        # driver/infra faults after the fact.
        nvidia-smi --query-gpu=index,memory.used,memory.total,ecc.errors.uncorrected.volatile.total,pstate \
            --format=csv >> "$task_log" 2>&1

        if is_gpu_init_failure "$task_log"; then
            if [ $attempt -lt $max_attempts ]; then
                echo "[$job_name] transient GPU-init signature (CUDA_ERROR_NO_DEVICE) — retrying after backoff."
                sleep $((attempt * 10))
                attempt=$((attempt + 1))
                continue
            fi
            echo "[$job_name] still hitting GPU-init failures after $max_attempts attempts — giving up on this individual, moving on to the next one."
            echo "gpu_unavailable" > "$out_dir/.failed"
            return 1
        fi

        # Non-GPU-init failure (e.g. real OOM, bad input) — don't retry,
        # record the real exit code so downstream fitness logic can see it.
        echo "$rc" > "$out_dir/.failed"
        return $rc
    done
}

slot_failed=0
while IFS=$'\t' read -r idx job_name json_path out_dir; do
    [ -z "$job_name" ] && continue
    run_one "$job_name" "$json_path" "$out_dir" || slot_failed=1
done <<< "$LINES"

echo "Slot $SLURM_ARRAY_TASK_ID: finished at $(date)."
[ "$slot_failed" -ne 0 ] && exit 1
exit 0
