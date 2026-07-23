#!/bin/bash
#SBATCH --job-name=af3_gen0009_b000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=aoraki_gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-9
#SBATCH --output=optimization_a7/alphafold3/optimization_test/gen_0009/af3_jobs/batch_0000/logs/af3_%A_%a.log

module purge
module load alphafold3/3.0.3

export AF3_DB=/opt/alphafold3/databases
export AF3_WD=optimization_a7/alphafold3/optimization_test/gen_0009
export AF3_MODELS=/projects/health_sciences/bms/biochemistry/gardner_group/alphafold3-weights

PER_GPU=10
SLOT_START=$((SLURM_ARRAY_TASK_ID * PER_GPU))
SLOT_END=$((SLOT_START + PER_GPU - 1))
LINES=$(awk -F'\t' -v s="$SLOT_START" -v e="$SLOT_END" '$1 >= s && $1 <= e' "optimization_a7/alphafold3/optimization_test/gen_0009/af3_jobs/batch_0000/jobs.tsv")

echo "Slot $SLURM_ARRAY_TASK_ID: $(echo "$LINES" | wc -l) job(s) running sequentially on this GPU."

run_one() {
    local job_name="$1" json_path="$2" out_dir="$3"
    mkdir -p "$out_dir"
    rm -f "$out_dir/.failed"
    echo "[$job_name] starting at $(date)"
    af3 "$json_path" "$out_dir" --run_data_pipeline=true --run_inference=true \
        > "optimization_a7/alphafold3/optimization_test/gen_0009/af3_jobs/batch_0000/logs/af3_task_${job_name}.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[$job_name] FAILED at $(date) (exit $rc) — see optimization_a7/alphafold3/optimization_test/gen_0009/af3_jobs/batch_0000/logs/af3_task_${job_name}.log"
        echo "$rc" > "$out_dir/.failed"
        return $rc
    fi
    local n_cifs
    n_cifs=$(find "$out_dir" -path "$out_dir/**/*_model.cif" 2>/dev/null | wc -l)
    if [ "$n_cifs" -eq 0 ]; then
        echo "[$job_name] exited 0 but no output matching '**/*_model.cif' found — treating as FAILED (likely OOM)."
        echo "no_output_despite_exit_0" > "$out_dir/.failed"
        return 1
    fi
    echo "[$job_name] finished at $(date), $n_cifs structure file(s) confirmed."
}

slot_failed=0
while IFS=$'\t' read -r idx job_name json_path out_dir; do
    [ -z "$job_name" ] && continue
    run_one "$job_name" "$json_path" "$out_dir" || slot_failed=1
done <<< "$LINES"

echo "Slot $SLURM_ARRAY_TASK_ID: all assigned jobs finished at $(date)."
[ "$slot_failed" -ne 0 ] && exit 1
exit 0
