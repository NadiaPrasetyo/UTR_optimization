#!/usr/bin/env python3
"""
prepare_af3_jobs.py

Read a FASTA file (one or more sequences) and prepare everything needed to
run AlphaFold3 on a SLURM cluster:

  1. One AlphaFold3 JSON input file per sequence.
       - default mode: sequences are grouped into batches of --batch-size
         (default 25) and submitted as a SLURM array job, with each array
         task calling AF3 ONCE per batch via --input_dir/--output_dir
         (AF3's documented multi-input mode) instead of once per sequence.
         This matters a lot in practice: AF3 pays its model/database load
         cost once per array task, not once per sequence, so batching
         sharply cuts both wall-clock time and the number of array tasks
         competing for GPU concurrency. Use --batch-size 1 for the old
         one-sequence-per-task behaviour.
       - --complex mode: all FASTA records are combined into ONE job as
         multiple chains of a single complex (e.g. for predicting a
         protein-protein / protein-ligand-free complex) — already a single
         af3 invocation, so batching doesn't apply here.
  2. A batch manifest (batches.tsv, or jobs.tsv in --complex mode) mapping
     SLURM array task IDs to their input/output directories.
  3. A ready-to-submit SLURM array script (run_af3_array.sh) that sources
     the manifest and calls `af3` for the task's batch.

This module also doubles as a small library: mutate_3utr_ga.py imports
build_af3_json / sanitize_name / clean_RNA_seq from it directly so both
scripts share exactly the same AF3 JSON conventions.

Usage examples
--------------
# Each sequence in multi.fasta becomes its own AF3 job, run as an array:
python3 prepare_af3_jobs.py \\
    --fasta multi.fasta \\
    --work-dir /projects/.../a7s \\
    --db /opt/alphafold3/databases \\
    --models /projects/.../alphafold3-weights

# Combine every sequence in the fasta into a single multi-chain complex job:
python3 prepare_af3_jobs.py \\
    --fasta complex.fasta \\
    --work-dir /projects/.../a7s \\
    --complex \\
    --job-name my_complex

After running, submit with:
    sbatch <work-dir>/af3_jobs/run_af3_array.sh
"""

import argparse
import json
import re
import sys
from pathlib import Path


# --------------------------------------------------------------------------
# FASTA parsing
# --------------------------------------------------------------------------
def parse_fasta(fasta_path: Path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_chunks = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)

    if header is None:
        raise ValueError(f"No FASTA records found in {fasta_path}")


def sanitize_name(name: str) -> str:
    """Turn a FASTA header (or any string) into a safe job/file name."""
    name = name.split()[0] if name.split() else name  # first whitespace token
    name = name.lower()
    name = re.sub(r"[^a-z0-9_\-]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_-")
    return name or "job"


def clean_RNA_seq(seq: str) -> str:
    return seq.replace("T", "U").upper()


# --------------------------------------------------------------------------
# AlphaFold3 JSON construction
# --------------------------------------------------------------------------
def build_af3_json(job_name: str, chains, model_seeds):
    """
    chains: list of dicts {"id": "A", "sequence": "...", "type": "protein"}
    """
    sequences = []
    for chain in chains:
        mol_type = chain["type"]
        sequences.append({
            mol_type: {
                "id": chain["id"],
                "sequence": chain["sequence"],
                "unpairedMsa": ""
            }
        })

    return {
        "name": job_name,
        "sequences": sequences,
        "modelSeeds": model_seeds,
        "dialect": "alphafold3",
        "version": 1,
    }


def chain_ids_generator():
    """A, B, C, ... Z, AA, AB, ... for chain IDs beyond 26 chains."""
    import string
    letters = string.ascii_uppercase
    n = 1
    while True:
        for combo in _product(letters, n):
            yield "".join(combo)
        n += 1


def _product(letters, n):
    import itertools
    return itertools.product(letters, repeat=n)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Prepare AlphaFold3 JSON inputs and a SLURM array script from a FASTA file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--fasta", required=True, type=Path,
                    help="Input FASTA file (one or more sequences).")
    p.add_argument("--work-dir", required=True, type=Path,
                    help="AlphaFold3 working directory (AF3_WD). JSON inputs, "
                         "job list, logs, and the SLURM script are written under "
                         "<work-dir>/af3_jobs/, and af3 output goes to "
                         "<work-dir>/af3_output/.")
    p.add_argument("--db", default="/opt/alphafold3/databases",
                    help="Path for AF3_DB (databases). Default: %(default)s")
    p.add_argument("--models", required=True,
                    help="Path for AF3_MODELS (model weights).")
    p.add_argument("--module", default="alphafold3/3.0.3",
                    help="Environment module to load. Default: %(default)s")
    p.add_argument("--complex", action="store_true",
                    help="Combine ALL sequences in the FASTA into a single "
                         "multi-chain job instead of one job per sequence.")
    p.add_argument("--job-name", default=None,
                    help="Job name to use in --complex mode (default: FASTA "
                         "filename stem). Ignored otherwise (each record's "
                         "header is used as its job name).")
    p.add_argument("--model-seed", type=int, action="append", default=None,
                    help="Model seed(s) to use, can repeat this flag for "
                         "multiple seeds. Default: single seed 1.")
    p.add_argument("--batch-size", type=int, default=25,
                    help="Number of sequences packed into a single AF3 "
                         "invocation (via --input_dir) per SLURM array task "
                         "(default: %(default)s). Ignored in --complex mode, "
                         "which is already one job. AF3 pays its "
                         "model/database load cost once per array task, not "
                         "once per sequence, so batching is the main lever "
                         "on wall-clock time for large sequence sets — set "
                         "to 1 to get the old one-sequence-per-task behaviour.")
    # SLURM options (mirrors the example script's defaults)
    p.add_argument("--partition", default="aoraki_gpu_L40", help="SLURM partition default: %(default)s")
    p.add_argument("--time", default="02:00:00", help="SLURM time default: %(default)s "
                    "(scale this up with --batch-size — a bigger batch takes longer per task)")
    p.add_argument("--mem", default="32G", help="SLURM mem default: %(default)s")
    p.add_argument("--cpus-per-task", type=int, default=8, help="SLURM cpus-per-task default: %(default)s")
    p.add_argument("--gres", default="gpu:1", help="SLURM gres default: %(default)s")
    p.add_argument("--exclude", default="aoraki18",
                    help="Nodes to exclude, comma separated. Pass '' to not exclude any. Default: %(default)s")
    p.add_argument("--array-max-parallel", type=int, default=None,
                    help="Cap on simultaneously running array tasks, e.g. 4 "
                         "for '--array=0-9%%4'. Default: no cap.")

    args = p.parse_args()

    if not args.fasta.exists():
        sys.exit(f"ERROR: FASTA file not found: {args.fasta}")

    records = list(parse_fasta(args.fasta))
    if not records:
        sys.exit(f"ERROR: no sequences found in {args.fasta}")

    model_seeds = args.model_seed if args.model_seed else [1]

    jobs_dir = args.work_dir / "af3_jobs"
    inputs_dir = jobs_dir / "inputs"
    logs_dir = jobs_dir / "logs"
    output_root = args.work_dir / "af3_output"
    for d in (inputs_dir, logs_dir, output_root):
        d.mkdir(parents=True, exist_ok=True)

    exclude_line = f"#SBATCH --exclude={args.exclude}\n" if args.exclude else ""

    if args.complex:
        # A single multi-chain job is already one af3 invocation — batching
        # doesn't apply, so this keeps the original single-JSON-per-task
        # (positional-args) script shape.
        job_name = sanitize_name(args.job_name or args.fasta.stem)
        chain_ids = chain_ids_generator()
        chains = []
        for header, seq in records:
            if not seq:
                sys.exit(f"ERROR: empty sequence for record '{header}'")
            chains.append({
                "id": next(chain_ids),
                "sequence": clean_RNA_seq(seq),
                "type": "rna",
            })
        af3_json = build_af3_json(job_name, chains, model_seeds)
        json_path = inputs_dir / f"{job_name}.json"
        json_path.write_text(json.dumps(af3_json, indent=2))
        print(f"[complex mode] Combined {len(chains)} chain(s) into job '{job_name}'")

        jobs_tsv = jobs_dir / "jobs.tsv"
        out_dir = output_root / job_name
        with open(jobs_tsv, "w") as fh:
            fh.write(f"0\t{job_name}\t{json_path}\t{out_dir}\n")

        script = f"""#!/bin/bash
#SBATCH --job-name=af3_{job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cpus_per_task}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --partition={args.partition}
#SBATCH --gres={args.gres}
#SBATCH --output={logs_dir}/af3_%j.log
{exclude_line}
echo "Job started on $(hostname) at $(date)"

module purge
module load {args.module}

export AF3_DB={args.db}
export AF3_WD={args.work_dir}
export AF3_MODELS={args.models}

JSON_PATH="{json_path}"
OUT_DIR="{out_dir}"
mkdir -p "$OUT_DIR"

echo "Running AlphaFold3 for job: {job_name}"
af3 "$JSON_PATH" "$OUT_DIR" --run_data_pipeline=true --run_inference=true

echo "AlphaFold3 job '{job_name}' finished at $(date)."
"""
        script_path = jobs_dir / "run_af3_array.sh"
        script_path.write_text(script)
        script_path.chmod(0o755)

        print(f"\nWrote 1 AF3 input JSON file to: {inputs_dir}")
        print(f"Wrote job list to:              {jobs_tsv}")
        print(f"Wrote SLURM script to:          {script_path}")
        print(f"\nSubmit with:\n  sbatch {script_path}")
        return

    # ---- per-sequence mode: batch multiple JSONs per array task ----------
    seen_names = set()
    named_records = []  # (job_name, seq)
    for header, seq in records:
        if not seq:
            sys.exit(f"ERROR: empty sequence for record '{header}'")
        job_name = sanitize_name(header)
        base_name = job_name
        suffix = 2
        while job_name in seen_names:
            job_name = f"{base_name}_{suffix}"
            suffix += 1
        seen_names.add(job_name)
        named_records.append((job_name, seq))

    batch_size = max(1, args.batch_size)
    batches = [named_records[i:i + batch_size] for i in range(0, len(named_records), batch_size)]

    batch_rows = []
    for batch_idx, batch in enumerate(batches):
        batch_input_dir = inputs_dir / f"batch_{batch_idx:04d}"
        batch_output_dir = output_root / f"batch_{batch_idx:04d}"
        batch_input_dir.mkdir(parents=True, exist_ok=True)
        batch_output_dir.mkdir(parents=True, exist_ok=True)
        for job_name, seq in batch:
            chains = [{"id": "A", "sequence": clean_RNA_seq(seq), "type": "rna"}]
            af3_json = build_af3_json(job_name, chains, model_seeds)
            (batch_input_dir / f"{job_name}.json").write_text(json.dumps(af3_json, indent=2))
        batch_rows.append((batch_idx, batch_input_dir, batch_output_dir, len(batch)))

    print(f"[per-sequence mode] Prepared {len(named_records)} sequence(s) from {args.fasta} "
          f"as {len(batches)} batch(es) of <= {batch_size} sequence(s)/array task")

    # Batch manifest (tab separated: array_index, batch_input_dir, batch_output_dir, n_seq)
    jobs_tsv = jobs_dir / "batches.tsv"
    with open(jobs_tsv, "w") as fh:
        for batch_idx, batch_input_dir, batch_output_dir, n_seq in batch_rows:
            fh.write(f"{batch_idx}\t{batch_input_dir}\t{batch_output_dir}\t{n_seq}\n")

    n_tasks = len(batches)
    array_range = f"0-{n_tasks - 1}"
    if args.array_max_parallel:
        array_range += f"%{args.array_max_parallel}"

    script = f"""#!/bin/bash
#SBATCH --job-name=af3_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cpus_per_task}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --partition={args.partition}
#SBATCH --gres={args.gres}
#SBATCH --array={array_range}
#SBATCH --output={logs_dir}/af3_%A_%a.log
{exclude_line}
echo "=== AF3 batch task starting on $(hostname) at $(date) ==="
echo "Array task ID: $SLURM_ARRAY_TASK_ID  (batch size <= {batch_size} sequences/task)"

# Load the AlphaFold3 module
module purge
module load {args.module}

# AlphaFold3 environment
export AF3_DB={args.db}
export AF3_WD={args.work_dir}
export AF3_MODELS={args.models}
echo "AF3_MODELS is: $AF3_MODELS"

# Look up this array task's batch in batches.tsv (columns: idx, input_dir, output_dir, n_seq)
BATCHES_TSV="{jobs_tsv}"
LINE=$(awk -F'\\t' -v idx="$SLURM_ARRAY_TASK_ID" '$1 == idx {{print}}' "$BATCHES_TSV")

if [ -z "$LINE" ]; then
    echo "ERROR: no batch found for array task ID $SLURM_ARRAY_TASK_ID in $BATCHES_TSV"
    exit 1
fi

BATCH_INPUT_DIR=$(echo "$LINE" | cut -f2)
BATCH_OUTPUT_DIR=$(echo "$LINE" | cut -f3)
N_SEQ=$(echo "$LINE" | cut -f4)
mkdir -p "$BATCH_OUTPUT_DIR"

echo "Batch input dir  : $BATCH_INPUT_DIR  ($N_SEQ sequences)"
echo "Batch output dir : $BATCH_OUTPUT_DIR"
echo "--- af3 command ---"
set -x
af3 --input_dir="$BATCH_INPUT_DIR" --output_dir="$BATCH_OUTPUT_DIR" \\
    --model_dir="$AF3_MODELS" --db_dir="$AF3_DB" \\
    --run_data_pipeline=true --run_inference=true
AF3_EXIT=$?
set +x

echo "=== AF3 batch task finished at $(date) with exit code $AF3_EXIT ==="
exit $AF3_EXIT
"""

    script_path = jobs_dir / "run_af3_array.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)

    print(f"\nWrote {len(named_records)} AF3 input JSON file(s) to: {inputs_dir}")
    print(f"  (grouped into {n_tasks} batch(es) of <= {batch_size} sequence(s) each)")
    print(f"Wrote batch manifest to:  {jobs_tsv}")
    print(f"Wrote SLURM array script: {script_path}")
    print(f"\nSubmit with:\n  sbatch {script_path}")
    print(f"\nEach individual sequence's output will land under:\n"
          f"  {output_root}/batch_NNNN/<sanitized_sequence_name>/")


if __name__ == "__main__":
    main()