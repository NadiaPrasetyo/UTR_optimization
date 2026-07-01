#!/usr/bin/env python3
"""
prepare_af3_jobs.py

Read a FASTA file (one or more sequences) and prepare everything needed to
run AlphaFold3 on a SLURM cluster:

  1. One AlphaFold3 JSON input file per job.
       - default mode: each FASTA record -> its own separate job (monomer
         prediction), submitted together as a SLURM array job.
       - --complex mode: all FASTA records are combined into ONE job as
         multiple chains of a single complex (e.g. for predicting a
         protein-protein / protein-ligand-free complex).
  2. A job list file (jobs.tsv) mapping SLURM array task IDs to job
     name / JSON path / output directory.
  3. A ready-to-submit SLURM array script (run_af3_array.sh) that sources
     the job list and calls `af3` for the task's job.

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


def detect_molecule_type(sequence: str) -> str:
    """
    Very simple heuristic to decide if a sequence is protein, DNA, or RNA.
    AlphaFold3 JSON needs the sequence tagged as protein/dna/rna.
    """
    s = sequence.upper()
    bases = set(s)
    if bases <= set("ACGT N"):
        return "dna"
    if bases <= set("ACGU N"):
        return "rna"
    return "protein"


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
                         "<work-dir>/af3_output/<job_name>/.")
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
    # SLURM options (mirrors the example script's defaults)
    p.add_argument("--partition", default="aoraki_gpu_L40")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--mem", default="32G")
    p.add_argument("--cpus-per-task", type=int, default=8)
    p.add_argument("--gres", default="gpu:1")
    p.add_argument("--exclude", default="aoraki18",
                    help="Nodes to exclude, comma separated. Pass '' to not exclude any.")
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

    jobs = []  # list of (job_name, json_path)

    if args.complex:
        job_name = sanitize_name(args.job_name or args.fasta.stem)
        chain_ids = chain_ids_generator()
        chains = []
        for header, seq in records:
            if not seq:
                sys.exit(f"ERROR: empty sequence for record '{header}'")
            chains.append({
                "id": next(chain_ids),
                "sequence": seq,
                "type": detect_molecule_type(seq),
            })
        af3_json = build_af3_json(job_name, chains, model_seeds)
        json_path = inputs_dir / f"{job_name}.json"
        json_path.write_text(json.dumps(af3_json, indent=2))
        jobs.append((job_name, json_path))
        print(f"[complex mode] Combined {len(chains)} chain(s) into job '{job_name}'")
    else:
        seen_names = set()
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

            chains = [{
                "id": "A",
                "sequence": seq,
                "type": detect_molecule_type(seq),
            }]
            af3_json = build_af3_json(job_name, chains, model_seeds)
            json_path = inputs_dir / f"{job_name}.json"
            json_path.write_text(json.dumps(af3_json, indent=2))
            jobs.append((job_name, json_path))
        print(f"[per-sequence mode] Prepared {len(jobs)} job(s) from {args.fasta}")

    # Write job list (tab separated: array_index, job_name, json_path, out_dir)
    jobs_tsv = jobs_dir / "jobs.tsv"
    with open(jobs_tsv, "w") as fh:
        for idx, (job_name, json_path) in enumerate(jobs):
            out_dir = output_root / job_name
            fh.write(f"{idx}\t{job_name}\t{json_path}\t{out_dir}\n")

    # Build the SLURM array script
    n_jobs = len(jobs)
    array_range = f"0-{n_jobs - 1}"
    if args.array_max_parallel:
        array_range += f"%{args.array_max_parallel}"

    exclude_line = f"#SBATCH --exclude={args.exclude}\n" if args.exclude else ""

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
echo "Job started on $(hostname) at $(date)"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

# Load the AlphaFold3 module
module purge
module load {args.module}

# AlphaFold3 environment
export AF3_DB={args.db}
export AF3_WD={args.work_dir}
export AF3_MODELS={args.models}
echo "AF3_MODELS is: $AF3_MODELS"

# Look up this array task's job in jobs.tsv (columns: idx, name, json_path, out_dir)
JOBS_TSV="{jobs_tsv}"
LINE=$(awk -F'\\t' -v idx="$SLURM_ARRAY_TASK_ID" '$1 == idx {{print}}' "$JOBS_TSV")

if [ -z "$LINE" ]; then
    echo "ERROR: no job found for array task ID $SLURM_ARRAY_TASK_ID in $JOBS_TSV"
    exit 1
fi

JOB_NAME=$(echo "$LINE" | cut -f2)
JSON_PATH=$(echo "$LINE" | cut -f3)
OUT_DIR=$(echo "$LINE" | cut -f4)

mkdir -p "$OUT_DIR"

echo "Running AlphaFold3 for job: $JOB_NAME"
echo "  input JSON: $JSON_PATH"
echo "  output dir: $OUT_DIR"

af3 "$JSON_PATH" "$OUT_DIR" --run_data_pipeline=true --run_inference=true

echo "AlphaFold3 job '$JOB_NAME' finished at $(date)."
"""

    script_path = jobs_dir / "run_af3_array.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)

    print(f"\nWrote {n_jobs} AF3 input JSON file(s) to: {inputs_dir}")
    print(f"Wrote job list to:                    {jobs_tsv}")
    print(f"Wrote SLURM array script to:           {script_path}")
    print(f"\nSubmit with:\n  sbatch {script_path}")


if __name__ == "__main__":
    main()