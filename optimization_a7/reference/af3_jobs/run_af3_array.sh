#!/bin/bash
#SBATCH --job-name=af3_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=aoraki_gpu_L40
#SBATCH --gres=gpu:1
#SBATCH --array=0-2
#SBATCH --output=optimization_a7/reference/af3_jobs/logs/af3_%A_%a.log
#SBATCH --exclude=aoraki08

echo "Job started on $(hostname) at $(date)"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

# Load the AlphaFold3 module
module purge
module load alphafold3/3.0.3

# AlphaFold3 environment
export AF3_DB=/opt/alphafold3/databases
export AF3_WD=optimization_a7/reference
export AF3_MODELS=/projects/health_sciences/bms/biochemistry/gardner_group/alphafold3-weights
echo "AF3_MODELS is: $AF3_MODELS"

# Look up this array task's job in jobs.tsv (columns: idx, name, json_path, out_dir)
JOBS_TSV="optimization_a7/reference/af3_jobs/jobs.tsv"
LINE=$(awk -F'\t' -v idx="$SLURM_ARRAY_TASK_ID" '$1 == idx {print}' "$JOBS_TSV")

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
