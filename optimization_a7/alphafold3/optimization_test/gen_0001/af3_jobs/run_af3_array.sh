#!/bin/bash
#SBATCH --job-name=af3_gen0001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=aoraki_gpu_L40
#SBATCH --gres=gpu:1
#SBATCH --array=0-999
#SBATCH --output=optimization_a7/alphafold3/optimization_test/gen_0001/af3_jobs/logs/af3_%A_%a.log

module purge
module load alphafold3/3.0.3

export AF3_DB=/opt/alphafold3/databases
export AF3_WD=optimization_a7/alphafold3/optimization_test/gen_0001
export AF3_MODELS=/projects/health_sciences/bms/biochemistry/gardner_group/alphafold3-weights

JOBS_TSV="optimization_a7/alphafold3/optimization_test/gen_0001/af3_jobs/jobs.tsv"
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
af3 "$JSON_PATH" "$OUT_DIR" --run_data_pipeline=true --run_inference=true
echo "AlphaFold3 job '$JOB_NAME' finished at $(date)."
