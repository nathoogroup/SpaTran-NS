#!/usr/bin/env bash
#SBATCH --job-name=spatran_cov_prepare
#SBATCH --account=def-nathoo
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cov_prepare_%j.out
#SBATCH --error=logs/cov_prepare_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
source hpc/covariance_comparison/environment.sh
: "${RUN_DIR:?Set RUN_DIR to a NEW output directory}"
Rscript hpc/covariance_comparison/prepare.R \
  "${DATA_DIR:-data/spatial_datasets}" "$RUN_DIR" "${GENES_PER_TASK:-100}" \
  "${DATASET_LIST:-hpc/analysis_datasets.txt}"
