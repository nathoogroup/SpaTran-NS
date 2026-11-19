#!/usr/bin/env bash
#SBATCH --job-name=spatran_cov_combine
#SBATCH --account=def-nathoo
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cov_combine_%j.out
#SBATCH --error=logs/cov_combine_%j.err
set -euo pipefail
source hpc/covariance_comparison/environment.sh
: "${RUN_DIR:?Set RUN_DIR}"
Rscript hpc/covariance_comparison/combine.R "$RUN_DIR"
