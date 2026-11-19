#!/usr/bin/env bash
#SBATCH --job-name=spatran_cov_gene
#SBATCH --account=def-nathoo
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cov_gene_%A_%a.out
#SBATCH --error=logs/cov_gene_%A_%a.err
set -euo pipefail
# submit.sh sets --chdir to the repository root, independent of the caller.
source hpc/covariance_comparison/environment.sh
: "${RUN_DIR:?Use submit.sh or export RUN_DIR}"
Rscript hpc/covariance_comparison/run_task.R "$RUN_DIR"
