#!/bin/bash
#SBATCH --job-name=spatran_ns
#SBATCH --output=logs/spatran_ns_%A_%a.out
#SBATCH --error=logs/spatran_ns_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=1-18
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Spatial Nonstationarity Analysis - Array Job Submission Script
# =============================================================================
# One SLURM array task per dataset.  Tasks run in parallel on the cluster,
# each processing one dataset with 8 cores for gene-level parallelism.
#
# Runtime estimate (per task):
#   ~16,000 genes / 8 cores × ~10 sec/gene ≈ 333 min → 6h walltime is safe
#
# Prerequisites:
#   1. bash hpc/install_inla.sh      (once, on login node)
#   2. sbatch hpc/submit_download.sbatch  (download datasets)
#
# Update --array upper bound to match number of unique datasets:
#   ls data/spatial_datasets/*.rds | grep -v "v3_13\|log\|manifest" | wc -l
#
# Submit with: sbatch hpc/submit_analysis.sbatch
# =============================================================================

echo "=========================================="
echo "Array job started at:   $(date)"
echo "Array job ID:           $SLURM_ARRAY_JOB_ID"
echo "Array task ID:          $SLURM_ARRAY_TASK_ID"
echo "Node:                   $SLURM_NODELIST"
echo "CPUs:                   $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

# Core modules (same as install_inla.sh) + Bioconductor bundle for system packages
module load StdEnv/2023 gcc/12.3 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

# Set R library path (must match install_inla.sh)
export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

DATA_DIR="data/spatial_datasets"
OUTPUT_DIR="results"
N_CORES=$SLURM_CPUS_PER_TASK

mkdir -p "$OUTPUT_DIR"

echo ""
echo "Configuration:"
echo "  Data directory:   $DATA_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Cores:            $N_CORES"
echo "  Dataset index:    $SLURM_ARRAY_TASK_ID"
echo ""

# SLURM_ARRAY_TASK_ID is read inside run_analysis.R to select the dataset
Rscript hpc/run_analysis.R "$DATA_DIR" "$OUTPUT_DIR" "$N_CORES"

EXIT_STATUS=$?
echo ""
echo "=========================================="
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully at: $(date)"
else
    echo "Task $SLURM_ARRAY_TASK_ID FAILED (exit $EXIT_STATUS) at: $(date)"
fi
echo "=========================================="
exit $EXIT_STATUS
