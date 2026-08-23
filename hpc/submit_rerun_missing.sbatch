#!/bin/bash
#SBATCH --job-name=spatran_rerun
#SBATCH --output=logs/spatran_rerun_%A_%a.out
#SBATCH --error=logs/spatran_rerun_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --array=11,12
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Targeted retry for datasets with explicit in_tissue annotations:
#   task 11 = Visium_humanDLPFC   (3,639 in-tissue spots)
#   task 12 = Visium_mouseCoronal (2,702 in-tissue spots)
#
# Uses 32 cores (vs 8 in original) to scale runtime from ~12h to ~3h.
# 5h walltime gives ~2h buffer.  OVERWRITE=1 forces replacement of stale files.
#
# Direct retry command (the dedicated launcher also schedules finalization):
#   sbatch hpc/submit_rerun_missing.sbatch
# =============================================================================

echo "=========================================="
echo "Array job started at:   $(date)"
echo "Array job ID:           $SLURM_ARRAY_JOB_ID"
echo "Array task ID:          $SLURM_ARRAY_TASK_ID"
echo "Node:                   $SLURM_NODELIST"
echo "CPUs:                   $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

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

export OVERWRITE=1
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
