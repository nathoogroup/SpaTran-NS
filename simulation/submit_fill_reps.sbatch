#!/bin/bash
#SBATCH --job-name=sim_fill
#SBATCH --output=logs/sim_%x_%A_%a.out
#SBATCH --error=logs/sim_%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --array=101,102
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Simulation Study: Fill Missing Replicates (ILC reps 14 and 100)
# =============================================================================
# Runs seeds 101 and 102 (different from the failed reps 14 and 100) and
# renames the output files to replicate_014.rds and replicate_100.rds so
# the aggregate script picks them up correctly.
#
# Usage:
#   sbatch simulation/submit_fill_reps.sbatch
# =============================================================================

DATASET="HumanBreastCancerILC"

# Map seed ID -> target replicate ID
if   [ "$SLURM_ARRAY_TASK_ID" -eq 101 ]; then TARGET_REP=14
elif [ "$SLURM_ARRAY_TASK_ID" -eq 102 ]; then TARGET_REP=100
else
    echo "ERROR: unexpected array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

TARGET_FILE=$(printf "simulation/output/%s/replicate_%03d.rds" "$DATASET" "$TARGET_REP")

echo "=========================================="
echo "Fill replicate at:        $(date)"
echo "Dataset:                  $DATASET"
echo "Seed replicate ID:        $SLURM_ARRAY_TASK_ID"
echo "Target replicate ID:      $TARGET_REP"
echo "Target file:              $TARGET_FILE"
echo "Node:                     $SLURM_NODELIST"
echo "CPUs:                     $SLURM_CPUS_PER_TASK"
echo "=========================================="

# Skip if target already exists
if [ -f "$TARGET_FILE" ]; then
    echo "Target file already exists, nothing to do."
    exit 0
fi

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

# Run with seed replicate ID (101 or 102) — output goes to replicate_101.rds / replicate_102.rds
Rscript simulation/run_replicate.R \
    "$DATASET" \
    "$SLURM_ARRAY_TASK_ID" \
    "$SLURM_CPUS_PER_TASK" \
    "simulation/configs" \
    "simulation/output"

EXIT_STATUS=$?
echo ""
echo "=========================================="

if [ $EXIT_STATUS -eq 0 ]; then
    echo "Simulation completed at: $(date)"

    SRC_RDS=$(printf "simulation/output/%s/replicate_%03d.rds" "$DATASET" "$SLURM_ARRAY_TASK_ID")
    SRC_LOG=$(printf "simulation/output/%s/replicate_%03d_gene_log.txt" "$DATASET" "$SLURM_ARRAY_TASK_ID")
    TGT_LOG=$(printf "simulation/output/%s/replicate_%03d_gene_log.txt" "$DATASET" "$TARGET_REP")

    mv "$SRC_RDS" "$TARGET_FILE"
    echo "Renamed: $SRC_RDS -> $TARGET_FILE"

    if [ -f "$SRC_LOG" ]; then
        mv "$SRC_LOG" "$TGT_LOG"
        echo "Renamed: $SRC_LOG -> $TGT_LOG"
    fi

    echo "Replicate $TARGET_REP filled at: $(date)"
else
    echo "Simulation FAILED (exit $EXIT_STATUS) at: $(date)"
fi
echo "=========================================="
exit $EXIT_STATUS
