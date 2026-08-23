#!/bin/bash
#SBATCH --job-name=sim_rep
#SBATCH --output=logs/sim_%x_%A_%a.out
#SBATCH --error=logs/sim_%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --array=1-100
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Simulation Study: Run Replicates (Array Job)
# =============================================================================
# Each array task runs one simulation replicate for a given dataset.
# SLURM_ARRAY_TASK_ID is used as the replicate ID.
#
# Usage:
#   sbatch --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch
#   sbatch --export=DATASET=HumanOvarianCancer   simulation/submit_simulation.sbatch
#
# To run fewer replicates for testing:
#   sbatch --array=1-5 --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch
#
# Runtime estimate per replicate (simulate + fit, 32 cores):
#   ILC:          ~14 sec/gene × 16564 genes / 32 cores ≈ 120 min
#   OvarianCancer:~120 sec/gene × 16684 genes / 32 cores ≈ 17h (larger dataset, more spots)
#   Per-gene 120s timeout caps stuck INLA calls.
# =============================================================================

if [ -z "$DATASET" ]; then
    echo "ERROR: DATASET not set. Use --export=DATASET=<name>"
    echo "  Example: sbatch --export=DATASET=HumanBreastCancerILC simulation/submit_simulation.sbatch"
    exit 1
fi

echo "=========================================="
echo "Simulation replicate at:  $(date)"
echo "Dataset:                  $DATASET"
echo "Array job ID:             $SLURM_ARRAY_JOB_ID"
echo "Replicate (array task):   $SLURM_ARRAY_TASK_ID"
echo "Node:                     $SLURM_NODELIST"
echo "CPUs:                     $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

N_CORES=$SLURM_CPUS_PER_TASK

Rscript simulation/run_replicate.R \
    "$DATASET" \
    "$SLURM_ARRAY_TASK_ID" \
    "$N_CORES" \
    "simulation/configs" \
    "simulation/output"

EXIT_STATUS=$?
echo ""
echo "=========================================="
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Replicate $SLURM_ARRAY_TASK_ID completed at: $(date)"
else
    echo "Replicate $SLURM_ARRAY_TASK_ID FAILED (exit $EXIT_STATUS) at: $(date)"
fi
echo "=========================================="
exit $EXIT_STATUS
