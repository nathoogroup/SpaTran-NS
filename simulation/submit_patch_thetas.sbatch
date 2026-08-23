#!/bin/bash
#SBATCH --job-name=patch_thetas
#SBATCH --output=logs/patch_thetas_%j.out
#SBATCH --error=logs/patch_thetas_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Fit nonstationary model for H1 genes and patch theta1_ns–theta4_ns into
# the existing results RDS, then regenerate the simulation config.
#
# Usage:
#   sbatch --export=DATASET=HumanBreastCancerILC simulation/submit_patch_thetas.sbatch
# =============================================================================

if [ -z "$DATASET" ]; then
    echo "ERROR: DATASET not set. Use --export=DATASET=<name>"
    exit 1
fi

echo "=========================================="
echo "Patch NS thetas started at: $(date)"
echo "Dataset:  $DATASET"
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

# Step 1: Fit NS thetas for H1 genes and patch results RDS
echo ""
echo "Step 1: Fitting NS thetas for H1 genes..."
Rscript simulation/patch_ns_thetas.R "$DATASET" "$SLURM_CPUS_PER_TASK"

if [ $? -ne 0 ]; then
    echo "ERROR: patch_ns_thetas.R failed"
    exit 1
fi

# Step 2: Regenerate simulation config with fitted thetas
echo ""
echo "Step 2: Regenerating simulation config..."
Rscript simulation/prepare_simulation.R "$DATASET"

if [ $? -ne 0 ]; then
    echo "ERROR: prepare_simulation.R failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "Patch complete at: $(date)"
echo "Submitting simulation replicates..."
sbatch --export=DATASET="$DATASET" simulation/submit_simulation.sbatch
echo "=========================================="
