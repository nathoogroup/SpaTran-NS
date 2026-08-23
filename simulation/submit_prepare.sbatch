#!/bin/bash
#SBATCH --job-name=sim_prep
#SBATCH --output=logs/sim_prepare_%j.out
#SBATCH --error=logs/sim_prepare_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Simulation Study: Prepare Configurations for Both Datasets
# =============================================================================
# Run once before launching simulation replicates.
# Creates simulation config RDS files in simulation/configs/.
#
# Submit: sbatch simulation/submit_prepare.sbatch
# =============================================================================

echo "=========================================="
echo "Prepare simulation configs at: $(date)"
echo "Node: $SLURM_NODELIST"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

echo ""
echo "--- HumanBreastCancerILC ---"
Rscript simulation/prepare_simulation.R HumanBreastCancerILC

echo ""
echo "--- HumanOvarianCancer ---"
Rscript simulation/prepare_simulation.R HumanOvarianCancer

echo ""
echo "=========================================="
echo "Preparation complete at: $(date)"
echo "=========================================="
