#!/bin/bash

# =============================================================================
# Simulation Study: Prepare Configurations for Both Datasets
# =============================================================================
# Run once before launching simulation replicates.
# Creates simulation config RDS files in simulation/configs/.
#
# Submit: bash simulation/prepare.sh
# =============================================================================

mkdir -p logs

module load StdEnv/2023 gcc/12.3 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 boost/1.85.0 r-bundle-bioconductor/3.20 

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

echo ""
echo "--- HumanBreastCancerILC ---"
cd $HOME/projects/def-nathoo/pvelidi/SpaTran-NS

Rscript simulation/prepare_simulation.R HumanBreastCancerILC

echo ""
echo "--- HumanOvarianCancer ---"
Rscript simulation/prepare_simulation.R HumanOvarianCancer

echo ""
echo "=========================================="
echo "Preparation complete at: $(date)"
echo "=========================================="
