#!/bin/bash
# =============================================================================
# Simulation Study: Aggregate Results
# =============================================================================
# Run after all replicates complete for both datasets.
# Plotting is done locally, not on the cluster.
#
# Submit: bash simulation/aggregate.sh
# =============================================================================

echo "=========================================="
echo "Aggregation started at: $(date)"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 boost/1.85.0 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"


echo ""
echo "--- Aggregating HumanBreastCancerILC ---"
cd $HOME/projects/def-nathoo/pvelidi/SpaTran-NS
Rscript simulation/aggregate_results.R HumanBreastCancerILC 100

echo ""
echo "--- Aggregating HumanOvarianCancer ---"
cd $HOME/projects/def-nathoo/pvelidi/SpaTran-NS
Rscript simulation/aggregate_results.R HumanOvarianCancer 100

echo ""
echo "=========================================="
echo "Aggregation complete at: $(date)"
echo "Download results and run plotting locally."
echo "=========================================="
