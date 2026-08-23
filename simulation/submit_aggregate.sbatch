#!/bin/bash
#SBATCH --job-name=sim_agg
#SBATCH --output=logs/sim_aggregate_%j.out
#SBATCH --error=logs/sim_aggregate_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Simulation Study: Aggregate Results
# =============================================================================
# Run after all replicates complete for both datasets.
# Plotting is done locally, not on the cluster.
#
# Submit: sbatch simulation/submit_aggregate.sbatch
# =============================================================================

echo "=========================================="
echo "Aggregation started at: $(date)"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

# Ensure patchwork for plotting
Rscript -e 'if (!requireNamespace("patchwork", quietly=TRUE)) install.packages("patchwork", repos="https://cloud.r-project.org")' 2>/dev/null || true

echo ""
echo "--- Aggregating HumanBreastCancerILC ---"
Rscript simulation/aggregate_results.R HumanBreastCancerILC 100

echo ""
echo "--- Aggregating HumanOvarianCancer ---"
Rscript simulation/aggregate_results.R HumanOvarianCancer 100

echo ""
echo "=========================================="
echo "Aggregation complete at: $(date)"
echo "Download results and run plotting locally."
echo "=========================================="
