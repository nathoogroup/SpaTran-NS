#!/bin/bash
# =============================================================================
# Simulation Study: Aggregate Results (Login Node)
# =============================================================================
# Run on the login node after all replicates complete for both datasets.
# Computes FDR/power curves and creates figures.
#
# Usage: bash simulation/aggregate_simulation.sh
#    or: ./simulation/aggregate_simulation.sh
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

echo "=========================================="
echo "Aggregation started at: $(date)"
echo "Project root: $PROJECT_ROOT"
echo "=========================================="

# Load modules (Alliance Canada; no-op if module not available)
if command -v module &>/dev/null; then
  module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
              gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20 2>/dev/null || true
  export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
fi

# Ensure patchwork for plotting
Rscript -e 'if (!requireNamespace("patchwork", quietly=TRUE)) install.packages("patchwork", repos="https://cloud.r-project.org")' 2>/dev/null || true

echo ""
echo "--- Aggregating HumanBreastCancerILC ---"
Rscript simulation/aggregate_results.R HumanBreastCancerILC 100

echo ""
echo "--- Aggregating HumanOvarianCancer ---"
Rscript simulation/aggregate_results.R HumanOvarianCancer 100

echo ""
echo "--- Creating figure ---"
Rscript simulation/plot_fdr_power.R

echo ""
echo "=========================================="
echo "Aggregation complete at: $(date)"
echo "=========================================="
