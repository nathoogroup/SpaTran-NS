#!/bin/bash
# =============================================================================
# Simulation Study: Prepare Configurations (Login Node)
# =============================================================================
# Run on the login node before launching simulation replicates.
# Creates simulation config RDS files in simulation/configs/.
#
# Usage: bash simulation/prepare_simulation.sh
#    or: ./simulation/prepare_simulation.sh
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

echo "=========================================="
echo "Prepare simulation configs at: $(date)"
echo "Project root: $PROJECT_ROOT"
echo "=========================================="

# Load modules (Alliance Canada; no-op if module not available)
if command -v module &>/dev/null; then
  module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
              gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20 2>/dev/null || true
  export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
fi

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
