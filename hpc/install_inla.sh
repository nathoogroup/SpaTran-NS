#!/bin/bash
# =============================================================================
# R-INLA Installation Script for Alliance Canada Clusters
# =============================================================================
# Source: https://docs.alliancecan.ca/wiki/R-INLA
#
# Run this script ONCE on a login node (not in a job) to install INLA.
# After installation, the same modules must be loaded in every job script.
#
# Usage:
#   bash hpc/install_inla.sh
# =============================================================================

set -e

# (1) Load required modules
module load StdEnv/2023 gcc/12.3 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 gsl/2.7 jags/4.3.2

LOGFILE=r_INLA_install_${EBVERSIONR}_${CC_CLUSTER}_$(date --iso=min).log
echo "Installation log: $LOGFILE"

# (2) Set R library path and install INLA package with dependencies
export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"
echo "R_LIBS is $R_LIBS"
mkdir -p "$R_LIBS"

R -e 'install.packages("remotes", repos=c("https://mirror.csclub.uwaterloo.ca/CRAN/"))'

R -e 'install.packages("INLA",
      repos=c("https://mirror.csclub.uwaterloo.ca/CRAN/",
              INLA="https://inla.r-inla-download.org/R/stable"),
      dep=TRUE, Ncpus=2)' \
    |& tee "$LOGFILE"

# (3) Install pre-compiled INLA binaries for Rocky Linux 8
R -e 'library(INLA); inla.binary.install(os="Rocky Linux-8")' |& tee -a "$LOGFILE"

# (4) Patch binaries for Alliance Canada library paths
chmod u+x "$R_LIBS/INLA/bin/linux/64bit/"*.so.* \
          "$R_LIBS/INLA/bin/linux/64bit/"*.so \
          "$R_LIBS/INLA/bin/linux/64bit/first/"*.so

sed -i 's/\(^.*export LD_LIBRARY_PATH\)/echo "Skipping LD_LIBRARY_PATH." #\1/g' \
    "$R_LIBS/INLA/bin/linux/64bit/"*.run

setrpaths.sh --path "$R_LIBS/INLA/bin/linux/64bit/malloc" \
    --add_path "\$ORIGIN:$EBROOTGENTOO/lib/gcc/x86_64-pc-linux-gnu/${EBVERSIONGCC::2}"

setrpaths.sh --path "$R_LIBS/INLA/bin/linux/64bit" \
    --add_path '$ORIGIN/first:$ORIGIN:$ORIGIN/malloc'

echo ""
echo "=== INLA installation complete ==="
echo "Log saved to: $LOGFILE"
echo ""
echo "Test with:"
echo "  module load StdEnv/2023 gcc/12.3 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 gsl/2.7 jags/4.3.2"
echo "  export R_LIBS=\$HOME/R/x86_64-pc-linux-gnu-library/4.4"
echo "  R -e 'library(INLA); inla.mesh.1d(1:5); cat(\"INLA OK\\n\")'"
