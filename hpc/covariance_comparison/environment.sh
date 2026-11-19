#!/usr/bin/env bash
# Source from the Slurm jobs (or a cluster login shell for direct R commands).
module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 \
  gdal/3.9.1 udunits/2.2.28 gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20
export R_LIBS="${CC_R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export INLA_NUM_THREADS=1
