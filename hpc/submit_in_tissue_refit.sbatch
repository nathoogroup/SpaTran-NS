#!/bin/bash
#SBATCH --job-name=spatran_tissue
#SBATCH --output=logs/spatran_tissue_%A_%a.out
#SBATCH --error=logs/spatran_tissue_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --array=11,12
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

set -euo pipefail

# Refit only the datasets whose Bioconductor objects contain explicit
# off-tissue spots. Each task uses all in-tissue spots jointly and parallelizes
# the independent gene fits over 32 cores.
#
#   task 11: Visium_humanDLPFC   (3,639/4,992 spots retained)
#   task 12: Visium_mouseCoronal (2,702/4,992 spots retained)

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
mkdir -p logs results

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20
export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"
export OVERWRITE=1

echo "Refitting in-tissue dataset task ${SLURM_ARRAY_TASK_ID} with ${SLURM_CPUS_PER_TASK} gene workers"
Rscript hpc/run_analysis.R data/spatial_datasets results "$SLURM_CPUS_PER_TASK"
