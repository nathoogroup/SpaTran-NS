#!/bin/bash
#SBATCH --job-name=download_spe
#SBATCH --output=logs/download_spe_%j.out
#SBATCH --error=logs/download_spe_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

set -euo pipefail

# =============================================================================
# Download SpatialExperiment Datasets - SBATCH Submission Script
# =============================================================================
# Submit with:
#   mkdir -p logs
#   sbatch hpc/submit_download.sbatch
#
# Run this BEFORE the analysis job.
# No INLA needed here; only Bioconductor packages.
# =============================================================================

echo "=========================================="
echo "Download job started at: $(date)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Node:    ${SLURM_NODELIST}"
echo "=========================================="

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
mkdir -p logs

# Download only needs ExperimentHub + SpatialExperiment — covered by the
# Bioconductor bundle module; no user library or INLA modules required.
module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            r-bundle-bioconductor/3.20

OUTPUT_DIR="data/spatial_datasets"
mkdir -p "$OUTPUT_DIR"

echo "Downloading SpatialExperiment datasets to: $OUTPUT_DIR"
if Rscript hpc/download_datasets.R "$OUTPUT_DIR"; then
    echo "=========================================="
    echo "Download completed successfully at: $(date)"
    echo "Next step: sbatch hpc/submit_analysis.sbatch"
else
    EXIT_STATUS=$?
    echo "==========================================" >&2
    echo "Download FAILED (exit ${EXIT_STATUS}) at: $(date)" >&2
    echo "==========================================" >&2
    exit "${EXIT_STATUS}"
fi
echo "=========================================="
