#!/bin/bash
#SBATCH --job-name=download_spe
#SBATCH --output=logs/download_spe_%j.out
#SBATCH --error=logs/download_spe_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=standard
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com

# =============================================================================
# Download SpatialExperiment Datasets - SBATCH Submission Script
# =============================================================================
# Submit with: sbatch submit_download.sbatch
#
# This downloads all SpatialExperiment datasets from ExperimentHub
# Run this BEFORE running the main analysis
# =============================================================================

echo "=========================================="
echo "Download job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "=========================================="

# Create logs directory
mkdir -p logs

# Set up environment
# Load R module (modify version as needed for your HPC)
module load StdEnv/2023 r/4.5.0

# Set R library path
export R_LIBS=~/.local/R/$EBVERSIONR/

# Set working directory
cd $SLURM_SUBMIT_DIR

# Configuration
OUTPUT_DIR="data/spatial_datasets"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Downloading datasets to: $OUTPUT_DIR"
echo ""

# Run download script
Rscript hpc/download_datasets.R "$OUTPUT_DIR"

# Check exit status
EXIT_STATUS=$?

echo ""
echo "=========================================="
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Download completed successfully at: $(date)"
    echo "Datasets saved to: $OUTPUT_DIR"
    echo ""
    echo "To run the analysis, submit:"
    echo "  sbatch hpc/submit_analysis.sbatch"
else
    echo "Download failed with exit code $EXIT_STATUS at: $(date)"
fi
echo "=========================================="

exit $EXIT_STATUS
