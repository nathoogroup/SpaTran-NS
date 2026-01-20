#!/bin/bash
# =============================================================================
# Download SpatialExperiment Datasets from ExperimentHub
# =============================================================================
# Usage: ./download_datasets.sh [output_dir]
# =============================================================================

set -e  # Exit on error

# Set R library path
export R_LIBS=~/.local/R/$EBVERSIONR/

# Configuration
OUTPUT_DIR="${1:-data/spatial_datasets}"

echo "=========================================="
echo "SpatialExperiment Dataset Downloader"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Started at: $(date)"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Run download script
echo "Running download script..."
Rscript hpc/download_datasets.R "$OUTPUT_DIR"

echo ""
echo "=========================================="
echo "Download completed at: $(date)"
echo "Datasets saved to: $OUTPUT_DIR"
echo "=========================================="
