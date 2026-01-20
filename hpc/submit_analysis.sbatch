#!/bin/bash
#SBATCH --job-name=spatran_ns
#SBATCH --output=logs/spatran_ns_%j.out
#SBATCH --error=logs/spatran_ns_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com

# =============================================================================
# Spatial Nonstationarity Analysis - SBATCH Submission Script
# =============================================================================
# Submit with: sbatch submit_analysis.sbatch
#
# Modify the following as needed:
#   --time          : Wall clock time limit
#   --cpus-per-task : Number of cores (should match n_cores in R script)
#   --mem           : Memory allocation
#   --partition     : Your HPC partition name
#   --mail-user     : Your email for notifications
# =============================================================================

echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="

# Create logs directory if it doesn't exist
mkdir -p logs

# Set up environment
# Load R module (modify version as needed for your HPC)
module load StdEnv/2023 r/4.4.0

# Set R library path
export R_LIBS=~/.local/R/$EBVERSIONR/

# Set working directory
cd $SLURM_SUBMIT_DIR

# Configuration
DATA_DIR="data/spatial_datasets"
OUTPUT_DIR="results/$(date +%Y%m%d_%H%M%S)"
N_CORES=$SLURM_CPUS_PER_TASK
GENES_PER_DATASET="all"  # Use "all" or a number like "1000" for top N genes

# Create output directory
mkdir -p $OUTPUT_DIR

echo ""
echo "Configuration:"
echo "  Data directory: $DATA_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Number of cores: $N_CORES"
echo "  Genes per dataset: $GENES_PER_DATASET"
echo ""

# Run analysis
echo "Starting R analysis..."
Rscript hpc/run_analysis.R "$DATA_DIR" "$OUTPUT_DIR" "$N_CORES" "$GENES_PER_DATASET"

# Check exit status
EXIT_STATUS=$?

echo ""
echo "=========================================="
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Job completed successfully at: $(date)"
else
    echo "Job failed with exit code $EXIT_STATUS at: $(date)"
fi
echo "=========================================="

exit $EXIT_STATUS
