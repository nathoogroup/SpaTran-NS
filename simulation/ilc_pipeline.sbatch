#!/bin/bash
#SBATCH --job-name=ilc_pipeline
#SBATCH --output=logs/ilc_pipeline_%j.out
#SBATCH --error=logs/ilc_pipeline_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# ILC Pipeline: Fit → Check Thetas → Submit Simulation
# =============================================================================
# Step 1: Run INLA fit for HumanBreastCancerILC (dataset index 2)
# Step 2: Verify theta3_ns / theta4_ns are present in results
# Step 3: If thetas look good, submit the simulation array job
# =============================================================================

set -e

echo "=========================================="
echo "ILC Pipeline at: $(date)"
echo "Node:            $SLURM_NODELIST"
echo "CPUs:            $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

DATA_DIR="data/spatial_datasets"
OUTPUT_DIR="results"
N_CORES=$SLURM_CPUS_PER_TASK

# =============================================================================
# STEP 1: Run analysis fit for ILC (array task 2)
# =============================================================================
echo ""
echo "--- STEP 1: Running INLA fit for HumanBreastCancerILC ---"

OVERWRITE=1 SLURM_ARRAY_TASK_ID=2 \
    Rscript hpc/run_analysis.R "$DATA_DIR" "$OUTPUT_DIR" "$N_CORES"

echo ""
echo "Fit step completed at: $(date)"

# =============================================================================
# STEP 2: Check thetas in results
# =============================================================================
echo ""
echo "--- STEP 2: Checking theta coverage in results ---"

RESULTS_FILE="$OUTPUT_DIR/HumanBreastCancerILC_results.rds"

if [ ! -f "$RESULTS_FILE" ]; then
    echo "ERROR: Results file not found: $RESULTS_FILE"
    exit 1
fi

THETA_CHECK=$(Rscript --vanilla -e "
r <- readRDS('$RESULTS_FILE')
n_total  <- nrow(r)
n_theta3 <- sum(!is.na(r\$theta3_ns))
n_theta4 <- sum(!is.na(r\$theta4_ns))
pct <- round(100 * n_theta3 / max(n_total, 1), 1)
cat(sprintf('Genes: %d | theta3_ns non-NA: %d | theta4_ns non-NA: %d | Coverage: %.1f%%\n',
            n_total, n_theta3, n_theta4, pct))
cat(sprintf('PASS=%s\n', ifelse(n_theta3 > 0 && n_theta4 > 0, 'yes', 'no')))
")

echo "$THETA_CHECK"

PASS=$(echo "$THETA_CHECK" | grep "^PASS=" | cut -d= -f2)

if [ "$PASS" != "yes" ]; then
    echo "ERROR: theta3_ns / theta4_ns not found or all NA in results. Aborting simulation."
    exit 1
fi

# =============================================================================
# STEP 3: Submit simulation array job
# =============================================================================
echo ""
echo "--- STEP 3: Submitting simulation array job for HumanBreastCancerILC ---"

SIM_JOB=$(sbatch --parsable \
    --export=DATASET=HumanBreastCancerILC \
    simulation/submit_simulation.sbatch)

echo "Simulation job submitted: $SIM_JOB"
echo ""
echo "=========================================="
echo "Pipeline complete at: $(date)"
echo "Fit results: $RESULTS_FILE"
echo "Simulation job ID: $SIM_JOB"
echo "=========================================="
