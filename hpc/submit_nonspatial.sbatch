#!/bin/bash
#SBATCH --job-name=spatran_null
#SBATCH --output=logs/spatran_null_%A_%a.out
#SBATCH --error=logs/spatran_null_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-13
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

# =============================================================================
# Patch: Non-spatial (intercept + nugget) model — Array Job
# =============================================================================
# Appends log_ml_nonspatial, sigma_eps_sq_nonspatial, mu_nonspatial to each
# per-dataset results RDS.  Must be run after the main analysis is complete.
#
# Runtime estimate (per task):
#   ~16,000 genes / 8 cores × ~0.5 sec/gene ≈ 17 min → 3h is generous
#
# Submit: sbatch hpc/submit_nonspatial.sbatch
# Rerun:  sbatch --export=OVERWRITE=1 hpc/submit_nonspatial.sbatch
# =============================================================================

echo "=========================================="
echo "Non-spatial patch started at: $(date)"
echo "Array job ID:                 $SLURM_ARRAY_JOB_ID"
echo "Array task ID:                $SLURM_ARRAY_TASK_ID"
echo "Node:                         $SLURM_NODELIST"
echo "CPUs:                         $SLURM_CPUS_PER_TASK"
echo "=========================================="

mkdir -p logs

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

cd "$SLURM_SUBMIT_DIR"

export OVERWRITE=${OVERWRITE:-0}
Rscript hpc/patch_nonspatial.R data/spatial_datasets results $SLURM_CPUS_PER_TASK

EXIT_STATUS=$?
echo ""
echo "=========================================="
if [ $EXIT_STATUS -eq 0 ]; then
    echo "Task $SLURM_ARRAY_TASK_ID completed at: $(date)"
else
    echo "Task $SLURM_ARRAY_TASK_ID FAILED (exit $EXIT_STATUS) at: $(date)"
fi
echo "=========================================="
exit $EXIT_STATUS
