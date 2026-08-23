#!/bin/bash
set -euo pipefail

# Submit the two-dataset refit and finalize only if both array tasks succeed.
mkdir -p logs
refit_job=$(sbatch --parsable hpc/submit_in_tissue_refit.sbatch)
finalize_job=$(sbatch --parsable --dependency="afterok:${refit_job}" \
  hpc/finalize_in_tissue_refit.sbatch)

echo "In-tissue refit array job: ${refit_job}"
echo "Dependent migration/combine job: ${finalize_job}"
