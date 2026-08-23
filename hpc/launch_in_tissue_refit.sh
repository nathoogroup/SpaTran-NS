#!/bin/bash
set -euo pipefail

# Submit 16 gene chunks for each of two datasets and finalize only if all 32
# interleaved array tasks succeed.
mkdir -p logs
refit_job=$(sbatch --parsable hpc/submit_in_tissue_refit.sbatch)
finalize_job=$(sbatch --parsable --dependency="afterok:${refit_job}" \
  hpc/finalize_in_tissue_refit.sbatch)

echo "In-tissue gene-chunk array job: ${refit_job} (32 tasks)"
echo "Dependent migration/combine job: ${finalize_job}"
