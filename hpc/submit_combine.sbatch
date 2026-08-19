#!/bin/bash
#SBATCH --job-name=spatran_combine
#SBATCH --output=logs/spatran_combine_%j.out
#SBATCH --error=logs/spatran_combine_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

set -euo pipefail

# Submit this job with an afterok dependency on the full analysis array:
#   analysis_job=$(sbatch --parsable hpc/submit_analysis.sbatch)
#   sbatch --dependency="afterok:${analysis_job}" hpc/submit_combine.sbatch

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20

export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"
export N_EXPECTED_DATASETS="${N_EXPECTED_DATASETS:-12}"

echo "Combining ${N_EXPECTED_DATASETS} validated dataset results at: $(date)"
if Rscript hpc/combine_results.R results "${N_EXPECTED_DATASETS}"; then
    echo "Combine job completed successfully at: $(date)"
else
    EXIT_STATUS=$?
    echo "Combine job FAILED (exit ${EXIT_STATUS}) at: $(date)" >&2
    exit "${EXIT_STATUS}"
fi
