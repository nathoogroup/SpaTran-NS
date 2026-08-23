#!/bin/bash
#SBATCH --job-name=spatran_tissue_finalize
#SBATCH --output=logs/spatran_tissue_finalize_%j.out
#SBATCH --error=logs/spatran_tissue_finalize_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
mkdir -p logs results

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20
export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"

# First require all 16 gene chunks for each refitted dataset and merge them into
# canonical per-dataset files. The other ten source objects have no in_tissue
# annotation, so their numerical fits are unchanged; validate and migrate only
# their metadata. Finally, require all twelve files to pass the strict combiner.
Rscript hpc/merge_gene_chunks.R data/spatial_datasets results \
  results/gene_chunks Visium_humanDLPFC 16
Rscript hpc/merge_gene_chunks.R data/spatial_datasets results \
  results/gene_chunks Visium_mouseCoronal 16
Rscript hpc/migrate_unannotated_results_to_schema4.R data/spatial_datasets results
Rscript hpc/combine_results.R results 12
