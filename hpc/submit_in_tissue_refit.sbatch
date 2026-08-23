#!/bin/bash
#SBATCH --job-name=spatran_tissue
#SBATCH --output=logs/spatran_tissue_%A_%a.out
#SBATCH --error=logs/spatran_tissue_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-32%32
#SBATCH --account=def-nathoo
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=pvelidi@uvic.ca

set -euo pipefail

# Refit only the datasets whose Bioconductor objects contain explicit
# off-tissue spots. Each dataset is split into 16 disjoint round-robin gene
# chunks, giving 32 interleaved array tasks in total. Every individual gene fit
# still uses all in-tissue spots jointly.
#
#   odd tasks: Visium_humanDLPFC   (dataset index 11; 3,639 spots)
#   even tasks: Visium_mouseCoronal (dataset index 12; 2,702 spots)

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
mkdir -p logs results

module load StdEnv/2023 gcc/12.3 boost/1.85.0 r/4.4.0 geos/3.12.0 gdal/3.9.1 udunits/2.2.28 \
            gsl/2.7 jags/4.3.2 r-bundle-bioconductor/3.20
export R_LIBS="$HOME/R/x86_64-pc-linux-gnu-library/${EBVERSIONR:0:3}"
export OVERWRITE=1
export GENE_CHUNK_COUNT=16
export GENE_CHUNK_OUTPUT_DIR=results/gene_chunks

task_offset=$((SLURM_ARRAY_TASK_ID - 1))
dataset_offset=$((task_offset % 2))
export ANALYSIS_DATASET_INDEX=$((11 + dataset_offset))
export GENE_CHUNK_ID=$((task_offset / 2 + 1))

dataset_names=(Visium_humanDLPFC Visium_mouseCoronal)
dataset_name=${dataset_names[$dataset_offset]}

echo "Refitting ${dataset_name}, gene chunk ${GENE_CHUNK_ID}/${GENE_CHUNK_COUNT}, with ${SLURM_CPUS_PER_TASK} workers"
Rscript hpc/run_analysis.R data/spatial_datasets results "$SLURM_CPUS_PER_TASK"
