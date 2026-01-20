# HPC Scripts for Spatial Nonstationarity Analysis

This directory contains scripts for running the spatial nonstationarity analysis on an HPC cluster using SLURM.

## Overview

The analysis:
1. Downloads all SpatialExperiment datasets from ExperimentHub
2. For each dataset, analyzes all genes (or top N expressed genes)
3. Fits stationary and nonstationary Matérn models using INLA
4. Calculates Bayes factors comparing the models
5. Saves results for downstream analysis

## Files

| File | Description |
|------|-------------|
| `download_datasets.R` | Downloads SpatialExperiment datasets from ExperimentHub |
| `analysis_functions.R` | Helper functions for the nonstationarity analysis |
| `run_analysis.R` | Main analysis script with doParallel support |
| `submit_download.sbatch` | SBATCH script to submit the download job |
| `submit_analysis.sbatch` | SBATCH script to submit the analysis job |

## Prerequisites

### R Packages
```r
install.packages(c("doParallel", "foreach", "dplyr"))

# Bioconductor packages
BiocManager::install(c("SpatialExperiment", "ExperimentHub", "SummarizedExperiment"))

# INLA (from INLA repository)
install.packages("INLA", 
                 repos = c(getOption("repos"), 
                          INLA = "https://inla.r-inla-download.org/R/testing"), 
                 type = "binary")
```

## Usage

### Step 1: Download Datasets

First, download all SpatialExperiment datasets:

```bash
# Interactive
Rscript hpc/download_datasets.R data/spatial_datasets

# Or submit to SLURM
sbatch hpc/submit_download.sbatch
```

### Step 2: Run Analysis

After datasets are downloaded:

```bash
# Interactive (for testing)
Rscript hpc/run_analysis.R data/spatial_datasets results 8 100

# Submit to SLURM
sbatch hpc/submit_analysis.sbatch
```

### Command Line Arguments

**download_datasets.R:**
```
Rscript download_datasets.R [output_dir]

  output_dir : Directory to save datasets (default: data/spatial_datasets)
```

**run_analysis.R:**
```
Rscript run_analysis.R [data_dir] [output_dir] [n_cores] [genes_per_dataset]

  data_dir          : Directory with .rds SpatialExperiment files
  output_dir        : Directory for results
  n_cores           : Number of parallel cores (default: all available)
  genes_per_dataset : "all" or a number like "1000" for top N genes
```

## SBATCH Configuration

Edit the SBATCH scripts to match your HPC environment:

```bash
#SBATCH --time=48:00:00        # Wall clock time
#SBATCH --cpus-per-task=32     # Number of cores
#SBATCH --mem=128G             # Memory
#SBATCH --partition=standard   # Your partition name
#SBATCH --mail-user=you@email  # Your email
```

Also uncomment/modify the module loading section:
```bash
# module load R/4.4.0
# source activate your_conda_env
```

## Output

Results are saved to the output directory:

```
results/
├── Dataset1_results.rds       # Per-dataset results
├── Dataset2_results.rds
├── ...
├── all_results_combined.rds   # Combined results (R object)
└── all_results_combined.csv   # Combined results (CSV)
```

### Results Columns

| Column | Description |
|--------|-------------|
| `gene_name` | Gene symbol/name |
| `sigma_b_sq_stationary` | Spatial variance (stationary model) |
| `sigma_eps_sq_stationary` | Noise variance (stationary model) |
| `range_stationary` | Spatial range (stationary model) |
| `prop_spatial_stationary` | Proportion of spatial variance |
| `log_ml_stationary` | Log marginal likelihood (stationary) |
| `sigma_b_sq_nonstationary` | Spatial variance (nonstationary model) |
| `sigma_eps_sq_nonstationary` | Noise variance (nonstationary model) |
| `range_nonstationary` | Spatial range (nonstationary model) |
| `prop_spatial_nonstationary` | Proportion of spatial variance |
| `log_ml_nonstationary` | Log marginal likelihood (nonstationary) |
| `log_bayes_factor` | log(BF_10) = log(p(Y|M1)/p(Y|M0)) |
| `bayes_factor` | BF_10 on original scale |
| `bf_interpretation` | Wagenmakers/Jeffreys interpretation |
| `dataset_name` | Source dataset name |

### Bayes Factor Interpretation

| log(BF) Range | Evidence |
|---------------|----------|
| 0 - 1.1 | Anecdotal |
| 1.1 - 2.3 | Moderate |
| 2.3 - 3.4 | Strong |
| 3.4 - 4.6 | Very strong |
| > 4.6 | Extreme |

Positive log(BF) favors nonstationary; negative favors stationary.

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View output logs
tail -f logs/spatran_ns_JOBID.out

# Cancel a job
scancel JOBID
```

## Troubleshooting

### INLA not working
Make sure INLA binary is available. On some systems you may need:
```r
install.packages("INLA", 
                 repos = c(getOption("repos"), 
                          INLA = "https://inla.r-inla-download.org/R/testing"), 
                 type = "binary")
```

### Out of memory
- Reduce `genes_per_dataset` 
- Increase `--mem` in SBATCH
- Process fewer datasets per job

### Job timeout
- Increase `--time` in SBATCH
- Reduce `genes_per_dataset`
- Split into multiple jobs
