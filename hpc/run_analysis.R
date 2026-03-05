#!/usr/bin/env Rscript
# =============================================================================
# Spatial Nonstationarity Analysis - HPC Script with Parallel Processing
# =============================================================================
# For each 10x Visium SpatialExperiment dataset:
#   1. Filter lowly-expressed genes
#   2. Normalise (library-size + log2)
#   3. Fit stationary and non-stationary Matern models via INLA SPDE for ALL genes
#   4. Compute Bayes factor for non-stationary covariance kernel
#
# Usage:
#   Rscript hpc/run_analysis.R [data_dir] [output_dir] [n_cores]
#
#   data_dir   : directory with .rds SpatialExperiment files (from download_datasets.R)
#   output_dir : output directory for results
#   n_cores    : parallel cores (default: all available)
# =============================================================================

# When run as a SLURM array job, SLURM_ARRAY_TASK_ID selects the dataset index.
# When run directly, data_dir is processed sequentially (useful for testing).
args       <- commandArgs(trailingOnly = TRUE)
data_dir   <- if (length(args) >= 1) args[1] else "data/spatial_datasets"
output_dir <- if (length(args) >= 2) args[2] else "results"
n_cores    <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores()
task_id    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "hpc")

cat("=============================================================\n")
cat("Spatial Nonstationarity Analysis - HPC Parallel Processing\n")
cat("=============================================================\n")
cat("Data directory:   ", data_dir, "\n")
cat("Output directory: ", output_dir, "\n")
cat("Cores:            ", n_cores, "\n")
cat("Gene selection:    all expressed genes (after quality filter)\n")
cat("SLURM array task: ", if (task_id > 0) task_id else "N/A (sequential mode)", "\n")
cat("Start time:       ", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# Load packages
cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(INLA)
  library(dplyr)
  library(doParallel)
  library(foreach)
})

# Source analysis functions
cat("Loading analysis functions...\n")
source(file.path(script_dir, "analysis_functions.R"))

# Verify INLA
cat("Testing INLA...\n")
tryCatch({
  INLA::inla.mesh.1d(seq(0, 1, length.out = 5))
  cat("INLA OK\n\n")
}, error = function(e) {
  stop("INLA not working: ", e$message,
       "\nRun bash hpc/install_inla.sh first.")
})

# Find datasets — exclude log/manifest files and _v3_13 duplicates.
# The _v3_13 variants are re-releases of the same underlying tissue samples;
# keeping only one version avoids duplicated computation.
dataset_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
dataset_files <- dataset_files[!grepl("(download_log|manifest)", basename(dataset_files))]
dataset_files <- dataset_files[!grepl("_v3_13", basename(dataset_files))]
dataset_files <- sort(dataset_files)
cat("Found", length(dataset_files), "unique dataset files (v3.13 duplicates excluded)\n")

# In array mode, pick the single dataset for this task
if (task_id > 0) {
  if (task_id > length(dataset_files)) {
    stop("SLURM_ARRAY_TASK_ID (", task_id, ") exceeds number of datasets (",
         length(dataset_files), ")")
  }
  cat("Array mode: processing dataset", task_id, "of", length(dataset_files), "\n")
  dataset_files <- dataset_files[task_id]
}
cat("\n")

# Setup parallel backend
cat("Setting up", n_cores, "parallel workers...\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, suppressPackageStartupMessages({
  library(Matrix)
  library(INLA)
  library(dplyr)
  # Prevent INLA from spawning its own threads — we are already parallelising
  # externally, so each worker should use a single thread.
  INLA::inla.setOption(num.threads = "1:1")
}))

clusterExport(cl, c(
  "normalize_expression", "detrend_expression",
  "create_spde_mesh",
  "fit_stationary_matern", "fit_nonstationary_matern",
  "extract_variance_components", "calculate_bayes_factor",
  "analyze_gene_nonstationarity", "analyze_gene_safe"
))

# ---- Per-dataset loop -------------------------------------------------------

all_results <- list()

for (ds_idx in seq_along(dataset_files)) {
  ds_file <- dataset_files[ds_idx]
  ds_name <- tools::file_path_sans_ext(basename(ds_file))

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Dataset", ds_idx, "/", length(dataset_files), ":", ds_name, "\n")
  cat(rep("=", 60), "\n")

  # Skip if results already exist (resumable)
  ds_out <- file.path(output_dir, paste0(ds_name, "_results.rds"))
  if (file.exists(ds_out)) {
    cat("  Results already exist, skipping (delete to rerun)\n")
    all_results[[ds_name]] <- readRDS(ds_out)
    next
  }

  # ---- Load ----------------------------------------------------------------
  spe <- tryCatch(readRDS(ds_file), error = function(e) {
    cat("  ERROR loading:", e$message, "\n"); NULL
  })
  if (is.null(spe) || !inherits(spe, "SpatialExperiment")) {
    cat("  Skipping (not a SpatialExperiment)\n"); next
  }

  n_genes_total <- nrow(spe)
  n_spots       <- ncol(spe)
  cat(sprintf("  Genes: %d, Spots: %d\n", n_genes_total, n_spots))

  # ---- Extract counts (once, with error handling for HDF5-backed objects) --
  # Some datasets (e.g. Janesick) are stored as HDF5-backed DelayedArrays;
  # their .rds is tiny and the actual data lives in a separate HDF5 file that
  # may be absent.  We catch that failure here and skip the dataset.
  counts_mat <- tryCatch(counts(spe), error = function(e) {
    cat("  ERROR reading counts (HDF5-backed assay with missing backing file?):\n  ",
        e$message, "\n")
    NULL
  })
  if (is.null(counts_mat)) next

  # ---- Filter lowly-expressed genes ----------------------------------------
  min_spots  <- ceiling(n_spots * 0.005)  # detected in >= 0.5% of spots
  keep       <- (Matrix::rowSums(counts_mat > 0) >= min_spots) &
                (Matrix::rowSums(counts_mat) >= 3)
  counts_mat <- counts_mat[keep, ]   # subset the matrix directly
  spe        <- spe[keep, ]
  cat(sprintf("  After quality filter: %d genes (analyse all)\n", nrow(spe)))

  # All filtered genes
  gene_indices <- seq_len(nrow(spe))

  # ---- Coordinates and gene names ------------------------------------------
  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords)[1:2] <- c("x", "y")

  rd         <- rowData(spe)
  gene_names <- if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
                else if ("symbol"    %in% colnames(rd)) as.character(rd$symbol)
                else rownames(spe)

  # Pre-build the SPDE mesh once; all workers share the same mesh object
  # rather than rebuilding it for every gene.
  cat("  Building shared SPDE mesh...\n")
  mesh <- create_spde_mesh(coords)
  cat(sprintf("  Mesh: %d nodes\n", mesh$n))

  # Free the SPE before the cluster export — this is the main memory saving.
  # coords, gene_names, counts_mat, and mesh are all we need from here on.
  rm(spe); gc()

  # ---- Parallel analysis ---------------------------------------------------
  clusterExport(cl, c("counts_mat", "coords", "gene_names", "mesh"),
                envir = environment())

  cat("  Running parallel INLA fits for", length(gene_indices), "genes...\n")
  t0 <- proc.time()

  results_list <- foreach(
    i              = seq_along(gene_indices),
    .combine       = rbind,
    .errorhandling = "pass",
    .packages      = c("Matrix", "INLA", "dplyr")
  ) %dopar% {
    analyze_gene_safe(
      gene_idx       = gene_indices[i],
      counts_mat     = counts_mat,
      coords         = coords,
      gene_names     = gene_names,
      normalize      = TRUE,
      detrend        = TRUE,
      detrend_method = "polynomial",
      mesh           = mesh
    )
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  Completed in %.1f min\n", elapsed / 60))

  # ---- Save dataset results ------------------------------------------------
  if (is.data.frame(results_list) && nrow(results_list) > 0) {
    results_list$dataset_name  <- ds_name
    results_list$dataset_file  <- basename(ds_file)
    results_list$n_spots       <- n_spots
    results_list$n_genes_total <- n_genes_total

    saveRDS(results_list, ds_out)
    cat("  Saved:", ds_out, "\n")

    all_results[[ds_name]] <- results_list

    n_ok       <- sum(!is.na(results_list$log_bayes_factor))
    n_ns_favor <- sum(results_list$log_bayes_factor > 0, na.rm = TRUE)
    cat(sprintf("  Successful: %d/%d  |  Favour non-stationary: %d (%.1f%%)\n",
                n_ok, nrow(results_list), n_ns_favor,
                100 * n_ns_favor / max(n_ok, 1)))
  } else {
    cat("  WARNING: no results returned for this dataset\n")
  }

  rm(counts_mat, mesh); gc()
}

stopCluster(cl)

# ---- Combine and save all results ------------------------------------------
cat("\n", rep("=", 60), "\n", sep = "")
cat("Combining results across all datasets...\n")

combined <- do.call(rbind, all_results)

saveRDS(combined, file.path(output_dir, "all_results_combined.rds"))
write.csv(combined, file.path(output_dir, "all_results_combined.csv"),
          row.names = FALSE)
cat("Saved combined results to:", output_dir, "\n")

# ---- Summary ----------------------------------------------------------------
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n")
cat("Datasets processed:", length(all_results), "\n")
cat("Total genes:       ", nrow(combined), "\n")
cat("Successful fits:   ", sum(!is.na(combined$log_bayes_factor)), "\n")
cat("Failed fits:       ", sum(is.na(combined$log_bayes_factor)), "\n")

if (sum(!is.na(combined$log_bayes_factor)) > 0) {
  bf_sum <- combined %>%
    filter(!is.na(log_bayes_factor)) %>%
    summarise(
      mean_log_bf        = mean(log_bayes_factor),
      median_log_bf      = median(log_bayes_factor),
      prop_nonstationary = mean(log_bayes_factor > 0),
      prop_extreme_ns    = mean(log_bayes_factor >  log(100)),
      prop_extreme_stat  = mean(log_bayes_factor < -log(100))
    )
  cat("\nBayes Factor summary:\n")
  cat("  Mean log(BF):                ", round(bf_sum$mean_log_bf, 3), "\n")
  cat("  Median log(BF):              ", round(bf_sum$median_log_bf, 3), "\n")
  cat("  % favouring non-stationary:  ", round(100 * bf_sum$prop_nonstationary, 1), "%\n")
  cat("  % extreme for stationary:    ", round(100 * bf_sum$prop_extreme_stat,  1), "%\n")
  cat("  % extreme for non-stationary:", round(100 * bf_sum$prop_extreme_ns,    1), "%\n")
}

cat("\nEnd time:", format(Sys.time()), "\n")
cat("Done!\n")
