#!/usr/bin/env Rscript
# =============================================================================
# Spatial Nonstationarity Analysis - HPC Script with Parallel Processing
# =============================================================================
# Analyzes all genes across all SpatialExperiment datasets using doParallel
#
# Usage: Rscript run_analysis.R [data_dir] [output_dir] [n_cores] [genes_per_dataset]
#
# Arguments:
#   data_dir          - Directory containing .rds SpatialExperiment files
#   output_dir        - Directory for output results
#   n_cores           - Number of cores for parallel processing (default: all available)
#   genes_per_dataset - Number of top expressed genes to analyze per dataset 
#                       (default: "all" for all genes, or a number like 1000)
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data/spatial_datasets"
output_dir <- if (length(args) >= 2) args[2] else "results"
n_cores <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores()
genes_per_dataset <- if (length(args) >= 4) args[4] else "all"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get script directory for sourcing functions
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") {
  script_dir <- "."
}

cat("=============================================================\n")
cat("Spatial Nonstationarity Analysis - HPC Parallel Processing\n")
cat("=============================================================\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Number of cores:", n_cores, "\n")
cat("Genes per dataset:", genes_per_dataset, "\n")
cat("Start time:", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# Load required packages
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

# Test INLA
cat("Testing INLA installation...\n")
tryCatch({
  test_mesh <- INLA::inla.mesh.1d(seq(0, 1, length.out = 5))
  cat("✓ INLA is working\n\n")
}, error = function(e) {
  stop("INLA is not working properly: ", e$message)
})

# Find all dataset files
dataset_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(dataset_files), "dataset files\n\n")

if (length(dataset_files) == 0) {
  stop("No .rds files found in ", data_dir)
}

# Setup parallel backend
cat("Setting up parallel backend with", n_cores, "cores...\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export required packages and functions to cluster
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SummarizedExperiment)
    library(Matrix)
    library(INLA)
    library(dplyr)
  })
})

# Export functions to cluster
clusterExport(cl, c(
  "normalize_expression",
  "detrend_expression", 
  "create_spde_mesh",
  "fit_stationary_matern",
  "fit_nonstationary_matern_simple",
  "extract_variance_components",
  "calculate_bayes_factor",
  "analyze_gene_nonstationarity",
  "analyze_gene_safe"
))

# Process each dataset
all_results <- list()

for (ds_idx in seq_along(dataset_files)) {
  ds_file <- dataset_files[ds_idx]
  ds_name <- tools::file_path_sans_ext(basename(ds_file))
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Dataset", ds_idx, "/", length(dataset_files), ":", ds_name, "\n")
  cat(rep("=", 60), "\n")
  
  # Load dataset
  cat("Loading dataset...\n")
  tryCatch({
    spe <- readRDS(ds_file)
    
    if (!inherits(spe, "SpatialExperiment")) {
      cat("  WARNING: Not a SpatialExperiment, skipping\n")
      next
    }
    
    n_genes_total <- nrow(spe)
    n_spots <- ncol(spe)
    cat(sprintf("  Genes: %d, Spots: %d\n", n_genes_total, n_spots))
    
  }, error = function(e) {
    cat("  ERROR loading dataset:", e$message, "\n")
    next
  })
  
  # Get coordinates
  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords)[1:2] <- c("x", "y")
  
  # Get gene names
  rd <- rowData(spe)
  gene_names <- if ("symbol" %in% colnames(rd)) {
    as.character(rd$symbol)
  } else if ("gene_name" %in% colnames(rd)) {
    as.character(rd$gene_name)
  } else {
    rownames(spe)
  }
  
  # Select genes to analyze
  if (genes_per_dataset == "all") {
    gene_indices <- 1:n_genes_total
  } else {
    n_genes_analyze <- min(as.integer(genes_per_dataset), n_genes_total)
    # Select top expressed genes
    gene_totals <- Matrix::rowSums(assay(spe, "counts"))
    gene_indices <- order(gene_totals, decreasing = TRUE)[1:n_genes_analyze]
  }
  
  n_genes <- length(gene_indices)
  cat(sprintf("  Analyzing %d genes\n", n_genes))
  
  # Export dataset-specific variables to cluster
  clusterExport(cl, c("spe", "coords", "gene_names"), envir = environment())
  
  # Run parallel analysis
  cat("  Running parallel analysis...\n")
  start_time <- Sys.time()
  
  results_list <- foreach(
    i = seq_along(gene_indices),
    .combine = rbind,
    .errorhandling = "pass",
    .packages = c("SpatialExperiment", "SummarizedExperiment", "Matrix", "INLA", "dplyr")
  ) %dopar% {
    gene_idx <- gene_indices[i]
    analyze_gene_safe(
      gene_idx = gene_idx,
      spe = spe,
      coords = coords,
      gene_names = gene_names,
      normalize = TRUE,
      detrend = TRUE,
      detrend_method = "polynomial"
    )
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("  Completed in %.2f minutes\n", as.numeric(elapsed)))
  
  # Add dataset info to results
  if (is.data.frame(results_list) && nrow(results_list) > 0) {
    results_list$dataset_name <- ds_name
    results_list$dataset_file <- basename(ds_file)
    results_list$n_spots <- n_spots
    results_list$n_genes_total <- n_genes_total
    
    # Save individual dataset results
    ds_output_file <- file.path(output_dir, paste0(ds_name, "_results.rds"))
    saveRDS(results_list, ds_output_file)
    cat("  Saved results to:", ds_output_file, "\n")
    
    # Store in combined results
    all_results[[ds_name]] <- results_list
    
    # Quick summary
    n_successful <- sum(!is.na(results_list$log_bayes_factor))
    n_favor_nonstat <- sum(results_list$log_bayes_factor > 0, na.rm = TRUE)
    cat(sprintf("  Successful: %d/%d, Favor nonstationary: %d (%.1f%%)\n",
                n_successful, nrow(results_list), n_favor_nonstat,
                100 * n_favor_nonstat / max(n_successful, 1)))
  } else {
    cat("  WARNING: No results returned\n")
  }
  
  # Clean up memory
  rm(spe)
  gc()
}

# Stop cluster
stopCluster(cl)

# Combine all results
cat("\n", rep("=", 60), "\n", sep = "")
cat("Combining results...\n")

combined_results <- do.call(rbind, all_results)

# Save combined results
combined_file <- file.path(output_dir, "all_results_combined.rds")
saveRDS(combined_results, combined_file)
cat("Saved combined results to:", combined_file, "\n")

# Also save as CSV for easy viewing
csv_file <- file.path(output_dir, "all_results_combined.csv")
write.csv(combined_results, csv_file, row.names = FALSE)
cat("Saved CSV to:", csv_file, "\n")

# Summary statistics
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n")
cat("Total datasets processed:", length(all_results), "\n")
cat("Total genes analyzed:", nrow(combined_results), "\n")
cat("Successful analyses:", sum(!is.na(combined_results$log_bayes_factor)), "\n")
cat("Failed analyses:", sum(is.na(combined_results$log_bayes_factor)), "\n")

if (sum(!is.na(combined_results$log_bayes_factor)) > 0) {
  bf_summary <- combined_results %>%
    filter(!is.na(log_bayes_factor)) %>%
    summarise(
      mean_log_bf = mean(log_bayes_factor),
      median_log_bf = median(log_bayes_factor),
      prop_favor_nonstat = mean(log_bayes_factor > 0),
      prop_extreme_stat = mean(log_bayes_factor < -log(100)),
      prop_extreme_nonstat = mean(log_bayes_factor > log(100))
    )
  
  cat("\nBayes Factor Summary:\n")
  cat("  Mean log(BF):", round(bf_summary$mean_log_bf, 3), "\n")
  cat("  Median log(BF):", round(bf_summary$median_log_bf, 3), "\n")
  cat("  % favoring nonstationary:", round(100 * bf_summary$prop_favor_nonstat, 1), "%\n")
  cat("  % extreme for stationary:", round(100 * bf_summary$prop_extreme_stat, 1), "%\n")
  cat("  % extreme for nonstationary:", round(100 * bf_summary$prop_extreme_nonstat, 1), "%\n")
}

cat("\nEnd time:", format(Sys.time()), "\n")
cat("Done!\n")
