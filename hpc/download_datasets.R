#!/usr/bin/env Rscript
# =============================================================================
# Download SpatialExperiment Datasets from ExperimentHub
# =============================================================================
# This script downloads ALL available SpatialExperiment datasets
# from ExperimentHub and saves them locally for HPC analysis.
#
# Usage: Rscript download_datasets.R [output_dir]
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "data/spatial_datasets"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SpatialExperiment)
})

cat("=== SpatialExperiment Dataset Downloader ===\n")
cat("Output directory:", output_dir, "\n\n")

# Initialize ExperimentHub
eh <- ExperimentHub()

# Query for ALL SpatialExperiment datasets (not just Visium)
cat("Querying ExperimentHub for ALL SpatialExperiment datasets...\n")
spe_query <- query(eh, "SpatialExperiment")
cat("Found", length(spe_query), "SpatialExperiment datasets\n\n")

# Create manifest
manifest <- data.frame(
  eh_id = names(spe_query),
  title = mcols(spe_query)$title,
  species = mcols(spe_query)$species,
  rdataclass = mcols(spe_query)$rdataclass,
  stringsAsFactors = FALSE
)

# Save manifest
manifest_file <- file.path(output_dir, "dataset_manifest.csv")
write.csv(manifest, manifest_file, row.names = FALSE)
cat("Saved dataset manifest to:", manifest_file, "\n\n")

# Download each dataset
cat("Downloading datasets...\n")
cat(rep("=", 50), "\n", sep = "")

successful_downloads <- 0
failed_downloads <- c()

for (i in seq_along(spe_query)) {
  eh_id <- names(spe_query)[i]
  title <- mcols(spe_query)$title[i]
  
  cat(sprintf("\n[%d/%d] %s\n", i, length(spe_query), title))
  cat("  EH ID:", eh_id, "\n")
  
  # Create safe filename
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", title)
  safe_name <- gsub("_+", "_", safe_name)
  safe_name <- substr(safe_name, 1, 50)
  output_file <- file.path(output_dir, paste0(safe_name, ".rds"))
  
  # Skip if already downloaded
  if (file.exists(output_file)) {
    cat("  Already downloaded, skipping...\n")
    successful_downloads <- successful_downloads + 1
    next
  }
  
  # Download and save
  tryCatch({
    cat("  Downloading...\n")
    obj <- spe_query[[i]]
    
    # Verify it's a SpatialExperiment
    if (!inherits(obj, "SpatialExperiment")) {
      cat("  WARNING: Not a SpatialExperiment object (class:", class(obj)[1], "), skipping\n")
      failed_downloads <- c(failed_downloads, eh_id)
      next
    }
    
    # Get basic info
    n_genes <- nrow(obj)
    n_spots <- ncol(obj)
    cat(sprintf("  Genes: %d, Spots: %d\n", n_genes, n_spots))
    
    # Save to RDS
    cat("  Saving to:", output_file, "\n")
    saveRDS(obj, output_file)
    
    successful_downloads <- successful_downloads + 1
    cat("  SUCCESS\n")
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    failed_downloads <<- c(failed_downloads, eh_id)
  })
}

# Summary
cat("\n\n=== Download Summary ===\n")
cat("Successful downloads:", successful_downloads, "\n")
cat("Failed downloads:", length(failed_downloads), "\n")

if (length(failed_downloads) > 0) {
  cat("\nFailed dataset IDs:\n")
  for (id in failed_downloads) {
    cat("  -", id, "\n")
  }
}

# Save download log
log_data <- list(
  timestamp = Sys.time(),
  total_queried = length(spe_query),
  successful = successful_downloads,
  failed = failed_downloads
)
saveRDS(log_data, file.path(output_dir, "download_log.rds"))

cat("\nDone! Datasets saved to:", output_dir, "\n")
