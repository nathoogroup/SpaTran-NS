#!/usr/bin/env Rscript
# =============================================================================
# Download SpatialExperiment Datasets from ExperimentHub
# =============================================================================
# This script downloads available Visium resources from ExperimentHub, saves
# SpatialExperiment objects locally, and verifies the study dataset manifest.
#
# Usage: Rscript download_datasets.R [output_dir]
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "data/spatial_datasets"

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  normalizePath("hpc")
}
dataset_list_file <- Sys.getenv(
  "ANALYSIS_DATASET_LIST",
  unset = file.path(script_dir, "analysis_datasets.txt")
)
if (!file.exists(dataset_list_file)) {
  stop("Dataset list not found: ", dataset_list_file, call. = FALSE)
}
required_datasets <- trimws(readLines(dataset_list_file, warn = FALSE))
required_datasets <- required_datasets[
  nzchar(required_datasets) & !startsWith(required_datasets, "#")
]
if (!length(required_datasets) || anyDuplicated(required_datasets)) {
  stop("Dataset list must contain unique study dataset basenames", call. = FALSE)
}

read_spatial_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

save_spatial_rds_atomic <- function(object, path) {
  tmp_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
  saveRDS(object, tmp_path)
  if (!file.rename(tmp_path, path)) {
    stop("Could not atomically replace dataset file: ", path, call. = FALSE)
  }
  invisible(path)
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SpatialExperiment)
})

cat("=== 10x Visium SpatialExperiment Dataset Downloader ===\n")
cat("Output directory:", output_dir, "\n\n")

# Initialize ExperimentHub
eh <- ExperimentHub()

# Query for 10x Visium SpatialExperiment datasets
cat("Querying ExperimentHub for 10x Visium SpatialExperiment datasets...\n")
spe_query <- query(eh, c("SpatialExperiment", "Visium"))
cat("Found", length(spe_query), "Visium SpatialExperiment datasets\n\n")

# Show manifest before downloading
cat("Dataset manifest:\n")
manifest_meta <- mcols(spe_query)
for (i in seq_len(nrow(manifest_meta))) {
  cat(sprintf("  [%s] %s\n", names(spe_query)[i], manifest_meta$title[i]))
}
cat("\n")

# Create manifest
manifest <- data.frame(
  eh_id      = names(spe_query),
  title      = mcols(spe_query)$title,
  species    = mcols(spe_query)$species,
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
  
  # Reuse only a readable SpatialExperiment. A corrupt or wrong-class file is
  # downloaded again and atomically replaced.
  if (file.exists(output_file)) {
    existing <- read_spatial_rds(output_file)
    if (inherits(existing, "SpatialExperiment")) {
      cat("  Existing SpatialExperiment is readable, skipping...\n")
      successful_downloads <- successful_downloads + 1
      next
    }
    cat("  Existing file is invalid; downloading a replacement...\n")
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

    # Drop image data — not needed for analysis and causes saveRDS to fail
    # when newer BiocGenerics::containsOutOfMemoryData encounters LoadedSpatialImage
    if (nrow(SpatialExperiment::imgData(obj)) > 0) {
      SpatialExperiment::imgData(obj) <- S4Vectors::DataFrame()
    }

    # Save to RDS
    cat("  Saving to:", output_file, "\n")
    save_spatial_rds_atomic(obj, output_file)
    
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
  failed = failed_downloads,
  required_datasets = required_datasets
)
saveRDS(log_data, file.path(output_dir, "download_log.rds"))

required_paths <- file.path(output_dir, paste0(required_datasets, ".rds"))
required_valid <- vapply(required_paths, function(path) {
  object <- read_spatial_rds(path)
  inherits(object, "SpatialExperiment")
}, logical(1))
if (!all(required_valid)) {
  stop(
    "Required study datasets are missing or invalid: ",
    paste(basename(required_paths[!required_valid]), collapse = ", "),
    call. = FALSE
  )
}

if (length(failed_downloads) > 0L) {
  cat("\nOptional query hits failed or were not SpatialExperiment objects; ",
      "all required study datasets are present.\n", sep = "")
}
cat("\nDone! All", length(required_datasets),
    "required study datasets are available in:", output_dir, "\n")
