#!/usr/bin/env Rscript

# Strictly validate and merge every gene chunk for one dataset into the
# canonical Dataset_results.rds consumed by combine_results.R.

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1]] else "data/spatial_datasets"
results_dir <- if (length(args) >= 2L) args[[2]] else "results"
chunk_dir <- if (length(args) >= 3L) args[[3]] else file.path(results_dir, "gene_chunks")
dataset_name <- if (length(args) >= 4L) args[[4]] else ""
chunk_count_arg <- if (length(args) >= 5L) args[[5]] else ""

if (!nzchar(dataset_name) || !nzchar(chunk_count_arg)) {
  stop(
    paste(
      "Usage: Rscript hpc/merge_gene_chunks.R",
      "<data_dir> <results_dir> <chunk_dir> <dataset_name> <chunk_count>"
    ),
    call. = FALSE
  )
}

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]])))
} else {
  normalizePath("hpc")
}

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
})
source(file.path(script_dir, "analysis_functions.R"))
source(file.path(script_dir, "spot_selection.R"))
source(file.path(script_dir, "gene_chunks.R"))

chunk_count <- positive_integer(chunk_count_arg, "chunk_count")
implementation_md5 <- fingerprint_analysis_implementation(script_dir)
data_path <- file.path(data_dir, paste0(dataset_name, ".rds"))
if (!file.exists(data_path)) {
  stop("Dataset file not found: ", data_path, call. = FALSE)
}

spe <- readRDS(data_path)
if (!inherits(spe, "SpatialExperiment")) {
  stop(data_path, " is not a SpatialExperiment", call. = FALSE)
}
source_rds_md5 <- unname(as.character(tools::md5sum(data_path)))
n_genes_total <- nrow(spe)
spe <- select_in_tissue_spots(spe, dataset_name = dataset_name)
selection <- attr(spe, "spot_selection")
n_spots <- ncol(spe)

counts_mat <- counts(spe)
min_spots <- ceiling(n_spots * 0.005)
keep <- (Matrix::rowSums(counts_mat > 0) >= min_spots) &
  (Matrix::rowSums(counts_mat) >= 3)
spe <- spe[keep, ]
n_genes_filtered <- nrow(spe)
if (n_genes_filtered < 1L) {
  stop("No genes passed the expression filter", call. = FALSE)
}

gene_ids <- rownames(spe)
all_gene_indices <- seq_len(n_genes_filtered)
if (is.null(gene_ids)) gene_ids <- as.character(all_gene_indices)
missing_ids <- is.na(gene_ids) | !nzchar(gene_ids)
gene_ids[missing_ids] <- as.character(all_gene_indices[missing_ids])

metadata_integer <- function(metadata, field) {
  value <- if (is.list(metadata)) metadata[[field]] else NULL
  if (length(value) == 1L) suppressWarnings(as.integer(value)) else NA_integer_
}

require_constant_column <- function(result, field, expected, path) {
  if (!field %in% names(result)) {
    stop(basename(path), " is missing ", field, call. = FALSE)
  }
  observed <- unique(result[[field]])
  if (length(observed) != 1L || is.na(observed) ||
      !identical(as.character(observed), as.character(expected))) {
    stop(basename(path), " has inconsistent ", field, call. = FALSE)
  }
}

chunks <- lapply(seq_len(chunk_count), function(chunk_id) {
  path <- gene_chunk_result_path(
    chunk_dir, dataset_name, chunk_id, chunk_count
  )
  if (!file.exists(path)) {
    stop("Missing gene chunk result: ", path, call. = FALSE)
  }

  result <- readRDS(path)
  metadata <- attr(result, "analysis_metadata")
  expected_indices <- gene_chunk_indices(
    n_genes_filtered, chunk_id, chunk_count
  )
  expected_ids <- gene_ids[expected_indices]

  if (!is.list(metadata) ||
      !isTRUE(metadata$complete) ||
      !identical(
        metadata_integer(metadata, "schema_version"),
        analysis_schema_version
      ) ||
      !analysis_implementation_matches(
        metadata$analysis_implementation_md5,
        implementation_md5
      ) ||
      !identical(as.character(metadata$dataset_name), dataset_name) ||
      !identical(as.character(metadata$dataset_file), basename(data_path)) ||
      !identical(as.character(metadata$source_rds_md5), source_rds_md5) ||
      !identical(as.character(metadata$spot_selection), "in_tissue_spots") ||
      !identical(metadata_integer(metadata, "n_spots"), n_spots) ||
      !identical(
        metadata_integer(metadata, "n_spots_supplied"),
        as.integer(selection$n_spots_supplied)
      ) ||
      !identical(
        metadata_integer(metadata, "n_spots_removed"),
        as.integer(selection$n_spots_removed)
      ) ||
      !identical(
        isTRUE(metadata$in_tissue_annotation),
        isTRUE(selection$annotation_present)
      ) ||
      !identical(metadata_integer(metadata, "n_genes_total"), n_genes_total) ||
      !identical(
        metadata_integer(metadata, "n_genes_filtered"),
        n_genes_filtered
      ) ||
      !identical(
        metadata_integer(metadata, "n_genes_analyzed"),
        length(expected_indices)
      ) ||
      !identical(metadata_integer(metadata, "gene_chunk_id"), chunk_id) ||
      !identical(
        metadata_integer(metadata, "gene_chunk_count"),
        chunk_count
      )) {
    stop(
      basename(path),
      " has missing or inconsistent gene-chunk metadata",
      call. = FALSE
    )
  }

  check <- validate_gene_results(result, expected_indices, expected_ids)
  if (!check$valid) {
    stop(basename(path), " failed validation: ", check$message, call. = FALSE)
  }

  require_constant_column(result, "dataset_name", dataset_name, path)
  require_constant_column(result, "dataset_file", basename(data_path), path)
  require_constant_column(result, "n_spots", n_spots, path)
  require_constant_column(
    result, "n_spots_supplied", selection$n_spots_supplied, path
  )
  require_constant_column(
    result, "n_spots_removed", selection$n_spots_removed, path
  )
  require_constant_column(
    result, "in_tissue_annotation", selection$annotation_present, path
  )
  require_constant_column(result, "n_genes_total", n_genes_total, path)

  cat(sprintf(
    "  [OK] chunk %d/%d: %d genes (%d failed fits retained)\n",
    chunk_id, chunk_count, nrow(result),
    sum(!is.na(result$error_message))
  ))
  result
})

merged <- dplyr::bind_rows(chunks)
merged <- merged[order(merged$gene_index), , drop = FALSE]
assert_complete_gene_results(merged, all_gene_indices, gene_ids)
if (sum(!is.na(merged$log_bayes_factor)) == 0L) {
  stop("All merged gene fits failed", call. = FALSE)
}

attr(merged, "analysis_metadata") <- list(
  schema_version = analysis_schema_version,
  complete = TRUE,
  analysis_implementation_md5 = implementation_md5,
  dataset_name = dataset_name,
  dataset_file = basename(data_path),
  source_rds_md5 = source_rds_md5,
  n_genes_analyzed = n_genes_filtered,
  n_genes_filtered = n_genes_filtered,
  n_genes_total = n_genes_total,
  gene_chunk_count = chunk_count,
  n_spots = n_spots,
  n_spots_supplied = as.integer(selection$n_spots_supplied),
  n_spots_removed = as.integer(selection$n_spots_removed),
  in_tissue_annotation = isTRUE(selection$annotation_present),
  spot_selection = "in_tissue_spots",
  merged_from_gene_chunks = TRUE,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

output_path <- file.path(results_dir, paste0(dataset_name, "_results.rds"))
save_rds_atomic(merged, output_path)
cat(sprintf(
  "Merged %d validated genes from %d chunks into %s\n",
  nrow(merged), chunk_count, output_path
))
