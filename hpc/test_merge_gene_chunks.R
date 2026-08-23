#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]])))
} else {
  normalizePath("hpc")
}

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
})
source(file.path(script_dir, "analysis_functions.R"))
source(file.path(script_dir, "gene_chunks.R"))

test_root <- tempfile("gene_chunk_merge_test_")
data_dir <- file.path(test_root, "data")
results_dir <- file.path(test_root, "results")
chunk_dir <- file.path(results_dir, "gene_chunks")
dir.create(data_dir, recursive = TRUE)
dir.create(chunk_dir, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE), add = TRUE)

dataset_name <- "SyntheticDataset"
gene_ids <- paste0("gene", seq_len(5L))
spot_ids <- paste0("spot", seq_len(4L))
spe <- SpatialExperiment(
  assays = list(counts = matrix(
    2L,
    nrow = length(gene_ids),
    ncol = length(spot_ids),
    dimnames = list(gene_ids, spot_ids)
  )),
  colData = S4Vectors::DataFrame(
    in_tissue = c(1L, 1L, 0L, 0L),
    row.names = spot_ids
  ),
  spatialCoords = cbind(x = seq_along(spot_ids), y = seq_along(spot_ids))
)
data_path <- file.path(data_dir, paste0(dataset_name, ".rds"))
saveRDS(spe, data_path)

implementation_md5 <- fingerprint_analysis_implementation(script_dir)
chunk_count <- 2L
for (chunk_id in seq_len(chunk_count)) {
  indices <- gene_chunk_indices(length(gene_ids), chunk_id, chunk_count)
  numeric_df <- as.data.frame(matrix(
    1,
    nrow = length(indices),
    ncol = length(analysis_result_numeric_columns)
  ))
  names(numeric_df) <- analysis_result_numeric_columns
  numeric_df$prop_spatial_stationary <- 0.5
  numeric_df$prop_spatial_nonstationary <- 0.5
  numeric_df$log_ml_stationary <- 1
  numeric_df$log_ml_nonstationary <- 2
  numeric_df$log_bayes_factor <- 1
  numeric_df$bayes_factor <- exp(1)

  result <- data.frame(
    gene_index = indices,
    gene_id = gene_ids[indices],
    gene_name = gene_ids[indices],
    numeric_df,
    bf_interpretation = "Positive",
    error_message = NA_character_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  result$dataset_name <- dataset_name
  result$dataset_file <- basename(data_path)
  result$n_spots <- 2L
  result$n_spots_supplied <- 4L
  result$n_spots_removed <- 2L
  result$in_tissue_annotation <- TRUE
  result$n_genes_total <- 5L
  attr(result, "analysis_metadata") <- list(
    schema_version = analysis_schema_version,
    complete = TRUE,
    analysis_implementation_md5 = implementation_md5,
    dataset_name = dataset_name,
    dataset_file = basename(data_path),
    source_rds_md5 = unname(as.character(tools::md5sum(data_path))),
    n_genes_analyzed = length(indices),
    n_genes_filtered = 5L,
    n_genes_total = 5L,
    gene_chunk_id = chunk_id,
    gene_chunk_count = chunk_count,
    n_spots = 2L,
    n_spots_supplied = 4L,
    n_spots_removed = 2L,
    in_tissue_annotation = TRUE,
    spot_selection = "in_tissue_spots"
  )
  saveRDS(
    result,
    gene_chunk_result_path(
      chunk_dir, dataset_name, chunk_id, chunk_count
    )
  )
}

merge_output <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    shQuote(file.path(script_dir, "merge_gene_chunks.R")),
    shQuote(data_dir),
    shQuote(results_dir),
    shQuote(chunk_dir),
    dataset_name,
    as.character(chunk_count)
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(merge_output, "status")) && attr(merge_output, "status") != 0L) {
  stop(paste(merge_output, collapse = "\n"), call. = FALSE)
}

merged <- readRDS(file.path(results_dir, paste0(dataset_name, "_results.rds")))
metadata <- attr(merged, "analysis_metadata")
stopifnot(
  nrow(merged) == 5L,
  identical(as.integer(merged$gene_index), seq_len(5L)),
  identical(as.character(merged$gene_id), gene_ids),
  isTRUE(metadata$complete),
  identical(metadata$n_genes_analyzed, 5L),
  identical(metadata$gene_chunk_count, 2L),
  isTRUE(metadata$merged_from_gene_chunks)
)

cat("Gene-chunk merge test passed.\n")
