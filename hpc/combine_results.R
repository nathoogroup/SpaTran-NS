#!/usr/bin/env Rscript
# Validate and combine completed per-dataset result files after a SLURM array.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[1] else "results"
expected_datasets <- if (length(args) >= 2L) {
  suppressWarnings(as.integer(args[2]))
} else {
  suppressWarnings(as.integer(Sys.getenv("N_EXPECTED_DATASETS", unset = "12")))
}

if (is.na(expected_datasets) || expected_datasets < 1L) {
  stop("expected_datasets must be a positive integer", call. = FALSE)
}

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  normalizePath("hpc")
}

suppressPackageStartupMessages(library(dplyr))
source(file.path(script_dir, "analysis_functions.R"))
spot_selection_policy <- "all_bioconductor_supplied_spots"
analysis_implementation_md5 <- fingerprint_analysis_implementation(script_dir)

dataset_list_file <- Sys.getenv(
  "ANALYSIS_DATASET_LIST",
  unset = file.path(script_dir, "analysis_datasets.txt")
)
if (!file.exists(dataset_list_file)) {
  stop("Dataset list not found: ", dataset_list_file, call. = FALSE)
}

dataset_names_expected <- trimws(readLines(dataset_list_file, warn = FALSE))
dataset_names_expected <- dataset_names_expected[
  nzchar(dataset_names_expected) & !startsWith(dataset_names_expected, "#")
]
if (anyDuplicated(dataset_names_expected)) {
  stop("Dataset list contains duplicate basenames: ", dataset_list_file, call. = FALSE)
}
if (length(dataset_names_expected) != expected_datasets) {
  stop(
    sprintf(
      "%s lists %d datasets; expected %d",
      dataset_list_file, length(dataset_names_expected), expected_datasets
    ),
    call. = FALSE
  )
}

result_files <- file.path(
  output_dir,
  paste0(dataset_names_expected, "_results.rds")
)
missing_result_files <- result_files[!file.exists(result_files)]
if (length(missing_result_files) > 0L) {
  stop(
    "Missing per-dataset result files: ",
    paste(basename(missing_result_files), collapse = ", "),
    call. = FALSE
  )
}

if (length(result_files) != expected_datasets) {
  stop(
    sprintf(
      "Found %d per-dataset result files in %s; expected %d",
      length(result_files), output_dir, expected_datasets
    ),
    call. = FALSE
  )
}

cat(sprintf("Validating %d per-dataset result files...\n", length(result_files)))
all_results <- lapply(result_files, function(path) {
  result <- readRDS(path)
  metadata <- attr(result, "analysis_metadata")
  metadata_integer <- function(field) {
    value <- if (is.list(metadata)) metadata[[field]] else NULL
    if (length(value) == 1L) suppressWarnings(as.integer(value)) else NA_integer_
  }
  schema_version <- if (is.list(metadata) && length(metadata$schema_version) == 1L) {
    suppressWarnings(as.integer(metadata$schema_version))
  } else {
    NA_integer_
  }
  n_genes_analyzed <- metadata_integer("n_genes_analyzed")
  n_genes_total <- metadata_integer("n_genes_total")
  n_spots <- metadata_integer("n_spots")

  if (!is.list(metadata) ||
      !isTRUE(metadata$complete) ||
      is.na(schema_version) ||
      schema_version != analysis_schema_version ||
      !analysis_implementation_matches(
        metadata$analysis_implementation_md5,
        analysis_implementation_md5
      ) ||
      is.na(n_genes_analyzed) ||
      n_genes_analyzed < 1L ||
      is.na(n_genes_total) ||
      n_genes_total < n_genes_analyzed ||
      is.na(n_spots) ||
      n_spots < 1L ||
      length(metadata$source_rds_md5) != 1L ||
      !grepl("^[0-9a-fA-F]{32}$", as.character(metadata$source_rds_md5)) ||
      !identical(as.character(metadata$spot_selection), spot_selection_policy)) {
    stop(
      basename(path),
      " is a legacy or unverified result; rerun its analysis task first",
      call. = FALSE
    )
  }

  expected_indices <- seq_len(n_genes_analyzed)
  check <- validate_gene_results(result, expected_indices)
  if (!check$valid) {
    stop(basename(path), " failed validation: ", check$message, call. = FALSE)
  }

  required_dataset_columns <- c(
    "dataset_name", "dataset_file", "n_spots", "n_genes_total"
  )
  missing_dataset_columns <- setdiff(required_dataset_columns, names(result))
  if (length(missing_dataset_columns) > 0L) {
    stop(
      basename(path),
      " is missing dataset columns: ",
      paste(missing_dataset_columns, collapse = ", "),
      call. = FALSE
    )
  }

  dataset_names <- unique(as.character(result$dataset_name))
  dataset_files <- unique(as.character(result$dataset_file))
  result_spot_counts <- unique(suppressWarnings(as.integer(result$n_spots)))
  result_gene_counts <- unique(suppressWarnings(as.integer(result$n_genes_total)))
  if (length(dataset_names) != 1L ||
      !identical(dataset_names, as.character(metadata$dataset_name)) ||
      !identical(dataset_names, dataset_names_expected[match(path, result_files)]) ||
      length(dataset_files) != 1L ||
      !identical(dataset_files, as.character(metadata$dataset_file)) ||
      length(result_spot_counts) != 1L ||
      is.na(result_spot_counts) ||
      !identical(result_spot_counts, n_spots) ||
      length(result_gene_counts) != 1L ||
      is.na(result_gene_counts) ||
      !identical(result_gene_counts, n_genes_total)) {
    stop(basename(path), " has inconsistent dataset metadata", call. = FALSE)
  }

  cat(sprintf(
    "  [OK] %-45s %d rows (%d failed fits retained)\n",
    basename(path),
    nrow(result),
    sum(!is.na(result$error_message))
  ))
  result
})

combined <- dplyr::bind_rows(all_results)
dataset_gene_key <- paste(combined$dataset_name, combined$gene_index, sep = "::")
if (anyDuplicated(dataset_gene_key)) {
  stop("Combined results contain duplicate dataset_name/gene_index keys", call. = FALSE)
}

attr(combined, "analysis_metadata") <- list(
  schema_version = analysis_schema_version,
  complete       = TRUE,
  analysis_implementation_md5 = analysis_implementation_md5,
  dataset_names  = dataset_names_expected,
  n_datasets     = length(all_results),
  n_rows         = nrow(combined),
  spot_selection = spot_selection_policy,
  source_rds_md5 = stats::setNames(
    vapply(
      all_results,
      function(result) as.character(attr(result, "analysis_metadata")$source_rds_md5),
      character(1)
    ),
    dataset_names_expected
  ),
  completed_at   = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

rds_path <- file.path(output_dir, "all_results_combined.rds")
csv_path <- file.path(output_dir, "all_results_combined.csv")
save_rds_atomic(combined, rds_path)
write_csv_atomic(combined, csv_path, row.names = FALSE)

cat(sprintf("Saved %d validated rows to %s and %s\n", nrow(combined), rds_path, csv_path))
