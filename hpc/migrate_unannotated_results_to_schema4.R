#!/usr/bin/env Rscript

# Migrate validated schema-v3 results to schema v4 only when the source
# SpatialExperiment has no in_tissue annotation. In that case the new spot
# policy retains the exact same spots, so model estimates are not refitted or
# modified. Datasets with an in_tissue annotation must be refitted instead.

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1]] else "data/spatial_datasets"
results_dir <- if (length(args) >= 2L) args[[2]] else "results"
dry_run <- identical(Sys.getenv("MIGRATION_DRY_RUN", unset = "0"), "1")

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
source(file.path(script_dir, "spot_selection.R"))

legacy_implementation_md5 <- c(
  "analysis_functions.R" = "72c1620d5b29f26f8fcd93060f6a13e8",
  "run_analysis.R" = "8353bd7dab128333f7af131275067dfb"
)
prior_schema4_implementation_md5 <- c(
  "analysis_functions.R" = "f8d3a5b7a37adaffa49bfacb92757fe8",
  "spot_selection.R" = "01f0a1b19dfdb1cac4a7e0dd6ea27b3a",
  "run_analysis.R" = "8541bdb3c357719d0163fbfb29fe34c5"
)
current_implementation_md5 <- fingerprint_analysis_implementation(script_dir)

dataset_list_file <- Sys.getenv(
  "ANALYSIS_DATASET_LIST",
  unset = file.path(script_dir, "analysis_datasets.txt")
)
dataset_names <- trimws(readLines(dataset_list_file, warn = FALSE))
dataset_names <- dataset_names[
  nzchar(dataset_names) & !startsWith(dataset_names, "#")
]

for (dataset_name in dataset_names) {
  data_path <- file.path(data_dir, paste0(dataset_name, ".rds"))
  result_path <- file.path(results_dir, paste0(dataset_name, "_results.rds"))
  if (!file.exists(data_path) || !file.exists(result_path)) {
    stop("Missing data or result file for ", dataset_name, call. = FALSE)
  }

  spe <- readRDS(data_path)
  has_annotation <- "in_tissue" %in%
    colnames(SummarizedExperiment::colData(spe))
  if (has_annotation) {
    cat("Skipping ", dataset_name, ": requires numerical refit\n", sep = "")
    next
  }

  result <- readRDS(result_path)
  metadata <- attr(result, "analysis_metadata")
  schema_version <- if (is.list(metadata) &&
                        length(metadata$schema_version) == 1L) {
    suppressWarnings(as.integer(metadata$schema_version))
  } else {
    NA_integer_
  }
  is_current <- is.list(metadata) &&
    isTRUE(metadata$complete) &&
    identical(schema_version, analysis_schema_version) &&
    identical(as.character(metadata$spot_selection), "in_tissue_spots") &&
    analysis_implementation_matches(
      metadata$analysis_implementation_md5,
      current_implementation_md5
    )
  is_approved_input <- is.list(metadata) && isTRUE(metadata$complete) && (
    (
      identical(schema_version, 3L) &&
        identical(
          as.character(metadata$spot_selection),
          "all_bioconductor_supplied_spots"
        ) &&
        analysis_implementation_matches(
          metadata$analysis_implementation_md5,
          legacy_implementation_md5
        )
    ) ||
      (
        identical(schema_version, 4L) &&
          identical(as.character(metadata$spot_selection), "in_tissue_spots") &&
          analysis_implementation_matches(
            metadata$analysis_implementation_md5,
            prior_schema4_implementation_md5
          )
      )
  )
  if (!is_current && !is_approved_input) {
    stop(dataset_name, " is not an approved migration input", call. = FALSE)
  }

  source_md5 <- unname(as.character(tools::md5sum(data_path)))
  if (!identical(as.character(metadata$source_rds_md5), source_md5)) {
    stop(dataset_name, " source dataset fingerprint changed", call. = FALSE)
  }
  if (!identical(as.character(metadata$dataset_name), dataset_name) ||
      !identical(as.character(metadata$dataset_file), basename(data_path))) {
    stop(dataset_name, " dataset metadata is inconsistent", call. = FALSE)
  }

  selected <- select_in_tissue_spots(spe, dataset_name, verbose = FALSE)
  selection <- attr(selected, "spot_selection")
  if (isTRUE(selection$annotation_present) ||
      selection$n_spots_removed != 0L ||
      selection$n_spots_selected != ncol(spe)) {
    stop(dataset_name, " would change under the in-tissue policy", call. = FALSE)
  }
  if (!identical(
    suppressWarnings(as.integer(metadata$n_spots)),
    as.integer(selection$n_spots_selected)
  )) {
    stop(dataset_name, " result spot count does not match its source", call. = FALSE)
  }

  check <- validate_gene_results(result, seq_len(nrow(result)))
  if (!check$valid) {
    stop(dataset_name, " failed result validation: ", check$message, call. = FALSE)
  }

  if (is_current) {
    expected_supplied <- as.integer(selection$n_spots_supplied)
    if (!identical(
          suppressWarnings(as.integer(metadata$n_spots_supplied)),
          expected_supplied
        ) ||
        !identical(suppressWarnings(as.integer(metadata$n_spots_removed)), 0L) ||
        !identical(isTRUE(metadata$in_tissue_annotation), FALSE) ||
        !all(c(
          "n_spots_supplied", "n_spots_removed", "in_tissue_annotation"
        ) %in% names(result))) {
      stop(dataset_name, " has incomplete current metadata", call. = FALSE)
    }
    cat("Already current: ", dataset_name, "\n", sep = "")
    next
  }

  result$n_spots_supplied <- as.integer(selection$n_spots_supplied)
  result$n_spots_removed <- 0L
  result$in_tissue_annotation <- FALSE
  metadata$schema_version <- analysis_schema_version
  metadata$analysis_implementation_md5 <- current_implementation_md5
  metadata$n_spots_supplied <- as.integer(selection$n_spots_supplied)
  metadata$n_spots_removed <- 0L
  metadata$in_tissue_annotation <- FALSE
  metadata$spot_selection <- "in_tissue_spots"
  metadata$migrated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  metadata$migration <- paste(
    "Schema-v3 numerical results retained unchanged because the source object",
    "has no in_tissue annotation and the schema-v4 policy selects every supplied spot."
  )
  attr(result, "analysis_metadata") <- metadata

  if (dry_run) {
    cat("Validated migration: ", dataset_name, "\n", sep = "")
  } else {
    save_rds_atomic(result, result_path)
    cat("Migrated atomically: ", dataset_name, "\n", sep = "")
  }
}

cat(if (dry_run) "Dry-run migration validation complete.\n" else "Migration complete.\n")
