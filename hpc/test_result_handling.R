#!/usr/bin/env Rscript
# Lightweight regression test for per-gene failure handling and completeness.
# This test does not fit an INLA model.

args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  normalizePath("hpc")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(foreach)
})
source(file.path(script_dir, "analysis_functions.R"))

original_analyzer <- analyze_gene_nonstationarity
analyze_gene_nonstationarity <- function(coords, expr, gene_name = NULL, ...) {
  if (identical(gene_name, "fails")) {
    stop("synthetic Newton-Raphson non-convergence")
  }

  result <- empty_gene_analysis_result(gene_name, NA_character_)
  result$sigma_b_sq_stationary <- 1
  result$sigma_eps_sq_stationary <- 1
  result$range_stationary <- 2
  result$prop_spatial_stationary <- 0.5
  result$log_ml_stationary <- -10
  result$sigma_b_sq_nonstationary <- 1.1
  result$sigma_eps_sq_nonstationary <- 0.9
  result$range_nonstationary <- 2.1
  result$prop_spatial_nonstationary <- 0.55
  result$log_ml_nonstationary <- -9.5
  result$theta1_ns <- 0
  result$theta2_ns <- 0
  result$theta3_ns <- 0
  result$theta4_ns <- 0
  result$log_ml_nonspatial <- -11
  result$sigma_eps_sq_nonspatial <- 1
  result$mu_nonspatial <- 0
  result$log_bayes_factor <- 0.5
  result$bayes_factor <- exp(0.5)
  result$bf_interpretation <- "Anecdotal for nonstationary"
  if (identical(gene_name, "nonfinite")) {
    result$sigma_b_sq_stationary <- Inf
  }
  result
}

counts_mat <- matrix(seq_len(16), nrow = 4L)
coords <- data.frame(x = seq_len(4), y = seq_len(4))
gene_names <- c("before", "fails", "nonfinite", "after")
gene_ids <- paste0("ENSG_TEST_", seq_len(4))

rows <- lapply(seq_len(4), function(i) {
  analyze_gene_safe(
    gene_idx   = i,
    counts_mat = counts_mat,
    coords     = coords,
    gene_names = gene_names,
    gene_ids   = gene_ids,
    mesh       = structure(list(), class = "synthetic_mesh")
  )
})

stopifnot(all(vapply(rows[-1], function(x) identical(names(x), names(rows[[1]])), logical(1))))

results <- dplyr::bind_rows(rows)
stopifnot(
  nrow(results) == 4L,
  identical(results$gene_index, 1:4),
  identical(results$gene_id, gene_ids),
  is.na(results$error_message[1]),
  grepl("Newton-Raphson", results$error_message[2], fixed = TRUE),
  grepl("non-finite", results$error_message[3], fixed = TRUE),
  is.na(results$error_message[4]),
  !is.na(results$log_bayes_factor[4])
)

failed_metric_columns <- setdiff(
  names(results)[vapply(results, is.numeric, logical(1))],
  "gene_index"
)
stopifnot(all(is.na(results[2, failed_metric_columns, drop = TRUE])))

good_check <- validate_gene_results(results, 1:4, gene_ids)
truncated_check <- validate_gene_results(results[1:2, ], 1:4, gene_ids)
duplicated_results <- results
duplicated_results$gene_index[4] <- 2L
duplicate_check <- validate_gene_results(duplicated_results, 1:4, gene_ids)
missing_column_check <- validate_gene_results(
  results[, setdiff(names(results), "log_ml_nonspatial")],
  1:4,
  gene_ids
)
partial_failure <- results
partial_failure$log_ml_stationary[2] <- -10
partial_failure_check <- validate_gene_results(partial_failure, 1:4, gene_ids)
extreme_bf_results <- results
extreme_bf_results$log_ml_nonstationary[c(1, 4)] <- c(-1010, 990)
extreme_bf_results$log_bayes_factor[c(1, 4)] <- c(-1000, 1000)
extreme_bf_results$bayes_factor[c(1, 4)] <- c(0, Inf)
extreme_bf_check <- validate_gene_results(extreme_bf_results, 1:4, gene_ids)
stopifnot(
  good_check$valid,
  extreme_bf_check$valid,
  !truncated_check$valid,
  !duplicate_check$valid,
  !missing_column_check$valid,
  !partial_failure_check$valid
)

# An incomplete collector result must fail before it can replace an old file.
sentinel_path <- tempfile(fileext = ".rds")
saveRDS("sentinel", sentinel_path)
validation_error <- tryCatch({
  assert_complete_gene_results(results[1:2, ], 1:4, gene_ids)
  save_rds_atomic(results[1:2, ], sentinel_path)
  NULL
}, error = identity)
stopifnot(inherits(validation_error, "error"), identical(readRDS(sentinel_path), "sentinel"))
unlink(sentinel_path)

atomic_rds_path <- tempfile(fileext = ".rds")
atomic_csv_path <- tempfile(fileext = ".csv")
saveRDS("old", atomic_rds_path)
writeLines("old", atomic_csv_path)
save_rds_atomic(results, atomic_rds_path)
write_csv_atomic(results, atomic_csv_path, row.names = FALSE)
stopifnot(
  validate_gene_results(readRDS(atomic_rds_path), 1:4, gene_ids)$valid,
  nrow(read.csv(atomic_csv_path, stringsAsFactors = FALSE)) == 4L
)
unlink(c(atomic_rds_path, atomic_csv_path))

# Unexpected collector errors must propagate instead of returning a prefix.
registerDoSEQ()
collector_error <- tryCatch(
  foreach(i = 1:3, .errorhandling = "stop") %do% {
    if (i == 2L) stop("unexpected worker error")
    i
  },
  error = identity
)
stopifnot(inherits(collector_error, "error"))

# Verify the functions required by a PSOCK worker are explicitly exportable.
cl <- parallel::makeCluster(2L)
tryCatch({
  parallel::clusterExport(
    cl,
    c(
      "analyze_gene_nonstationarity", "empty_gene_analysis_result",
      "validate_gene_analysis_payload", "analysis_result_numeric_columns",
      "analysis_result_required_columns", "analyze_gene_safe",
      "counts_mat", "coords", "gene_names", "gene_ids"
    ),
    envir = environment()
  )
  parallel_rows <- parallel::parLapply(cl, seq_len(4), function(i) {
    analyze_gene_safe(
      gene_idx   = i,
      counts_mat = counts_mat,
      coords     = coords,
      gene_names = gene_names,
      gene_ids   = gene_ids,
      mesh       = structure(list(), class = "synthetic_mesh")
    )
  })
  parallel_results <- dplyr::bind_rows(parallel_rows)
  stopifnot(validate_gene_results(parallel_results, 1:4, gene_ids)$valid)
}, finally = {
  parallel::stopCluster(cl)
})

analyze_gene_nonstationarity <- original_analyzer
cat("PASS: per-gene failures remain explicit rows and completeness checks reject truncation.\n")
