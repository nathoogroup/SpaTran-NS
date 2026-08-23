#!/usr/bin/env Rscript
# =============================================================================
# Patch: fit non-spatial (intercept + nugget) model for all genes
# =============================================================================
# Loads existing per-dataset results, re-preprocesses each gene's expression
# identically to the main analysis (log1p -> quantile normalise -> polynomial
# detrend), fits the non-spatial INLA model, and appends three columns:
#
#   log_ml_nonspatial     : INLA marginal log-likelihood
#   sigma_eps_sq_nonspatial : Gaussian noise variance (1 / precision)
#   mu_nonspatial         : posterior mean intercept
#
# Then updates all_results_combined.rds/csv.
#
# Usage (SLURM array, one task per dataset):
#   sbatch hpc/submit_nonspatial.sbatch
#
# Usage (local, all datasets sequentially):
#   Rscript hpc/patch_nonspatial.R [data_dir] [results_dir] [n_cores]
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
data_dir    <- if (length(args) >= 1) args[1] else "data/spatial_datasets"
results_dir <- if (length(args) >= 2) args[2] else "results"
n_cores     <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores()
task_id     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
overwrite   <- isTRUE(Sys.getenv("OVERWRITE", "0") == "1")

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "hpc")

cat("=============================================================\n")
cat("Patch: Non-spatial model (intercept + nugget)\n")
cat("=============================================================\n")
cat("Data directory:   ", data_dir, "\n")
cat("Results directory:", results_dir, "\n")
cat("Cores:            ", n_cores, "\n")
cat("SLURM array task: ", if (task_id > 0) task_id else "N/A (sequential)", "\n")
cat("Start time:       ", format(Sys.time()), "\n")
cat("=============================================================\n\n")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(INLA)
  library(dplyr)
  library(doParallel)
  library(foreach)
})

source(file.path(script_dir, "analysis_functions.R"))

# Build dataset list from existing results files — ensures we only process
# datasets that were actually analysed and skips any data files without results.
results_files <- sort(list.files(results_dir, pattern = "_results\\.rds$", full.names = TRUE))
results_files <- results_files[!grepl("combined|original", basename(results_files))]

dataset_files <- vapply(results_files, function(rf) {
  ds_name <- sub("_results\\.rds$", "", basename(rf))
  file.path(data_dir, paste0(ds_name, ".rds"))
}, character(1))

exist <- file.exists(dataset_files)
if (any(!exist)) {
  cat("WARNING: missing data files for:", paste(basename(dataset_files[!exist]), collapse=", "), "\n")
}
dataset_files <- dataset_files[exist]
cat("Found", length(dataset_files), "datasets with both results and data files\n")

if (task_id > 0) {
  if (task_id > length(dataset_files))
    stop("SLURM_ARRAY_TASK_ID (", task_id, ") exceeds number of datasets (",
         length(dataset_files), ")")
  cat("Array mode: processing dataset", task_id, "of", length(dataset_files), "\n")
  dataset_files <- dataset_files[task_id]
}
cat("\n")

# --- Safe per-gene wrapper ---------------------------------------------------

fit_nonspatial_safe <- function(gene_idx, counts_mat, coords, gene_names,
                                 timeout_secs = 60L) {
  gene_name <- gene_names[gene_idx]
  tryCatch({
    setTimeLimit(elapsed = timeout_secs, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)

    expr <- log1p(counts_mat[gene_idx, ])
    expr <- normalize_expression(expr, method = "quantile")
    expr <- detrend_expression(coords, expr, method = "polynomial")$detrended

    fit  <- fit_nonspatial(expr, verbose = FALSE)
    comp <- extract_nonspatial_components(fit)

    data.frame(
      gene_name             = gene_name,
      log_ml_nonspatial     = comp$log_ml,
      sigma_eps_sq_nonspatial = comp$sigma_eps_sq,
      mu_nonspatial         = comp$mu,
      stringsAsFactors      = FALSE
    )
  }, error = function(e) {
    data.frame(
      gene_name             = gene_name,
      log_ml_nonspatial     = NA_real_,
      sigma_eps_sq_nonspatial = NA_real_,
      mu_nonspatial         = NA_real_,
      stringsAsFactors      = FALSE
    )
  })
}

# --- Parallel setup ----------------------------------------------------------

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, suppressPackageStartupMessages({
  library(Matrix)
  library(INLA)
  library(dplyr)
  INLA::inla.setOption(num.threads = "1:1")
}))

clusterExport(cl, c(
  "normalize_expression", "detrend_expression",
  "fit_nonspatial", "extract_nonspatial_components",
  "fit_nonspatial_safe"
))

# --- Per-dataset loop --------------------------------------------------------

for (ds_file in dataset_files) {
  ds_name <- tools::file_path_sans_ext(basename(ds_file))
  ds_out  <- file.path(results_dir, paste0(ds_name, "_results.rds"))

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Dataset:", ds_name, "\n")
  cat(rep("=", 60), "\n")

  if (!file.exists(ds_out)) {
    cat("  No results file found, skipping\n"); next
  }

  existing <- readRDS(ds_out)

  if ("log_ml_nonspatial" %in% names(existing) && !overwrite) {
    cat("  log_ml_nonspatial already present, skipping (OVERWRITE=1 to rerun)\n"); next
  }

  spe <- tryCatch(readRDS(ds_file), error = function(e) {
    cat("  ERROR loading SPE:", e$message, "\n"); NULL
  })
  if (is.null(spe) || !inherits(spe, "SpatialExperiment")) {
    cat("  Skipping (not a SpatialExperiment)\n"); next
  }

  counts_mat <- tryCatch(counts(spe), error = function(e) {
    cat("  ERROR reading counts:", e$message, "\n"); NULL
  })
  if (is.null(counts_mat)) next

  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords)[1:2] <- c("x", "y")

  rd         <- rowData(spe)
  gene_names <- if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
                else if ("symbol" %in% colnames(rd)) as.character(rd$symbol)
                else rownames(spe)

  # Apply the same quality filter as the main analysis
  n_spots   <- ncol(spe)
  min_spots <- ceiling(n_spots * 0.005)
  keep      <- (Matrix::rowSums(counts_mat > 0) >= min_spots) &
               (Matrix::rowSums(counts_mat) >= 3)
  counts_mat <- counts_mat[keep, ]
  gene_names <- gene_names[keep]

  cat(sprintf("  Genes to fit: %d\n", length(gene_names)))

  rm(spe); gc()

  clusterExport(cl, c("counts_mat", "coords", "gene_names"), envir = environment())

  cat("  Running non-spatial fits...\n")
  t0 <- proc.time()

  ns_raw <- foreach(
    i              = seq_along(gene_names),
    .combine       = c,
    .multicombine  = TRUE,
    .errorhandling = "pass",
    .packages      = c("Matrix", "INLA", "dplyr")
  ) %dopar% {
    list(fit_nonspatial_safe(i, counts_mat, coords, gene_names))
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  Completed in %.1f min\n", elapsed / 60))

  ns_results <- dplyr::bind_rows(Filter(is.data.frame, ns_raw))

  # Merge into existing results by gene_name
  existing <- existing %>%
    select(-any_of(c("log_ml_nonspatial", "sigma_eps_sq_nonspatial", "mu_nonspatial"))) %>%
    left_join(ns_results, by = "gene_name")

  saveRDS(existing, ds_out)
  cat("  Saved:", ds_out, "\n")

  n_ok <- sum(!is.na(existing$log_ml_nonspatial))
  cat(sprintf("  Non-spatial fits OK: %d / %d\n", n_ok, nrow(existing)))

  rm(counts_mat); gc()
}

stopCluster(cl)

cat("\nEnd time:", format(Sys.time()), "\n")
cat("Done!\n")
