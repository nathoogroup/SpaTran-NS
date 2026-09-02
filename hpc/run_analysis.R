#!/usr/bin/env Rscript
# =============================================================================
# Spatial Nonstationarity Analysis - HPC Script with Parallel Processing
# =============================================================================
# For each 10x Visium SpatialExperiment dataset:
#   1. Filter lowly-expressed genes
#   2. Apply log1p, rank-based quantile normalization, and mean detrending
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
task_id    <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0")))
overwrite  <- identical(Sys.getenv("OVERWRITE", unset = "0"), "1")
spot_selection_policy <- "all_bioconductor_supplied_spots"

if (is.na(n_cores) || n_cores < 1L) {
  stop("n_cores must be a positive integer", call. = FALSE)
}
if (is.na(task_id) || task_id < 0L) {
  stop("SLURM_ARRAY_TASK_ID must be a non-negative integer", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  normalizePath("hpc")
}

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
analysis_implementation_md5 <- fingerprint_analysis_implementation(script_dir)

# Verify INLA
cat("Testing INLA...\n")
tryCatch({
  INLA::inla.mesh.1d(seq(0, 1, length.out = 5))
  cat("INLA OK\n\n")
}, error = function(e) {
  stop("INLA not working: ", e$message,
       "\nRun bash hpc/install_inla.sh first.")
})

# Use a version-controlled dataset list so SLURM indices cannot shift when
# ExperimentHub adds resources or unrelated .rds files appear in data_dir.
dataset_list_file <- Sys.getenv(
  "ANALYSIS_DATASET_LIST",
  unset = file.path(script_dir, "analysis_datasets.txt")
)
if (!file.exists(dataset_list_file)) {
  stop("Dataset list not found: ", dataset_list_file, call. = FALSE)
}

dataset_names <- trimws(readLines(dataset_list_file, warn = FALSE))
dataset_names <- dataset_names[nzchar(dataset_names) & !startsWith(dataset_names, "#")]
if (anyDuplicated(dataset_names)) {
  stop("Dataset list contains duplicate basenames: ", dataset_list_file, call. = FALSE)
}
dataset_files <- file.path(data_dir, paste0(dataset_names, ".rds"))
missing_dataset_files <- dataset_files[!file.exists(dataset_files)]
if (length(missing_dataset_files) > 0L) {
  stop(
    "The following listed datasets are missing: ",
    paste(basename(missing_dataset_files), collapse = ", "),
    call. = FALSE
  )
}
cat("Loaded", length(dataset_files), "study datasets from", dataset_list_file, "\n")

if (length(dataset_files) == 0L) {
  stop("No dataset .rds files found in: ", data_dir, call. = FALSE)
}

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

tryCatch({
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
  "fit_stationary_matern", "fit_nonstationary_matern", "fit_nonspatial_gaussian",
  "extract_variance_components", "extract_nonstationary_thetas",
  "extract_nonspatial_components",
  "calculate_bayes_factor", "empty_gene_analysis_result",
  "validate_gene_analysis_payload", "analysis_result_numeric_columns",
  "analysis_result_required_columns",
  "analyze_gene_nonstationarity", "analyze_gene_safe"
))

worker_exports_ok <- unlist(clusterCall(cl, function() {
  all(vapply(
    c("extract_nonstationary_thetas", "extract_nonspatial_components",
      "fit_nonspatial_gaussian", "empty_gene_analysis_result",
      "validate_gene_analysis_payload",
      "analyze_gene_nonstationarity", "analyze_gene_safe"),
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  ))
}))
if (!all(worker_exports_ok)) {
  stop("One or more analysis functions were not exported to every worker", call. = FALSE)
}

# ---- Per-dataset loop -------------------------------------------------------

all_results <- list()

for (ds_idx in seq_along(dataset_files)) {
  ds_file <- dataset_files[ds_idx]
  ds_name <- tools::file_path_sans_ext(basename(ds_file))
  source_rds_md5 <- unname(as.character(tools::md5sum(ds_file)))
  if (length(source_rds_md5) != 1L || is.na(source_rds_md5)) {
    stop("Could not fingerprint source dataset: ", ds_file, call. = FALSE)
  }

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Dataset", ds_idx, "/", length(dataset_files), ":", ds_name, "\n")
  cat(rep("=", 60), "\n")

  # Existing outputs are validated after filtering the source dataset.  A file
  # is never skipped merely because it exists.
  ds_out <- file.path(output_dir, paste0(ds_name, "_results.rds"))

  # ---- Load ----------------------------------------------------------------
  load_error <- NULL
  spe <- tryCatch(readRDS(ds_file), error = function(e) {
    load_error <<- conditionMessage(e)
    NULL
  })
  if (is.null(spe) || !inherits(spe, "SpatialExperiment")) {
    problem <- if (!is.null(load_error)) {
      paste0("could not load dataset: ", load_error)
    } else {
      "object is not a SpatialExperiment"
    }
    stop(problem, call. = FALSE)
  }

  n_genes_total <- nrow(spe)
  n_spots       <- ncol(spe)
  cat(sprintf("  Genes: %d, Spots: %d\n", n_genes_total, n_spots))
  cat("  Spot selection: all spots supplied in the Bioconductor object\n")

  # ---- Extract counts (once, with error handling for HDF5-backed objects) --
  # Some datasets (e.g. Janesick) are stored as HDF5-backed DelayedArrays;
  # their .rds is tiny and the actual data lives in a separate HDF5 file that
  # may be absent.  We surface that failure explicitly and stop rather than
  # creating a partial study result.
  counts_error <- NULL
  counts_mat <- tryCatch(counts(spe), error = function(e) {
    counts_error <<- conditionMessage(e)
    NULL
  })
  if (is.null(counts_mat)) {
    problem <- paste0(
      "could not read counts (possibly a missing HDF5 backing file): ",
      counts_error
    )
    stop(problem, call. = FALSE)
  }

  # ---- Filter lowly-expressed genes ----------------------------------------
  min_spots  <- ceiling(n_spots * 0.005)  # detected in >= 0.5% of spots
  keep       <- (Matrix::rowSums(counts_mat > 0) >= min_spots) &
                (Matrix::rowSums(counts_mat) >= 3)
  counts_mat <- counts_mat[keep, ]   # subset the matrix directly
  spe        <- spe[keep, ]
  cat(sprintf("  After quality filter: %d genes (analyse all)\n", nrow(spe)))

  if (nrow(spe) == 0L) {
    problem <- "no genes passed the expression filter"
    stop(problem, call. = FALSE)
  }

  # All filtered genes
  gene_indices <- seq_len(nrow(spe))

  # ---- Coordinates and gene names ------------------------------------------
  coords <- as.data.frame(spatialCoords(spe))
  if (ncol(coords) < 2L) {
    problem <- "spatialCoords does not contain at least two coordinate columns"
    stop(problem, call. = FALSE)
  }
  colnames(coords)[1:2] <- c("x", "y")

  rd         <- rowData(spe)
  gene_names <- if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
                else if ("symbol"    %in% colnames(rd)) as.character(rd$symbol)
                else rownames(spe)
  gene_ids   <- rownames(spe)
  if (is.null(gene_ids)) gene_ids <- as.character(gene_indices)
  missing_ids <- is.na(gene_ids) | !nzchar(gene_ids)
  gene_ids[missing_ids] <- as.character(gene_indices[missing_ids])
  if (is.null(gene_names)) gene_names <- gene_ids
  missing_names <- is.na(gene_names) | !nzchar(gene_names)
  gene_names[missing_names] <- gene_ids[missing_names]

  # A cached output is reusable only if it has the strict one-row-per-gene
  # schema and exactly matches this filtered dataset.  Legacy/partial outputs
  # are rerun automatically instead of being silently accepted.
  if (file.exists(ds_out) && !overwrite) {
    cached_error <- NULL
    cached <- tryCatch(readRDS(ds_out), error = function(e) {
      cached_error <<- conditionMessage(e)
      NULL
    })
    cached_check <- if (is.null(cached)) {
      list(valid = FALSE, message = paste0("could not read cached file: ", cached_error))
    } else {
      cached_metadata <- attr(cached, "analysis_metadata")
      cached_schema <- if (is.list(cached_metadata) &&
                           length(cached_metadata$schema_version) == 1L) {
        suppressWarnings(as.integer(cached_metadata$schema_version))
      } else {
        NA_integer_
      }
      if (!is.list(cached_metadata) ||
          !isTRUE(cached_metadata$complete) ||
          is.na(cached_schema) ||
          cached_schema != analysis_schema_version ||
          !analysis_implementation_matches(
            cached_metadata$analysis_implementation_md5,
            analysis_implementation_md5
          ) ||
          !identical(as.character(cached_metadata$dataset_name), ds_name) ||
          !identical(as.character(cached_metadata$dataset_file), basename(ds_file)) ||
          !identical(as.character(cached_metadata$source_rds_md5), source_rds_md5) ||
          !identical(
            as.character(cached_metadata$spot_selection),
            spot_selection_policy
          ) ||
          !identical(
            suppressWarnings(as.integer(cached_metadata$n_spots)),
            n_spots
          ) ||
          !identical(
            suppressWarnings(as.integer(cached_metadata$n_genes_total)),
            n_genes_total
          ) ||
          !identical(
            suppressWarnings(as.integer(cached_metadata$n_genes_analyzed)),
            length(gene_indices)
          )) {
        list(valid = FALSE, message = "missing or inconsistent schema-v3 completion metadata")
      } else {
        validate_gene_results(cached, gene_indices, gene_ids)
      }
    }

    if (cached_check$valid) {
      cat("  Existing results are complete; skipping (set OVERWRITE=1 to rerun)\n")
      all_results[[ds_name]] <- cached
      rm(spe, counts_mat, cached); gc()
      next
    }

    cat("  Existing results are not reusable:", cached_check$message, "\n")
    cat("  Rerunning dataset; the old file remains until a complete replacement is ready.\n")
    rm(cached); gc()
  } else if (file.exists(ds_out) && overwrite) {
    cat("  OVERWRITE=1: rerunning existing result file\n")
  }

  # Pre-build the SPDE mesh once; all workers share the same mesh object
  # rather than rebuilding it for every gene.
  cat("  Building shared SPDE mesh...\n")
  mesh <- create_spde_mesh(coords)
  cat(sprintf("  Mesh: %d nodes\n", mesh$n))

  # Free the SPE before the cluster export — this is the main memory saving.
  # coords, gene identifiers, counts_mat, and mesh are all we need from here on.
  rm(spe); gc()

  # ---- Parallel analysis ---------------------------------------------------
  clusterExport(cl, c("counts_mat", "coords", "gene_names", "gene_ids", "mesh"),
                envir = environment())

  cat("  Running parallel INLA fits for", length(gene_indices), "genes...\n")
  t0 <- proc.time()

  # Collect a list first.  Combining inside foreach with base rbind can return
  # a silently truncated prefix when one worker row has an unexpected schema.
  results_list <- foreach(
    i              = seq_along(gene_indices),
    .errorhandling = "stop",
    .inorder       = TRUE,
    .packages      = c("Matrix", "INLA", "dplyr")
  ) %dopar% {
    analyze_gene_safe(
      gene_idx       = gene_indices[i],
      counts_mat     = counts_mat,
      coords         = coords,
      gene_names     = gene_names,
      gene_ids       = gene_ids,
      normalize      = TRUE,
      detrend        = TRUE,
      detrend_method = "polynomial",
      mesh           = mesh
    )
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  Completed in %.1f min\n", elapsed / 60))

  # ---- Save dataset results ------------------------------------------------
  results_df <- dplyr::bind_rows(results_list)
  assert_complete_gene_results(results_df, gene_indices, gene_ids)

  n_ok       <- sum(!is.na(results_df$log_bayes_factor))
  n_failed   <- sum(!is.na(results_df$error_message))
  n_ns_favor <- sum(results_df$log_bayes_factor > 0, na.rm = TRUE)
  if (n_ok == 0L) {
    stop(
      "All gene fits failed; refusing to save a systemically invalid result",
      call. = FALSE
    )
  }

  results_df$dataset_name  <- ds_name
  results_df$dataset_file  <- basename(ds_file)
  results_df$n_spots       <- n_spots
  results_df$n_genes_total <- n_genes_total
  attr(results_df, "analysis_metadata") <- list(
    schema_version     = analysis_schema_version,
    complete           = TRUE,
    analysis_implementation_md5 = analysis_implementation_md5,
    dataset_name       = ds_name,
    dataset_file       = basename(ds_file),
    source_rds_md5     = source_rds_md5,
    n_genes_analyzed   = length(gene_indices),
    n_genes_total      = n_genes_total,
    n_spots            = n_spots,
    spot_selection     = spot_selection_policy,
    completed_at       = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  # The previous file is replaced only after the complete object has been
  # serialized successfully in the same directory.
  save_rds_atomic(results_df, ds_out)
  cat("  Saved complete result:", ds_out, "\n")

  all_results[[ds_name]] <- results_df

  cat(sprintf(
    "  Successful fits: %d/%d  |  Failed fits retained as NA rows: %d\n",
    n_ok, nrow(results_df), n_failed
  ))
  cat(sprintf(
    "  Favour non-stationary: %d (%.1f%% of successful fits)\n",
    n_ns_favor, 100 * n_ns_favor / max(n_ok, 1)
  ))

  rm(counts_mat, mesh); gc()
}

if (length(all_results) == 0L) {
  stop("No datasets produced valid results", call. = FALSE)
}
if (task_id == 0L && !identical(names(all_results), dataset_names)) {
  stop("Sequential analysis did not produce every listed dataset", call. = FALSE)
}

if (task_id > 0L) {
  # Array tasks share output_dir, so they must never race to overwrite a common
  # all_results_combined file.  Combine only after the entire array succeeds.
  cat("\nArray mode: per-dataset result complete.\n")
  cat("Shared combined outputs were not written by this array task.\n")
} else {
  # ---- Combine and save all results ----------------------------------------
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Combining results across all datasets...\n")

  combined <- dplyr::bind_rows(all_results)
  attr(combined, "analysis_metadata") <- list(
    schema_version = analysis_schema_version,
    complete       = TRUE,
    analysis_implementation_md5 = analysis_implementation_md5,
    dataset_names  = names(all_results),
    n_datasets     = length(all_results),
    n_rows         = nrow(combined),
    spot_selection = spot_selection_policy,
    source_rds_md5 = stats::setNames(
      vapply(
        all_results,
        function(result) as.character(
          attr(result, "analysis_metadata")$source_rds_md5
        ),
        character(1)
      ),
      names(all_results)
    ),
    completed_at   = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  save_rds_atomic(combined, file.path(output_dir, "all_results_combined.rds"))
  write_csv_atomic(
    combined,
    file.path(output_dir, "all_results_combined.csv"),
    row.names = FALSE
  )
  cat("Saved combined results to:", output_dir, "\n")

  # ---- Summary --------------------------------------------------------------
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("ANALYSIS SUMMARY\n")
  cat(rep("=", 60), "\n")
  cat("Datasets processed:", length(all_results), "\n")
  cat("Total genes:       ", nrow(combined), "\n")
  cat("Successful fits:   ", sum(!is.na(combined$log_bayes_factor)), "\n")
  cat("Failed fits:       ", sum(!is.na(combined$error_message)), "\n")

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
}

cat("\nEnd time:", format(Sys.time()), "\n")
cat("Done!\n")
}, finally = {
  try(stopCluster(cl), silent = TRUE)
})
