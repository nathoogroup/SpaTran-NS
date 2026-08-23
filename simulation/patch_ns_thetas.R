#!/usr/bin/env Rscript
# =============================================================================
# Patch NS Thetas: Fit nonstationary model for H1 genes and add theta1_ns
# through theta4_ns to an existing results RDS file.
#
# Only the nonstationary model is fitted (stationary thetas are already stored
# as sigma_b_sq / range in the results and are not re-fitted).
#
# Usage:
#   Rscript simulation/patch_ns_thetas.R <dataset_name> [n_cores] [data_dir] [results_dir]
#
# Example:
#   Rscript simulation/patch_ns_thetas.R HumanBreastCancerILC 32
# =============================================================================

args         <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript patch_ns_thetas.R <dataset_name> [n_cores] [data_dir] [results_dir]")

dataset_name <- args[1]
n_cores      <- if (length(args) >= 2) as.integer(args[2]) else parallel::detectCores()
data_dir     <- if (length(args) >= 3) args[3] else "data/spatial_datasets"
results_dir  <- if (length(args) >= 4) args[4] else "results"

cat("=============================================================\n")
cat("Patch NS Thetas\n")
cat("=============================================================\n")
cat("Dataset:   ", dataset_name, "\n")
cat("Cores:     ", n_cores, "\n")
cat("Start:     ", format(Sys.time()), "\n")
cat("=============================================================\n\n")

suppressPackageStartupMessages({
  library(INLA)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(doParallel)
  library(foreach)
})

INLA::inla.setOption(num.threads = "1:1")

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "simulation")
hpc_dir    <- file.path(dirname(script_dir), "hpc")
if (!file.exists(file.path(hpc_dir, "analysis_functions.R"))) hpc_dir <- "hpc"
source(file.path(hpc_dir, "analysis_functions.R"))

# --- Load existing results ---------------------------------------------------

results_file <- file.path(results_dir, paste0(dataset_name, "_results.rds"))
if (!file.exists(results_file)) stop("Results file not found: ", results_file)
results <- readRDS(results_file)
cat(sprintf("Loaded results: %d genes\n", nrow(results)))

# Identify H1 genes needing NS theta fits
h1_genes <- results %>%
  filter(!is.na(log_bayes_factor), log_bayes_factor > 0) %>%
  pull(gene_name)
cat(sprintf("H1 genes to fit: %d\n", length(h1_genes)))

if (length(h1_genes) == 0) stop("No H1 genes found.")

# --- Load spatial dataset ----------------------------------------------------

spe_file <- file.path(data_dir, paste0(dataset_name, ".rds"))
if (!file.exists(spe_file)) stop("Dataset not found: ", spe_file)
cat("Loading spatial dataset...\n")
spe <- readRDS(spe_file)

coords     <- as.data.frame(spatialCoords(spe))
colnames(coords)[1:2] <- c("x", "y")
rd             <- rowData(spe)
gene_names_spe <- if ("gene_name" %in% colnames(rd)) {
                    as.character(rd$gene_name)
                  } else if ("symbol" %in% colnames(rd)) {
                    as.character(rd$symbol)
                  } else {
                    rownames(spe)
                  }

# Filter to H1 genes present in SPE
h1_in_spe  <- intersect(h1_genes, gene_names_spe)
n_missing  <- length(h1_genes) - length(h1_in_spe)
if (n_missing > 0) cat(sprintf("  WARNING: %d H1 genes not found in SPE\n", n_missing))
cat(sprintf("  Fitting %d H1 genes\n", length(h1_in_spe)))

h1_idx_spe <- match(h1_in_spe, gene_names_spe)
counts_mat  <- as.matrix(assay(spe, "counts")[h1_idx_spe, , drop = FALSE])
rm(spe); gc()

# --- Build mesh (same procedure as main analysis) ----------------------------

cat("Building SPDE mesh...\n")
mesh <- create_spde_mesh(coords)
cat(sprintf("  Mesh nodes: %d\n", mesh$n))

# --- Parallel NS fitting for H1 genes ----------------------------------------

cat(sprintf("Setting up %d parallel workers...\n", n_cores))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(dplyr)
  INLA::inla.setOption(num.threads = "1:1")
}))

clusterExport(cl, c(
  "fit_nonstationary_matern", "extract_nonstationary_thetas",
  "normalize_expression", "detrend_expression",
  "counts_mat", "coords", "mesh", "h1_in_spe"
), envir = environment())

# --- Checkpoint: skip genes already fitted in a previous run -----------------

checkpoint_file <- file.path(results_dir,
                             paste0(dataset_name, "_patch_checkpoint.rds"))
if (file.exists(checkpoint_file)) {
  ckpt <- readRDS(checkpoint_file)
  cat(sprintf("Checkpoint found: %d genes already fitted — skipping them\n",
              nrow(ckpt)))
  already_done <- ckpt$gene_name
  todo_idx     <- which(!(h1_in_spe %in% already_done))
  counts_mat   <- counts_mat[todo_idx, , drop = FALSE]
  h1_in_spe    <- h1_in_spe[todo_idx]
} else {
  ckpt <- NULL
}

cat(sprintf("Fitting NS model for %d H1 genes...\n", length(h1_in_spe)))
t0 <- proc.time()

batch_size <- n_cores * 4L
n_genes    <- length(h1_in_spe)
batches    <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / batch_size))
batch_results <- vector("list", length(batches))

for (b in seq_along(batches)) {
  idx_b   <- batches[[b]]
  cat(sprintf("  Batch %d/%d (genes %d–%d)...\n",
              b, length(batches), idx_b[1], idx_b[length(idx_b)]))

  raw_b <- foreach(
    i              = idx_b,
    .errorhandling = "pass",
    .packages      = c("INLA", "Matrix", "dplyr")
  ) %dopar% {
    gene_name <- h1_in_spe[i]
    tryCatch({
      setTimeLimit(elapsed = 120L, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
      expr      <- log1p(counts_mat[i, ])
      expr_proc <- normalize_expression(expr)
      expr_proc <- detrend_expression(coords, expr_proc, "polynomial")$detrended
      fit       <- fit_nonstationary_matern(coords, expr_proc, mesh = mesh, verbose = FALSE)
      thetas    <- extract_nonstationary_thetas(fit$result)
      data.frame(gene_name = gene_name,
                 theta1_ns = thetas$theta1,
                 theta2_ns = thetas$theta2,
                 theta3_ns = thetas$theta3,
                 theta4_ns = thetas$theta4,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(gene_name = gene_name,
                 theta1_ns = NA_real_,
                 theta2_ns = NA_real_,
                 theta3_ns = NA_real_,
                 theta4_ns = NA_real_,
                 stringsAsFactors = FALSE)
    })
  }

  batch_results[[b]] <- do.call(rbind, Filter(is.data.frame, raw_b))

  # Incremental checkpoint after every batch
  ckpt_so_far <- do.call(rbind, Filter(Negate(is.null), batch_results))
  if (!is.null(ckpt)) ckpt_so_far <- rbind(ckpt, ckpt_so_far)
  saveRDS(ckpt_so_far, checkpoint_file)
}

stopCluster(cl)

thetas_df <- do.call(rbind, Filter(Negate(is.null), batch_results))
elapsed   <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("  Done in %.1f min\n", elapsed / 60))
cat(sprintf("  Successful fits: %d / %d\n",
            sum(!is.na(thetas_df$theta3_ns)), nrow(thetas_df)))

# Merge with checkpoint results
if (!is.null(ckpt)) {
  thetas_df <- rbind(ckpt, thetas_df)
}
saveRDS(thetas_df, checkpoint_file)
cat(sprintf("Checkpoint saved: %d total fitted genes\n", nrow(thetas_df)))

# --- Patch results and save --------------------------------------------------

# Add theta columns (NA for H0 genes, fitted values for H1 genes)
results$theta1_ns <- NA_real_
results$theta2_ns <- NA_real_
results$theta3_ns <- NA_real_
results$theta4_ns <- NA_real_

idx <- match(thetas_df$gene_name, results$gene_name)
results$theta1_ns[idx] <- thetas_df$theta1_ns
results$theta2_ns[idx] <- thetas_df$theta2_ns
results$theta3_ns[idx] <- thetas_df$theta3_ns
results$theta4_ns[idx] <- thetas_df$theta4_ns

n_patched <- sum(!is.na(results$theta3_ns))
cat(sprintf("Patched %d genes with fitted NS thetas\n", n_patched))

saveRDS(results, results_file)
cat("Saved:", results_file, "\n")

if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
cat("Done!\n")
