#!/usr/bin/env Rscript
# Quick login-node test of the batched checkpoint approach.
# Runs 10 H1 genes with 2 cores, batch_size=3.
# Usage:  Rscript simulation/test_patch_batches.R

n_cores    <- 2L
n_test     <- 10L
batch_size <- 3L

dataset_name <- "HumanBreastCancerILC"
data_dir     <- "data/spatial_datasets"
results_dir  <- "results"
checkpoint_file <- file.path(results_dir, paste0(dataset_name, "_test_checkpoint.rds"))

cat("=== Batch checkpoint test ===\n")
cat(sprintf("n_test=%d  n_cores=%d  batch_size=%d\n\n", n_test, n_cores, batch_size))

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

source("hpc/analysis_functions.R")

# Load a small slice of H1 genes
results_file <- file.path(results_dir, paste0(dataset_name, "_results.rds"))
results  <- readRDS(results_file)
h1_genes <- results %>%
  filter(!is.na(log_bayes_factor), log_bayes_factor > 0) %>%
  pull(gene_name)
h1_genes <- head(h1_genes, n_test)
cat(sprintf("H1 genes selected: %d\n", length(h1_genes)))

# Load SPE
spe        <- readRDS(file.path(data_dir, paste0(dataset_name, ".rds")))
coords     <- as.data.frame(spatialCoords(spe))
colnames(coords)[1:2] <- c("x", "y")
rd         <- rowData(spe)
gene_names_spe <- if ("gene_name" %in% colnames(rd)) {
  as.character(rd$gene_name)
} else if ("symbol" %in% colnames(rd)) {
  as.character(rd$symbol)
} else {
  rownames(spe)
}

cat(sprintf("gene_names_spe sample: %s\n", paste(head(gene_names_spe, 3), collapse=", ")))
cat(sprintf("h1_genes sample: %s\n", paste(head(h1_genes, 3), collapse=", ")))
h1_genes   <- intersect(h1_genes, gene_names_spe)
cat(sprintf("After intersect: %d genes remain\n", length(h1_genes)))
if (length(h1_genes) == 0) stop("No H1 genes found in SPE — check gene name format mismatch")
h1_idx_spe <- match(h1_genes, gene_names_spe)
counts_mat  <- as.matrix(assay(spe, "counts")[h1_idx_spe, , drop = FALSE])
rm(spe); gc(verbose = FALSE)

mesh <- create_spde_mesh(coords)
cat(sprintf("Mesh nodes: %d\n\n", mesh$n))

# Set up workers
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, suppressPackageStartupMessages({
  library(INLA); library(Matrix); library(dplyr)
  INLA::inla.setOption(num.threads = "1:1")
}))
clusterExport(cl, c(
  "fit_nonstationary_matern", "extract_nonstationary_thetas",
  "normalize_expression", "detrend_expression",
  "counts_mat", "coords", "mesh", "h1_genes"
), envir = environment())

# --- Batched loop with incremental checkpoint --------------------------------
n_genes       <- length(h1_genes)
batches       <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / batch_size))
batch_results <- vector("list", length(batches))
t0            <- proc.time()

for (b in seq_along(batches)) {
  idx_b <- batches[[b]]
  cat(sprintf("Batch %d/%d (genes %d-%d)...\n",
              b, length(batches), idx_b[1], idx_b[length(idx_b)]))

  raw_b <- foreach(
    i              = idx_b,
    .errorhandling = "pass",
    .packages      = c("INLA", "Matrix", "dplyr")
  ) %dopar% {
    gene_name <- h1_genes[i]
    tryCatch({
      setTimeLimit(elapsed = 120L, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
      expr      <- log1p(counts_mat[i, ])
      expr_proc <- normalize_expression(expr)
      expr_proc <- detrend_expression(coords, expr_proc, "polynomial")$detrended
      fit       <- fit_nonstationary_matern(coords, expr_proc, mesh = mesh, verbose = FALSE)
      thetas    <- extract_nonstationary_thetas(fit$result)
      data.frame(gene_name = gene_name,
                 theta1_ns = thetas$theta1, theta2_ns = thetas$theta2,
                 theta3_ns = thetas$theta3, theta4_ns = thetas$theta4,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(gene_name = gene_name,
                 theta1_ns = NA_real_, theta2_ns = NA_real_,
                 theta3_ns = NA_real_, theta4_ns = NA_real_,
                 stringsAsFactors = FALSE)
    })
  }

  # Report any worker errors
  errs <- Filter(function(x) inherits(x, "error") || inherits(x, "condition"), raw_b)
  if (length(errs) > 0)
    cat(sprintf("  Worker errors (%d): %s\n", length(errs), conditionMessage(errs[[1]])))

  batch_results[[b]] <- do.call(rbind, Filter(is.data.frame, raw_b))

  # Incremental checkpoint
  ckpt_so_far <- do.call(rbind, Filter(Negate(is.null), batch_results))
  if (!is.null(ckpt_so_far))
    saveRDS(ckpt_so_far, checkpoint_file)
  cat(sprintf("  -> checkpoint saved (%d genes so far)\n", if (is.null(ckpt_so_far)) 0L else nrow(ckpt_so_far)))
}

stopCluster(cl)

thetas_df <- do.call(rbind, Filter(Negate(is.null), batch_results))
elapsed   <- (proc.time() - t0)[["elapsed"]]

cat(sprintf("\nDone in %.1f sec\n", elapsed))
cat(sprintf("Successful fits: %d / %d\n", sum(!is.na(thetas_df$theta3_ns)), nrow(thetas_df)))
cat(sprintf("Sec per gene (wall, 2 cores): %.1f\n", elapsed / ceiling(n_genes / n_cores)))
cat("\nThetas:\n")
print(thetas_df)

# Clean up test checkpoint
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
cat("\nTest checkpoint removed. DONE.\n")
