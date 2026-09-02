#!/usr/bin/env Rscript
# =============================================================================
# Login-node smoke test: run BF analysis on 5 genes from the first dataset
# Usage: Rscript hpc/test_analysis.R [n_genes]
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
n_test <- if (length(args) >= 1) as.integer(args[1]) else 5

args_all    <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  SCRIPT_DIR <- getwd()
}
PROJECT <- dirname(SCRIPT_DIR)

cat("=== SpaTran-NS Login-Node Test ===\n")
cat("Project dir:", PROJECT, "\n")
cat("Genes to test:", n_test, "\n\n")

# ---- Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(INLA)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
})

source(file.path(PROJECT, "hpc", "analysis_functions.R"))

# ---- Verify INLA ------------------------------------------------------------
cat("Checking INLA... ")
INLA::inla.mesh.1d(seq(0, 1, length.out = 5))
cat("OK\n\n")

# ---- Pick first available dataset -------------------------------------------
data_dir <- file.path(PROJECT, "data", "spatial_datasets")
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files <- rds_files[!grepl("download_log|manifest", rds_files)]

if (length(rds_files) == 0) stop("No datasets in ", data_dir)

ds_file <- rds_files[1]
cat("Dataset:", basename(ds_file), "\n")
spe <- readRDS(ds_file)
cat(sprintf("Dimensions: %d genes x %d spots\n\n", nrow(spe), ncol(spe)))

# ---- Filter -----------------------------------------------------------------
counts_mat <- counts(spe)
n_spots    <- ncol(spe)
keep <- (Matrix::rowSums(counts_mat > 0) >= ceiling(n_spots * 0.005)) &
        (Matrix::rowSums(counts_mat) >= 3)
spe  <- spe[keep, ]
cat(sprintf("After filter: %d genes — testing top %d by total count\n\n",
            nrow(spe), n_test))

# Select top expressed genes for a fast, representative test
top_idx  <- order(Matrix::rowSums(counts(spe)), decreasing = TRUE)[seq_len(n_test)]
rd       <- rowData(spe)
if ("gene_name" %in% colnames(rd)) {
  gnames <- as.character(rd$gene_name[top_idx])
} else if ("symbol" %in% colnames(rd)) {
  gnames <- as.character(rd$symbol[top_idx])
} else {
  gnames <- rownames(spe)[top_idx]
}

coords <- as.data.frame(spatialCoords(spe))
colnames(coords)[1:2] <- c("x", "y")

# ---- Run analysis -----------------------------------------------------------
cat(rep("-", 60), "\n", sep = "")
results <- vector("list", n_test)
t0      <- proc.time()

for (i in seq_len(n_test)) {
  cat(sprintf("[%d/%d] %s ... ", i, n_test, gnames[i]))
  flush.console()
  expr <- log1p(counts(spe)[top_idx[i], ])

  results[[i]] <- tryCatch(
    analyze_gene_nonstationarity(
      coords         = coords,
      expr           = expr,
      gene_name      = gnames[i],
      normalize      = TRUE,
      detrend        = TRUE,
      detrend_method = "polynomial",
      verbose        = FALSE
    ),
    error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
  )

  if (!is.null(results[[i]])) {
    r <- results[[i]]
    cat(sprintf("log(BF)=%+.2f  [%s]\n",
                r$log_bayes_factor, r$bf_interpretation))
  }
}

elapsed <- (proc.time() - t0)[["elapsed"]]

# ---- Summary ----------------------------------------------------------------
df <- do.call(rbind, Filter(Negate(is.null), results))
rownames(df) <- NULL

cat("\n", rep("=", 60), "\n", sep = "")
cat("RESULTS\n")
cat(rep("=", 60), "\n")
print(df[, c("gene_name", "log_bayes_factor", "bayes_factor", "bf_interpretation")],
      row.names = FALSE)

cat(sprintf("\nTotal time: %.1f min (%.1f sec/gene)\n",
            elapsed / 60, elapsed / n_test))
cat("Done.\n")
