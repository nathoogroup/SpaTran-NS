#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: Run One Replicate
# =============================================================================
# Simulates spatial expression for all genes under their assigned hypothesis
# (H0 = stationary, H1 = nonstationary), fits both models, and computes
# BF_{10} for each gene. Designed for SLURM array jobs.
#
# Usage:
#   Rscript simulation/run_replicate.R <dataset_name> <replicate_id> [n_cores] [config_dir] [output_dir] [n_genes_subset]
#
# SLURM array mode:
#   SLURM_ARRAY_TASK_ID overrides <replicate_id> when set.
#
# Examples:
#   Rscript simulation/run_replicate.R HumanBreastCancerILC 1 8
#   Rscript simulation/run_replicate.R HumanOvarianCancer 1 4 simulation/configs simulation/output 100
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_replicate.R <dataset_name> <replicate_id> ",
       "[n_cores] [config_dir] [output_dir] [n_genes_subset]")
}

dataset_name   <- args[1]
replicate_id   <- as.integer(args[2])
n_cores        <- if (length(args) >= 3) as.integer(args[3]) else parallel::detectCores()
config_dir     <- if (length(args) >= 4) args[4] else "simulation/configs"
output_dir     <- if (length(args) >= 5) args[5] else "simulation/output"
n_genes_subset <- if (length(args) >= 6) as.integer(args[6]) else 0L

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
if (task_id > 0) replicate_id <- task_id

dir.create(file.path(output_dir, dataset_name), recursive = TRUE, showWarnings = FALSE)

# Skip if already completed
out_file <- file.path(output_dir, dataset_name,
                      sprintf("replicate_%03d.rds", replicate_id))
if (file.exists(out_file)) {
  cat("Replicate", replicate_id, "already exists, skipping.\n")
  quit(save = "no", status = 0)
}

cat("=============================================================\n")
cat("Simulation Replicate\n")
cat("=============================================================\n")
cat("Dataset:       ", dataset_name, "\n")
cat("Replicate:     ", replicate_id, "\n")
cat("Cores:         ", n_cores, "\n")
cat("Config dir:    ", config_dir, "\n")
cat("Output dir:    ", output_dir, "\n")
if (n_genes_subset > 0) cat("Gene subset:   ", n_genes_subset, "\n")
cat("Start time:    ", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# --- Load packages -----------------------------------------------------------

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(dplyr)
  library(doParallel)
  library(foreach)
})

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "simulation")
source(file.path(script_dir, "simulation_functions.R"))

INLA::inla.setOption(num.threads = "1:1")

# --- Load simulation config --------------------------------------------------

config_file <- file.path(config_dir, paste0(dataset_name, "_sim_config.rds"))
if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file,
       "\nRun prepare_simulation.R first.")
}

cat("Loading configuration...\n")
config <- readRDS(config_file)

coords     <- config$coords
mesh       <- config$mesh
sim_params <- config$sim_params

if (n_genes_subset > 0 && n_genes_subset < nrow(sim_params)) {
  n_h0_sub <- min(
    sum(sim_params$hypothesis == "H0"),
    round(n_genes_subset * config$n_h0 / config$n_genes)
  )
  n_h1_sub <- n_genes_subset - n_h0_sub

  set.seed(replicate_id * 1000 + 999)
  h0_idx <- sample(which(sim_params$hypothesis == "H0"), n_h0_sub)
  h1_idx <- sample(which(sim_params$hypothesis == "H1"),
                   min(n_h1_sub, sum(sim_params$hypothesis == "H1")))
  sim_params <- sim_params[sort(c(h0_idx, h1_idx)), ]
  cat(sprintf("  Subset: %d genes (H0=%d, H1=%d)\n",
              nrow(sim_params), length(h0_idx), length(h1_idx)))
}

n_genes <- nrow(sim_params)
n_h0    <- sum(sim_params$hypothesis == "H0")
n_h1    <- sum(sim_params$hypothesis == "H1")

cat(sprintf("  Genes: %d (H0=%d, H1=%d)  |  Spots: %d  |  Mesh: %d nodes\n",
            n_genes, n_h0, n_h1, nrow(coords), mesh$n))

# --- Build SPDE models -------------------------------------------------------

cat("Building SPDE models and projection matrix...\n")
spde_stat <- build_spde_stationary(mesh)
spde_ns   <- build_spde_nonstationary(mesh)
A         <- INLA::inla.spde.make.A(mesh = mesh,
                                     loc = as.matrix(coords[, c("x", "y")]))
s.index   <- INLA::inla.spde.make.index(name = "spatial", n.spde = mesh$n)

# --- Set reproducible seed base ----------------------------------------------

seed_base <- switch(dataset_name,
                    HumanBreastCancerILC = 42000,
                    HumanOvarianCancer   = 43000,
                    99000)

# --- Parallel simulate + fit -------------------------------------------------
# Simulation and model fitting are combined in one parallel loop.
# Each worker simulates gene expression and fits both models for its assigned
# genes. Per-gene seeds (seed_base + replicate_id * 100000 + i) ensure
# reproducibility regardless of worker scheduling order.

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
  "build_spde_stationary", "build_spde_nonstationary",
  "simulate_expression", "fit_and_get_logbf",
  "mesh", "A", "s.index", "sim_params",
  "seed_base", "replicate_id"
), envir = environment())

# Build SPDE models on each worker once at startup — avoids serialising
# large SPDE/sparse-matrix objects across the socket connection.
clusterEvalQ(cl, {
  spde_stat_w <- build_spde_stationary(mesh)
  spde_ns_w   <- build_spde_nonstationary(mesh)
})

log_file <- file.path(output_dir, dataset_name,
                      sprintf("replicate_%03d_gene_log.txt", replicate_id))
cat(sprintf("Simulating and fitting models for %d genes... (logging to %s)\n",
            n_genes, log_file))
t0_fit <- proc.time()

clusterExport(cl, c("n_genes", "log_file"), envir = environment())

# Collect as list (O(n) combine) then rbind once — avoids O(n²) sequential
# rbind that stalls workers on large gene sets.
raw_list <- foreach(
  i              = seq_len(n_genes),
  .combine       = c,
  .multicombine  = TRUE,
  .errorhandling = "pass",
  .packages      = c("INLA", "Matrix", "dplyr")
) %dopar% {
  gene_name <- sim_params$gene_name[i]
  t0_gene   <- proc.time()
  res <- tryCatch({
    # Hard wall-clock cap: prevents stuck INLA optimisers blocking a worker
    setTimeLimit(elapsed = 120L, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
    # Per-gene seed: reproducible regardless of worker scheduling order
    set.seed(seed_base + replicate_id * 100000L + i)
    gene <- sim_params[i, ]
    if (gene$hypothesis == "H0") {
      theta <- c(gene$log_tau, gene$log_kappa)
      y <- simulate_expression(spde_stat_w, theta, A, gene$sigma_eps_sq)
    } else {
      theta <- c(gene$log_tau, gene$log_kappa, gene$theta3, gene$theta4)
      y <- simulate_expression(spde_ns_w, theta, A, gene$sigma_eps_sq)
    }
    log_bf <- fit_and_get_logbf(y, mesh, spde_stat_w, spde_ns_w, A, s.index)
    data.frame(gene_name  = gene_name,
               hypothesis = sim_params$hypothesis[i],
               log_bf_sim = log_bf,
               bf_sim     = exp(log_bf),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(gene_name  = gene_name,
               hypothesis = sim_params$hypothesis[i],
               log_bf_sim = NA_real_,
               bf_sim     = NA_real_,
               stringsAsFactors = FALSE)
  })
  elapsed <- (proc.time() - t0_gene)["elapsed"]
  cat(sprintf("%s | Gene %d/%d %s [%s] | %.1f sec\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), i, n_genes, gene_name,
              sim_params$hypothesis[i], elapsed), file = log_file, append = TRUE)
  list(res)
}

stopCluster(cl)

results_list <- do.call(rbind, Filter(is.data.frame, raw_list))

elapsed_fit <- (proc.time() - t0_fit)[["elapsed"]]
cat(sprintf("  Model fitting complete in %.1f min\n", elapsed_fit / 60))

# --- Save replicate results --------------------------------------------------

if (!is.data.frame(results_list) || nrow(results_list) == 0) {
  stop("No results returned - check for errors.")
}

n_dropped <- n_genes - nrow(results_list)
if (n_dropped > 0) {
  cat(sprintf("  WARNING: %d genes dropped due to errors\n", n_dropped))
}

results_list$replicate_id <- replicate_id
results_list$dataset_name <- dataset_name

saveRDS(results_list, out_file)
cat("Saved:", out_file, "\n")

# --- Quick summary -----------------------------------------------------------

n_ok       <- sum(!is.na(results_list$log_bf_sim))
n_reject_1 <- sum(results_list$bf_sim > 1, na.rm = TRUE)
n_fp_1     <- sum(results_list$hypothesis == "H0" &
                    results_list$bf_sim > 1, na.rm = TRUE)
n_tp_1     <- sum(results_list$hypothesis == "H1" &
                    results_list$bf_sim > 1, na.rm = TRUE)

cat("\n=============================================================\n")
cat("REPLICATE SUMMARY (c = 1)\n")
cat("=============================================================\n")
cat(sprintf("Successful fits:   %d / %d\n", n_ok, n_genes))
cat(sprintf("Rejections (BF>1): %d\n", n_reject_1))
cat(sprintf("  True positives:  %d\n", n_tp_1))
cat(sprintf("  False positives: %d\n", n_fp_1))
if (n_reject_1 > 0) {
  cat(sprintf("  FDR at c=1:      %.3f\n", n_fp_1 / n_reject_1))
}
if (n_h1 > 0) {
  cat(sprintf("  Power at c=1:    %.3f\n", n_tp_1 / n_h1))
}
cat(sprintf("Total time:        %.1f min\n", elapsed_fit / 60))
cat("=============================================================\n")
cat("Done!\n")
