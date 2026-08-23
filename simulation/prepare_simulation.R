#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: Prepare Simulation Configuration
# =============================================================================
# One-time setup: loads the actual data and analysis results, classifies genes
# as H0 (stationary, BF <= 1) or H1 (nonstationary, BF > 1), extracts model
# parameters for simulation, and saves a configuration RDS file.
#
# Usage:
#   Rscript simulation/prepare_simulation.R <dataset_name> [data_dir] [results_dir] [output_dir]
#
# Example:
#   Rscript simulation/prepare_simulation.R HumanBreastCancerILC
#   Rscript simulation/prepare_simulation.R HumanOvarianCancer
# =============================================================================

args       <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript prepare_simulation.R <dataset_name>")

dataset_name <- args[1]
data_dir     <- if (length(args) >= 2) args[2] else "data/spatial_datasets"
results_dir  <- if (length(args) >= 3) args[3] else "results"
output_dir   <- if (length(args) >= 4) args[4] else "simulation/configs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat("Simulation Study: Preparing Configuration\n")
cat("=============================================================\n")
cat("Dataset:         ", dataset_name, "\n")
cat("Data directory:  ", data_dir, "\n")
cat("Results directory:", results_dir, "\n")
cat("Output directory: ", output_dir, "\n")
cat("=============================================================\n\n")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(INLA)
  library(dplyr)
})

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "simulation")
source(file.path(script_dir, "simulation_functions.R"))

# Also need create_spde_mesh from analysis_functions.R
hpc_dir <- file.path(dirname(script_dir), "hpc")
if (!file.exists(file.path(hpc_dir, "analysis_functions.R"))) {
  hpc_dir <- "hpc"
}
source(file.path(hpc_dir, "analysis_functions.R"))

# --- Load analysis results ---------------------------------------------------

results_file <- file.path(results_dir, paste0(dataset_name, "_results.rds"))
if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}
cat("Loading results:", results_file, "\n")
results <- readRDS(results_file)

cat(sprintf("  Total genes in results: %d\n", nrow(results)))
cat(sprintf("  Genes with valid BF:    %d\n", sum(!is.na(results$log_bayes_factor))))

# --- Load spatial dataset (for coordinates) ----------------------------------

spe_file <- file.path(data_dir, paste0(dataset_name, ".rds"))
if (!file.exists(spe_file)) {
  stop("Dataset file not found: ", spe_file)
}
cat("Loading spatial data:", spe_file, "\n")
spe <- readRDS(spe_file)

coords <- as.data.frame(spatialCoords(spe))
colnames(coords)[1:2] <- c("x", "y")
n_spots <- nrow(coords)
cat(sprintf("  Spots: %d\n", n_spots))

rm(spe); gc()

# --- Build mesh (same procedure as main analysis) ----------------------------

cat("Building SPDE mesh...\n")
mesh <- create_spde_mesh(coords)
cat(sprintf("  Mesh nodes: %d\n", mesh$n))

A <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))

# --- Filter to genes with valid parameters -----------------------------------

results_valid <- results %>%
  filter(
    !is.na(log_bayes_factor),
    !is.na(sigma_b_sq_stationary),
    !is.na(sigma_eps_sq_stationary),
    !is.na(range_stationary),
    sigma_b_sq_stationary > 0,
    sigma_eps_sq_stationary > 0,
    range_stationary > 0,
    is.finite(sigma_b_sq_stationary),
    is.finite(sigma_eps_sq_stationary),
    is.finite(range_stationary)
  )

cat(sprintf("  Valid genes after filtering: %d\n", nrow(results_valid)))

# --- Classify genes as H0 or H1 ---------------------------------------------

results_valid$hypothesis <- ifelse(
  results_valid$log_bayes_factor > 0, "H1", "H0"
)

n_h0 <- sum(results_valid$hypothesis == "H0")
n_h1 <- sum(results_valid$hypothesis == "H1")
cat(sprintf("\n  H0 (stationary, BF <= 1): %d genes (%.1f%%)\n",
            n_h0, 100 * n_h0 / nrow(results_valid)))
cat(sprintf("  H1 (nonstationary, BF > 1): %d genes (%.1f%%)\n",
            n_h1, 100 * n_h1 / nrow(results_valid)))

# --- Derive simulation parameters for each gene -----------------------------

cat("\nDeriving simulation parameters...\n")

theta_cols <- c("theta1_ns", "theta2_ns", "theta3_ns", "theta4_ns")
if (!all(theta_cols %in% colnames(results_valid))) {
  stop("Results file is missing fitted nonstationary theta columns (",
       paste(theta_cols, collapse = ", "), "). ",
       "Re-run the analysis pipeline to produce these columns before preparing the simulation config.")
}

sim_params <- results_valid %>%
  mutate(
    log_kappa    = log(sqrt(8) / range_stationary),
    log_tau      = log(1 / sqrt(4 * pi * (sqrt(8) / range_stationary)^2 *
                                sigma_b_sq_stationary)),
    sigma_eps_sq = sigma_eps_sq_stationary,
    theta3       = 0,
    theta4       = 0
  )

# For H1 genes: use fitted thetas. Drop any H1 gene missing finite thetas.
h1_idx <- which(sim_params$hypothesis == "H1")
n_use_fitted  <- 0L
n_dropped_h1  <- 0L

for (i in h1_idx) {
  t1 <- sim_params$theta1_ns[i]
  t2 <- sim_params$theta2_ns[i]
  t3 <- sim_params$theta3_ns[i]
  t4 <- sim_params$theta4_ns[i]
  if (all(is.finite(c(t1, t2, t3, t4)))) {
    sim_params$log_tau[i]   <- t1
    sim_params$log_kappa[i] <- t2
    sim_params$theta3[i]    <- t3
    sim_params$theta4[i]    <- t4
    n_use_fitted <- n_use_fitted + 1L
  } else {
    sim_params$hypothesis[i] <- "excluded"
    n_dropped_h1 <- n_dropped_h1 + 1L
  }
}

cat(sprintf("  Using fitted nonstationary thetas for %d H1 genes\n", n_use_fitted))
if (n_dropped_h1 > 0) {
  dropped_names <- sim_params$gene_name[sim_params$hypothesis == "excluded"]
  cat(sprintf("  Excluded %d H1 gene(s) with missing fitted thetas: %s\n",
              n_dropped_h1, paste(dropped_names, collapse = ", ")))
}

sim_params <- sim_params %>% filter(hypothesis != "excluded")
n_h0 <- sum(sim_params$hypothesis == "H0")
n_h1 <- sum(sim_params$hypothesis == "H1")

sim_params <- sim_params %>%
  select(gene_name, hypothesis, log_bayes_factor, bayes_factor,
         sigma_b_sq_stationary, sigma_eps_sq_stationary, range_stationary,
         log_tau, log_kappa, sigma_eps_sq, theta3, theta4)

# --- Build and save configuration --------------------------------------------

cat("\nSaving simulation configuration...\n")

config <- list(
  dataset_name        = dataset_name,
  n_spots             = n_spots,
  n_genes             = nrow(sim_params),
  n_h0                = n_h0,
  n_h1                = n_h1,
  n_use_fitted_thetas = n_use_fitted,
  n_dropped_h1        = n_dropped_h1,
  coords              = coords,
  mesh                = mesh,
  sim_params          = sim_params,
  created             = Sys.time()
)

config_file <- file.path(output_dir, paste0(dataset_name, "_sim_config.rds"))
saveRDS(config, config_file)
cat("  Saved:", config_file, "\n")

cat(sprintf("\n  Config size: %.1f MB\n",
            file.info(config_file)$size / 1e6))

# --- Print summary -----------------------------------------------------------

cat("\n=============================================================\n")
cat("CONFIGURATION SUMMARY\n")
cat("=============================================================\n")
cat("Dataset:             ", dataset_name, "\n")
cat("Total genes:         ", config$n_genes, "\n")
cat("H0 (stationary):     ", config$n_h0, "\n")
cat("H1 (nonstationary):  ", config$n_h1, "\n")
cat("Spots per tissue:    ", config$n_spots, "\n")
cat("Mesh nodes:          ", config$mesh$n, "\n")
cat("Config file:         ", config_file, "\n")
cat("=============================================================\n")
cat("Done!\n")
