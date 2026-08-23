#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: Aggregate Results and Compute FDR/Power Curves
# =============================================================================
# Loads all replicate results for a dataset, computes FDR(c) and average
# power(c) for a grid of Bayes factor thresholds c.
#
# Usage:
#   Rscript simulation/aggregate_results.R <dataset_name> [n_replicates] [output_dir]
#
# Examples:
#   Rscript simulation/aggregate_results.R HumanBreastCancerILC 100
#   Rscript simulation/aggregate_results.R HumanOvarianCancer 100
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript aggregate_results.R <dataset_name> [n_replicates]")

dataset_name <- args[1]
n_replicates <- if (length(args) >= 2) as.integer(args[2]) else 100L
output_dir   <- if (length(args) >= 3) args[3] else "simulation/output"

cat("=============================================================\n")
cat("Simulation Study: Aggregating Results\n")
cat("=============================================================\n")
cat("Dataset:     ", dataset_name, "\n")
cat("Replicates:  ", n_replicates, "\n")
cat("Output dir:  ", output_dir, "\n")
cat("=============================================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
})

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "simulation")
source(file.path(script_dir, "simulation_functions.R"))

# --- Load replicate results --------------------------------------------------

rep_dir <- file.path(output_dir, dataset_name)
cat("Loading replicates from:", rep_dir, "\n")

replicate_results <- list()
missing_reps <- c()

for (r in seq_len(n_replicates)) {
  rep_file <- file.path(rep_dir, sprintf("replicate_%03d.rds", r))
  if (file.exists(rep_file)) {
    replicate_results[[length(replicate_results) + 1]] <- readRDS(rep_file)
  } else {
    missing_reps <- c(missing_reps, r)
  }
}

n_loaded <- length(replicate_results)
cat(sprintf("  Loaded: %d / %d replicates\n", n_loaded, n_replicates))

if (length(missing_reps) > 0) {
  if (length(missing_reps) <= 20) {
    cat("  Missing:", paste(missing_reps, collapse = ", "), "\n")
  } else {
    cat("  Missing:", length(missing_reps), "replicates\n")
  }
}

if (n_loaded < 2) stop("Need at least 2 completed replicates to aggregate.")

# --- Basic info from first replicate -----------------------------------------

first_rep <- replicate_results[[1]]
n_genes <- nrow(first_rep)
n_h0    <- sum(first_rep$hypothesis == "H0")
n_h1    <- sum(first_rep$hypothesis == "H1")

cat(sprintf("\n  Genes per replicate: %d (H0=%d, H1=%d)\n", n_genes, n_h0, n_h1))

# --- Define threshold grid ---------------------------------------------------

c_grid <- sort(unique(c(
  0,
  seq(0.1, 1.0, by = 0.1),
  seq(1.5, 5.0, by = 0.5),
  seq(6, 10, by = 1),
  seq(15, 50, by = 5),
  seq(60, 100, by = 10),
  150, 200, 300, 500, 1000, 5000, 10000
)))

# --- Compute FDR and power ---------------------------------------------------

cat("\nComputing FDR and power curves...\n")
t0 <- proc.time()

curves <- compute_fdr_power(replicate_results, c_grid)

elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("  Completed in %.1f sec\n", elapsed))

# --- Report key thresholds ---------------------------------------------------

cat("\n=============================================================\n")
cat("FDR AND POWER AT KEY THRESHOLDS\n")
cat("=============================================================\n")

key_thresholds <- c(1, 3, 10, 30, 100)
for (ct in key_thresholds) {
  row <- curves[which.min(abs(curves$c_threshold - ct)), ]
  cat(sprintf("  c = %5.0f : FDR = %.4f (sd=%.4f)  |  Power = %.4f (sd=%.4f)\n",
              ct, row$fdr, row$fdr_sd, row$avg_power, row$avg_power_sd))
}

# --- Also compute per-replicate FDR at c=1 for reporting ---------------------

fdr_at_1 <- numeric(n_loaded)
power_at_1 <- numeric(n_loaded)
for (r in seq_len(n_loaded)) {
  rep_data <- replicate_results[[r]] %>% filter(!is.na(bf_sim))
  fp <- sum(rep_data$hypothesis == "H0" & rep_data$bf_sim > 1)
  tp <- sum(rep_data$hypothesis == "H1" & rep_data$bf_sim > 1)
  disc <- fp + tp
  fdr_at_1[r] <- if (disc > 0) fp / disc else 0
  power_at_1[r] <- tp / sum(rep_data$hypothesis == "H1")
}

cat(sprintf("\n  FDR at c=1:   mean=%.4f, median=%.4f, sd=%.4f, [%.4f, %.4f]\n",
            mean(fdr_at_1), median(fdr_at_1), sd(fdr_at_1),
            quantile(fdr_at_1, 0.025), quantile(fdr_at_1, 0.975)))
cat(sprintf("  Power at c=1: mean=%.4f, median=%.4f, sd=%.4f, [%.4f, %.4f]\n",
            mean(power_at_1), median(power_at_1), sd(power_at_1),
            quantile(power_at_1, 0.025), quantile(power_at_1, 0.975)))

# --- Save curves -------------------------------------------------------------

curves$dataset_name <- dataset_name
curves$n_replicates <- n_loaded
curves$n_genes      <- n_genes
curves$n_h0         <- n_h0
curves$n_h1         <- n_h1

curves_file <- file.path(output_dir, paste0(dataset_name, "_fdr_power_curves.rds"))
saveRDS(curves, curves_file)
cat("\nSaved curves:", curves_file, "\n")

curves_csv <- file.path(output_dir, paste0(dataset_name, "_fdr_power_curves.csv"))
write.csv(curves, curves_csv, row.names = FALSE)
cat("Saved CSV:  ", curves_csv, "\n")

# --- Save per-replicate summaries for diagnostics ----------------------------

rep_summary <- data.frame(
  replicate = seq_len(n_loaded),
  fdr_c1    = fdr_at_1,
  power_c1  = power_at_1
)
rep_summary_file <- file.path(output_dir,
                              paste0(dataset_name, "_replicate_summaries.csv"))
write.csv(rep_summary, rep_summary_file, row.names = FALSE)
cat("Saved replicate summaries:", rep_summary_file, "\n")

cat("\nDone!\n")
