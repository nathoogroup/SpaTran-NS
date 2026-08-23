#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: Local Test (CLI entry point)
# =============================================================================
# Runs the local simulation test. When executed as a script, parses
# command-line args and calls run_local_simulation_test().
#
# Usage:
#   Rscript simulation/test_simulation_local.R [n_genes] [n_reps] [dataset]
#
# Defaults: 20 genes, 3 replicates, HumanBreastCancerILC
#
# To use as a function from R:
#   source("simulation/run_local_test.R")
#   run_local_simulation_test(n_genes = 10, n_reps = 2)
#   run_local_simulation_test(10, 2, "HumanOvarianCancer")
# =============================================================================

# When run via Rscript, switch to project root (parent of simulation/)
script_path <- tryCatch(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  error = function(e) ""
)
if (nzchar(script_path)) {
  proj_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
  if (file.exists(file.path(proj_root, "hpc", "run_analysis.R"))) {
    setwd(proj_root)
  }
}

source(file.path("simulation", "run_local_test.R"))

args <- commandArgs(trailingOnly = TRUE)
n_genes <- if (length(args) >= 1) as.integer(args[1]) else 20L
n_reps  <- if (length(args) >= 2) as.integer(args[2]) else 3L
dataset <- if (length(args) >= 3) args[3] else "HumanBreastCancerILC"

run_local_simulation_test(
  n_genes = n_genes,
  n_reps = n_reps,
  dataset = dataset,
  project_root = getwd(),
  install_pkgs = TRUE,
  verbose = TRUE
)
