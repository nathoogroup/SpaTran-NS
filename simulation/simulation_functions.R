# =============================================================================
# Simulation Study: Core Functions
# =============================================================================
# Functions for simulating spatial gene expression under stationary (H0) and
# nonstationary (H1) Matern models, fitting models, and computing FDR/power.
#
# Uses the INLA SPDE framework (same as the main analysis) for both
# simulation and inference.
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(dplyr)
})

# --- Derive SPDE theta from variance components -----------------------------
#
# For inla.spde2.matern with alpha=2 (nu=1, d=2):
#   sigma_b^2 = 1 / (4*pi * kappa^2 * tau^2)
#   range     = sqrt(8) / kappa
#
# Returns theta = c(log(tau), log(kappa))

derive_theta_stationary <- function(sigma_b_sq, spatial_range) {
  kappa <- sqrt(8) / spatial_range
  tau   <- 1 / sqrt(4 * pi * kappa^2 * sigma_b_sq)
  c(log(tau), log(kappa))
}

# --- Build SPDE models -------------------------------------------------------
# Pre-build both stationary and nonstationary SPDE models for a given mesh.
# These are reused across all genes.

build_spde_stationary <- function(mesh) {
  INLA::inla.spde2.matern(
    mesh    = mesh,
    alpha   = 2,
    B.tau   = matrix(c(0, 1, 0), nrow = 1, ncol = 3),
    B.kappa = matrix(c(0, 0, 1), nrow = 1, ncol = 3)
  )
}

build_spde_nonstationary <- function(mesh) {
  n_mesh <- mesh$n
  mesh_x <- as.vector(scale(mesh$loc[, 1]))
  mesh_y <- as.vector(scale(mesh$loc[, 2]))

  B.tau <- matrix(c(0, 1, 0, 0, 0), nrow = n_mesh, ncol = 5, byrow = TRUE)
  B.kappa <- cbind(
    rep(0, n_mesh),
    rep(0, n_mesh),
    rep(1, n_mesh),
    mesh_x,
    mesh_y
  )

  INLA::inla.spde2.matern(
    mesh    = mesh,
    alpha   = 2,
    B.tau   = B.tau,
    B.kappa = B.kappa
  )
}

# --- Simulate gene expression ------------------------------------------------
#
# Simulate from a Gaussian spatial process on the mesh, project to observation
# locations, and add independent Gaussian noise.
#
# spde  : INLA SPDE model object (stationary or nonstationary)
# theta : hyperparameter vector (length 2 for stationary, 4 for nonstationary)
# A     : projection matrix (n_obs x n_mesh)
# sigma_eps_sq : noise variance

simulate_expression <- function(spde, theta, A, sigma_eps_sq) {
  Q <- INLA::inla.spde2.precision(spde, theta = theta)
  z <- as.vector(INLA::inla.qsample(n = 1, Q = Q))
  u <- as.vector(A %*% z)
  y <- u + rnorm(nrow(A), mean = 0, sd = sqrt(sigma_eps_sq))
  return(y)
}

# --- Fit both models and return log BF ---------------------------------------
#
# Fit stationary and nonstationary Matern models to simulated expression data.
# Returns log(BF_10). Uses pre-built SPDE models, mesh, projection matrix A,
# and spatial index to avoid redundant recomputation across genes.

fit_and_get_logbf <- function(expr, mesh, spde_stat, spde_ns, A, s.index) {
  stack <- INLA::inla.stack(
    data    = list(y = expr),
    A       = list(A, 1),
    effects = list(s.index, list(intercept = rep(1, length(expr)))),
    tag     = "est"
  )

  stack_data_stat <- INLA::inla.stack.data(stack, spde = spde_stat)
  stack_data_ns   <- INLA::inla.stack.data(stack, spde = spde_ns)
  stack_A         <- INLA::inla.stack.A(stack)

  fit_stat <- tryCatch({
    INLA::inla(
      formula           = y ~ -1 + intercept + f(spatial, model = spde_stat),
      data              = stack_data_stat,
      family            = "gaussian",
      control.predictor = list(A = stack_A, compute = FALSE),
      control.compute   = list(mlik = TRUE),
      verbose           = FALSE
    )
  }, error = function(e) NULL)

  fit_ns <- tryCatch({
    INLA::inla(
      formula           = y ~ -1 + intercept + f(spatial, model = spde_ns),
      data              = stack_data_ns,
      family            = "gaussian",
      control.predictor = list(A = stack_A, compute = FALSE),
      control.compute   = list(mlik = TRUE),
      verbose           = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(fit_stat) || is.null(fit_ns)) return(NA_real_)

  log_bf <- fit_ns$mlik[1, 1] - fit_stat$mlik[1, 1]
  return(log_bf)
}

# --- Safe wrapper for parallel execution -------------------------------------

simulate_and_fit_gene <- function(gene_idx, sim_config, mesh,
                                   spde_stat, spde_ns, A, s.index) {
  gene <- sim_config[gene_idx, ]
  gene_name <- gene$gene_name

  tryCatch({
    if (gene$hypothesis == "H0") {
      theta <- c(gene$log_tau, gene$log_kappa)
      y <- simulate_expression(spde_stat, theta, A, gene$sigma_eps_sq)
    } else {
      theta <- c(gene$log_tau, gene$log_kappa, gene$theta3, gene$theta4)
      y <- simulate_expression(spde_ns, theta, A, gene$sigma_eps_sq)
    }

    log_bf <- fit_and_get_logbf(y, mesh, spde_stat, spde_ns, A, s.index)

    data.frame(
      gene_name  = gene_name,
      hypothesis = gene$hypothesis,
      log_bf_sim = log_bf,
      bf_sim     = exp(log_bf),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      gene_name  = gene_name,
      hypothesis = gene$hypothesis,
      log_bf_sim = NA_real_,
      bf_sim     = NA_real_,
      error_msg  = e$message,
      stringsAsFactors = FALSE
    )
  })
}

# --- Compute FDR and power curves -------------------------------------------
#
# Given BF results from all replicates, compute FDR(c) and average power(c)
# for a grid of decision thresholds c.
#
# replicate_results : list of data frames (one per replicate), each with
#                     columns: gene_name, hypothesis, bf_sim
# c_grid            : vector of threshold values

compute_fdr_power <- function(replicate_results, c_grid = NULL) {
  if (is.null(c_grid)) {
    c_grid <- sort(unique(c(
      0, seq(0.1, 1.0, by = 0.1),
      seq(1.5, 5, by = 0.5),
      seq(6, 10, by = 1),
      seq(15, 50, by = 5),
      seq(60, 100, by = 10),
      150, 200, 300, 500, 1000, 5000, 10000
    )))
  }

  n_reps <- length(replicate_results)
  n_c    <- length(c_grid)

  fdr_mat   <- matrix(NA, nrow = n_reps, ncol = n_c)
  power_mat <- matrix(NA, nrow = n_reps, ncol = n_c)

  for (r in seq_len(n_reps)) {
    rep_data <- replicate_results[[r]] %>% filter(!is.na(bf_sim))

    h0_bfs <- rep_data$bf_sim[rep_data$hypothesis == "H0"]
    h1_bfs <- rep_data$bf_sim[rep_data$hypothesis == "H1"]
    n_h1   <- length(h1_bfs)

    for (k in seq_len(n_c)) {
      fp   <- sum(h0_bfs > c_grid[k])
      tp   <- sum(h1_bfs > c_grid[k])
      disc <- fp + tp

      fdr_mat[r, k]   <- if (disc > 0) fp / disc else 0
      power_mat[r, k] <- if (n_h1 > 0) tp / n_h1 else 0
    }
  }

  fdr_curve   <- colMeans(fdr_mat, na.rm = TRUE)
  power_curve <- colMeans(power_mat, na.rm = TRUE)

  fdr_sd   <- apply(fdr_mat, 2, sd, na.rm = TRUE)
  power_sd <- apply(power_mat, 2, sd, na.rm = TRUE)

  data.frame(
    c_threshold  = c_grid,
    fdr          = fdr_curve,
    fdr_sd       = fdr_sd,
    avg_power    = power_curve,
    avg_power_sd = power_sd
  )
}

cat("Simulation functions loaded.\n")
