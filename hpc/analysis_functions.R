# =============================================================================
# Spatial Nonstationarity Analysis Functions
# =============================================================================
# Implements Bayes factor comparison between stationary and non-stationary
# Matern covariance kernels via the SPDE approach (Lindgren et al. 2011).
#
# Non-stationary model: Ingebrigtsen et al. (2014) spatially-varying
# log(kappa(s)) = theta2 + theta3*x(s) + theta4*y(s) via B.kappa matrix.
#
# Both models use inla.spde2.matern with the same prior family, ensuring
# a valid and consistent Bayes factor from the INLA marginal log-likelihoods.
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(dplyr)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

# --- Normalization -----------------------------------------------------------

normalize_expression <- function(expr, method = "quantile") {
  if (method == "quantile") {
    expr_rank <- rank(expr, ties.method = "average")
    return(qnorm((expr_rank - 0.5) / length(expr)))
  } else if (method == "zscore") {
    return((expr - mean(expr)) / sd(expr))
  } else {
    return(expr)
  }
}

# --- Detrending (mean trend only) -------------------------------------------

detrend_expression <- function(coords, expr, method = "polynomial") {
  if (method == "loess") {
    loess_fit <- loess(expr ~ x + y, data = coords, span = 0.5)
    trend <- predict(loess_fit, newdata = coords)
    return(list(detrended = expr - trend, trend = trend))
  } else if (method == "polynomial") {
    x_c <- coords$x - mean(coords$x)
    y_c <- coords$y - mean(coords$y)
    poly_fit <- lm(expr ~ x_c + y_c + I(x_c^2) + I(y_c^2) + I(x_c * y_c),
                   data = data.frame(expr = expr, x_c = x_c, y_c = y_c))
    trend <- predict(poly_fit)
    return(list(detrended = expr - trend, trend = trend))
  } else {
    return(list(detrended = expr, trend = rep(0, length(expr))))
  }
}

# --- Mesh creation -----------------------------------------------------------

create_spde_mesh <- function(coords, max.edge = NULL, cutoff = NULL) {
  x_range <- diff(range(coords$x))
  y_range <- diff(range(coords$y))
  max_range <- max(x_range, y_range)
  min_range <- min(x_range, y_range)

  if (is.null(max.edge)) max.edge <- c(max_range / 5, max_range / 3)
  if (is.null(cutoff))  cutoff  <- min_range / 20

  INLA::inla.mesh.2d(
    loc      = as.matrix(coords[, c("x", "y")]),
    max.edge = max.edge,
    cutoff   = cutoff
  )
}

# --- Stationary Matern model (M0) -------------------------------------------
#
# SPDE: (kappa^2 - Delta)^(alpha/2) (tau * u(s)) = W(s)
#   log(tau(s))   = theta1  [constant]
#   log(kappa(s)) = theta2  [constant]
#
# B.tau   = [0, 1, 0]   (columns: W, theta1, theta2)
# B.kappa = [0, 0, 1]
# Hyperparameters: Theta1 for spatial (log tau), Theta2 for spatial (log kappa)

fit_stationary_matern <- function(coords, expr, mesh = NULL, verbose = FALSE) {
  if (is.null(mesh)) mesh <- create_spde_mesh(coords)

  spde <- INLA::inla.spde2.matern(
    mesh    = mesh,
    alpha   = 2,
    B.tau   = matrix(c(0, 1, 0), nrow = 1, ncol = 3),
    B.kappa = matrix(c(0, 0, 1), nrow = 1, ncol = 3)
  )

  A       <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))
  s.index <- INLA::inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)

  stack <- INLA::inla.stack(
    data    = list(y = expr),
    A       = list(A, 1),
    effects = list(s.index, list(intercept = rep(1, length(expr)))),
    tag     = "est"
  )

  result <- INLA::inla(
    formula           = y ~ -1 + intercept + f(spatial, model = spde),
    data              = INLA::inla.stack.data(stack, spde = spde),
    family            = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute   = list(config = TRUE, dic = TRUE, cpo = TRUE),
    verbose           = verbose
  )

  list(result = result, mesh = mesh, spde = spde, stack = stack)
}

# --- Non-stationary Matern model (M1) ---------------------------------------
#
# Ingebrigtsen et al. (2014): spatially-varying range parameter.
#   log(tau(s))   = theta1                              [constant]
#   log(kappa(s)) = theta2 + theta3*x(s) + theta4*y(s) [linear in space]
#
# B.tau  [n_mesh x 5]: each row = c(0, 1, 0, 0, 0)
# B.kappa[n_mesh x 5]: row i  = c(0, 0, 1, x_i, y_i)
#   where x_i, y_i are scaled coordinates of mesh node i
# Hyperparameters: Theta1..Theta4 for spatial

fit_nonstationary_matern <- function(coords, expr, mesh = NULL, verbose = FALSE) {
  if (is.null(mesh)) mesh <- create_spde_mesh(coords)

  n_mesh <- mesh$n
  mesh_x <- as.vector(scale(mesh$loc[, 1]))
  mesh_y <- as.vector(scale(mesh$loc[, 2]))

  # Columns: [W-coeff, theta1-coeff, theta2-coeff, theta3-coeff, theta4-coeff]
  B.tau   <- matrix(c(0, 1, 0, 0, 0), nrow = n_mesh, ncol = 5, byrow = TRUE)
  B.kappa <- cbind(
    rep(0, n_mesh),  # W coefficient (unused)
    rep(0, n_mesh),  # theta1 (applies to tau, not kappa)
    rep(1, n_mesh),  # theta2: log(kappa) intercept
    mesh_x,          # theta3: log(kappa) x-slope
    mesh_y           # theta4: log(kappa) y-slope
  )

  spde <- INLA::inla.spde2.matern(
    mesh    = mesh,
    alpha   = 2,
    B.tau   = B.tau,
    B.kappa = B.kappa
  )

  A       <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))
  s.index <- INLA::inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)

  stack <- INLA::inla.stack(
    data    = list(y = expr),
    A       = list(A, 1),
    effects = list(s.index, list(intercept = rep(1, length(expr)))),
    tag     = "est"
  )

  result <- INLA::inla(
    formula           = y ~ -1 + intercept + f(spatial, model = spde),
    data              = INLA::inla.stack.data(stack, spde = spde),
    family            = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute   = list(config = TRUE, dic = TRUE, cpo = TRUE),
    verbose           = verbose
  )

  list(result = result, mesh = mesh, spde = spde, stack = stack)
}

# --- Extract variance components --------------------------------------------
#
# For inla.spde2.matern with alpha=2 (nu=1, d=2):
#   sigma_b^2 = 1 / (4*pi * kappa^2 * tau^2)
#   range     = sqrt(8) / kappa  (practical range where correlation ~ 0.13)
#
# Stationary : kappa = exp(Theta2), tau = exp(Theta1)
# Non-stationary: reports values at domain centre (x_scaled=0, y_scaled=0),
#   where log(kappa) = theta2 (the intercept), so kappa = exp(Theta2)

extract_variance_components <- function(fit_result) {
  hyper  <- fit_result$summary.hyperpar
  log_ml <- fit_result$mlik[1, 1]

  theta1_idx <- grep("Theta1 for spatial", rownames(hyper), ignore.case = TRUE)
  theta2_idx <- grep("Theta2 for spatial", rownames(hyper), ignore.case = TRUE)

  sigma_b_sq    <- NA_real_
  spatial_range <- NA_real_

  if (length(theta1_idx) > 0 && length(theta2_idx) > 0) {
    tau   <- exp(hyper$mean[theta1_idx[1]])
    kappa <- exp(hyper$mean[theta2_idx[1]])
    sigma_b_sq    <- 1 / (4 * pi * kappa^2 * tau^2)
    spatial_range <- sqrt(8) / kappa
  }

  noise_idx    <- grep("Precision for the Gaussian", rownames(hyper), ignore.case = TRUE)
  sigma_eps_sq <- if (length(noise_idx) > 0) 1 / hyper$mean[noise_idx[1]] else NA_real_

  prop_spatial <- if (!is.na(sigma_b_sq) && !is.na(sigma_eps_sq))
    sigma_b_sq / (sigma_b_sq + sigma_eps_sq) else NA_real_

  list(
    sigma_b_sq              = sigma_b_sq,
    sigma_eps_sq            = sigma_eps_sq,
    spatial_range           = spatial_range,
    prop_spatial            = prop_spatial,
    log_marginal_likelihood = log_ml,
    hyperparameters         = hyper
  )
}

# --- Bayes factor (Wagenmakers / Jeffreys scale) ----------------------------
#
# BF_10 = p(Y | M1) / p(Y | M0)
# log(BF_10) = log_ml_nonstationary - log_ml_stationary
# Positive => evidence for non-stationary covariance kernel

calculate_bayes_factor <- function(log_ml_stationary, log_ml_nonstationary) {
  log_bf <- log_ml_nonstationary - log_ml_stationary
  bf     <- exp(log_bf)

  abs_lbf   <- abs(log_bf)
  direction <- if (log_bf >= 0) "nonstationary" else "stationary"

  strength <- dplyr::case_when(
    abs_lbf < log(3)   ~ "Anecdotal",
    abs_lbf < log(10)  ~ "Moderate",
    abs_lbf < log(30)  ~ "Strong",
    abs_lbf < log(100) ~ "Very strong",
    TRUE               ~ "Extreme"
  )

  list(
    log_bayes_factor = log_bf,
    bayes_factor     = bf,
    interpretation   = paste(strength, "for", direction)
  )
}

# --- Single-gene analysis ---------------------------------------------------

analyze_gene_nonstationarity <- function(coords, expr, gene_name = NULL,
                                          normalize      = TRUE,
                                          detrend        = TRUE,
                                          detrend_method = "polynomial",
                                          mesh           = NULL,
                                          verbose        = FALSE) {
  expr_proc <- if (normalize) normalize_expression(expr) else expr
  if (detrend) expr_proc <- detrend_expression(coords, expr_proc, detrend_method)$detrended

  if (is.null(mesh)) mesh <- create_spde_mesh(coords)

  fit_stat <- fit_stationary_matern(coords, expr_proc, mesh = mesh, verbose = verbose)
  fit_ns   <- fit_nonstationary_matern(coords, expr_proc, mesh = mesh, verbose = verbose)

  res_stat <- extract_variance_components(fit_stat$result)
  res_ns   <- extract_variance_components(fit_ns$result)
  bf       <- calculate_bayes_factor(res_stat$log_marginal_likelihood,
                                     res_ns$log_marginal_likelihood)

  data.frame(
    gene_name                  = ifelse(is.null(gene_name), "unknown", gene_name),
    sigma_b_sq_stationary      = res_stat$sigma_b_sq,
    sigma_eps_sq_stationary    = res_stat$sigma_eps_sq,
    range_stationary           = res_stat$spatial_range,
    prop_spatial_stationary    = res_stat$prop_spatial,
    log_ml_stationary          = res_stat$log_marginal_likelihood,
    sigma_b_sq_nonstationary   = res_ns$sigma_b_sq,
    sigma_eps_sq_nonstationary = res_ns$sigma_eps_sq,
    range_nonstationary        = res_ns$spatial_range,
    prop_spatial_nonstationary = res_ns$prop_spatial,
    log_ml_nonstationary       = res_ns$log_marginal_likelihood,
    log_bayes_factor           = bf$log_bayes_factor,
    bayes_factor               = bf$bayes_factor,
    bf_interpretation          = bf$interpretation,
    stringsAsFactors           = FALSE
  )
}

# --- Safe wrapper for parallel execution ------------------------------------

analyze_gene_safe <- function(gene_idx, counts_mat, coords, gene_names,
                               normalize      = TRUE,
                               detrend        = TRUE,
                               detrend_method = "polynomial",
                               mesh           = NULL) {
  gene_name <- gene_names[gene_idx]
  tryCatch({
    expr <- log1p(counts_mat[gene_idx, ])
    analyze_gene_nonstationarity(
      coords         = coords,
      expr           = expr,
      gene_name      = gene_name,
      normalize      = normalize,
      detrend        = detrend,
      detrend_method = detrend_method,
      mesh           = mesh,
      verbose        = FALSE
    )
  }, error = function(e) {
    data.frame(
      gene_name                  = gene_name,
      sigma_b_sq_stationary      = NA_real_,
      sigma_eps_sq_stationary    = NA_real_,
      range_stationary           = NA_real_,
      prop_spatial_stationary    = NA_real_,
      log_ml_stationary          = NA_real_,
      sigma_b_sq_nonstationary   = NA_real_,
      sigma_eps_sq_nonstationary = NA_real_,
      range_nonstationary        = NA_real_,
      prop_spatial_nonstationary = NA_real_,
      log_ml_nonstationary       = NA_real_,
      log_bayes_factor           = NA_real_,
      bayes_factor               = NA_real_,
      bf_interpretation          = NA_character_,
      error_message              = e$message,
      stringsAsFactors           = FALSE
    )
  })
}

cat("Analysis functions loaded.\n")
