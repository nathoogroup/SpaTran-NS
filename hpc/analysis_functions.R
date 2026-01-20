# =============================================================================
# Spatial Nonstationarity Analysis Functions
# =============================================================================
# Helper functions for INLA-based spatial nonstationarity analysis
# Source this file in your analysis scripts
# =============================================================================

# Required packages
required_packages <- c("INLA", "Matrix", "dplyr", "SpatialExperiment")

# --- Normalization ---
normalize_expression <- function(expr, method = "quantile") {
  if (method == "quantile") {
    expr_rank <- rank(expr, ties.method = "average")
    expr_norm <- qnorm((expr_rank - 0.5) / length(expr))
    return(expr_norm)
  } else if (method == "zscore") {
    return((expr - mean(expr)) / sd(expr))
  } else {
    return(expr)
  }
}

# --- Detrending ---
detrend_expression <- function(coords, expr, method = "polynomial") {
  if (method == "loess") {
    loess_fit <- loess(expr ~ x + y, data = coords, span = 0.5)
    trend <- predict(loess_fit, newdata = coords)
    detrended <- expr - trend
    return(list(detrended = detrended, trend = trend))
  } else if (method == "polynomial") {
    x_centered <- coords$x - mean(coords$x)
    y_centered <- coords$y - mean(coords$y)
    poly_fit <- lm(expr ~ x_centered + y_centered + 
                   I(x_centered^2) + I(y_centered^2) + I(x_centered * y_centered),
                   data = data.frame(expr = expr, x_centered = x_centered, y_centered = y_centered))
    trend <- predict(poly_fit)
    detrended <- expr - trend
    return(list(detrended = detrended, trend = trend))
  } else {
    return(list(detrended = expr, trend = rep(0, length(expr))))
  }
}

# --- INLA Mesh Creation ---
create_spde_mesh <- function(coords, max.edge = NULL, cutoff = NULL) {
  if (is.null(max.edge)) {
    x_range <- diff(range(coords$x))
    y_range <- diff(range(coords$y))
    max_range <- max(x_range, y_range)
    max.edge <- c(max_range / 5, max_range / 3)
  }
  
  if (is.null(cutoff)) {
    x_range <- diff(range(coords$x))
    y_range <- diff(range(coords$y))
    min_range <- min(x_range, y_range)
    cutoff <- min_range / 20
  }
  
  mesh <- INLA::inla.mesh.2d(
    loc = as.matrix(coords[, c("x", "y")]),
    max.edge = max.edge,
    cutoff = cutoff
  )
  
  return(mesh)
}

# --- Stationary Matern Model ---
fit_stationary_matern <- function(coords, expr, mesh = NULL, verbose = FALSE) {
  if (is.null(mesh)) {
    mesh <- create_spde_mesh(coords)
  }
  
  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,
    prior.range = c(0.5, 0.5),
    prior.sigma = c(1, 0.1)
  )
  
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))
  s.index <- INLA::inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)
  
  stack <- INLA::inla.stack(
    data = list(y = expr),
    A = list(A, 1),
    effects = list(
      s.index,
      list(intercept = rep(1, length(expr)))
    ),
    tag = "est"
  )
  
  formula <- y ~ -1 + intercept + f(spatial, model = spde)
  
  result <- INLA::inla(
    formula = formula,
    data = INLA::inla.stack.data(stack, spde = spde),
    family = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(config = TRUE, dic = TRUE, cpo = TRUE),
    verbose = verbose
  )
  
  return(list(result = result, mesh = mesh, spde = spde, stack = stack))
}

# --- Nonstationary Matern Model (with polynomial covariates) ---
fit_nonstationary_matern_simple <- function(coords, expr, mesh = NULL, verbose = FALSE) {
  if (is.null(mesh)) {
    mesh <- create_spde_mesh(coords)
  }
  
  coords_scaled <- coords
  coords_scaled$x_scaled <- scale(coords$x)[, 1]
  coords_scaled$y_scaled <- scale(coords$y)[, 1]
  
  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,
    prior.range = c(0.5, 0.5),
    prior.sigma = c(1, 0.1)
  )
  
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))
  s.index <- INLA::inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)
  
  stack <- INLA::inla.stack(
    data = list(y = expr),
    A = list(A, 1),
    effects = list(
      s.index,
      list(
        intercept = rep(1, length(expr)),
        x_scaled = coords_scaled$x_scaled,
        y_scaled = coords_scaled$y_scaled,
        x2 = coords_scaled$x_scaled^2,
        y2 = coords_scaled$y_scaled^2,
        xy = coords_scaled$x_scaled * coords_scaled$y_scaled
      )
    ),
    tag = "est"
  )
  
  formula <- y ~ -1 + intercept + 
    x_scaled + y_scaled + x2 + y2 + xy +
    f(spatial, model = spde)
  
  result <- INLA::inla(
    formula = formula,
    data = INLA::inla.stack.data(stack, spde = spde),
    family = "gaussian",
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.compute = list(config = TRUE, dic = TRUE, cpo = TRUE),
    verbose = verbose
  )
  
  return(list(result = result, mesh = mesh, spde = spde, stack = stack))
}

# --- Extract Variance Components ---
extract_variance_components <- function(fit_result) {
  hyper <- fit_result$summary.hyperpar
  
  sigma_b_sq <- NA
  sigma_eps_sq <- NA
  spatial_range <- NA
  
  # Extract spatial Stdev
  stdev_idx <- grep("Stdev for spatial", rownames(hyper), ignore.case = TRUE)
  if (length(stdev_idx) > 0) {
    spatial_stdev <- hyper$mean[stdev_idx[1]]
    sigma_b_sq <- spatial_stdev^2
  }
  
  # Extract range
  range_idx <- grep("Range for spatial", rownames(hyper), ignore.case = TRUE)
  if (length(range_idx) > 0) {
    spatial_range <- hyper$mean[range_idx[1]]
  }
  
  # Extract noise variance
  noise_prec_idx <- grep("Precision for the Gaussian observations", rownames(hyper))
  if (length(noise_prec_idx) > 0) {
    noise_prec <- hyper$mean[noise_prec_idx]
    sigma_eps_sq <- 1 / noise_prec
  }
  
  # Calculate proportion
  prop_spatial <- NA
  if (!is.na(sigma_b_sq) && !is.na(sigma_eps_sq)) {
    prop_spatial <- sigma_b_sq / (sigma_b_sq + sigma_eps_sq)
  }
  
  log_ml <- fit_result$mlik[1, 1]
  
  return(list(
    sigma_b_sq = sigma_b_sq,
    sigma_eps_sq = sigma_eps_sq,
    spatial_range = spatial_range,
    prop_spatial = prop_spatial,
    log_marginal_likelihood = log_ml,
    hyperparameters = hyper
  ))
}

# --- Calculate Bayes Factor (Wagenmakers/Jeffreys scale) ---
calculate_bayes_factor <- function(log_ml_stationary, log_ml_nonstationary) {
  log_bf <- log_ml_nonstationary - log_ml_stationary
  bf <- exp(log_bf)
  
  abs_log_bf <- abs(log_bf)
  direction <- ifelse(log_bf >= 0, "nonstationary", "stationary")
  
  # Wagenmakers/Jeffreys thresholds
  strength <- dplyr::case_when(
    abs_log_bf < log(3)   ~ "Anecdotal",
    abs_log_bf < log(10)  ~ "Moderate",
    abs_log_bf < log(30)  ~ "Strong",
    abs_log_bf < log(100) ~ "Very strong",
    TRUE                  ~ "Extreme"
  )
  
  interpretation <- paste(strength, "for", direction)
  
  return(list(
    log_bayes_factor = log_bf,
    bayes_factor = bf,
    interpretation = interpretation
  ))
}

# --- Analyze Single Gene ---
analyze_gene_nonstationarity <- function(coords, expr, gene_name = NULL, 
                                         normalize = TRUE, detrend = TRUE,
                                         detrend_method = "polynomial",
                                         verbose = FALSE) {
  
  # Normalize
  if (normalize) {
    expr_processed <- normalize_expression(expr, method = "quantile")
  } else {
    expr_processed <- expr
  }
  
  # Detrend
  if (detrend) {
    detrend_result <- detrend_expression(coords, expr_processed, method = detrend_method)
    expr_processed <- detrend_result$detrended
  }
  
  # Create mesh
  mesh <- create_spde_mesh(coords)
  
  # Fit models
  fit_stationary <- fit_stationary_matern(coords, expr_processed, mesh = mesh, verbose = verbose)
  fit_nonstationary <- fit_nonstationary_matern_simple(coords, expr_processed, mesh = mesh, verbose = verbose)
  
  # Extract results
  results_stationary <- extract_variance_components(fit_stationary$result)
  results_nonstationary <- extract_variance_components(fit_nonstationary$result)
  
  # Calculate Bayes factor
  bf_result <- calculate_bayes_factor(
    results_stationary$log_marginal_likelihood,
    results_nonstationary$log_marginal_likelihood
  )
  
  # Compile results
  results <- data.frame(
    gene_name = ifelse(is.null(gene_name), "unknown", gene_name),
    sigma_b_sq_stationary = results_stationary$sigma_b_sq,
    sigma_eps_sq_stationary = results_stationary$sigma_eps_sq,
    range_stationary = results_stationary$spatial_range,
    prop_spatial_stationary = results_stationary$prop_spatial,
    log_ml_stationary = results_stationary$log_marginal_likelihood,
    sigma_b_sq_nonstationary = results_nonstationary$sigma_b_sq,
    sigma_eps_sq_nonstationary = results_nonstationary$sigma_eps_sq,
    range_nonstationary = results_nonstationary$spatial_range,
    prop_spatial_nonstationary = results_nonstationary$prop_spatial,
    log_ml_nonstationary = results_nonstationary$log_marginal_likelihood,
    log_bayes_factor = bf_result$log_bayes_factor,
    bayes_factor = bf_result$bayes_factor,
    bf_interpretation = bf_result$interpretation,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# --- Analyze Single Gene (wrapper for parallel execution) ---
analyze_gene_safe <- function(gene_idx, spe, coords, gene_names, 
                               normalize = TRUE, detrend = TRUE,
                               detrend_method = "polynomial") {
  
  gene_name <- gene_names[gene_idx]
  
  tryCatch({
    # Get expression
    expr <- log1p(SummarizedExperiment::assay(spe, "counts")[gene_idx, ])
    
    # Analyze
    result <- analyze_gene_nonstationarity(
      coords = coords,
      expr = expr,
      gene_name = gene_name,
      normalize = normalize,
      detrend = detrend,
      detrend_method = detrend_method,
      verbose = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    # Return NA row on error
    return(data.frame(
      gene_name = gene_name,
      sigma_b_sq_stationary = NA_real_,
      sigma_eps_sq_stationary = NA_real_,
      range_stationary = NA_real_,
      prop_spatial_stationary = NA_real_,
      log_ml_stationary = NA_real_,
      sigma_b_sq_nonstationary = NA_real_,
      sigma_eps_sq_nonstationary = NA_real_,
      range_nonstationary = NA_real_,
      prop_spatial_nonstationary = NA_real_,
      log_ml_nonstationary = NA_real_,
      log_bayes_factor = NA_real_,
      bayes_factor = NA_real_,
      bf_interpretation = NA_character_,
      error_message = e$message,
      stringsAsFactors = FALSE
    ))
  })
}

cat("Analysis functions loaded successfully.\n")
