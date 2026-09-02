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

# --- Non-spatial Gaussian model ---------------------------------------------

fit_nonspatial_gaussian <- function(expr, verbose = FALSE) {
  INLA::inla(
    formula         = y ~ 1,
    data            = data.frame(y = as.numeric(expr)),
    family          = "gaussian",
    control.compute = list(config = TRUE),
    verbose         = verbose
  )
}

extract_nonspatial_components <- function(fit_result) {
  hyper <- fit_result$summary.hyperpar
  noise_idx <- grep(
    "Precision for the Gaussian",
    rownames(hyper),
    ignore.case = TRUE
  )
  intercept_idx <- match("(Intercept)", rownames(fit_result$summary.fixed))
  if (is.na(intercept_idx)) intercept_idx <- 1L

  list(
    log_marginal_likelihood = fit_result$mlik[1, 1],
    sigma_eps_sq = if (length(noise_idx)) 1 / hyper$mean[noise_idx[1]] else NA_real_,
    mu = fit_result$summary.fixed$mean[intercept_idx]
  )
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

# --- Extract nonstationary SPDE thetas ---------------------------------------
#
# For the nonstationary Matern (Ingebrigtsen et al. 2014):
#   log(tau(s))   = theta1
#   log(kappa(s)) = theta2 + theta3*x(s) + theta4*y(s)
#
# Returns theta1, theta2, theta3, theta4 from the INLA fit (for simulation).
# Used only when the fit has 4 spatial thetas (nonstationary model).

extract_nonstationary_thetas <- function(fit_result) {
  hyper <- fit_result$summary.hyperpar
  rn   <- rownames(hyper)

  theta1_idx <- grep("^Theta1 for spatial", rn, ignore.case = TRUE)
  theta2_idx <- grep("^Theta2 for spatial", rn, ignore.case = TRUE)
  theta3_idx <- grep("^Theta3 for spatial", rn, ignore.case = TRUE)
  theta4_idx <- grep("^Theta4 for spatial", rn, ignore.case = TRUE)

  theta1 <- if (length(theta1_idx) > 0) hyper$mean[theta1_idx[1]] else NA_real_
  theta2 <- if (length(theta2_idx) > 0) hyper$mean[theta2_idx[1]] else NA_real_
  theta3 <- if (length(theta3_idx) > 0) hyper$mean[theta3_idx[1]] else NA_real_
  theta4 <- if (length(theta4_idx) > 0) hyper$mean[theta4_idx[1]] else NA_real_

  list(theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4)
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
  fit_nonspatial <- fit_nonspatial_gaussian(expr_proc, verbose = verbose)

  res_stat <- extract_variance_components(fit_stat$result)
  res_ns   <- extract_variance_components(fit_ns$result)
  res_nonspatial <- extract_nonspatial_components(fit_nonspatial)
  bf       <- calculate_bayes_factor(res_stat$log_marginal_likelihood,
                                     res_ns$log_marginal_likelihood)

  thetas_ns <- extract_nonstationary_thetas(fit_ns$result)

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
    theta1_ns                  = thetas_ns$theta1,
    theta2_ns                  = thetas_ns$theta2,
    theta3_ns                  = thetas_ns$theta3,
    theta4_ns                  = thetas_ns$theta4,
    log_ml_nonspatial          = res_nonspatial$log_marginal_likelihood,
    sigma_eps_sq_nonspatial    = res_nonspatial$sigma_eps_sq,
    mu_nonspatial              = res_nonspatial$mu,
    log_bayes_factor           = bf$log_bayes_factor,
    bayes_factor               = bf$bayes_factor,
    bf_interpretation          = bf$interpretation,
    error_message              = NA_character_,
    stringsAsFactors           = FALSE
  )
}

# --- Safe wrapper for parallel execution ------------------------------------

empty_gene_analysis_result <- function(gene_name, error_message) {
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
    theta1_ns                  = NA_real_,
    theta2_ns                  = NA_real_,
    theta3_ns                  = NA_real_,
    theta4_ns                  = NA_real_,
    log_ml_nonspatial          = NA_real_,
    sigma_eps_sq_nonspatial    = NA_real_,
    mu_nonspatial              = NA_real_,
    log_bayes_factor           = NA_real_,
    bayes_factor               = NA_real_,
    bf_interpretation          = NA_character_,
    error_message              = as.character(error_message),
    stringsAsFactors           = FALSE
  )
}

analyze_gene_safe <- function(gene_idx, counts_mat, coords, gene_names,
                               gene_ids       = NULL,
                               normalize      = TRUE,
                               detrend        = TRUE,
                               detrend_method = "polynomial",
                               mesh           = NULL) {
  gene_name <- gene_names[gene_idx]
  gene_id <- if (is.null(gene_ids)) as.character(gene_idx) else gene_ids[gene_idx]

  result <- tryCatch({
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
    ) |>
      validate_gene_analysis_payload()
  }, error = function(e) {
    empty_gene_analysis_result(gene_name, conditionMessage(e))
  })

  data.frame(
    gene_index = as.integer(gene_idx),
    gene_id    = as.character(gene_id),
    result,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# --- Result integrity and atomic persistence ---------------------------------

analysis_schema_version <- 3L
analysis_implementation_files <- c("analysis_functions.R", "run_analysis.R")

analysis_result_numeric_columns <- c(
  "sigma_b_sq_stationary", "sigma_eps_sq_stationary", "range_stationary",
  "prop_spatial_stationary", "log_ml_stationary",
  "sigma_b_sq_nonstationary", "sigma_eps_sq_nonstationary",
  "range_nonstationary", "prop_spatial_nonstationary",
  "log_ml_nonstationary", "theta1_ns", "theta2_ns", "theta3_ns",
  "theta4_ns", "log_ml_nonspatial", "sigma_eps_sq_nonspatial",
  "mu_nonspatial", "log_bayes_factor", "bayes_factor"
)

analysis_result_required_columns <- c(
  "gene_index", "gene_id", "gene_name",
  analysis_result_numeric_columns,
  "bf_interpretation", "error_message"
)

fingerprint_analysis_implementation <- function(script_dir) {
  paths <- file.path(script_dir, analysis_implementation_files)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Cannot fingerprint missing analysis files: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  hashes <- unname(as.character(tools::md5sum(paths)))
  if (anyNA(hashes) || any(!grepl("^[0-9a-fA-F]{32}$", hashes))) {
    stop("Could not fingerprint the analysis implementation", call. = FALSE)
  }
  stats::setNames(hashes, analysis_implementation_files)
}

analysis_implementation_matches <- function(observed, expected) {
  is.character(observed) &&
    identical(names(observed), names(expected)) &&
    identical(tolower(unname(observed)), tolower(unname(expected)))
}

validate_gene_analysis_payload <- function(result) {
  required <- setdiff(
    analysis_result_required_columns,
    c("gene_index", "gene_id", "error_message")
  )
  if (!is.data.frame(result) || nrow(result) != 1L) {
    stop("model fit did not return exactly one result row", call. = FALSE)
  }

  missing <- setdiff(required, names(result))
  if (length(missing) > 0L) {
    stop(
      "model fit omitted outputs: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric <- analysis_result_numeric_columns[
    !vapply(result[analysis_result_numeric_columns], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "model fit returned non-numeric outputs: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  finite_columns <- setdiff(analysis_result_numeric_columns, "bayes_factor")
  nonfinite <- finite_columns[!vapply(
    result[finite_columns],
    function(x) length(x) == 1L && is.finite(x),
    logical(1)
  )]
  if (length(nonfinite) > 0L) {
    stop(
      "model fit returned non-finite outputs: ",
      paste(nonfinite, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(result$bayes_factor) || result$bayes_factor < 0) {
    stop("model fit returned an invalid bayes_factor", call. = FALSE)
  }
  if (!is.character(result$gene_name) || is.na(result$gene_name) ||
      !nzchar(result$gene_name)) {
    stop("model fit returned an invalid gene_name", call. = FALSE)
  }
  if (!is.character(result$bf_interpretation) ||
      is.na(result$bf_interpretation) ||
      !nzchar(result$bf_interpretation)) {
    stop("model fit returned an invalid bf_interpretation", call. = FALSE)
  }

  expected_log_bf <- result$log_ml_nonstationary - result$log_ml_stationary
  if (abs(result$log_bayes_factor - expected_log_bf) >
      1e-8 * (1 + abs(expected_log_bf))) {
    stop(
      "model fit returned an inconsistent log_bayes_factor",
      call. = FALSE
    )
  }
  if (result$prop_spatial_stationary < 0 ||
      result$prop_spatial_stationary > 1 ||
      result$prop_spatial_nonstationary < 0 ||
      result$prop_spatial_nonstationary > 1) {
    stop("model fit returned a spatial proportion outside [0, 1]", call. = FALSE)
  }

  result
}

validate_gene_results <- function(results, expected_gene_indices,
                                  expected_gene_ids = NULL) {
  problems <- character()

  if (!is.data.frame(results)) {
    return(list(valid = FALSE, message = "result object is not a data.frame"))
  }

  missing_columns <- setdiff(analysis_result_required_columns, names(results))
  if (length(missing_columns) > 0L) {
    problems <- c(
      problems,
      paste0("missing required columns: ", paste(missing_columns, collapse = ", "))
    )
  }

  numeric_columns <- intersect(analysis_result_numeric_columns, names(results))
  non_numeric_columns <- numeric_columns[
    !vapply(results[numeric_columns], is.numeric, logical(1))
  ]
  if (length(non_numeric_columns) > 0L) {
    problems <- c(
      problems,
      paste0(
        "model columns are not numeric: ",
        paste(non_numeric_columns, collapse = ", ")
      )
    )
  }

  character_columns <- intersect(
    c("gene_id", "gene_name", "bf_interpretation", "error_message"),
    names(results)
  )
  non_character_columns <- character_columns[
    !vapply(results[character_columns], is.character, logical(1))
  ]
  if (length(non_character_columns) > 0L) {
    problems <- c(
      problems,
      paste0(
        "identifier/status columns are not character: ",
        paste(non_character_columns, collapse = ", ")
      )
    )
  }

  expected_gene_indices <- as.integer(expected_gene_indices)
  if (nrow(results) != length(expected_gene_indices)) {
    problems <- c(
      problems,
      sprintf(
        "row count is %d but %d scheduled genes were expected",
        nrow(results), length(expected_gene_indices)
      )
    )
  }

  if ("gene_index" %in% names(results)) {
    observed_indices <- as.integer(results$gene_index)
    if (anyNA(observed_indices)) {
      problems <- c(problems, "gene_index contains missing values")
    }
    if (anyDuplicated(observed_indices)) {
      problems <- c(problems, "gene_index contains duplicates")
    }
    if (!identical(observed_indices, expected_gene_indices)) {
      problems <- c(problems, "gene_index does not exactly match the scheduled order")
    }
  }

  if (!is.null(expected_gene_ids) && "gene_id" %in% names(results)) {
    if (!identical(as.character(results$gene_id), as.character(expected_gene_ids))) {
      problems <- c(problems, "gene_id does not exactly match the filtered dataset")
    }
  }

  for (field in intersect(c("gene_id", "gene_name"), names(results))) {
    values <- as.character(results[[field]])
    if (anyNA(values) || any(!nzchar(values))) {
      problems <- c(problems, paste0(field, " contains missing or empty values"))
    }
  }

  if (length(missing_columns) == 0L && length(non_numeric_columns) == 0L) {
    failed_fit <- is.na(results$log_bayes_factor)
    has_error <- !is.na(results$error_message) & nzchar(results$error_message)
    if (any(failed_fit & !has_error)) {
      problems <- c(problems, "failed fits are missing an error_message")
    }
    if (any(!failed_fit & has_error)) {
      problems <- c(problems, "successful fits contain an error_message")
    }

    successful_fit <- !failed_fit
    finite_success_columns <- setdiff(
      analysis_result_numeric_columns,
      "bayes_factor"
    )
    success_values <- as.matrix(results[successful_fit, finite_success_columns, drop = FALSE])
    if (length(success_values) > 0L && any(!is.finite(success_values))) {
      problems <- c(problems, "successful fits contain non-finite model outputs")
    }
    if (any(
      successful_fit &
        (is.na(results$bayes_factor) | results$bayes_factor < 0)
    )) {
      problems <- c(problems, "successful fits contain an invalid bayes_factor")
    }
    if (any(
      successful_fit &
        (is.na(results$bf_interpretation) | !nzchar(results$bf_interpretation))
    )) {
      problems <- c(problems, "successful fits are missing bf_interpretation")
    }

    failed_values <- as.matrix(results[failed_fit, analysis_result_numeric_columns, drop = FALSE])
    if (length(failed_values) > 0L && any(!is.na(failed_values))) {
      problems <- c(problems, "failed fits contain partial numeric model outputs")
    }
    if (any(failed_fit & !is.na(results$bf_interpretation))) {
      problems <- c(problems, "failed fits contain bf_interpretation")
    }

    expected_log_bf <- results$log_ml_nonstationary - results$log_ml_stationary
    inconsistent_log_bf <- successful_fit & (
      abs(results$log_bayes_factor - expected_log_bf) >
        1e-8 * (1 + abs(expected_log_bf))
    )
    if (any(inconsistent_log_bf)) {
      problems <- c(
        problems,
        "log_bayes_factor is inconsistent with the two spatial marginal likelihoods"
      )
    }

    if (any(
      successful_fit &
        (results$prop_spatial_stationary < 0 |
           results$prop_spatial_stationary > 1 |
           results$prop_spatial_nonstationary < 0 |
           results$prop_spatial_nonstationary > 1)
    )) {
      problems <- c(problems, "successful fits contain spatial proportions outside [0, 1]")
    }
  }

  list(
    valid   = length(problems) == 0L,
    message = if (length(problems) == 0L) "complete" else paste(unique(problems), collapse = "; ")
  )
}

assert_complete_gene_results <- function(results, expected_gene_indices,
                                         expected_gene_ids = NULL) {
  check <- validate_gene_results(results, expected_gene_indices, expected_gene_ids)
  if (!check$valid) {
    stop("Incomplete gene results: ", check$message, call. = FALSE)
  }
  invisible(results)
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir  = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)

  saveRDS(object, tmp_path)
  if (!file.rename(tmp_path, path)) {
    stop("Could not atomically replace result file: ", path, call. = FALSE)
  }

  invisible(path)
}

write_csv_atomic <- function(object, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir  = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)

  write.csv(object, tmp_path, row.names = row.names)
  if (!file.rename(tmp_path, path)) {
    stop("Could not atomically replace CSV file: ", path, call. = FALSE)
  }

  invisible(path)
}

cat("Analysis functions loaded.\n")
