# =============================================================================
# run_local_simulation_test: Local Simulation Test Function
# =============================================================================
# Callable function for running the local simulation test. Source this file
# and call with desired parameters.
#
# Usage from R:
#   source("simulation/run_local_test.R")
#   run_local_simulation_test(n_genes = 10, n_reps = 2)
#   run_local_simulation_test(10, 2, "HumanOvarianCancer")
#
# Or run the script (which calls the function):
#   Rscript simulation/test_simulation_local.R 10 2 HumanBreastCancerILC
# =============================================================================

run_local_simulation_test <- function(
  n_genes = 10,
  n_reps = 2,
  dataset = "HumanBreastCancerILC",
  project_root = getwd(),
  output_dir = NULL,
  install_pkgs = TRUE,
  verbose = TRUE
) {
  msg <- function(...) if (verbose) cat(..., "\n", sep = "")

  msg("=============================================================")
  msg("LOCAL SIMULATION TEST")
  msg("=============================================================")
  msg("Dataset:    ", dataset)
  msg("Genes:      ", n_genes)
  msg("Replicates: ", n_reps)
  msg("=============================================================")
  msg("")

  out_dir <- output_dir %||% file.path(project_root, "simulation/output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Install and load packages --------------------------------------------
  cran_mirror <- "https://cloud.r-project.org"
  if (install_pkgs) {
    for (pkg in c("Matrix", "dplyr", "tidyr", "ggplot2", "patchwork", "scales")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        msg("Installing ", pkg, "...")
        install.packages(pkg, repos = cran_mirror)
      }
    }
    if (!requireNamespace("INLA", quietly = TRUE)) {
      msg("Installing INLA...")
      install.packages("INLA",
        repos = c(CRAN = cran_mirror, INLA = "https://inla.r-inla-download.org/R/stable"),
        dep = TRUE
      )
    }
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
      msg("Installing SpatialExperiment...")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = cran_mirror)
      BiocManager::install("SpatialExperiment", ask = FALSE, update = FALSE)
    }
  }

  suppressPackageStartupMessages({
    library(INLA)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
  })
  INLA::inla.setOption(num.threads = "1:1")
  msg("All packages loaded.")
  msg("")

  # --- Source simulation functions ------------------------------------------
  sim_dir <- file.path(project_root, "simulation")
  source(file.path(sim_dir, "simulation_functions.R"))

  # --- Load data ------------------------------------------------------------
  results_file <- file.path(project_root, "results", paste0(dataset, "_results.rds"))
  spe_file <- file.path(project_root, "data/spatial_datasets", paste0(dataset, ".rds"))
  if (!file.exists(results_file)) stop("Results not found: ", results_file)
  if (!file.exists(spe_file)) stop("Dataset not found: ", spe_file)

  results <- readRDS(results_file)
  n_spots <- unique(results$n_spots)
  msg("Loaded results: ", nrow(results), " genes, ", n_spots, " spots")

  suppressPackageStartupMessages(library(SpatialExperiment))
  spe <- readRDS(spe_file)
  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords)[1:2] <- c("x", "y")
  rm(spe); gc(verbose = FALSE)
  msg("  Loaded ", nrow(coords), " spots")

  # --- Build mesh and SPDE --------------------------------------------------
  msg("Building mesh...")
  xr <- diff(range(coords$x)); yr <- diff(range(coords$y))
  mesh <- INLA::inla.mesh.2d(
    loc = as.matrix(coords[, c("x", "y")]),
    max.edge = c(max(xr, yr) / 5, max(xr, yr) / 3),
    cutoff = min(xr, yr) / 20
  )
  msg("  Mesh nodes: ", mesh$n)

  spde_stat <- build_spde_stationary(mesh)
  spde_ns <- build_spde_nonstationary(mesh)
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(coords[, c("x", "y")]))
  s.index <- INLA::inla.spde.make.index(name = "spatial", n.spde = mesh$n)

  # --- Select genes ---------------------------------------------------------
  results_valid <- results %>%
    filter(!is.na(log_bayes_factor), sigma_b_sq_stationary > 0,
           sigma_eps_sq_stationary > 0, range_stationary > 0,
           is.finite(sigma_b_sq_stationary), is.finite(sigma_eps_sq_stationary),
           is.finite(range_stationary))

  n_h0 <- sum(results_valid$log_bayes_factor <= 0)
  n_h1 <- sum(results_valid$log_bayes_factor > 0)
  h0_frac <- n_h0 / nrow(results_valid)
  n_h0_sub <- max(2, round(n_genes * h0_frac))
  n_h1_sub <- max(2, n_genes - n_h0_sub)

  set.seed(777)
  h0_rows <- results_valid %>% filter(log_bayes_factor <= 0) %>% slice_sample(n = n_h0_sub)
  h1_rows <- results_valid %>% filter(log_bayes_factor > 0) %>% slice_sample(n = n_h1_sub)
  sub_results <- bind_rows(h0_rows, h1_rows)
  msg("Selected ", nrow(sub_results), " genes: H0=", nrow(h0_rows), ", H1=", nrow(h1_rows))

  # --- Derive sim params ----------------------------------------------------
  sub_results$hypothesis <- ifelse(sub_results$log_bayes_factor > 0, "H1", "H0")
  sim_params <- sub_results %>%
    mutate(
      log_kappa = log(sqrt(8) / range_stationary),
      log_tau = log(1 / sqrt(4 * pi * (sqrt(8) / range_stationary)^2 * sigma_b_sq_stationary)),
      sigma_eps_sq = sigma_eps_sq_stationary, theta3 = 0, theta4 = 0
    )

  has_fitted <- all(c("theta1_ns", "theta2_ns", "theta3_ns", "theta4_ns") %in% colnames(sub_results))
  h1_idx <- which(sim_params$hypothesis == "H1")
  h1_use_fitted <- rep(FALSE, nrow(sim_params))

  if (has_fitted && length(h1_idx) > 0) {
    for (i in h1_idx) {
      t1 <- sub_results$theta1_ns[match(sim_params$gene_name[i], sub_results$gene_name)]
      t2 <- sub_results$theta2_ns[match(sim_params$gene_name[i], sub_results$gene_name)]
      t3 <- sub_results$theta3_ns[match(sim_params$gene_name[i], sub_results$gene_name)]
      t4 <- sub_results$theta4_ns[match(sim_params$gene_name[i], sub_results$gene_name)]
      if (all(is.finite(c(t1, t2, t3, t4)))) {
        sim_params$log_tau[i] <- t1; sim_params$log_kappa[i] <- t2
        sim_params$theta3[i] <- t3; sim_params$theta4[i] <- t4
        h1_use_fitted[i] <- TRUE
      }
    }
    if (sum(h1_use_fitted) > 0) msg("  Using fitted thetas for ", sum(h1_use_fitted), " H1 genes")
  }
  if (any(!h1_use_fitted[h1_idx])) {
    ns_params <- assign_nonstationarity_params(sub_results, seed = 42)
    for (i in h1_idx) {
      if (!h1_use_fitted[i]) {
        j <- match(sim_params$gene_name[i], ns_params$gene_name)
        if (!is.na(j)) {
          sim_params$theta3[i] <- ns_params$theta3[j]
          sim_params$theta4[i] <- ns_params$theta4[j]
        }
      }
    }
  }

  param_table <- sim_params %>%
    mutate(magnitude = sqrt(theta3^2 + theta4^2),
           label = ifelse(hypothesis == "H0", sprintf("range=%.0f", range_stationary),
                         sprintf("range=%.0f, mag=%.2f", range_stationary, magnitude)))

  n_genes_act <- nrow(sim_params)

  # --- Spatial maps (replicate 1) --------------------------------------------
  msg("Simulating replicate 1 for spatial maps...")
  set.seed(42001)
  sim_expr_list <- vector("list", n_genes_act)
  for (i in seq_len(n_genes_act)) {
    g <- sim_params[i, ]
    th <- if (g$hypothesis == "H0") c(g$log_tau, g$log_kappa)
         else c(g$log_tau, g$log_kappa, g$theta3, g$theta4)
    spde <- if (g$hypothesis == "H0") spde_stat else spde_ns
    sim_expr_list[[i]] <- simulate_expression(spde, th, A, g$sigma_eps_sq)
  }

  spatial_df <- do.call(rbind, lapply(seq_len(n_genes_act), function(i) {
    g <- sim_params[i, ]; mag <- sqrt(g$theta3^2 + g$theta4^2)
    sub <- if (g$hypothesis == "H0") sprintf("range=%.0f", g$range_stationary)
           else sprintf("range=%.0f, mag=%.2f", g$range_stationary, mag)
    data.frame(x = coords$x, y = coords$y, expression = sim_expr_list[[i]],
               gene_name = g$gene_name, hypothesis = g$hypothesis,
               panel_label = sprintf("%s [%s]\n%s", g$gene_name, g$hypothesis, sub),
               stringsAsFactors = FALSE)
  }))
  spatial_df$panel_label <- factor(spatial_df$panel_label, levels = unique(spatial_df$panel_label))

  nc <- ceiling(sqrt(n_genes_act)); nr <- ceiling(n_genes_act / nc)
  p_spatial <- ggplot(spatial_df, aes(x = x, y = y, color = expression)) +
    geom_point(size = 0.3, alpha = 0.8) + scale_color_viridis_c(option = "magma", name = "Expression") +
    facet_wrap(~ panel_label, ncol = nc) + coord_fixed() +
    labs(title = paste0("Simulated Spatial Expression  |  ", dataset),
         subtitle = paste0("Replicate 1  |  ", nrow(coords), " spots  |  H0: ", sum(sim_params$hypothesis == "H0"),
                          " genes  |  H1: ", sum(sim_params$hypothesis == "H1"), " genes")) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, color = "grey30"),
          strip.text = element_text(size = 7, face = "bold"), axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right")

  spatial_file <- file.path(out_dir, "test_spatial_expression.png")
  ggsave(spatial_file, p_spatial, width = 3 * nc + 1.5, height = 3 * nr + 1, dpi = 200, bg = "white")
  msg("  Saved: ", spatial_file)

  # --- Run replicates -------------------------------------------------------
  msg("")
  msg("Running ", n_reps, " replicates...")
  msg("")
  replicate_results <- list()
  all_gene_bfs <- list()
  t0 <- proc.time()
  seed_base <- switch(dataset, HumanBreastCancerILC = 42000, HumanOvarianCancer = 43000, 99000)

  for (r in seq_len(n_reps)) {
    t0r <- proc.time()
    msg("--- Replicate ", r, "/", n_reps, " ---")
    set.seed(seed_base + r)
    rep_results <- data.frame(gene_name = character(n_genes_act), hypothesis = character(n_genes_act),
                              log_bf_sim = numeric(n_genes_act), bf_sim = numeric(n_genes_act),
                              stringsAsFactors = FALSE)

    for (i in seq_len(n_genes_act)) {
      t0_gene <- proc.time()
      g <- sim_params[i, ]
      th <- if (g$hypothesis == "H0") c(g$log_tau, g$log_kappa)
           else c(g$log_tau, g$log_kappa, g$theta3, g$theta4)
      spde <- if (g$hypothesis == "H0") spde_stat else spde_ns
      y <- simulate_expression(spde, th, A, g$sigma_eps_sq)
      log_bf <- tryCatch(fit_and_get_logbf(y, mesh, spde_stat, spde_ns, A, s.index),
                         error = function(e) { if (verbose) cat("    Error:", e$message, "\n"); NA_real_ })
      rep_results$gene_name[i] <- g$gene_name
      rep_results$hypothesis[i] <- g$hypothesis
      rep_results$log_bf_sim[i] <- log_bf
      rep_results$bf_sim[i] <- exp(log_bf)
      elapsed_gene <- (proc.time() - t0_gene)["elapsed"]
      if (verbose) cat(sprintf("  Gene %2d/%d [%s] %-15s %s  [%.1f sec]\n", i, n_genes_act, g$hypothesis, g$gene_name,
                               if (is.na(log_bf)) "FAIL" else sprintf("logBF=%.2f", log_bf), elapsed_gene))
    }

    replicate_results[[r]] <- rep_results
    rep_results$replicate <- r
    all_gene_bfs[[r]] <- rep_results
    fp <- sum(rep_results$hypothesis == "H0" & rep_results$bf_sim > 1, na.rm = TRUE)
    tp <- sum(rep_results$hypothesis == "H1" & rep_results$bf_sim > 1, na.rm = TRUE)
    disc <- fp + tp
    msg(sprintf("  -> FP=%d, TP=%d, FDR=%.3f, Power=%.3f  [%.0f sec]", fp, tp,
                if (disc > 0) fp / disc else 0, tp / sum(rep_results$hypothesis == "H1"),
                (proc.time() - t0r)["elapsed"]))
    msg("")
  }

  elapsed <- (proc.time() - t0)["elapsed"]
  msg("All replicates completed in ", round(elapsed / 60, 1), " min")
  msg("")

  # --- FDR/power curves -----------------------------------------------------
  curves <- compute_fdr_power(replicate_results)
  fdr_1 <- sapply(replicate_results, function(rd) {
    rd <- rd %>% filter(!is.na(bf_sim))
    fp <- sum(rd$hypothesis == "H0" & rd$bf_sim > 1)
    tp <- sum(rd$hypothesis == "H1" & rd$bf_sim > 1)
    if (fp + tp > 0) fp / (fp + tp) else 0
  })
  power_1 <- sapply(replicate_results, function(rd) {
    rd <- rd %>% filter(!is.na(bf_sim))
    tp <- sum(rd$hypothesis == "H1" & rd$bf_sim > 1)
    tp / sum(rd$hypothesis == "H1")
  })

  msg("RESULTS AT KEY THRESHOLDS")
  msg("  c=1:   FDR=", round(mean(fdr_1), 4), " +/- ", round(sd(fdr_1), 4),
      "  Power=", round(mean(power_1), 4), " +/- ", round(sd(power_1), 4))
  for (ct in c(3, 10, 30, 100)) {
    row <- curves[which.min(abs(curves$c_threshold - ct)), ]
    msg("  c=", ct, ": FDR=", round(row$fdr, 4), "  Power=", round(row$avg_power, 4))
  }
  msg("")

  # --- FDR/Power figure -----------------------------------------------------
  df_plot <- curves %>% filter(c_threshold > 0) %>% mutate(log10_c = log10(c_threshold))
  c1_fdr <- curves$fdr[which.min(abs(curves$c_threshold - 1))]
  c1_pow <- curves$avg_power[which.min(abs(curves$c_threshold - 1))]

  p_fdr <- ggplot(df_plot, aes(x = log10_c, y = fdr)) +
    geom_ribbon(aes(ymin = pmax(fdr - fdr_sd, 0), ymax = pmin(fdr + fdr_sd, 1)), fill = "#D55E00", alpha = 0.15) +
    geom_line(color = "#D55E00", linewidth = 1) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(data = df_plot %>% filter(abs(c_threshold - 1) < 0.01), color = "#D55E00", size = 3) +
    annotate("text", x = 0.15, y = c1_fdr + 0.05, label = sprintf("c=1: FDR=%.3f", c1_fdr),
             hjust = 0, size = 3, fontface = "bold", color = "#D55E00") +
    scale_x_continuous(name = "Threshold c (log scale)", breaks = c(-1, 0, 0.48, 1, 1.48, 2, 3, 4),
                      labels = c("0.1", "1", "3", "10", "30", "100", "1000", "10000")) +
    scale_y_continuous(limits = c(0, 1)) + labs(y = "FDR", title = "A") + theme_bw(base_size = 11)

  p_pow <- ggplot(df_plot, aes(x = log10_c, y = avg_power)) +
    geom_ribbon(aes(ymin = pmax(avg_power - avg_power_sd, 0), ymax = pmin(avg_power + avg_power_sd, 1)),
                fill = "#0072B2", alpha = 0.15) + geom_line(color = "#0072B2", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(data = df_plot %>% filter(abs(c_threshold - 1) < 0.01), color = "#0072B2", size = 3) +
    annotate("text", x = 0.15, y = c1_pow - 0.05, label = sprintf("c=1: Power=%.3f", c1_pow),
             hjust = 0, size = 3, fontface = "bold", color = "#0072B2") +
    scale_x_continuous(name = "Threshold c (log scale)", breaks = c(-1, 0, 0.48, 1, 1.48, 2, 3, 4),
                      labels = c("0.1", "1", "3", "10", "30", "100", "1000", "10000")) +
    scale_y_continuous(limits = c(0, 1)) + labs(y = "Avg Power", title = "B") + theme_bw(base_size = 11)

  fig_curves <- p_fdr + p_pow + plot_annotation(
    title = paste0("FDR and Power Curves  |  ", dataset),
    subtitle = paste0("Genes: ", n_genes_act, " (H0=", sum(sim_params$hypothesis == "H0"), ", H1=",
                      sum(sim_params$hypothesis == "H1"), ")  |  Replicates: ", n_reps,
                      "  |  Spots: ", nrow(coords), "  |  Mesh: ", mesh$n, " nodes")
  )
  curves_file <- file.path(out_dir, "test_fdr_power.png")
  ggsave(curves_file, fig_curves, width = 10, height = 4.5, dpi = 200, bg = "white")
  ggsave(sub("\\.png$", ".pdf", curves_file), fig_curves, width = 10, height = 4.5)
  msg("  Saved: ", curves_file)

  # --- Per-gene BF figure ---------------------------------------------------
  bf_df <- bind_rows(all_gene_bfs) %>% filter(!is.na(log_bf_sim)) %>%
    left_join(param_table %>% select(gene_name, hypothesis, magnitude), by = c("gene_name", "hypothesis")) %>%
    mutate(gene_label = ifelse(hypothesis == "H0", paste0(gene_name, " [H0]"),
                               sprintf("%s [H1] mag=%.2f", gene_name, magnitude)))
  gene_order <- bf_df %>% group_by(gene_label, hypothesis) %>% summarise(med_bf = median(log_bf_sim), .groups = "drop") %>% arrange(hypothesis, med_bf)
  bf_df$gene_label <- factor(bf_df$gene_label, levels = gene_order$gene_label)

  p_bf <- ggplot(bf_df, aes(x = gene_label, y = log_bf_sim, fill = hypothesis)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(alpha = 0.7, outlier.size = 1) + geom_jitter(aes(color = hypothesis), width = 0.15, size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = c(H0 = "#56B4E9", H1 = "#E69F00")) +
    scale_color_manual(values = c(H0 = "#0072B2", H1 = "#D55E00")) +
    labs(title = paste0("Per-Gene log(BF) Across ", n_reps, " Replicates  |  ", dataset),
         x = NULL, y = expression(log(BF[10])), fill = "Truth", color = "Truth") +
    coord_flip() + theme_bw(base_size = 11) + theme(legend.position = "bottom")

  bf_file <- file.path(out_dir, "test_gene_bf_distributions.png")
  ggsave(bf_file, p_bf, width = 9, height = max(4, 0.5 * n_genes_act + 1.5), dpi = 200, bg = "white")
  ggsave(sub("\\.png$", ".pdf", bf_file), p_bf, width = 9, height = max(4, 0.5 * n_genes_act + 1.5))
  msg("  Saved: ", bf_file)

  msg("")
  msg("TEST COMPLETE  |  Time: ", round(elapsed / 60, 1), " min")
  msg("Output: ", spatial_file, ", ", curves_file, ", ", bf_file)

  invisible(list(
    curves = curves,
    replicate_results = replicate_results,
    sim_params = sim_params,
    fdr_at_c1 = mean(fdr_1),
    power_at_c1 = mean(power_1),
    output_files = c(spatial = spatial_file, curves = curves_file, bf = bf_file),
    elapsed_min = elapsed / 60
  ))
}

# Null coalescing for optional args
`%||%` <- function(x, y) if (is.null(x)) y else x
