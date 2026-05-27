#!/usr/bin/env Rscript
# =============================================================================
# Spatial expression maps: spreadsheet-cell layout
# Page split: Breast Cancer (ILC) | Ovarian Cancer
# Columns A/B: A = no QN no detrend, B = QN + polynomial detrend
# Rows (i/ii/iii): non-stationary spatial / stationary spatial / stationary non-spatial
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_dir    <- if (length(args) >= 1) args[1] else "data/spatial_datasets"
results_dir <- if (length(args) >= 2) args[2] else "results"
output_dir  <- if (length(args) >= 3) args[3] else results_dir

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Matrix)
  library(scales)
})

# --- Helpers from analysis pipeline ---
normalize_expression <- function(expr) {
  expr_rank <- rank(expr, ties.method = "average")
  qnorm((expr_rank - 0.5) / length(expr))
}

detrend_expression <- function(coords, expr) {
  x_c <- coords$x - mean(coords$x)
  y_c <- coords$y - mean(coords$y)
  poly_fit <- lm(expr ~ x_c + y_c + I(x_c^2) + I(y_c^2) + I(x_c * y_c),
                 data = data.frame(expr = expr, x_c = x_c, y_c = y_c))
  list(detrended = expr - predict(poly_fit))
}

# --- Pick representative gene (excluding already-used symbols) ---
pick_gene <- function(res, category, exclude = character(0)) {
  if (category == "nonstationary") {
    sub <- res[res$log_bayes_factor > log(10) & !res$gene_name %in% exclude, ]
    sub[which.max(sub$log_bayes_factor), "gene_name"]
  } else if (category == "stationary_spatial") {
    sub <- res[res$log_bayes_factor < -log(10) &
               res$log_ml_stationary > res$log_ml_nonspatial &
               !res$gene_name %in% exclude, ]
    sub[which.max(sub$log_ml_stationary - sub$log_ml_nonspatial), "gene_name"]
  } else {
    sub <- res[res$log_bayes_factor < -log(10) &
               res$log_ml_nonspatial > res$log_ml_stationary &
               !res$gene_name %in% exclude, ]
    sub[which.max(sub$log_ml_nonspatial - sub$log_ml_stationary), "gene_name"]
  }
}

# --- Spatial map ---
# bf_label: a plotmath string for the BF annotation, e.g. "log~BF[NS*','*S]==354.0"
# Residuals are normal-score residuals; a shared +/-2 scale keeps colors
# comparable across genes while clipping only the most extreme values.
make_map <- function(coords, values, gene, tissue, legend_title, bf_label = NULL,
                     residual_limit = 2) {
  df  <- data.frame(x = coords$x, y = coords$y, v = values)
  lim <- residual_limit
  if (lim == 0) lim <- max(abs(values), na.rm = TRUE)

  p <- ggplot(df, aes(x = x, y = -y, color = v)) +
    geom_point(shape = 16, size = 0.95) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-lim, lim),
                          oob = squish,
                          name = legend_title,
                          breaks = c(-lim, 0, lim),
                          labels = c(sprintf("%.1f", -lim), "0", sprintf("%.1f", lim))) +
    labs(title = gene) +
    coord_fixed() +
    theme_void(base_size = 9) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 9, face = "bold.italic"),
      legend.title      = element_text(size = 6.5, hjust = 0.5),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width  = unit(0.18, "cm"),
      legend.text       = element_text(size = 5.5),
      plot.background   = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      plot.margin       = margin(4, 4, 4, 4)
    )
  if (!is.null(bf_label)) {
    p <- p + labs(caption = bf_label) +
      theme(plot.caption = element_text(size = 6.5, hjust = 0.5, color = "grey20",
                                        margin = margin(t = 2, b = 2)))
  }
  p
}

# --- Text strip (header / row label) ---
make_strip <- function(label, size = 4.5, bold = TRUE, fill = "#E8E8E8", angle = 0) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             fontface = if (bold) "bold" else "plain",
             size = size, hjust = 0.5, vjust = 0.5, angle = angle) +
    theme_void() +
    theme(plot.background = element_rect(fill = fill, color = "grey70", linewidth = 0.4))
}

# --- Process one dataset, return list of 6 maps (col_A/col_B x 3 rows) ---
process_dataset <- function(ds_name, spe_file, results_file) {
  cat("Processing", ds_name, "...\n")

  res <- readRDS(results_file)
  res <- res[!is.na(res$log_bayes_factor) & !is.na(res$log_ml_nonspatial), ]

  g_ns  <- pick_gene(res, "nonstationary")
  g_ss  <- pick_gene(res, "stationary_spatial",    exclude = g_ns)
  g_snp <- pick_gene(res, "stationary_nonspatial", exclude = c(g_ns, g_ss))
  genes <- c(nonstationary_spatial = g_ns,
             stationary_spatial    = g_ss,
             stationary_nonspatial = g_snp)
  cat("  Genes:", paste(names(genes), genes, sep = "=", collapse = ", "), "\n")

  spe <- readRDS(spe_file)
  rd  <- rowData(spe)
  g_names <- if ("symbol"    %in% colnames(rd)) as.character(rd$symbol)
             else if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
             else rownames(spe)

  coords   <- as.data.frame(spatialCoords(spe)); colnames(coords)[1:2] <- c("x", "y")
  lib_size <- Matrix::colSums(counts(spe))

  # BF label per row category
  bf_label <- function(cat_name, gene) {
    row <- res[res$gene_name == gene, ][1, ]
    if (cat_name == "nonstationary_spatial") {
      lbf <- round(row$log_bayes_factor, 1)
      sprintf("log BF(NS, S) = %s", lbf)
    } else if (cat_name == "stationary_spatial") {
      lbf <- round(row$log_ml_stationary - row$log_ml_nonspatial, 1)
      sprintf("log BF(S, NP) = %s", lbf)
    } else {
      lbf <- round(row$log_ml_nonspatial - row$log_ml_stationary, 1)
      sprintf("log BF(NP, S) = %s", lbf)
    }
  }

  plots <- list()
  for (cat_name in names(genes)) {
    gene <- genes[[cat_name]]
    idx  <- which(g_names == gene)[1]
    if (is.na(idx)) { cat("  WARNING: gene", gene, "not found\n"); next }

    raw      <- as.numeric(counts(spe)[idx, ])
    log_norm <- log2(raw / (lib_size / 1e6) + 1)
    norm_q   <- normalize_expression(log_norm)
    dt       <- detrend_expression(coords, norm_q)
    bfl      <- bf_label(cat_name, gene)

    plots[[cat_name]] <- make_map(coords, dt$detrended, gene, ds_name, "Residual",
                                  bf_label = bfl)
  }
  plots
}

# --- Build combined figure ---
#
#  cols: [row lbl] [Human1] [Human2] [Murine1] [Murine2]
#  row1: [       ] [  A: Human (spans 2)  ] [  B: Murine (spans 2) ]
#  row2: [       ] [tissue] [tissue]         [tissue]  [tissue]
#  row3: [ (i)   ] [ map  ] [  map ]         [ map  ]  [  map ]
#  row4: [ (ii)  ] [ map  ] [  map ]         [ map  ]  [  map ]
#  row5: [ (iii) ] [ map  ] [  map ]         [ map  ]  [  map ]
#
build_combined <- function(human_plots, murine_plots, human_names, murine_names) {
  row_keys <- c("nonstationary_spatial", "stationary_spatial", "stationary_nonspatial")
  row_labels <- c("(i)\nNon-stationary\nspatial",
                  "(ii)\nStationary\nspatial",
                  "(iii)\nStationary\nnon-spatial")

  h_empty  <- make_strip("",       size = 4, fill = "white")
  h_human  <- make_strip("A\nHuman",  size = 4, fill = "#D0D8E8")
  h_murine <- make_strip("B\nMurine", size = 4, fill = "#D0D8E8")

  sub_strips <- lapply(c(human_names, murine_names),
                       function(n) make_strip(n, size = 4, fill = "#E8EEF4"))

  all_pl <- c(list(h_empty, h_human, h_murine,
                   h_empty), sub_strips)   # h_empty = row-label col in row 2

  for (i in seq_along(row_keys)) {
    key <- row_keys[i]
    all_pl <- c(all_pl, list(
      make_strip(row_labels[i], size = 3.5, fill = "#F4F4F4"),
      human_plots[[1]][[key]],
      human_plots[[2]][[key]],
      murine_plots[[1]][[key]],
      murine_plots[[2]][[key]]
    ))
  }

  design <- c(
    area(1, 1),
    area(1, 2, 1, 3),   # A: Human spans cols 2-3
    area(1, 4, 1, 5),   # B: Murine spans cols 4-5
    area(2, 1), area(2, 2), area(2, 3), area(2, 4), area(2, 5),
    area(3, 1), area(3, 2), area(3, 3), area(3, 4), area(3, 5),
    area(4, 1), area(4, 2), area(4, 3), area(4, 4), area(4, 5),
    area(5, 1), area(5, 2), area(5, 3), area(5, 4), area(5, 5)
  )

  wrap_plots(all_pl) +
    plot_layout(design  = design,
                heights = c(0.18, 0.08, 1, 1, 1),
                widths  = c(0.28, 1, 1, 1, 1))
}

# --- Run ---
human_datasets <- list(
  list(name = "Breast Cancer (ILC)",
       spe  = file.path(data_dir, "HumanBreastCancerILC.rds"),
       res  = file.path(results_dir, "HumanBreastCancerILC_results.rds")),
  list(name = "Ovarian Cancer",
       spe  = file.path(data_dir, "HumanOvarianCancer.rds"),
       res  = file.path(results_dir, "HumanOvarianCancer_results.rds"))
)

murine_datasets <- list(
  list(name = "Brain (Coronal I)",
       spe  = file.path(data_dir, "MouseBrainCoronal.rds"),
       res  = file.path(results_dir, "MouseBrainCoronal_results.rds")),
  list(name = "Brain (Coronal II)",
       spe  = file.path(data_dir, "Visium_mouseCoronal.rds"),
       res  = file.path(results_dir, "Visium_mouseCoronal_results.rds"))
)

human_plots  <- lapply(human_datasets,  function(d) process_dataset(d$name, d$spe, d$res))
murine_plots <- lapply(murine_datasets, function(d) process_dataset(d$name, d$spe, d$res))

p_combined <- build_combined(
  human_plots, murine_plots,
  sapply(human_datasets,  `[[`, "name"),
  sapply(murine_datasets, `[[`, "name")
)

out_png <- file.path(output_dir, "spatial_maps_breast_ovarian.png")
out_pdf <- file.path(output_dir, "spatial_maps_breast_ovarian.pdf")

ggsave(out_png, p_combined, width = 18, height = 13, dpi = 300, bg = "white")
cat("Saved:", out_png, "\n")
ggsave(out_pdf, p_combined, width = 18, height = 13)
cat("Saved:", out_pdf, "\n")
cat("Done.\n")
