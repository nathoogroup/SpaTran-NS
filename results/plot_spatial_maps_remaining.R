#!/usr/bin/env Rscript
# =============================================================================
# Detrended spatial expression maps for remaining 11 tissues
# Two figures: Human (6 tissues) and Murine (5 tissues)
# Rows (i/ii/iii): non-stationary spatial / stationary spatial / stationary non-spatial
# Columns: one per tissue
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

# --- Helpers ---
normalize_expression <- function(expr) {
  qnorm((rank(expr, ties.method = "average") - 0.5) / length(expr))
}

detrend_expression <- function(coords, expr) {
  x_c <- coords$x - mean(coords$x)
  y_c <- coords$y - mean(coords$y)
  fit <- lm(expr ~ x_c + y_c + I(x_c^2) + I(y_c^2) + I(x_c * y_c),
            data = data.frame(expr = expr, x_c = x_c, y_c = y_c))
  expr - predict(fit)
}

pick_gene <- function(res, category, exclude = character(0)) {
  if (category == "nonstationary_spatial") {
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

make_map <- function(coords, values, gene, legend_title = "Residual",
                     residual_limit = 2) {
  df  <- data.frame(x = coords$x, y = coords$y, v = values)
  lim <- residual_limit
  if (lim == 0) lim <- max(abs(values), na.rm = TRUE)

  ggplot(df, aes(x = x, y = -y, color = v)) +
    geom_point(shape = 16, size = 0.85) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-lim, lim),
                          oob = squish,
                          name = legend_title,
                          breaks = c(-lim, 0, lim),
                          labels = c(sprintf("%.1f", -lim), "0", sprintf("%.1f", lim))) +
    labs(title = gene) +
    coord_fixed() +
    theme_void(base_size = 8) +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 8, face = "bold.italic"),
      legend.title      = element_text(size = 6, hjust = 0.5),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width  = unit(0.15, "cm"),
      legend.text       = element_text(size = 5),
      plot.background   = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      plot.margin       = margin(3, 3, 3, 3)
    )
}

make_strip <- function(label, size = 4, bold = TRUE, fill = "#E8E8E8", angle = 0) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             fontface = if (bold) "bold" else "plain",
             size = size, hjust = 0.5, vjust = 0.5, angle = angle) +
    theme_void() +
    theme(plot.background = element_rect(fill = fill, color = "grey70", linewidth = 0.4))
}

process_dataset <- function(ds_name, tissue_label, spe_file, results_file) {
  cat("  Processing", tissue_label, "...\n")
  res <- readRDS(results_file)
  res <- res[!is.na(res$log_bayes_factor) & !is.na(res$log_ml_nonspatial), ]

  g_ns  <- pick_gene(res, "nonstationary_spatial")
  g_ss  <- pick_gene(res, "stationary_spatial",    exclude = g_ns)
  g_snp <- pick_gene(res, "stationary_nonspatial", exclude = c(g_ns, g_ss))
  genes <- c(nonstationary_spatial = g_ns, stationary_spatial = g_ss, stationary_nonspatial = g_snp)
  cat("    Genes:", paste(names(genes), genes, sep = "=", collapse = ", "), "\n")

  spe <- readRDS(spe_file)
  rd  <- rowData(spe)
  g_names <- if ("symbol"    %in% colnames(rd)) as.character(rd$symbol)
             else if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
             else rownames(spe)

  coords   <- as.data.frame(spatialCoords(spe)); colnames(coords)[1:2] <- c("x", "y")
  lib_size <- Matrix::colSums(counts(spe))

  maps <- list()
  for (cat_name in names(genes)) {
    gene <- genes[[cat_name]]
    idx  <- which(g_names == gene)[1]
    if (is.na(idx)) { cat("    WARNING: gene", gene, "not found\n"); next }
    raw      <- as.numeric(counts(spe)[idx, ])
    log_norm <- log2(raw / (lib_size / 1e6) + 1)
    norm_q   <- normalize_expression(log_norm)
    detrended <- detrend_expression(coords, norm_q)
    maps[[cat_name]] <- make_map(coords, detrended, gene)
  }
  maps
}

# --- Build figure for one species ---
#
#  cols: [row lbl] [tissue1] [tissue2] ... [tissueN]
#  row1: [       ] [tissue1] [tissue2] ... [tissueN]   <- tissue headers
#  row2: [ (i)   ] [ map   ] [ map   ] ... [ map   ]
#  row3: [ (ii)  ] [ map   ] [ map   ] ... [ map   ]
#  row4: [ (iii) ] [ map   ] [ map   ] ... [ map   ]
#
build_figure <- function(all_maps, tissue_labels) {
  n <- length(tissue_labels)
  row_keys   <- c("nonstationary_spatial", "stationary_spatial", "stationary_nonspatial")
  row_labels <- c("(i)\nNon-stationary\nspatial",
                  "(ii)\nStationary\nspatial",
                  "(iii)\nStationary\nnon-spatial")

  # Header row: empty corner + tissue name strips
  all_pl <- list(make_strip("", fill = "white"))
  for (tl in tissue_labels)
    all_pl <- c(all_pl, list(make_strip(tl, size = 3.8, fill = "#D0D8E8")))

  # Data rows
  for (i in seq_along(row_keys)) {
    all_pl <- c(all_pl, list(make_strip(row_labels[i], size = 3.5, fill = "#F4F4F4")))
    for (maps in all_maps)
      all_pl <- c(all_pl, list(maps[[row_keys[i]]]))
  }

  # Design: (n+1) cols x 4 rows
  design <- list(area(1, 1))
  for (j in 2:(n + 1)) design <- c(design, list(area(1, j)))
  for (i in 2:4)
    for (j in 1:(n + 1)) design <- c(design, list(area(i, j)))

  wrap_plots(all_pl) +
    plot_layout(design  = do.call(c, design),
                heights = c(0.15, 1, 1, 1),
                widths  = c(0.45, rep(1, n)))
}

# --- Dataset definitions ---
human_datasets <- list(
  list(ds = "HumanSpinalCord",       tissue = "Spinal Cord",       spe = "HumanSpinalCord.rds"),
  list(ds = "Visium_humanDLPFC",     tissue = "DLPFC",             spe = "Visium_humanDLPFC.rds"),
  list(ds = "HumanGlioblastoma",     tissue = "Glioblastoma",      spe = "HumanGlioblastoma.rds"),
  list(ds = "HumanColorectalCancer", tissue = "Colorectal Cancer", spe = "HumanColorectalCancer.rds"),
  list(ds = "HumanLymphNode",        tissue = "Lymph Node",        spe = "HumanLymphNode.rds"),
  list(ds = "HumanHeart",            tissue = "Heart",             spe = "HumanHeart.rds")
)

murine_datasets <- list(
  list(ds = "MouseBrainCoronal",            tissue = "Brain (Coronal I)", spe = "MouseBrainCoronal.rds"),
  list(ds = "MouseBrainSagittalPosterior",  tissue = "Brain (Posterior)", spe = "MouseBrainSagittalPosterior.rds"),
  list(ds = "Visium_mouseCoronal",          tissue = "Brain (Coronal II)", spe = "Visium_mouseCoronal.rds"),
  list(ds = "MouseBrainSagittalAnterior",   tissue = "Brain (Anterior)",  spe = "MouseBrainSagittalAnterior.rds"),
  list(ds = "MouseKidneyCoronal",           tissue = "Kidney",            spe = "MouseKidneyCoronal.rds")
)

process_group <- function(datasets) {
  lapply(datasets, function(d)
    process_dataset(d$ds, d$tissue,
                    file.path(data_dir, d$spe),
                    file.path(results_dir, paste0(d$ds, "_results.rds"))))
}

cat("Processing Human datasets...\n")
human_maps <- process_group(human_datasets)
cat("Processing Murine datasets...\n")
murine_maps <- process_group(murine_datasets)

cat("Building Human figure...\n")
p_human <- build_figure(human_maps, sapply(human_datasets, `[[`, "tissue"))
ggsave(file.path(output_dir, "spatial_maps_human.png"),
       p_human, width = 18, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "spatial_maps_human.pdf"),
       p_human, width = 18, height = 11)
cat("Saved: spatial_maps_human.png\n")

cat("Building Murine figure...\n")
p_murine <- build_figure(murine_maps, sapply(murine_datasets, `[[`, "tissue"))
ggsave(file.path(output_dir, "spatial_maps_murine.png"),
       p_murine, width = 15, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "spatial_maps_murine.pdf"),
       p_murine, width = 15, height = 11)
cat("Saved: spatial_maps_murine.png\n")

cat("Done.\n")
