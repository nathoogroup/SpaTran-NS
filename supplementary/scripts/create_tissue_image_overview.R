#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "supplementary/scripts/create_tissue_image_overview.R"
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "data"))) root <- normalizePath(getwd(), mustWork = TRUE)

out_dir <- file.path(root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SpatialExperiment)
  library(ggplot2)
  library(patchwork)
})

datasets <- data.frame(
  dataset_name = c(
    "HumanBreastCancerILC", "HumanSpinalCord", "Visium_humanDLPFC",
    "HumanGlioblastoma", "HumanColorectalCancer", "HumanLymphNode",
    "HumanOvarianCancer", "HumanHeart", "MouseBrainCoronal",
    "MouseBrainSagittalPosterior", "Visium_mouseCoronal",
    "MouseBrainSagittalAnterior", "MouseKidneyCoronal"
  ),
  label = c(
    "Breast cancer (ILC)", "Spinal cord", "DLPFC", "Glioblastoma",
    "Colorectal cancer", "Lymph node", "Ovarian cancer", "Heart",
    "Brain (coronal I)", "Brain (posterior)", "Brain (coronal II)",
    "Brain (anterior)", "Kidney"
  ),
  species = c(rep("Human", 8), rep("Murine", 5)),
  stringsAsFactors = FALSE
)

manifest <- read.csv(file.path(root, "data", "spatial_datasets", "dataset_manifest.csv"),
  stringsAsFactors = FALSE
)
datasets <- merge(datasets, manifest[, c("eh_id", "title")],
  by.x = "dataset_name", by.y = "title", all.x = TRUE, sort = FALSE
)

# The EH9516/EH9517 records used when the original manifest was created have
# since been removed from ExperimentHub. Current replacement records provide the
# same named SpatialExperiment datasets and include image data.
datasets$eh_id[datasets$dataset_name == "Visium_humanDLPFC"] <- "EH9628"
datasets$eh_id[datasets$dataset_name == "Visium_mouseCoronal"] <- "EH9629"

eh <- ExperimentHub()

make_image_plot <- function(eh_id, label) {
  spe <- eh[[eh_id]]
  if (nrow(imgData(spe)) == 0) stop("No image data available for ", eh_id)
  img <- imgRaster(imgData(spe)$data[[1]])
  ggplot() +
    annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_fixed(expand = FALSE) +
    labs(title = label) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 2)),
      plot.background = element_rect(fill = "white", color = "grey75", linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2)
    )
}

plots <- lapply(seq_len(nrow(datasets)), function(i) {
  cat("Loading image:", datasets$dataset_name[i], datasets$eh_id[i], "\n")
  make_image_plot(datasets$eh_id[i], datasets$label[i])
})

human_idx <- which(datasets$species == "Human")
murine_idx <- which(datasets$species == "Murine")

human_row <- wrap_plots(plots[human_idx], nrow = 1) +
  plot_annotation(title = "Human",
    theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)))
murine_row <- wrap_plots(plots[murine_idx], nrow = 1) +
  plot_annotation(title = "Murine",
    theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)))

overview <- human_row / murine_row + plot_layout(heights = c(1, 1))

pdf_file <- file.path(out_dir, "tissue_image_overview.pdf")
png_file <- file.path(out_dir, "tissue_image_overview.png")
ggsave(pdf_file, overview, width = 11, height = 4.6, bg = "white")
ggsave(png_file, overview, width = 11, height = 4.6, dpi = 300, bg = "white")
cat("Saved:", pdf_file, "\n")
cat("Saved:", png_file, "\n")
