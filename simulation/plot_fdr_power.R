#!/usr/bin/env Rscript
# =============================================================================
# Simulation Study: FDR and Power Curves Figure
# =============================================================================
# Creates a publication-quality figure showing FDR and average power curves
# as a function of the BF threshold c for both simulation settings.
#
# Usage:
#   Rscript simulation/plot_fdr_power.R [output_dir] [fig_file]
# =============================================================================

args       <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "simulation/output"
fig_file   <- if (length(args) >= 2) args[2] else "simulation/output/fdr_power_figure.pdf"

cat("Creating FDR/Power figure...\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# --- Load curves for both datasets -------------------------------------------

datasets <- c("HumanBreastCancerILC", "HumanOvarianCancer")
labels   <- c("Human Breast Cancer (ILC)", "Human Ovarian Cancer")

all_curves <- list()
for (i in seq_along(datasets)) {
  ds <- datasets[i]
  f  <- file.path(output_dir, paste0(ds, "_fdr_power_curves.rds"))
  if (!file.exists(f)) {
    cat("  WARNING: missing curves for", ds, "- skipping\n")
    next
  }
  curves <- readRDS(f)
  curves$dataset_label <- labels[i]
  all_curves[[ds]] <- curves
}

if (length(all_curves) == 0) stop("No curve files found.")
df <- bind_rows(all_curves)

# --- Use log10(c) for x-axis, filter out c=0 --------------------------------

df_plot <- df %>%
  filter(c_threshold > 0) %>%
  mutate(log10_c = log10(c_threshold))

# --- Color palette -----------------------------------------------------------

cols <- c("Human Breast Cancer (ILC)" = "#D55E00",
          "Human Ovarian Cancer"      = "#0072B2")

# --- FDR panel ---------------------------------------------------------------

p_fdr <- ggplot(df_plot, aes(x = log10_c, y = fdr, color = dataset_label)) +
  geom_ribbon(aes(ymin = pmax(fdr - fdr_sd, 0),
                  ymax = pmin(fdr + fdr_sd, 1),
                  fill = dataset_label),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = 0.05, y = 0.98, label = "c = 1",
           hjust = 0, vjust = 1, size = 3, color = "grey40") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(
    breaks = c(-1, 0, log10(3), 1, log10(30), 2, 3, 4),
    labels = c("0.1", "1", "3", "10", "30", "100", "1000", "10000")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    x     = expression(paste("Threshold ", italic(c))),
    y     = "False Discovery Rate",
    color = NULL, fill = NULL,
    title = "A"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
    legend.key.size  = unit(0.8, "lines"),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# --- Power panel -------------------------------------------------------------

p_pow <- ggplot(df_plot, aes(x = log10_c, y = avg_power, color = dataset_label)) +
  geom_ribbon(aes(ymin = pmax(avg_power - avg_power_sd, 0),
                  ymax = pmin(avg_power + avg_power_sd, 1),
                  fill = dataset_label),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = 0.05, y = 0.98, label = "c = 1",
           hjust = 0, vjust = 1, size = 3, color = "grey40") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(
    breaks = c(-1, 0, log10(3), 1, log10(30), 2, 3, 4),
    labels = c("0.1", "1", "3", "10", "30", "100", "1000", "10000")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    x     = expression(paste("Threshold ", italic(c))),
    y     = "Average Power",
    color = NULL, fill = NULL,
    title = "B"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# --- Combine panels ----------------------------------------------------------

fig <- p_fdr + p_pow + plot_layout(ncol = 2)

# --- Annotate c=1 values in subtitle ----------------------------------------

c1_summary <- df_plot %>%
  filter(abs(c_threshold - 1) < 0.01) %>%
  mutate(label = sprintf("%s: FDR=%.3f, Power=%.3f",
                         dataset_label, fdr, avg_power))

if (nrow(c1_summary) > 0) {
  subtitle_text <- paste("At c = 1:", paste(c1_summary$label, collapse = "  |  "))
  fig <- fig + plot_annotation(
    subtitle = subtitle_text,
    theme = theme(plot.subtitle = element_text(size = 9, color = "grey30"))
  )
}

# --- Save --------------------------------------------------------------------

ggsave(fig_file, fig, width = 10, height = 4.5, dpi = 300)
cat("Saved figure:", fig_file, "\n")

fig_png <- sub("\\.pdf$", ".png", fig_file)
ggsave(fig_png, fig, width = 10, height = 4.5, dpi = 300)
cat("Saved PNG:  ", fig_png, "\n")

# --- Print table of key thresholds -------------------------------------------

cat("\n=============================================================\n")
cat("KEY THRESHOLD VALUES\n")
cat("=============================================================\n")

key_c <- c(1, 3, 10, 30, 100)
for (ds in unique(df$dataset_label)) {
  cat("\n", ds, ":\n", sep = "")
  cat(sprintf("  %8s  %8s  %12s\n", "c", "FDR", "Avg Power"))
  cat("  ", strrep("-", 32), "\n", sep = "")
  for (ct in key_c) {
    row <- df %>%
      filter(dataset_label == ds, abs(c_threshold - ct) < 0.01)
    if (nrow(row) > 0) {
      cat(sprintf("  %8.0f  %8.4f  %12.4f\n", ct, row$fdr[1], row$avg_power[1]))
    }
  }
}

cat("\nDone!\n")
