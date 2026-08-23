#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1]] else "data/spatial_datasets"

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]])))
} else {
  normalizePath("hpc")
}

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
})
source(file.path(script_dir, "spot_selection.R"))

make_spe <- function(in_tissue = NULL) {
  n <- if (is.null(in_tissue)) 4L else length(in_tissue)
  cd <- if (is.null(in_tissue)) {
    S4Vectors::DataFrame(row.names = paste0("spot", seq_len(n)))
  } else {
    S4Vectors::DataFrame(
      in_tissue = in_tissue,
      row.names = paste0("spot", seq_len(n))
    )
  }
  SpatialExperiment(
    assays = list(counts = matrix(seq_len(2L * n), nrow = 2L)),
    colData = cd,
    spatialCoords = cbind(x = seq_len(n), y = seq_len(n))
  )
}

selected <- select_in_tissue_spots(
  make_spe(c(1L, 0L, 1L, 0L)),
  verbose = FALSE
)
info <- attr(selected, "spot_selection")
stopifnot(
  ncol(selected) == 2L,
  identical(info$policy, "in_tissue_spots"),
  isTRUE(info$annotation_present),
  info$n_spots_supplied == 4L,
  info$n_spots_selected == 2L,
  info$n_spots_removed == 2L
)

unannotated <- select_in_tissue_spots(make_spe(), verbose = FALSE)
unannotated_info <- attr(unannotated, "spot_selection")
stopifnot(
  ncol(unannotated) == 4L,
  !isTRUE(unannotated_info$annotation_present),
  unannotated_info$n_spots_removed == 0L
)

invalid <- try(
  select_in_tissue_spots(make_spe(c("inside", "outside")), verbose = FALSE),
  silent = TRUE
)
stopifnot(inherits(invalid, "try-error"))

expected_counts <- c(
  Visium_humanDLPFC = 3639L,
  Visium_mouseCoronal = 2702L
)
for (dataset_name in names(expected_counts)) {
  path <- file.path(data_dir, paste0(dataset_name, ".rds"))
  if (!file.exists(path)) next
  spe <- select_in_tissue_spots(readRDS(path), dataset_name, verbose = FALSE)
  stopifnot(ncol(spe) == expected_counts[[dataset_name]])
}

cat("Spot-selection tests passed.\n")
