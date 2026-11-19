# Select the tissue spots used by the analysis and plotting workflows.
# Bioconductor objects that expose an in_tissue annotation are filtered to its
# true entries. Objects without that annotation are treated as already
# tissue-restricted and retain all supplied spots.
select_in_tissue_spots <- function(spe, dataset_name = "dataset") {
  if (!inherits(spe, "SpatialExperiment")) {
    stop(dataset_name, " is not a SpatialExperiment", call. = FALSE)
  }

  annotations <- SummarizedExperiment::colData(spe)
  if (!"in_tissue" %in% colnames(annotations)) {
    return(spe)
  }

  raw_values <- annotations$in_tissue
  values <- tolower(trimws(as.character(raw_values)))
  true_values <- c("1", "true", "t", "yes", "y")
  false_values <- c("0", "false", "f", "no", "n")
  recognized <- is.na(raw_values) | values %in% c(true_values, false_values)
  if (!all(recognized)) {
    stop(dataset_name, " has unrecognized in_tissue values", call. = FALSE)
  }

  keep <- !is.na(raw_values) & values %in% true_values
  if (!any(keep)) {
    stop(dataset_name, " has no spots marked in_tissue", call. = FALSE)
  }
  spe[, keep, drop = FALSE]
}
