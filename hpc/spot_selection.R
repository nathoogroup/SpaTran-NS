# =============================================================================
# Study-wide Visium spot selection
# =============================================================================
# Retain spots explicitly marked as in tissue. SpatialExperiment objects that
# do not provide an in_tissue field are treated as already tissue-restricted
# and are left unchanged.

in_tissue_spot_mask <- function(spe) {
  if (!inherits(spe, "SpatialExperiment")) {
    stop("spot selection requires a SpatialExperiment object", call. = FALSE)
  }

  if (!"in_tissue" %in% colnames(SummarizedExperiment::colData(spe))) {
    return(rep(TRUE, ncol(spe)))
  }

  raw_values <- SummarizedExperiment::colData(spe)$in_tissue
  values <- tolower(trimws(as.character(raw_values)))
  true_values <- c("1", "true", "t", "yes", "y")
  false_values <- c("0", "false", "f", "no", "n")
  recognized <- is.na(raw_values) | values %in% c(true_values, false_values)

  if (!all(recognized)) {
    stop(
      "in_tissue contains unrecognized values: ",
      paste(sort(unique(values[!recognized])), collapse = ", "),
      call. = FALSE
    )
  }

  keep <- !is.na(raw_values) & values %in% true_values
  if (!any(keep)) {
    stop("in_tissue annotation contains no tissue spots", call. = FALSE)
  }
  keep
}

select_in_tissue_spots <- function(spe, dataset_name = NULL, verbose = TRUE) {
  n_supplied <- ncol(spe)
  has_annotation <- "in_tissue" %in%
    colnames(SummarizedExperiment::colData(spe))
  keep <- in_tissue_spot_mask(spe)
  selected <- spe[, keep, drop = FALSE]
  n_selected <- ncol(selected)
  n_removed <- n_supplied - n_selected

  if (isTRUE(verbose)) {
    prefix <- if (!is.null(dataset_name) && nzchar(dataset_name)) {
      paste0("  ", dataset_name, ": ")
    } else {
      "  "
    }
    if (has_annotation) {
      cat(sprintf(
        "%sretained %d/%d in-tissue spots; removed %d off-tissue spots\n",
        prefix, n_selected, n_supplied, n_removed
      ))
    } else {
      cat(sprintf(
        "%sno in_tissue annotation; retained all %d supplied spots as tissue spots\n",
        prefix, n_selected
      ))
    }
  }

  attr(selected, "spot_selection") <- list(
    policy = "in_tissue_spots",
    annotation_present = has_annotation,
    n_spots_supplied = n_supplied,
    n_spots_selected = n_selected,
    n_spots_removed = n_removed
  )
  selected
}
