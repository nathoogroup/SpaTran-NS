#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "supplementary/scripts/create_supplementary_materials.R"
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "results"))) root <- normalizePath(getwd(), mustWork = TRUE)

out_dir <- file.path(root, "supplementary")
add_dir <- file.path(out_dir, "additional_files")
tab_dir <- file.path(out_dir, "generated_tables")
dir.create(add_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

analysis_implementation_files <- c("analysis_functions.R", "run_analysis.R")
analysis_implementation_paths <- file.path(
  root, "hpc", analysis_implementation_files
)
missing_implementation_files <- analysis_implementation_paths[
  !file.exists(analysis_implementation_paths)
]
if (length(missing_implementation_files) > 0L) {
  stop(
    "Missing analysis implementation files: ",
    paste(missing_implementation_files, collapse = ", ")
  )
}
expected_analysis_implementation_md5 <- stats::setNames(
  unname(as.character(tools::md5sum(analysis_implementation_paths))),
  analysis_implementation_files
)

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
  "bf_interpretation", "error_message",
  "dataset_name", "dataset_file", "n_spots", "n_genes_total"
)

read_rds_results <- function(results_dir, dataset_names) {
  files <- file.path(results_dir, paste0(dataset_names, "_results.rds"))
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop("Missing study result files: ", paste(basename(missing_files), collapse = ", "))
  }

  rows <- Map(function(f, expected_dataset) {
    x <- readRDS(f)
    if (!is.data.frame(x)) {
      stop(basename(f), " is not a data.frame")
    }
    metadata <- attr(x, "analysis_metadata")
    missing_columns <- setdiff(analysis_result_required_columns, names(x))
    schema_version <- if (is.list(metadata) && length(metadata$schema_version) == 1L) {
      suppressWarnings(as.integer(metadata$schema_version))
    } else {
      NA_integer_
    }
    n_expected <- if (is.list(metadata) && length(metadata$n_genes_analyzed) == 1L) {
      suppressWarnings(as.integer(metadata$n_genes_analyzed))
    } else {
      NA_integer_
    }
    n_spots_expected <- if (is.list(metadata) && length(metadata$n_spots) == 1L) {
      suppressWarnings(as.integer(metadata$n_spots))
    } else {
      NA_integer_
    }
    n_genes_total_expected <- if (is.list(metadata) && length(metadata$n_genes_total) == 1L) {
      suppressWarnings(as.integer(metadata$n_genes_total))
    } else {
      NA_integer_
    }
    gene_indices <- suppressWarnings(as.integer(x$gene_index))
    failed_fit <- is.na(x$log_bayes_factor)
    has_error <- !is.na(x$error_message) & nzchar(x$error_message)
    dataset_files <- unique(as.character(x$dataset_file))
    result_spot_counts <- unique(suppressWarnings(as.integer(x$n_spots)))
    result_gene_counts <- unique(suppressWarnings(as.integer(x$n_genes_total)))
    implementation_md5 <- if (is.list(metadata)) {
      metadata$analysis_implementation_md5
    } else {
      NULL
    }
    implementation_matches <- is.character(implementation_md5) &&
      identical(names(implementation_md5), names(expected_analysis_implementation_md5)) &&
      identical(
        tolower(unname(implementation_md5)),
        tolower(unname(expected_analysis_implementation_md5))
      )
    numeric_columns_ok <- all(vapply(
      x[intersect(analysis_result_numeric_columns, names(x))],
      is.numeric,
      logical(1)
    ))
    character_columns <- intersect(
      c("gene_id", "gene_name", "bf_interpretation", "error_message"),
      names(x)
    )
    character_columns_ok <- all(vapply(
      x[character_columns],
      is.character,
      logical(1)
    ))
    successful_fit <- !failed_fit
    finite_success_columns <- setdiff(
      analysis_result_numeric_columns,
      "bayes_factor"
    )
    successful_values <- if (!length(missing_columns) && numeric_columns_ok) {
      as.matrix(x[successful_fit, finite_success_columns, drop = FALSE])
    } else {
      matrix(NA_real_, nrow = 1L)
    }
    failed_values <- if (!length(missing_columns) && numeric_columns_ok) {
      as.matrix(x[failed_fit, analysis_result_numeric_columns, drop = FALSE])
    } else {
      matrix(0, nrow = 1L)
    }
    expected_log_bf <- if (!length(missing_columns) && numeric_columns_ok) {
      x$log_ml_nonstationary - x$log_ml_stationary
    } else {
      rep(NA_real_, nrow(x))
    }
    log_bf_consistent <- !any(
      successful_fit &
        abs(x$log_bayes_factor - expected_log_bf) >
          1e-8 * (1 + abs(expected_log_bf))
    )

    valid <- !length(missing_columns) &&
      is.list(metadata) &&
      isTRUE(metadata$complete) &&
      identical(schema_version, 3L) &&
      implementation_matches &&
      identical(as.character(metadata$dataset_name), expected_dataset) &&
      length(metadata$source_rds_md5) == 1L &&
      grepl("^[0-9a-fA-F]{32}$", as.character(metadata$source_rds_md5)) &&
      identical(
        as.character(metadata$spot_selection),
        "all_bioconductor_supplied_spots"
      ) &&
      !is.na(n_expected) &&
      identical(nrow(x), n_expected) &&
      identical(gene_indices, seq_len(n_expected)) &&
      identical(unique(as.character(x$dataset_name)), expected_dataset) &&
      length(dataset_files) == 1L &&
      identical(dataset_files, as.character(metadata$dataset_file)) &&
      !is.na(n_spots_expected) &&
      identical(result_spot_counts, n_spots_expected) &&
      !is.na(n_genes_total_expected) &&
      identical(result_gene_counts, n_genes_total_expected) &&
      numeric_columns_ok &&
      character_columns_ok &&
      (length(successful_values) == 0L || all(is.finite(successful_values))) &&
      !any(successful_fit & (is.na(x$bayes_factor) | x$bayes_factor < 0)) &&
      !any(successful_fit &
        (is.na(x$bf_interpretation) | !nzchar(x$bf_interpretation))) &&
      (length(failed_values) == 0L || all(is.na(failed_values))) &&
      !any(failed_fit & !is.na(x$bf_interpretation)) &&
      log_bf_consistent &&
      !any(failed_fit & !has_error) &&
      !any(!failed_fit & has_error)

    if (!valid) {
      details <- if (length(missing_columns)) {
        paste0("missing columns: ", paste(missing_columns, collapse = ", "))
      } else {
        "schema-v3 completeness or dataset metadata mismatch"
      }
      stop(basename(f), " is not a complete validated result (", details, ")")
    }
    x$source_result_file <- basename(f)
    x
  }, files, dataset_names)

  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (m in missing) x[[m]] <- NA
    x[all_names]
  })
  do.call(rbind, rows)
}

bf_strength <- function(log_bf) {
  out <- rep(NA_character_, length(log_bf))
  a <- abs(log_bf)
  out[!is.na(a) & a < log(3)] <- "Anecdotal"
  out[!is.na(a) & a >= log(3) & a < log(10)] <- "Moderate"
  out[!is.na(a) & a >= log(10) & a < log(30)] <- "Strong"
  out[!is.na(a) & a >= log(30) & a < log(100)] <- "Very strong"
  out[!is.na(a) & a >= log(100)] <- "Extreme"
  out
}

csv_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

latex_escape <- function(x) {
  x <- csv_escape(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

write_simple_table <- function(df, path, caption = NULL, label = NULL, digits = 3,
                               font_size = NULL, resize = FALSE) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (!is.null(caption)) {
    writeLines(sprintf("\\begin{table}[htbp]\n\\centering\n\\caption{%s}", latex_escape(caption)), con)
  } else {
    writeLines("\\begin{table}[htbp]\n\\centering", con)
  }
  if (!is.null(label)) writeLines(sprintf("\\label{%s}", label), con)
  if (!is.null(font_size)) writeLines(sprintf("\\%s", font_size), con)
  if (resize) writeLines("\\resizebox{\\linewidth}{!}{%", con)
  align <- paste(vapply(df, function(col) if (is.numeric(col)) "r" else "l", character(1)), collapse = "")
  writeLines(sprintf("\\begin{tabular}{%s}\n\\toprule", align), con)
  writeLines(paste(latex_escape(names(df)), collapse = " & "), con)
  writeLines(" \\\\\n\\midrule", con)
  whole_cols <- vapply(df, function(col) {
    if (!is.numeric(col)) return(FALSE)
    x <- col[!is.na(col)]
    length(x) > 0 && all(abs(x - round(x)) < 1e-9)
  }, logical(1))
  for (i in seq_len(nrow(df))) {
    vals <- vapply(seq_along(df), function(j) {
      v <- df[[j]][i]
      if (is.numeric(v)) {
        if (is.na(v)) {
          ""
        } else if (whole_cols[j]) {
          format(as.integer(round(v)), trim = TRUE, scientific = FALSE)
        } else {
          format(round(v, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
        }
      } else {
        latex_escape(v)
      }
    }, character(1))
    writeLines(paste(vals, collapse = " & "), con)
    writeLines(" \\\\", con)
  }
  writeLines("\\bottomrule\n\\end{tabular}%", con)
  if (resize) writeLines("}", con)
  writeLines("\\end{table}", con)
}

sample_type_map <- data.frame(
  dataset_name = c(
    "HumanBreastCancerILC", "HumanSpinalCord", "MouseBrainSagittalPosterior",
    "Visium_humanDLPFC", "Visium_mouseCoronal",
    "MouseBrainSagittalAnterior", "HumanGlioblastoma", "HumanColorectalCancer",
    "HumanLymphNode", "HumanOvarianCancer", "MouseKidneyCoronal", "HumanHeart",
    "HumanBreastCancerIDC", "HumanCerebellum"
  ),
  species_group = c(
    "Human", "Human", "Murine", "Human", "Murine", "Murine",
    "Human", "Human", "Human", "Human", "Murine", "Human", "Human", "Human"
  ),
  tissue = c(
    "Breast cancer (ILC)", "Spinal cord", "Brain (posterior)",
    "DLPFC", "Brain (coronal)", "Brain (anterior)", "Glioblastoma",
    "Colorectal cancer", "Lymph node", "Ovarian cancer", "Kidney", "Heart",
    "Breast cancer (IDC)", "Cerebellum"
  ),
  stringsAsFactors = FALSE
)

space_ranger_map <- data.frame(
  dataset_name = c(
    "HumanBreastCancerILC", "HumanSpinalCord", "MouseBrainSagittalPosterior",
    "Visium_humanDLPFC", "Visium_mouseCoronal",
    "MouseBrainSagittalAnterior", "HumanGlioblastoma", "HumanColorectalCancer",
    "HumanLymphNode", "HumanOvarianCancer", "MouseKidneyCoronal", "HumanHeart"
  ),
  space_ranger_version = c(
    "1.2.0", "1.2.0", "1.1.0", "1.0.0", "1.0.0", "1.1.0",
    "1.2.0", "1.2.0", "1.1.0", "1.2.0", "1.1.0", "1.1.0"
  ),
  space_ranger_version_source = c(
    rep("10x Genomics source dataset archived in ExperimentHub", 3),
    "spatialLIBD/Maynard et al. source documentation",
    "10x Genomics source dataset used by STexampleData",
    rep("10x Genomics source dataset archived in ExperimentHub", 7)
  ),
  space_ranger_source_url = c(
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_SpinalCord",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8095368/",
    "https://www.10xgenomics.com/datasets/mouse-brain-section-coronal-1-standard-1-0-0",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_OvarianCancer",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Kidney",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Heart"
  ),
  stringsAsFactors = FALSE
)

dataset_provenance_map <- data.frame(
  dataset_name = c(
    "HumanBreastCancerILC", "HumanSpinalCord", "MouseBrainSagittalPosterior",
    "Visium_humanDLPFC", "Visium_mouseCoronal",
    "MouseBrainSagittalAnterior", "HumanGlioblastoma", "HumanColorectalCancer",
    "HumanLymphNode", "HumanOvarianCancer", "MouseKidneyCoronal", "HumanHeart"
  ),
  repository = c(
    rep("ExperimentHub (TENxVisiumData)", 3),
    rep("ExperimentHub (STexampleData)", 2),
    rep("ExperimentHub (TENxVisiumData)", 7)
  ),
  accession = c(
    "EH6696", "EH6703", "EH6705", "EH9628", "EH9629", "EH6706",
    "EH6699", "EH6698", "EH6701", "EH6702", "EH6707", "EH6700"
  ),
  stringsAsFactors = FALSE
)

analysis_dataset_file <- file.path(root, "hpc", "analysis_datasets.txt")
if (!file.exists(analysis_dataset_file)) {
  stop("Study dataset list not found: ", analysis_dataset_file)
}
analysis_datasets <- trimws(readLines(analysis_dataset_file, warn = FALSE))
analysis_datasets <- analysis_datasets[
  nzchar(analysis_datasets) & !startsWith(analysis_datasets, "#")
]
if (length(analysis_datasets) != 12L || anyDuplicated(analysis_datasets)) {
  stop("Study dataset list must contain exactly 12 unique basenames")
}
results <- read_rds_results(file.path(root, "results"), analysis_datasets)
results$log_bayes_factor <- as.numeric(results$log_bayes_factor)
results$bayes_factor <- as.numeric(results$bayes_factor)
results$log_bf_nonstationary_vs_stationary <- results$log_bayes_factor
results$bf_strength_nonstationary_vs_stationary <- bf_strength(results$log_bayes_factor)
results$model_favored_nonstationary_vs_stationary <- ifelse(
  is.na(results$log_bayes_factor), NA,
  ifelse(results$log_bayes_factor > 0, "nonstationary_spatial", "stationary_spatial")
)

if ("log_ml_nonspatial" %in% names(results)) {
  results$log_bf_stationary_vs_nonspatial <- as.numeric(results$log_ml_stationary) -
    as.numeric(results$log_ml_nonspatial)
  results$bf_stationary_vs_nonspatial <- exp(results$log_bf_stationary_vs_nonspatial)
  results$bf_strength_stationary_vs_nonspatial <- bf_strength(results$log_bf_stationary_vs_nonspatial)
  results$final_three_model_category <- ifelse(
    is.na(results$log_bayes_factor), NA,
    ifelse(results$log_bayes_factor > 0, "nonstationary_spatial",
      ifelse(is.na(results$log_bf_stationary_vs_nonspatial), "stationary_unclassified",
        ifelse(results$log_bf_stationary_vs_nonspatial > 0, "stationary_spatial", "stationary_nonspatial")
      )
    )
  )
} else {
  results$log_bf_stationary_vs_nonspatial <- NA_real_
  results$bf_stationary_vs_nonspatial <- NA_real_
  results$bf_strength_stationary_vs_nonspatial <- NA_character_
  results$final_three_model_category <- NA_character_
}

per_gene_cols <- c(
  "dataset_name", "gene_index", "gene_id", "gene_name", "log_bayes_factor",
  "bf_strength_nonstationary_vs_stationary",
  "model_favored_nonstationary_vs_stationary",
  "log_bf_stationary_vs_nonspatial",
  "bf_strength_stationary_vs_nonspatial",
  "final_three_model_category", "error_message"
)
per_gene_cols <- intersect(per_gene_cols, names(results))
per_gene_results <- results[, per_gene_cols, drop = FALSE]
num_cols <- vapply(per_gene_results, is.numeric, logical(1))
per_gene_results[num_cols] <- lapply(per_gene_results[num_cols], function(x) signif(x, 6))
per_gene_path <- file.path(add_dir, "Additional_file_3_per_gene_bayes_factor_and_category_summary.csv")
unlink(file.path(add_dir, "Additional_file_3_full_gene_model_results.csv"))
write.csv(per_gene_results, per_gene_path, row.names = FALSE)

valid <- results[!is.na(results$log_bayes_factor), , drop = FALSE]
split_results <- split(results, factor(results$dataset_name, levels = analysis_datasets))
if (length(split_results) != length(analysis_datasets) || any(!lengths(split_results))) {
  stop("Validated results do not contain all 12 study datasets")
}
dataset_summary <- do.call(rbind, lapply(split_results, function(d_all) {
  d <- d_all[!is.na(d_all$log_bayes_factor), , drop = FALSE]
  if (!nrow(d)) stop("A dataset has no successful gene fits")
  stationary_candidates <- d[d$log_bayes_factor < 0, , drop = FALSE]
  data.frame(
    dataset_name = d$dataset_name[1],
    dataset_file = if ("dataset_file" %in% names(d)) d$dataset_file[1] else NA,
    n_spots = suppressWarnings(as.integer(d$n_spots[1])),
    n_genes_total = suppressWarnings(as.integer(d$n_genes_total[1])),
    n_genes_analyzed = nrow(d_all),
    n_gene_fits_failed = nrow(d_all) - nrow(d),
    n_nonstationary = sum(d$log_bayes_factor > 0, na.rm = TRUE),
    pct_nonstationary = 100 * mean(d$log_bayes_factor > 0, na.rm = TRUE),
    n_stationary_favored = sum(d$log_bayes_factor < 0, na.rm = TRUE),
    n_stationary_spatial = sum(d$final_three_model_category == "stationary_spatial", na.rm = TRUE),
    n_stationary_nonspatial = sum(d$final_three_model_category == "stationary_nonspatial", na.rm = TRUE),
    n_stationary_unclassified = sum(d$final_three_model_category == "stationary_unclassified", na.rm = TRUE),
    pct_stationary_spatial_among_stationary_favored =
      ifelse(nrow(stationary_candidates) > 0,
        100 * mean(stationary_candidates$final_three_model_category == "stationary_spatial", na.rm = TRUE),
        NA_real_),
    median_log_bf_ns_vs_s = median(d$log_bayes_factor, na.rm = TRUE),
    min_log_bf_ns_vs_s = min(d$log_bayes_factor, na.rm = TRUE),
    max_log_bf_ns_vs_s = max(d$log_bayes_factor, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
row.names(dataset_summary) <- NULL
dataset_summary <- merge(dataset_summary, sample_type_map, by = "dataset_name", all.x = TRUE)
dataset_summary <- merge(dataset_summary, space_ranger_map, by = "dataset_name", all.x = TRUE)
dataset_summary$species_group[is.na(dataset_summary$species_group)] <-
  ifelse(grepl("Mouse|mouse", dataset_summary$dataset_name[is.na(dataset_summary$species_group)]),
    "Murine", "Human")
dataset_summary$tissue[is.na(dataset_summary$tissue)] <-
  dataset_summary$dataset_name[is.na(dataset_summary$tissue)]

manifest_path <- file.path(root, "data", "spatial_datasets", "dataset_manifest.csv")
if (file.exists(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  names(manifest)[names(manifest) == "title"] <- "dataset_name"
  dataset_summary <- merge(dataset_summary, manifest, by = "dataset_name", all.x = TRUE)
}
dataset_summary$repository <- dataset_provenance_map$repository[
  match(dataset_summary$dataset_name, dataset_provenance_map$dataset_name)
]
dataset_summary$accession <- dataset_provenance_map$accession[
  match(dataset_summary$dataset_name, dataset_provenance_map$dataset_name)
]
dataset_summary <- dataset_summary[order(-dataset_summary$pct_nonstationary), ]
write.csv(dataset_summary, file.path(add_dir, "Additional_file_2_dataset_manifest_and_analysis_summary.csv"), row.names = FALSE)

build_mesh_summary <- function(dataset_summary) {
  data_dir <- file.path(root, "data", "spatial_datasets")
  rows <- lapply(seq_len(nrow(dataset_summary)), function(i) {
    ds <- dataset_summary$dataset_name[i]
    ds_file <- if ("dataset_file" %in% names(dataset_summary)) dataset_summary$dataset_file[i] else NA_character_
    candidates <- unique(file.path(data_dir, c(ds_file, paste0(ds, ".rds"))))
    candidates <- candidates[!is.na(candidates) & file.exists(candidates)]

    out <- data.frame(
      dataset_name = ds,
      n_spots = dataset_summary$n_spots[i],
      x_range = NA_real_,
      y_range = NA_real_,
      max_edge_inner = NA_real_,
      max_edge_outer = NA_real_,
      cutoff = NA_real_,
      mesh_nodes = NA_integer_,
      stringsAsFactors = FALSE
    )
    if (!length(candidates) || !requireNamespace("SpatialExperiment", quietly = TRUE)) return(out)

    spe <- tryCatch(readRDS(candidates[1]), error = function(e) NULL)
    if (is.null(spe)) return(out)
    coords <- tryCatch(as.data.frame(SpatialExperiment::spatialCoords(spe)), error = function(e) NULL)
    if (is.null(coords) || ncol(coords) < 2) return(out)
    colnames(coords)[1:2] <- c("x", "y")

    x_range <- diff(range(coords$x, na.rm = TRUE))
    y_range <- diff(range(coords$y, na.rm = TRUE))
    max_range <- max(x_range, y_range)
    min_range <- min(x_range, y_range)

    out$x_range <- x_range
    out$y_range <- y_range
    out$max_edge_inner <- max_range / 5
    out$max_edge_outer <- max_range / 3
    out$cutoff <- min_range / 20

    if (requireNamespace("INLA", quietly = TRUE)) {
      mesh <- tryCatch(
        INLA::inla.mesh.2d(
          loc = as.matrix(coords[, c("x", "y")]),
          max.edge = c(out$max_edge_inner, out$max_edge_outer),
          cutoff = out$cutoff
        ),
        error = function(e) NULL
      )
      if (!is.null(mesh)) out$mesh_nodes <- mesh$n
    }
    out
  })
  do.call(rbind, rows)
}

mesh_summary <- build_mesh_summary(dataset_summary)
write.csv(mesh_summary, file.path(tab_dir, "inla_mesh_parameters.csv"), row.names = FALSE)

strength_levels <- c("Extreme", "Very strong", "Strong", "Moderate", "Anecdotal")
valid$bf_strength_nonstationary_vs_stationary <- factor(valid$bf_strength_nonstationary_vs_stationary, levels = strength_levels)
evidence_summary <- aggregate(
  gene_name ~ dataset_name + model_favored_nonstationary_vs_stationary + bf_strength_nonstationary_vs_stationary,
  data = valid,
  FUN = length
)
names(evidence_summary)[names(evidence_summary) == "gene_name"] <- "n_genes"
totals <- aggregate(gene_name ~ dataset_name, data = valid, FUN = length)
names(totals)[2] <- "dataset_total_genes"
evidence_summary <- merge(evidence_summary, totals, by = "dataset_name")
evidence_summary$pct_of_dataset <- 100 * evidence_summary$n_genes / evidence_summary$dataset_total_genes
evidence_summary <- evidence_summary[order(evidence_summary$dataset_name,
  evidence_summary$model_favored_nonstationary_vs_stationary,
  evidence_summary$bf_strength_nonstationary_vs_stationary), ]
write.csv(evidence_summary, file.path(add_dir, "Additional_file_4_bayes_factor_evidence_summary.csv"), row.names = FALSE)

stationary <- valid[valid$log_bayes_factor < 0 & !is.na(valid$log_bf_stationary_vs_nonspatial), , drop = FALSE]
if (nrow(stationary)) {
  stationary$model_favored_stationary_vs_nonspatial <- ifelse(
    stationary$log_bf_stationary_vs_nonspatial > 0,
    "stationary_spatial", "stationary_nonspatial"
  )
  stationary_summary <- aggregate(
    gene_name ~ dataset_name + model_favored_stationary_vs_nonspatial + bf_strength_stationary_vs_nonspatial,
    data = stationary,
    FUN = length
  )
  names(stationary_summary)[names(stationary_summary) == "gene_name"] <- "n_genes"
  st_totals <- aggregate(gene_name ~ dataset_name, data = stationary, FUN = length)
  names(st_totals)[2] <- "stationary_favored_with_nonspatial_fit"
  stationary_summary <- merge(stationary_summary, st_totals, by = "dataset_name")
  stationary_summary$pct_of_stationary_favored <- 100 *
    stationary_summary$n_genes / stationary_summary$stationary_favored_with_nonspatial_fit
  stationary_summary <- stationary_summary[order(stationary_summary$dataset_name,
    stationary_summary$model_favored_stationary_vs_nonspatial,
    stationary_summary$bf_strength_stationary_vs_nonspatial), ]
} else {
  stationary_summary <- data.frame()
}
write.csv(stationary_summary, file.path(add_dir, "Additional_file_5_stationary_vs_nonspatial_summary.csv"), row.names = FALSE)

copy_study_csv_if_exists <- function(from, to) {
  if (file.exists(from)) {
    x <- read.csv(from, stringsAsFactors = FALSE)
    dataset_column <- intersect(c("dataset_name", "dataset"), names(x))
    if (length(dataset_column)) {
      x <- x[x[[dataset_column[1]]] %in% analysis_datasets, , drop = FALSE]
    }
    write.csv(x, to, row.names = FALSE)
    TRUE
  } else {
    FALSE
  }
}
copy_study_csv_if_exists(file.path(root, "results", "gsea_clusterprofiler.csv"),
  file.path(add_dir, "Additional_file_6_clusterprofiler_gsea_results.csv"))
copy_study_csv_if_exists(file.path(root, "results", "gsea_stationarity_enrichment.csv"),
  file.path(add_dir, "Additional_file_7_msigdb_fgsea_results.csv"))

curve_files <- file.path(root, "simulation", "output",
  c("HumanBreastCancerILC_fdr_power_curves.csv", "HumanOvarianCancer_fdr_power_curves.csv"))
existing_curve_files <- curve_files[file.exists(curve_files)]
curves <- if (length(existing_curve_files)) {
  do.call(rbind, lapply(existing_curve_files, read.csv, stringsAsFactors = FALSE))
} else {
  data.frame()
}
write.csv(curves, file.path(add_dir, "Additional_file_8_simulation_fdr_power_curves.csv"), row.names = FALSE)

rep_files <- file.path(root, "simulation", "output",
  c("HumanBreastCancerILC_replicate_summaries.csv", "HumanOvarianCancer_replicate_summaries.csv"))
existing_rep_files <- rep_files[file.exists(rep_files)]
rep_summaries <- if (length(existing_rep_files)) {
  do.call(rbind, lapply(existing_rep_files, function(f) {
    x <- read.csv(f, stringsAsFactors = FALSE)
    x$dataset_name <- sub("_replicate_summaries[.]csv$", "", basename(f))
    x
  }))
} else {
  data.frame()
}
write.csv(rep_summaries, file.path(add_dir, "Additional_file_9_simulation_replicate_summaries.csv"), row.names = FALSE)

script_manifest <- data.frame(
  path = c(
    "hpc/download_datasets.R",
    "hpc/run_analysis.R",
    "hpc/analysis_functions.R",
    "hpc/analysis_datasets.txt",
    "hpc/combine_results.R",
    "hpc/install_inla.sh",
    "hpc/preflight_check.sh",
    "hpc/submit_download.sbatch",
    "hpc/submit_analysis.sbatch",
    "hpc/submit_combine.sbatch",
    "hpc/test_result_handling.R",
    "supplementary/scripts/create_supplementary_materials.R",
    "results/plot_diverging_stacked_bar.R",
    "results/plot_nonstationarity_by_sample_type.R",
    "results/plot_stationary_vs_nonspatial.R",
    "results/plot_spatial_maps.R",
    "results/plot_spatial_maps_remaining.R",
    "results/gsea_clusterprofiler.R",
    "results/gsea_stationarity_enrichment.R",
    "simulation/prepare_simulation.R",
    "simulation/run_replicate.R",
    "simulation/aggregate_results.R",
    "simulation/plot_fdr_power.R"
  ),
  purpose = c(
    "Download SpatialExperiment datasets from ExperimentHub.",
    "Run per-gene stationary and nonstationary INLA model fits.",
    "Model-fitting, preprocessing, mesh, Bayes factor, and helper functions.",
    "Define the stable set and SLURM order of study datasets.",
    "Validate and combine completed per-dataset results.",
    "Install INLA and its cluster dependencies.",
    "Check cluster dependencies, syntax, datasets, and array mapping.",
    "Submit the ExperimentHub download job.",
    "Submit the per-dataset analysis array.",
    "Submit validated post-array result combination.",
    "Test failed-fit retention and result completeness checks.",
    "Generate supplementary tables and machine-readable files.",
    "Plot evidence strength for stationary versus nonstationary spatial covariance.",
    "Plot dataset-level proportion of nonstationary genes by sample type.",
    "Plot stationary spatial versus stationary non-spatial model evidence.",
    "Generate representative spatial residual maps for selected human and murine tissues.",
    "Generate representative spatial residual maps for remaining tissues.",
    "Run clusterProfiler GO biological process and KEGG GSEA.",
    "Run MSigDB/fgsea stationarity enrichment analysis.",
    "Create simulation configurations from fitted model parameters.",
    "Run one simulation replicate.",
    "Aggregate simulation replicate outputs into FDR and power curves.",
    "Plot simulation FDR and power curves."
  ),
  stringsAsFactors = FALSE
)
script_manifest$exists <- file.exists(file.path(root, script_manifest$path))
script_manifest$file_size_bytes <- ifelse(script_manifest$exists,
  file.info(file.path(root, script_manifest$path))$size, NA)
script_manifest$md5 <- ifelse(script_manifest$exists,
  as.character(tools::md5sum(file.path(root, script_manifest$path))), NA)
write.csv(script_manifest, file.path(add_dir, "Additional_file_10_software_reproducibility_manifest.csv"), row.names = FALSE)

additional_files <- data.frame(
  file_name = c(
    "Additional_file_1_Supplementary_Information.pdf",
    "Additional_file_2_dataset_manifest_and_analysis_summary.csv",
    "Additional_file_3_per_gene_bayes_factor_and_category_summary.csv",
    "Additional_file_4_bayes_factor_evidence_summary.csv",
    "Additional_file_5_stationary_vs_nonspatial_summary.csv",
    "Additional_file_6_clusterprofiler_gsea_results.csv",
    "Additional_file_7_msigdb_fgsea_results.csv",
    "Additional_file_8_simulation_fdr_power_curves.csv",
    "Additional_file_9_simulation_replicate_summaries.csv",
    "Additional_file_10_software_reproducibility_manifest.csv"
  ),
  format = c("PDF", rep("CSV", 9)),
  title = c(
    "Supplementary methods, figures, and summary tables",
    "Dataset manifest and analysis summary",
    "Per-gene Bayes factor and category summary",
    "Evidence-strength counts for nonstationary versus stationary models",
    "Stationary spatial versus stationary non-spatial summary",
    "clusterProfiler GSEA results",
    "MSigDB fgsea stationarity enrichment results",
    "Simulation FDR and power curves",
    "Simulation replicate summaries",
    "Software and reproducibility manifest"
  ),
  description = c(
    "Supplementary methods, supplemental figure legends, concise summary tables, and Genome Biology-oriented reproducibility notes.",
    "One row per analyzed dataset, including sample labels, Space Ranger version/source URL, number of spots and genes, ExperimentHub metadata when available, percentage of genes favoring nonstationary covariance, and stationary/non-spatial summaries.",
    "One row per analyzed gene and dataset, including log Bayes factors, Bayes factor evidence strengths, model favored by the nonstationary-versus-stationary comparison, stationary-versus-non-spatial comparison when available, and the final three-model category.",
    "Counts and percentages of genes by dataset, favored model, and Bayes factor evidence strength for the nonstationary versus stationary comparison.",
    "Counts and percentages of stationary-favored genes that favor either the stationary spatial model or the stationary non-spatial model.",
    "Full GO biological process and KEGG GSEA output from clusterProfiler, with genes ranked by log Bayes factor.",
    "Full Hallmark, KEGG, Reactome, and GO biological process fgsea output from MSigDB gene sets, with genes ranked by log Bayes factor.",
    "Aggregated threshold-specific false discovery rate and average power for the two simulation settings.",
    "Per-replicate FDR and power at the Bayes factor threshold c = 1 for the two simulation settings.",
    "List of analysis scripts, their role in the workflow, file sizes, and MD5 checksums."
  ),
  stringsAsFactors = FALSE
)
write.csv(additional_files, file.path(out_dir, "additional_file_metadata.csv"), row.names = FALSE)

ds_tex <- dataset_summary
ds_tex$repository <- dataset_provenance_map$repository[
  match(ds_tex$dataset_name, dataset_provenance_map$dataset_name)
]
ds_tex$accession <- dataset_provenance_map$accession[
  match(ds_tex$dataset_name, dataset_provenance_map$dataset_name)
]
ds_tex <- ds_tex[, c("dataset_name", "repository", "accession", "species_group", "tissue",
  "n_spots", "n_genes_analyzed", "space_ranger_version", "pct_nonstationary")]
names(ds_tex) <- c("Dataset", "Repository", "Accession", "Species", "Tissue", "Spots",
  "Genes", "Space Ranger", "% NS")
write_simple_table(ds_tex, file.path(tab_dir, "dataset_summary_table.tex"),
  paste(
    "Analyzed datasets, Bioconductor repository resources, ExperimentHub accessions,",
    "Space Ranger versions, and proportion of genes favoring the nonstationary spatial model.",
    "EH9628 and EH9629 are the current replacements for the accessions used during analysis,",
    "EH9516 and EH9517, respectively."
  ),
  "tab:dataset-summary", digits = 1, font_size = "scriptsize", resize = TRUE)

mesh_tex <- mesh_summary[, c("dataset_name", "n_spots", "x_range", "y_range",
  "max_edge_inner", "max_edge_outer", "cutoff", "mesh_nodes")]
names(mesh_tex) <- c("Dataset", "Spots", "x range", "y range",
  "max.edge[1]", "max.edge[2]", "cutoff", "Mesh nodes")
write_simple_table(mesh_tex, file.path(tab_dir, "inla_mesh_parameters_table.tex"),
  "INLA SPDE mesh settings used for each dataset. The mesh was generated with max.edge = c(max range / 5, max range / 3) and cutoff = min range / 20.",
  "tab:mesh-parameters", digits = 1)

if (nrow(curves) > 0L) {
  key <- curves[curves$c_threshold %in% c(1, 3, 10, 30, 100), ]
  key <- key[order(key$dataset_name, key$c_threshold), c("dataset_name", "c_threshold", "fdr", "avg_power", "n_replicates", "n_genes", "n_h0", "n_h1")]
  names(key) <- c("Dataset", "c", "FDR", "Power", "Replicates", "Genes", "H0 genes", "H1 genes")
  write_simple_table(key, file.path(tab_dir, "simulation_threshold_table.tex"),
    "Simulation false discovery rate and average power at representative Bayes factor thresholds.",
    "tab:simulation-thresholds", digits = 3)
}

gsea_path <- file.path(root, "results", "gsea_clusterprofiler.csv")
if (file.exists(gsea_path)) {
  gsea <- read.csv(gsea_path, stringsAsFactors = FALSE)
  dataset_column <- intersect(c("dataset_name", "dataset"), names(gsea))
  if (length(dataset_column)) {
    gsea <- gsea[gsea[[dataset_column[1]]] %in% analysis_datasets, , drop = FALSE]
  }
  gsea <- gsea[!is.na(gsea$padj), ]
  gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
  gsea_top <- head(gsea[gsea$padj < 0.05, ], 20)
  if (nrow(gsea_top)) {
    gsea_tex <- gsea_top[, intersect(c("dataset", "source", "pathway", "Description", "NES", "padj"), names(gsea_top))]
    names(gsea_tex) <- c("Dataset", "Source", "ID", "Description", "NES", "FDR")
    write_simple_table(gsea_tex, file.path(tab_dir, "top_gsea_table.tex"),
      "Top significant gene-set enrichment results from clusterProfiler.",
      "tab:top-gsea", digits = 3)
  }

  sig_gsea <- gsea[gsea$padj < 0.05 & gsea$source %in% c("GO_BP", "KEGG"), ]
  if (nrow(sig_gsea)) {
    sig_gsea <- sig_gsea[order(sig_gsea$source, sig_gsea$padj, -abs(sig_gsea$NES)), ]
    sig_gsea <- sig_gsea[!duplicated(sig_gsea$pathway), ]
    code_rows <- do.call(rbind, lapply(split(sig_gsea, sig_gsea$source), function(x) head(x, 12)))
    code_rows <- code_rows[order(code_rows$source, code_rows$padj), ]
    code_tex <- code_rows[, intersect(c("source", "pathway", "Description", "dataset", "NES", "padj"), names(code_rows))]
    names(code_tex) <- c("Source", "Identifier", "Full term or pathway name", "Example dataset", "NES", "FDR")

    con <- file(file.path(tab_dir, "go_kegg_identifier_expansions_table.tex"), open = "wt")
    writeLines("\\begin{table}[htbp]\n\\centering\n\\small", con)
    writeLines("\\caption{GO and KEGG identifier expansions for representative significant gene-set enrichment results. The complete identifier-to-description mapping is provided in Additional file 6.}", con)
    writeLines("\\label{tab:go-kegg-expansions}", con)
    writeLines("\\begin{tabular}{llp{0.38\\linewidth}lrr}\n\\toprule", con)
    writeLines("Source & Identifier & Full term or pathway name & Example dataset & NES & FDR \\\\\n\\midrule", con)
    for (i in seq_len(nrow(code_tex))) {
      vals <- c(
        latex_escape(code_tex$Source[i]),
        latex_escape(code_tex$Identifier[i]),
        latex_escape(code_tex$`Full term or pathway name`[i]),
        latex_escape(code_tex$`Example dataset`[i]),
        format(round(code_tex$NES[i], 3), nsmall = 3, trim = TRUE),
        format(signif(code_tex$FDR[i], 3), trim = TRUE, scientific = TRUE)
      )
      writeLines(paste(vals, collapse = " & "), con)
      writeLines(" \\\\", con)
    }
    writeLines("\\bottomrule\n\\end{tabular}\n\\end{table}", con)
    close(con)
  }
}

desc_entries <- paste0("\\textbf{Additional file ", seq_len(nrow(additional_files)), ".} ",
  additional_files$title, ". Format: ", additional_files$format, ". ",
  additional_files$description)
desc_lines <- c("\\section*{Additional files}", as.vector(rbind(desc_entries, "")))
desc_lines <- desc_lines[-length(desc_lines)]
writeLines(desc_lines, file.path(out_dir, "additional_file_descriptions_for_manuscript.tex"))

availability <- c(
  "\\section*{Availability of data and materials}",
  "All public spatial transcriptomics datasets analyzed in this study are listed in Additional file 2. Per-gene model-fitting results, enrichment outputs, simulation summaries, and reproducibility metadata are provided in Additional files 3--10. Source code for the analyses is available at \\url{https://github.com/nathoogroup/SpaTran-NS}; scripts used to generate the reported results are enumerated with checksums in Additional file 10. Reviewer-accessible or accession-specific links for any controlled or third-party datasets should be added here before submission.",
  "",
  "\\section*{Code availability}",
  "Analysis code is available at \\url{https://github.com/nathoogroup/SpaTran-NS}. The supplementary reproducibility manifest in Additional file 10 lists the principal scripts and MD5 checksums."
)
writeLines(availability, file.path(out_dir, "availability_and_code_declarations_draft.tex"))

checklist <- c(
  "# Genome Biology Supplement Readiness Checklist",
  "",
  "- [x] Supplementary Information file has a table of contents and mirrors the main analysis structure.",
  "- [x] Per-gene Bayes factor and model-category results are supplied as machine-readable CSV.",
  "- [x] Dataset-level summary and public dataset metadata are supplied as machine-readable CSV.",
  "- [x] Enrichment and simulation outputs are supplied as machine-readable CSV.",
  "- [x] Principal analysis scripts are listed with file size and MD5 checksum.",
  "- [ ] Add final persistent data accessions or reviewer-private links for every third-party dataset before submission.",
  "- [ ] Add ethics approval/waiver details for any human tissue, human data, or controlled-access material.",
  "- [ ] Confirm the repository URL, release tag, and license before submission.",
  "- [ ] Confirm each additional file remains below the journal upload limit after final PDF export.",
  "- [ ] Replace manuscript placeholder figures/captions and cite all additional files in sequence."
)
writeLines(checklist, file.path(out_dir, "Genome_Biology_submission_checklist.md"))

cat("Wrote supplementary files to ", out_dir, "\n", sep = "")
