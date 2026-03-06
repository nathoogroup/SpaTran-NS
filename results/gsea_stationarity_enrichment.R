#!/usr/bin/env Rscript
# =============================================================================
# GSEA enrichment: Gene sets associated with stationary vs non-stationary genes
# =============================================================================
# Uses MSigDB gene sets (Hallmark, KEGG) and fgsea to find pathways enriched
# in genes favoring non-stationary (log BF > 0) vs stationary (log BF < 0) models.
#
# Genes are ranked by log_bayes_factor; positive = non-stationary.
# GSEA finds gene sets enriched at top (non-stationary) or bottom (stationary).
#
# Requires: msigdbr, fgsea, dplyr
#   install.packages(c("msigdbr", "fgsea", "dplyr"))
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
  library(fgsea)
})

# ---- Load results (same as preliminary analysis) ----
rds_files <- c(
  list.files(".", pattern = "_results\\.rds$", full.names = TRUE),
  list.files("results", pattern = "_results\\.rds$", full.names = TRUE)
)
rds_files <- unique(rds_files)
rds_files <- rds_files[!grepl("all_results_combined", rds_files)]

if (length(rds_files) == 0) stop("No *_results.rds files found")

cat("Loading", length(rds_files), "result files...\n")
df_list <- lapply(rds_files, function(f) {
  x <- readRDS(f)
  if (!"dataset_name" %in% names(x)) {
    x$dataset_name <- sub("_results\\.rds$", "", basename(f))
  }
  x
})
df <- do.call(rbind, df_list)
df <- df %>% filter(!is.na(log_bayes_factor))

# ---- Species mapping ----
species_map <- c(
  HumanBreastCancerILC = "Homo sapiens", HumanBreastCancerIDC = "Homo sapiens",
  HumanColorectalCancer = "Homo sapiens", HumanGlioblastoma = "Homo sapiens",
  HumanHeart = "Homo sapiens", HumanLymphNode = "Homo sapiens",
  HumanOvarianCancer = "Homo sapiens", HumanSpinalCord = "Homo sapiens",
  HumanCerebellum = "Homo sapiens",
  MouseBrainCoronal = "Mus musculus", MouseBrainSagittalAnterior = "Mus musculus",
  MouseBrainSagittalPosterior = "Mus musculus", MouseKidneyCoronal = "Mus musculus",
  Visium_humanDLPFC = "Homo sapiens", Visium_mouseCoronal = "Mus musculus",
  WeberDivechaLCdata_Visium = "Homo sapiens", Janesick_breastCancer_Visium = "Homo sapiens"
)

# ---- Prepare gene sets ----
# Hallmark (H) + KEGG + Reactome + GO Biological Process (larger sets)
get_gene_sets <- function(species) {
  hallmark <- msigdbr(species = species, collection = "H") %>%
    select(gs_name, gene_symbol)
  kegg <- msigdbr(species = species, collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
    select(gs_name, gene_symbol)
  reactome <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
    select(gs_name, gene_symbol)
  go_bp <- msigdbr(species = species, collection = "C5", subcollection = "GO:BP") %>%
    select(gs_name, gene_symbol)
  gs <- bind_rows(hallmark, kegg, reactome, go_bp)
  split(gs$gene_symbol, gs$gs_name)
}

# ---- Run GSEA per dataset ----
run_gsea_dataset <- function(d, dataset_name, gene_sets, min_size = 15, max_size = 2000) {
  # Rank by log_bayes_factor (positive = non-stationary)
  ranks <- setNames(d$log_bayes_factor, d$gene_name)
  ranks <- ranks[!is.na(ranks) & !duplicated(names(ranks))]

  if (length(ranks) < 50) return(NULL)

  gs_filtered <- gene_sets[sapply(gene_sets, function(x) {
    n <- sum(unique(x) %in% names(ranks))
    n >= min_size && n <= max_size
  })]

  if (length(gs_filtered) == 0) return(NULL)

  set.seed(42)
  fgsea_res <- fgsea(
    pathways  = gs_filtered,
    stats     = ranks,
    minSize   = min_size,
    maxSize   = max_size,  # allow larger gene sets (Reactome/GO BP can be 500-2000 genes)
    nperm     = 10000,
    scoreType = "std"
  )

  fgsea_res %>%
    as.data.frame() %>%
    mutate(dataset = dataset_name) %>%
    arrange(padj)
}

# ---- Main ----
cat("\n=============================================================================\n")
cat("GSEA: Gene Set Enrichment for Stationary vs Non-Stationary Genes\n")
cat("=============================================================================\n\n")

datasets <- unique(df$dataset_name)
all_res <- list()

for (ds in datasets) {
  species <- species_map[ds]
  if (is.na(species)) species <- "Homo sapiens"  # default
  cat("  ", ds, " (", species, ") ... ", sep = "")

  d <- df %>% filter(dataset_name == ds)
  gene_sets <- get_gene_sets(species)

  res <- run_gsea_dataset(d, ds, gene_sets)
  if (is.null(res)) {
    cat("skipped (too few genes or gene sets)\n")
    next
  }

  all_res[[ds]] <- res
  n_sig <- sum(res$padj < 0.05)
  cat(n_sig, "gene sets padj < 0.05\n")
}

if (length(all_res) == 0) stop("No datasets produced GSEA results")

# ---- Report: Top gene sets for non-stationary (positive NES) ----
cat("\n--- Top gene sets enriched in NON-STATIONARY genes (positive NES) ---\n\n")

for (ds in names(all_res)) {
  r <- all_res[[ds]] %>% filter(NES > 0, padj < 0.05) %>% head(10)
  if (nrow(r) == 0) next
  cat("## ", ds, "\n", sep = "")
  print(r %>% select(pathway, NES, pval, padj, size) %>% as.data.frame(), row.names = FALSE)
  cat("\n")
}

# ---- Report: Top gene sets for stationary (negative NES) ----
cat("\n--- Top gene sets enriched in STATIONARY genes (negative NES) ---\n\n")

for (ds in names(all_res)) {
  r <- all_res[[ds]] %>% filter(NES < 0, padj < 0.05) %>% head(10)
  if (nrow(r) == 0) next
  cat("## ", ds, "\n", sep = "")
  print(r %>% select(pathway, NES, pval, padj, size) %>% as.data.frame(), row.names = FALSE)
  cat("\n")
}

# ---- Save full results ----
gsea_combined <- do.call(rbind, all_res)
out_dir <- if (dir.exists("results")) "results" else "."
out_csv <- file.path(out_dir, "gsea_stationarity_enrichment.csv")
out_rds <- file.path(out_dir, "gsea_stationarity_enrichment.rds")
# Remove leadingEdge (list column) for CSV
gsea_csv <- gsea_combined
if ("leadingEdge" %in% names(gsea_csv)) gsea_csv$leadingEdge <- NULL
write.csv(gsea_csv, out_csv, row.names = FALSE)
saveRDS(gsea_combined, out_rds)
cat("Saved:", out_csv, "\n")
cat("Saved:", out_rds, "\n")

cat("\nDone.\n")
