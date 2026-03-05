#!/usr/bin/env Rscript
# Preliminary analysis: what % of genes show evidence for non-stationarity?
# Evidence = positive log Bayes factor (BF favors non-stationary over stationary model)
# Loads all individual *_results.rds files (combined file may be incomplete)

suppressPackageStartupMessages(library(dplyr))

# Look in cwd and results/ subdir
rds_files <- c(
  list.files(".", pattern = "_results\\.rds$", full.names = TRUE),
  list.files("results", pattern = "_results\\.rds$", full.names = TRUE)
)
rds_files <- unique(rds_files)
rds_files <- rds_files[!grepl("all_results_combined", rds_files)]

if (length(rds_files) == 0) {
  stop("No *_results.rds files found in . or results/")
}

cat("Loading", length(rds_files), "result files...\n")
df_list <- lapply(rds_files, function(f) {
  x <- readRDS(f)
  if (!"dataset_name" %in% names(x)) {
    x$dataset_name <- sub("_results\\.rds$", "", basename(f))
  }
  x
})
df <- do.call(rbind, df_list)

# Exclude genes with NA log_bayes_factor (failed fits)
df_valid <- df %>% filter(!is.na(log_bayes_factor))

n_total <- nrow(df_valid)
n_failed <- nrow(df) - n_total

# --- Evidence for non-stationarity ---
# log_bf > 0  => BF > 1 => non-stationary model preferred
n_nonstationary <- sum(df_valid$log_bayes_factor > 0)
pct_nonstationary <- 100 * n_nonstationary / n_total

# By evidence strength (Jeffreys scale)
# Anecdotal: |log(BF)| < log(3)
# Moderate:  log(3) <= |log(BF)| < log(10)
# Strong:    log(10) <= |log(BF)| < log(30)
# Very strong: log(30) <= |log(BF)| < log(100)
# Extreme:   |log(BF)| >= log(100)

df_valid <- df_valid %>%
  mutate(
    favors_nonstationary = log_bayes_factor > 0,
    strength = case_when(
      abs(log_bayes_factor) < log(3)   ~ "Anecdotal",
      abs(log_bayes_factor) < log(10)  ~ "Moderate",
      abs(log_bayes_factor) < log(30)  ~ "Strong",
      abs(log_bayes_factor) < log(100) ~ "Very strong",
      TRUE                             ~ "Extreme"
    )
  )

# Non-stationary by strength
ns_by_strength <- df_valid %>%
  filter(favors_nonstationary) %>%
  count(strength) %>%
  arrange(match(strength, c("Anecdotal", "Moderate", "Strong", "Very strong", "Extreme")))

# At least moderate evidence (exclude anecdotal)
n_at_least_moderate <- sum(df_valid$favors_nonstationary & 
  df_valid$strength %in% c("Moderate", "Strong", "Very strong", "Extreme"))
pct_at_least_moderate <- 100 * n_at_least_moderate / n_total

# At least strong evidence
n_at_least_strong <- sum(df_valid$favors_nonstationary & 
  df_valid$strength %in% c("Strong", "Very strong", "Extreme"))
pct_at_least_strong <- 100 * n_at_least_strong / n_total

# --- Print report ---
cat("\n")
cat("=============================================================================\n")
cat("PRELIMINARY ANALYSIS: Evidence for Non-Stationarity\n")
cat("=============================================================================\n\n")

cat(sprintf("Total genes analyzed: %d\n", n_total))
if (n_failed > 0) cat(sprintf("Genes with failed fits (excluded): %d\n\n", n_failed))

cat("--- Overall ---\n")
cat(sprintf("Genes favoring NON-STATIONARY (log BF > 0): %d (%.2f%%)\n",
    n_nonstationary, pct_nonstationary))
cat(sprintf("Genes favoring stationary (log BF <= 0):    %d (%.2f%%)\n\n",
    n_total - n_nonstationary, 100 - pct_nonstationary))

cat("--- By evidence strength (non-stationary only) ---\n")
print(ns_by_strength)
cat("\n")

cat("--- Stringent thresholds ---\n")
cat(sprintf("At least MODERATE evidence for non-stationarity: %d (%.2f%%)\n",
    n_at_least_moderate, pct_at_least_moderate))
cat(sprintf("At least STRONG evidence for non-stationarity:   %d (%.2f%%)\n\n",
    n_at_least_strong, pct_at_least_strong))

cat("--- By dataset (summary) ---\n")
by_dataset <- df_valid %>%
  group_by(dataset_name) %>%
  summarise(
    n_genes = n(),
    n_nonstationary = sum(favors_nonstationary),
    pct_nonstationary = 100 * mean(favors_nonstationary),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_nonstationary))
print(by_dataset, n = Inf)

cat("\n")
cat("=============================================================================\n")
cat("INDIVIDUAL DATASET RESULTS\n")
cat("=============================================================================\n")

strength_order <- c("Anecdotal", "Moderate", "Strong", "Very strong", "Extreme")
datasets <- unique(df_valid$dataset_name)
by_dataset_sorted <- by_dataset %>% arrange(desc(pct_nonstationary))
datasets <- by_dataset_sorted$dataset_name

for (ds in datasets) {
  d <- df_valid %>% filter(dataset_name == ds)
  n <- nrow(d)
  n_ns <- sum(d$favors_nonstationary)
  pct <- 100 * n_ns / n
  n_mod <- sum(d$favors_nonstationary & d$strength %in% c("Moderate", "Strong", "Very strong", "Extreme"))
  n_str <- sum(d$favors_nonstationary & d$strength %in% c("Strong", "Very strong", "Extreme"))

  cat("\n")
  cat("--- ", ds, " ---\n", sep = "")
  cat(sprintf("  Genes: %d  |  Non-stationary: %d (%.2f%%)\n", n, n_ns, pct))
  cat(sprintf("  At least moderate: %d (%.2f%%)  |  At least strong: %d (%.2f%%)\n",
      n_mod, 100 * n_mod / n, n_str, 100 * n_str / n))

  ns_strength <- d %>%
    filter(favors_nonstationary) %>%
    count(strength) %>%
    arrange(match(strength, strength_order))
  cat("  By strength: ")
  cat(paste0(ns_strength$strength, "=", ns_strength$n, collapse = " | "))
  cat("\n")
}

cat("\n=============================================================================\n")
cat("Done.\n")
