#!/usr/bin/env Rscript

script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_flag) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]])))
} else {
  normalizePath("hpc")
}
source(file.path(script_dir, "gene_chunks.R"))

n_genes <- 101L
chunk_count <- 16L
chunks <- lapply(
  seq_len(chunk_count),
  function(chunk_id) gene_chunk_indices(n_genes, chunk_id, chunk_count)
)

stopifnot(
  identical(chunks[[1]], as.integer(c(1, 17, 33, 49, 65, 81, 97))),
  identical(chunks[[16]], as.integer(c(16, 32, 48, 64, 80, 96))),
  !anyDuplicated(unlist(chunks, use.names = FALSE)),
  identical(sort(unlist(chunks, use.names = FALSE)), seq_len(n_genes)),
  identical(
    basename(gene_chunk_result_path("chunks", "Dataset", 2L, 16L)),
    "Dataset_chunk_002_of_016_results.rds"
  )
)

invalid <- try(gene_chunk_indices(4L, 5L, 5L), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))

cat("Gene-chunk tests passed.\n")
