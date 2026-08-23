# Deterministic round-robin gene partitioning for SLURM array tasks.

positive_integer <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  parsed
}

gene_chunk_indices <- function(n_genes, chunk_id, chunk_count) {
  n_genes <- positive_integer(n_genes, "n_genes")
  chunk_id <- positive_integer(chunk_id, "chunk_id")
  chunk_count <- positive_integer(chunk_count, "chunk_count")
  if (chunk_id > chunk_count) {
    stop("chunk_id cannot exceed chunk_count", call. = FALSE)
  }

  indices <- seq.int(from = chunk_id, to = n_genes, by = chunk_count)
  if (length(indices) == 0L) {
    stop(
      sprintf(
        "gene chunk %d/%d is empty for %d filtered genes",
        chunk_id, chunk_count, n_genes
      ),
      call. = FALSE
    )
  }
  as.integer(indices)
}

gene_chunk_result_path <- function(chunk_dir, dataset_name, chunk_id,
                                   chunk_count) {
  chunk_id <- positive_integer(chunk_id, "chunk_id")
  chunk_count <- positive_integer(chunk_count, "chunk_count")
  file.path(
    chunk_dir,
    sprintf(
      "%s_chunk_%03d_of_%03d_results.rds",
      dataset_name, chunk_id, chunk_count
    )
  )
}
