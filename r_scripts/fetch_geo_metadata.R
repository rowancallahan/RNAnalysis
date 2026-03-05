#!/usr/bin/env Rscript
# Fetch sample metadata from GEO using GEOquery.
# Saves metadata.csv (and counts.csv if expression data available).
# Outputs JSON to stdout with file paths.
#
# Usage: Rscript fetch_geo_metadata.R <gse_id> <library_path> <output_dir>

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat('{"error": "Usage: Rscript fetch_geo_metadata.R <gse_id> <library_path> <output_dir>"}\n')
    quit(status = 0)
  }

  gse_id <- args[1]
  lib_path <- args[2]
  output_dir <- args[3]

  .libPaths(c(lib_path, .libPaths()))

  library(GEOquery, quietly = TRUE)
  library(jsonlite, quietly = TRUE)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  message("Fetching GEO data for: ", gse_id)

  # Download GEO series matrix (faster than full SOFT)
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = output_dir)

  if (is.list(gse) && length(gse) == 0) {
    cat(toJSON(list(error = paste0("No data found for ", gse_id)), auto_unbox = TRUE))
    cat("\n")
    quit(status = 0)
  }

  # If multiple platforms, use the first
  if (is.list(gse)) {
    eset <- gse[[1]]
  } else {
    eset <- gse
  }

  # Extract sample metadata (phenoData)
  meta_df <- pData(eset)

  # Clean up metadata — remove overly verbose/redundant columns
  drop_patterns <- c("contact", "supplementary_file", "relation",
                     "data_row_count", "status", "submission_date",
                     "last_update_date", "channel_count",
                     "extract_protocol", "treatment_protocol",
                     "growth_protocol", "data_processing",
                     "taxid_ch", "molecule_ch", "label_ch",
                     "label_protocol", "hyb_protocol", "scan_protocol",
                     "description")
  for (pat in drop_patterns) {
    drop_cols <- grep(pat, colnames(meta_df), ignore.case = TRUE, value = TRUE)
    meta_df <- meta_df[, !colnames(meta_df) %in% drop_cols, drop = FALSE]
  }

  # If parsed characteristic columns exist (e.g. "agent:ch1"), drop the raw
  # "characteristics_ch1", "characteristics_ch1.1" etc. to avoid duplication
  parsed_chars <- grep(":ch1$", colnames(meta_df), value = TRUE)
  if (length(parsed_chars) > 0) {
    raw_chars <- grep("^characteristics_ch", colnames(meta_df), value = TRUE)
    meta_df <- meta_df[, !colnames(meta_df) %in% raw_chars, drop = FALSE]
  }

  # Clean up parsed characteristic column names: "agent:ch1" -> "agent"
  colnames(meta_df) <- gsub(":ch\\d+$", "", colnames(meta_df))

  # Save metadata
  metadata_path <- file.path(output_dir, "metadata.csv")
  write.csv(meta_df, metadata_path)

  message("Saved metadata: ", nrow(meta_df), " samples x ", ncol(meta_df), " columns")

  # Check if expression data is available
  expr_mat <- tryCatch(exprs(eset), error = function(e) NULL)
  has_expression <- !is.null(expr_mat) && nrow(expr_mat) > 0 && ncol(expr_mat) > 0

  result <- list(
    metadata_path = metadata_path,
    metadata_shape = c(nrow(meta_df), ncol(meta_df)),
    expression_available = has_expression,
    source = "geo"
  )

  if (has_expression) {
    counts_path <- file.path(output_dir, "counts.csv")
    write.csv(expr_mat, counts_path)
    result$counts_path <- counts_path
    result$counts_shape <- c(nrow(expr_mat), ncol(expr_mat))
    message("Saved expression data: ", nrow(expr_mat), " features x ", ncol(expr_mat), " samples")
  }

  cat(toJSON(result, auto_unbox = TRUE))
  cat("\n")

  # Keep downloaded GEO series matrix files in the cache directory so
  # subsequent fetches for the same study skip the download.

}, error = function(e) {
  cat(toJSON(list(error = paste0("GEO fetch failed: ", e$message)), auto_unbox = TRUE))
  cat("\n")
})
