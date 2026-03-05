#!/usr/bin/env Rscript
# Fetch raw RNA-seq count matrix from GEO.
# Uses the GEO rnaseq_counts download endpoint (GEO2R-style),
# falling back to supplementary files and then series matrix.
# Saves counts.csv (+ metadata.csv if available).
# Outputs JSON to stdout with file paths.
#
# Usage: Rscript fetch_geo_counts.R <gse_id> <library_path> <output_dir>

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat('{"error": "Usage: Rscript fetch_geo_counts.R <gse_id> <library_path> <output_dir>"}\n')
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

# Helper: clean up GEO metadata columns
clean_geo_metadata <- function(meta_df) {
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
  # If parsed characteristic columns exist (e.g. "agent:ch1"), drop raw ones
  parsed_chars <- grep(":ch1$", colnames(meta_df), value = TRUE)
  if (length(parsed_chars) > 0) {
    raw_chars <- grep("^characteristics_ch", colnames(meta_df), value = TRUE)
    meta_df <- meta_df[, !colnames(meta_df) %in% raw_chars, drop = FALSE]
  }
  colnames(meta_df) <- gsub(":ch\\d+$", "", colnames(meta_df))
  meta_df
}

tryCatch({
  counts_mat <- NULL
  counts_source <- NULL
  metadata_path <- NULL
  meta_df <- NULL

  # ---------------------------------------------------------------
  # Strategy 1: GEO rnaseq_counts endpoint (GEO2R-style download)
  # ---------------------------------------------------------------
  message("Trying GEO rnaseq_counts endpoint for: ", gse_id)

  urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
  counts_file_name <- paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  url <- paste(urld, paste0("acc=", gse_id), paste0("file=", counts_file_name), sep = "&")

  tryCatch({
    # data.table::fread handles gzipped URLs efficiently
    if (requireNamespace("data.table", quietly = TRUE)) {
      tbl <- as.matrix(data.table::fread(url, header = TRUE, colClasses = "integer"), rownames = 1)
    } else {
      # Fallback: download to temp file and read
      tmp <- file.path(output_dir, counts_file_name)
      download.file(url, tmp, mode = "wb", quiet = TRUE)
      tbl <- as.matrix(read.table(gzfile(tmp), header = TRUE, row.names = 1,
                                  sep = "\t", check.names = FALSE))
      tryCatch(file.remove(tmp), error = function(e) NULL)
    }

    if (nrow(tbl) > 0 && ncol(tbl) > 0) {
      counts_mat <- tbl
      counts_source <- "geo_rnaseq_counts"
      message("Downloaded raw counts from GEO: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples")
    }
  }, error = function(e) {
    message("GEO rnaseq_counts endpoint not available: ", e$message)
  })

  # ---------------------------------------------------------------
  # Strategy 2: Supplementary count files
  # ---------------------------------------------------------------
  if (is.null(counts_mat)) {
    message("Trying GEO supplementary files...")

    suppl_files <- tryCatch({
      getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = output_dir)
    }, error = function(e) NULL)

    if (!is.null(suppl_files) && nrow(suppl_files) > 0) {
      fnames <- rownames(suppl_files)
      message("Found ", length(fnames), " supplementary file(s)")

      # Prioritize files that look like count matrices
      count_patterns <- c("count", "raw", "gene_expression", "readcount",
                          "featurecount", "htseq")
      count_file <- NULL

      for (pat in count_patterns) {
        matches <- grep(pat, basename(fnames), ignore.case = TRUE, value = TRUE)
        if (length(matches) > 0) {
          count_file <- matches[1]
          break
        }
      }

      # If no obvious count file, try the first tabular file
      if (is.null(count_file)) {
        tabular <- grep("\\.(txt|csv|tsv|tab)(\\.gz)?$", basename(fnames),
                        ignore.case = TRUE, value = TRUE)
        if (length(tabular) > 0) {
          count_file <- tabular[1]
        }
      }

      if (!is.null(count_file)) {
        message("Reading count file: ", basename(count_file))

        tryCatch({
          # Handle gzipped files
          if (grepl("\\.gz$", count_file)) {
            con <- gzfile(count_file)
          } else {
            con <- file(count_file)
          }

          # Read first lines to detect separator
          first_lines <- readLines(con, n = 5)
          close(con)

          if (any(grepl("\t", first_lines))) {
            sep <- "\t"
          } else if (any(grepl(",", first_lines))) {
            sep <- ","
          } else {
            sep <- ""
          }

          df <- read.table(count_file, header = TRUE, sep = sep, row.names = 1,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           comment.char = "", fill = TRUE)

          numeric_cols <- sapply(df, is.numeric)
          if (sum(numeric_cols) >= 2) {
            counts_mat <- as.matrix(df[, numeric_cols, drop = FALSE])
            counts_source <- "supplementary"
            message("Parsed count matrix: ", nrow(counts_mat), " genes x ",
                    ncol(counts_mat), " samples")
          }
        }, error = function(e) {
          message("Could not parse supplementary file: ", e$message)
        })
      }
    }
  }

  # ---------------------------------------------------------------
  # Strategy 3: Series matrix expression data
  # ---------------------------------------------------------------
  if (is.null(counts_mat)) {
    message("No supplementary count file found. Trying series matrix...")

    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = output_dir)

    if (is.list(gse) && length(gse) > 0) {
      eset <- gse[[1]]
    } else if (!is.list(gse)) {
      eset <- gse
    } else {
      cat(toJSON(list(error = paste0("No data found for ", gse_id)), auto_unbox = TRUE))
      cat("\n")
      quit(status = 0)
    }

    expr_mat <- tryCatch(exprs(eset), error = function(e) NULL)
    if (!is.null(expr_mat) && nrow(expr_mat) > 0 && ncol(expr_mat) > 0) {
      counts_mat <- expr_mat
      counts_source <- "series_matrix"
      message("Extracted expression data from series matrix: ",
              nrow(counts_mat), " features x ", ncol(counts_mat), " samples")
    }

    # Also extract metadata while we have the eset
    if (exists("eset")) {
      meta_df <- clean_geo_metadata(pData(eset))
      metadata_path <- file.path(output_dir, "metadata.csv")
      write.csv(meta_df, metadata_path)
      message("Saved metadata: ", nrow(meta_df), " samples x ", ncol(meta_df), " columns")
    }
  }

  # ---------------------------------------------------------------
  # Also fetch metadata if we got counts from strategy 1 or 2
  # ---------------------------------------------------------------
  if (!is.null(counts_mat) && is.null(metadata_path)) {
    message("Fetching metadata from GEO series matrix...")
    tryCatch({
      gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = output_dir)
      if (is.list(gse) && length(gse) > 0) {
        eset <- gse[[1]]
      } else if (!is.list(gse)) {
        eset <- gse
      }
      if (exists("eset")) {
        meta_df <- clean_geo_metadata(pData(eset))
        metadata_path <- file.path(output_dir, "metadata.csv")
        write.csv(meta_df, metadata_path)
        message("Saved metadata: ", nrow(meta_df), " samples x ", ncol(meta_df), " columns")
      }
    }, error = function(e) {
      message("Could not fetch metadata: ", e$message)
    })
  }

  # ---------------------------------------------------------------
  # No counts found at all
  # ---------------------------------------------------------------
  if (is.null(counts_mat)) {
    cat(toJSON(list(
      error = paste0("No count data found for ", gse_id,
                     ". Try recount3 for uniformly-processed RNA-seq counts.")
    ), auto_unbox = TRUE))
    cat("\n")
    quit(status = 0)
  }

  # Save counts
  counts_path <- file.path(output_dir, "counts.csv")
  write.csv(counts_mat, counts_path)

  result <- list(
    counts_path = counts_path,
    counts_shape = c(nrow(counts_mat), ncol(counts_mat)),
    counts_source = counts_source,
    source = "geo_counts"
  )

  if (!is.null(metadata_path) && file.exists(metadata_path)) {
    result$metadata_path <- metadata_path
    result$metadata_shape <- c(nrow(meta_df), ncol(meta_df))
  }

  cat(toJSON(result, auto_unbox = TRUE))
  cat("\n")

  # Keep downloaded GEO files in the cache directory so
  # subsequent fetches for the same study skip the download.

}, error = function(e) {
  cat(toJSON(list(error = paste0("GEO counts fetch failed: ", e$message)), auto_unbox = TRUE))
  cat("\n")
})
