#!/usr/bin/env Rscript
# Fetch RNA-seq counts + metadata from recount3 for a given study ID.
# Saves counts.csv and metadata.csv to a temp directory.
# Outputs JSON to stdout with file paths.
#
# Usage: Rscript fetch_recount3.R <study_id> <library_path> <output_dir>
#   study_id   : GEO (GSExxxxx) or SRA (SRPxxxxxx) accession
#   library_path: path to installed R packages
#   output_dir  : directory to write counts.csv and metadata.csv

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat('{"error": "Usage: Rscript fetch_recount3.R <study_id> <library_path> <output_dir>"}\n')
    quit(status = 0)
  }

  study_id <- args[1]
  lib_path <- args[2]
  output_dir <- args[3]

  .libPaths(c(lib_path, .libPaths()))

  library(recount3, quietly = TRUE)
  library(SummarizedExperiment, quietly = TRUE)
  library(jsonlite, quietly = TRUE)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  message("Searching for study: ", study_id)

  # Get available human projects
  human_projects <- available_projects()

  # Determine if input is GSE or SRP
  if (grepl("^GSE", study_id, ignore.case = TRUE)) {
    # GSE ID — look for matching project
    # recount3 stores projects by SRP ID, but some have GSE in the file_source
    # Try to find by matching the project field or external_id
    matching <- human_projects[grep(study_id, human_projects$project, ignore.case = TRUE), ]

    if (nrow(matching) == 0) {
      # Also search in file_source column
      matching <- human_projects[grep(study_id, human_projects$file_source, ignore.case = TRUE), ]
    }

    if (nrow(matching) == 0) {
      cat(toJSON(list(
        error = paste0("Study ", study_id, " not found in recount3. ",
                       "Try using the SRP ID instead, or use 'Fetch GEO Metadata' for metadata only.")
      ), auto_unbox = TRUE))
      cat("\n")
      quit(status = 0)
    }

    # Use the first match
    proj <- matching[1, ]
    message("Found matching project: ", proj$project, " (", proj$organism, ")")

  } else if (grepl("^SRP", study_id, ignore.case = TRUE)) {
    # SRP ID — direct lookup
    matching <- human_projects[human_projects$project == study_id, ]

    if (nrow(matching) == 0) {
      cat(toJSON(list(
        error = paste0("Study ", study_id, " not found in recount3.")
      ), auto_unbox = TRUE))
      cat("\n")
      quit(status = 0)
    }

    proj <- matching[1, ]
  } else {
    cat(toJSON(list(
      error = "Study ID must start with GSE (GEO) or SRP (SRA)."
    ), auto_unbox = TRUE))
    cat("\n")
    quit(status = 0)
  }

  # Create RangedSummarizedExperiment with gene-level counts
  message("Downloading counts from recount3...")
  rse <- create_rse(proj)

  # Extract raw counts matrix
  counts_mat <- assay(rse, "raw_counts")
  if (is.null(counts_mat)) {
    counts_mat <- assays(rse)[[1]]
  }

  # Use gene symbols if available, otherwise gene IDs
  if ("gene_name" %in% colnames(rowData(rse))) {
    gene_names <- rowData(rse)$gene_name
    # Handle duplicates by appending gene ID
    dups <- duplicated(gene_names) | gene_names == "" | is.na(gene_names)
    if (any(dups)) {
      gene_names[dups] <- paste0(gene_names[dups], "_", rownames(counts_mat)[dups])
    }
    rownames(counts_mat) <- gene_names
  }

  # Extract sample metadata
  meta_df <- as.data.frame(colData(rse))

  # Clean up metadata — keep useful columns
  useful_cols <- c("external_id", "study", "sra.sample_title",
                   "sra.experiment_title", "sra.sample_attributes",
                   "sra.library_layout", "sra.library_strategy",
                   "recount_qc.star.uniquely_mapped_reads_%_input",
                   "recount_qc.star.number_of_input_reads")
  keep_cols <- intersect(useful_cols, colnames(meta_df))
  if (length(keep_cols) > 0) {
    meta_clean <- meta_df[, keep_cols, drop = FALSE]
  } else {
    meta_clean <- meta_df
  }

  # Try to expand SRA attributes into separate columns
  if ("sra.sample_attributes" %in% colnames(meta_clean)) {
    tryCatch({
      attrs <- meta_clean$sra.sample_attributes
      # Parse "key;;value||key;;value" format
      parsed <- lapply(attrs, function(a) {
        pairs <- strsplit(as.character(a), "\\|\\|")[[1]]
        kv <- strsplit(pairs, ";;")
        setNames(sapply(kv, function(x) if (length(x) >= 2) x[2] else NA),
                 sapply(kv, function(x) x[1]))
      })
      attr_df <- do.call(rbind, lapply(parsed, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
      rownames(attr_df) <- rownames(meta_clean)
      meta_clean <- cbind(meta_clean, attr_df)
      meta_clean$sra.sample_attributes <- NULL
    }, error = function(e) {
      message("Could not parse SRA attributes: ", e$message)
    })
  }

  # Save to CSV
  counts_path <- file.path(output_dir, "counts.csv")
  metadata_path <- file.path(output_dir, "metadata.csv")

  write.csv(counts_mat, counts_path)
  write.csv(meta_clean, metadata_path)

  message("Saved counts: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples")
  message("Saved metadata: ", nrow(meta_clean), " samples x ", ncol(meta_clean), " columns")

  # Output paths as JSON to stdout
  result <- list(
    counts_path = counts_path,
    metadata_path = metadata_path,
    counts_shape = c(nrow(counts_mat), ncol(counts_mat)),
    metadata_shape = c(nrow(meta_clean), ncol(meta_clean)),
    source = "recount3",
    project_id = proj$project
  )
  cat(toJSON(result, auto_unbox = TRUE))
  cat("\n")

}, error = function(e) {
  cat(toJSON(list(error = paste0("recount3 fetch failed: ", e$message)), auto_unbox = TRUE))
  cat("\n")
})
