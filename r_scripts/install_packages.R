#!/usr/bin/env Rscript
# Install required R packages for GEO/recount3 data fetching.
# Called once during R environment setup.
#
# Usage: Rscript install_packages.R <library_path>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript install_packages.R <library_path>")
}

lib_path <- args[1]
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

message("Installing R packages to: ", lib_path)

# jsonlite for JSON output
message("Installing jsonlite...")
install.packages("jsonlite", repos = "https://cloud.r-project.org", lib = lib_path, quiet = FALSE)

# BiocManager for Bioconductor packages
message("Installing BiocManager...")
install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib_path, quiet = FALSE)

# recount3 for uniformly processed RNA-seq counts
message("Installing recount3 (this may take a few minutes)...")
BiocManager::install("recount3", lib = lib_path, update = FALSE, ask = FALSE)

# GEOquery for GEO SOFT metadata
message("Installing GEOquery...")
BiocManager::install("GEOquery", lib = lib_path, update = FALSE, ask = FALSE)

# SummarizedExperiment for working with recount3 objects
message("Installing SummarizedExperiment...")
BiocManager::install("SummarizedExperiment", lib = lib_path, update = FALSE, ask = FALSE)

message("All R packages installed successfully.")
