#!/usr/bin/env Rscript
# Download real PGC GWAS summary statistics for benchmarking.
# Uses: Anxiety (ANX), OCD, PTSD — 3 psychiatric traits commonly used with GenomicSEM.
#
# These files are freely available from PGC via figshare (CC BY 4.0).

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Trait 1: Anxiety Disorders (ANX) ---
anx_file <- file.path(data_dir, "anxiety.meta.full.cc.tbl.gz")
if (!file.exists(anx_file)) {
  cat("Downloading Anxiety GWAS (PGC, Otowa et al. 2016)...\n")
  download.file(
    "https://ndownloader.figshare.com/files/28570812",
    anx_file, method = "libcurl", mode = "wb"
  )
  cat("  -> ", anx_file, "\n")
} else {
  cat("Anxiety GWAS already present.\n")
}

# --- Trait 2: OCD ---
ocd_file <- file.path(data_dir, "ocd_aug2017.gz")
if (!file.exists(ocd_file)) {
  cat("Downloading OCD GWAS (PGC, IOCDF-GC 2018)...\n")
  download.file(
    "https://ndownloader.figshare.com/files/28169544",
    ocd_file, method = "libcurl", mode = "wb"
  )
  cat("  -> ", ocd_file, "\n")
} else {
  cat("OCD GWAS already present.\n")
}

# --- Trait 3: PTSD (European ancestry) ---
ptsd_zip <- file.path(data_dir, "ptsd_ea.zip")
ptsd_file <- file.path(data_dir, "pts_all_freeze1_ea.results.gz")
# Also check for unzipped variant
ptsd_alt <- Sys.glob(file.path(data_dir, "pts_*ea*"))
if (length(ptsd_alt) == 0 && !file.exists(ptsd_file)) {
  cat("Downloading PTSD GWAS (PGC, Duncan et al. 2018)...\n")
  download.file(
    "https://ndownloader.figshare.com/files/28169589",
    ptsd_zip, method = "libcurl", mode = "wb"
  )
  unzip(ptsd_zip, exdir = data_dir)
  file.remove(ptsd_zip)
  # Find the extracted file
  ptsd_candidates <- Sys.glob(file.path(data_dir, "pts_*ea*"))
  if (length(ptsd_candidates) > 0) {
    ptsd_file <- ptsd_candidates[1]
  }
  cat("  -> ", ptsd_file, "\n")
} else {
  cat("PTSD GWAS already present.\n")
}

cat("\nReal GWAS data ready.\n")
cat("Files:\n")
for (f in list.files(data_dir, pattern = "\\.gz$")) {
  cat("  ", f, "\n")
}
