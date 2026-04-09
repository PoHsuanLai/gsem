#!/usr/bin/env Rscript
# Download reference data for benchmarking (replaces setup.sh).
# Idempotent: skips downloads if files already exist.

ensure_bench_data <- function(bench_dir = NULL) {
  if (is.null(bench_dir)) {
    # Detect script directory when run via Rscript
    args <- commandArgs(trailingOnly = FALSE)
    script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
    if (!is.na(script_dir) && nzchar(script_dir)) {
      bench_dir <- script_dir
    } else {
      bench_dir <- "."
    }
  }

  data_dir <- file.path(bench_dir, "data")
  dir.create(file.path(data_dir, "eur_w_ld_chr"), showWarnings = FALSE, recursive = TRUE)

  # -- LD scores (European ancestry, ~50MB compressed) --
  ld_check <- file.path(data_dir, "eur_w_ld_chr", "1.l2.ldscore.gz")
  if (!file.exists(ld_check)) {
    cat("Downloading LD scores from Zenodo...\n")
    dest <- file.path(data_dir, "eur_w_ld_chr.tar.gz")
    download.file(
      "https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1",
      dest, method = "libcurl", mode = "wb"
    )
    untar(dest, exdir = data_dir)
    file.remove(dest)
    cat("LD scores extracted.\n")
  } else {
    cat("LD scores already present.\n")
  }

  # -- HapMap3 SNP list (~8MB compressed) --
  hm3_file <- file.path(data_dir, "w_hm3.snplist")
  if (!file.exists(hm3_file)) {
    cat("Downloading HapMap3 SNP list from Zenodo...\n")
    dest <- file.path(data_dir, "w_hm3.snplist.gz")
    download.file(
      "https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1",
      dest, method = "libcurl", mode = "wb"
    )
    system2("gunzip", dest)
    cat("HapMap3 SNP list extracted.\n")
  } else {
    cat("HapMap3 SNP list already present.\n")
  }

  # -- Simulated GWAS data (3 traits, N=50k) --
  gwas_check <- file.path(data_dir, "iter1GWAS1.sumstats.gz")
  # Also check the dot-separated variant from simLDSC
  gwas_check2 <- file.path(data_dir, "iter1.GWAS1.sumstats.gz")
  if (!file.exists(gwas_check) && !file.exists(gwas_check2)) {
    cat("Generating simulated GWAS data with simLDSC...\n")
    source(file.path(bench_dir, "generate_bench_data.R"))
  } else {
    cat("Simulated GWAS data already present.\n")
  }

  cat("\nSetup complete.\n")
  invisible(data_dir)
}

# Run if called directly
if (sys.nframe() == 0) {
  ensure_bench_data()
}
