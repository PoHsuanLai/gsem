#!/usr/bin/env Rscript
# Benchmark: R GenomicSEM vs gsemr (Rust) on real PGC GWAS data.
#
# Uses 3 PGC psychiatric traits: Anxiety, OCD, PTSD.
# Runs munge → LDSC → commonfactor/usermodel and compares timing + correctness.
#
# Prerequisites:
#   Rscript setup_data.R            # downloads LD scores + GWAS data
#   R CMD INSTALL ../bindings/r/    # installs gsemr
#
# Usage: cd bench && Rscript benchmark.R

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

library(bench)
library(GenomicSEM)

if (!requireNamespace("gsemr", quietly = TRUE)) {
  stop("gsemr not installed. Run: R CMD INSTALL ../bindings/r/")
}
library(gsemr)

# -- Configuration --
ptsd_file <- Sys.glob("data/*PTSD*")[1]
if (is.na(ptsd_file)) ptsd_file <- Sys.glob("data/pts_*ea*")[1]
raw_traits <- c(
  "data/anxiety.meta.full.cc.tbl.gz",
  "data/ocd_aug2017.gz",
  ptsd_file
)
if (any(is.na(raw_traits))) stop("Missing GWAS files. Run: Rscript setup_data.R")

ld     <- "data/eur_w_ld_chr/"
hm3    <- "data/w_hm3.snplist"
trait_names <- c("ANX", "OCD", "PTSD")

# Sample/population prevalences for binary traits
sample_prev <- c(0.5, 0.5, 0.5)  # case-control ~50/50
pop_prev    <- c(0.16, 0.02, 0.07)  # population prevalences

# OCD has no N column; PTSD Neff = 5831
n_overrides <- c(NA, 9725, 5831)

cat("\n========================================\n")
cat("  GenomicSEM Benchmark: R vs gsemr\n")
cat("  Data: PGC Anxiety, OCD, PTSD\n")
cat("========================================\n\n")

all_results <- list()

# ==========================================================================
# 1. Munge
# ==========================================================================
cat("--- [1/5] Munge ---\n")

# Check if already munged
munged_traits <- paste0(trait_names, ".sumstats.gz")
if (!all(file.exists(munged_traits))) {
  # Munge with R
  t_r_munge <- system.time({
    GenomicSEM::munge(
      files = raw_traits,
      hm3 = hm3,
      trait.names = trait_names,
      N = n_overrides,
      info.filter = 0.9,
      maf.filter = 0.01
    )
  })["elapsed"]
  cat(sprintf("  R munge: %.1fs\n", t_r_munge))

  # Munge with gsemr (files already exist from R, skip to avoid overwrite race)
  # gsemr::munge would produce identical output
} else {
  cat("  Already munged. Skipping.\n")
  t_r_munge <- NA
}

# ==========================================================================
# 2. LDSC
# ==========================================================================
cat("--- [2/5] LDSC ---\n")
all_results$ldsc <- bench::mark(
  R = GenomicSEM::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                       trait.names = trait_names, n.blocks = 200),
  gsemr = gsemr::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                       trait.names = trait_names, n.blocks = 200L),
  min_iterations = 3, check = FALSE, filter_gc = FALSE
)

# Pre-compute for downstream
r_cov <- GenomicSEM::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                           trait.names = trait_names, n.blocks = 200)
rust_cov <- gsemr::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                          trait.names = trait_names, n.blocks = 200L)

# ==========================================================================
# 3. Common Factor
# ==========================================================================
cat("--- [3/5] commonfactor ---\n")
all_results$commonfactor <- bench::mark(
  R = GenomicSEM::commonfactor(r_cov, estimation = "DWLS"),
  gsemr = gsemr::commonfactor(rust_cov, estimation = "DWLS"),
  min_iterations = 5, check = FALSE, filter_gc = FALSE
)

# ==========================================================================
# 4. User Model
# ==========================================================================
cat("--- [4/5] usermodel ---\n")
model_str <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~~ 1*F1\n",
  paste(trait_names, "~~", trait_names, collapse = "\n")
)

all_results$usermodel <- bench::mark(
  R = GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str),
  gsemr = gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str),
  min_iterations = 5, check = FALSE, filter_gc = FALSE
)

# ==========================================================================
# 5. rgmodel
# ==========================================================================
cat("--- [5/5] rgmodel ---\n")
tryCatch({
  all_results$rgmodel <- bench::mark(
    R = GenomicSEM::rgmodel(r_cov),
    gsemr = gsemr::rgmodel(rust_cov),
    min_iterations = 5, check = FALSE, filter_gc = FALSE
  )
}, error = function(e) {
  cat("  SKIPPED:", conditionMessage(e), "\n")
})

# ==========================================================================
# Summary
# ==========================================================================
cat("\n========================================\n")
cat("  Results\n")
cat("========================================\n\n")

for (name in names(all_results)) {
  res <- all_results[[name]]
  t_r <- as.numeric(res$median[1])
  t_rust <- as.numeric(res$median[2])
  speedup <- t_r / t_rust
  cat(sprintf("%-15s  R: %8s  gsemr: %8s  speedup: %.1fx\n",
    name, format(res$median[1]), format(res$median[2]), speedup
  ))
}

# ==========================================================================
# Correctness
# ==========================================================================
cat("\n--- Correctness (LDSC matrices) ---\n")
s_diff <- max(abs(as.matrix(r_cov$S) - as.matrix(rust_cov$S)))
v_diff <- max(abs(as.matrix(r_cov$V) - as.matrix(rust_cov$V)))
i_diff <- max(abs(as.matrix(r_cov$I) - as.matrix(rust_cov$I)))
cat(sprintf("  S max diff: %.2e\n  V max diff: %.2e\n  I max diff: %.2e\n", s_diff, v_diff, i_diff))

# ==========================================================================
# Write results.md
# ==========================================================================
fmt <- function(x) sprintf("%.2f", as.numeric(x))

lines <- c(
  "# Benchmark: gsemr (Rust) vs R GenomicSEM",
  "",
  "Data: PGC Anxiety, OCD, PTSD (real GWAS summary statistics)",
  "",
  "## Timing (median)",
  "",
  "| Function | R GenomicSEM | gsemr (Rust) | Speedup |",
  "|----------|-------------|-------------|---------|"
)
for (name in names(all_results)) {
  res <- all_results[[name]]
  speedup <- as.numeric(res$median[1]) / as.numeric(res$median[2])
  lines <- c(lines, sprintf("| %s | %ss | %ss | **%.1fx** |",
    name, fmt(res$median[1]), fmt(res$median[2]), speedup))
}
lines <- c(lines, "",
  "## Numerical Accuracy",
  "",
  "| Matrix | Max Diff |",
  "|--------|---------|",
  sprintf("| S (genetic cov) | %.2e |", s_diff),
  sprintf("| V (sampling cov) | %.2e |", v_diff),
  sprintf("| I (intercepts) | %.2e |", i_diff)
)

writeLines(lines, "results.md")
cat("\nResults written to bench/results.md\n")
