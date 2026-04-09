#!/usr/bin/env Rscript
# Benchmark: original R GenomicSEM vs Rust-backed gsemr package.
# Uses bench::mark() for statistically rigorous timing.
#
# Prerequisites:
#   install.packages(c("bench", "jsonlite"))
#   R CMD INSTALL bindings/r/
#
# Usage:
#   cd bench && Rscript benchmark.R

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

# -- Setup data --
source("setup_data.R")
ensure_bench_data(bench_dir = ".")

# -- Load packages --
library(bench)
library(jsonlite)
library(GenomicSEM)

if (!requireNamespace("gsemr", quietly = TRUE)) {
  stop(
    "gsemr package not installed.\n",
    "Install with: R CMD INSTALL ../bindings/r/\n",
    "(requires Rust toolchain)"
  )
}
library(gsemr)

# -- Configuration --
traits <- paste0("data/iter1GWAS", 1:3, ".sumstats.gz")
ld     <- "data/eur_w_ld_chr/"
trait_names <- c("V1", "V2", "V3")

cat("\n========================================\n")
cat("  GenomicSEM Benchmark: R vs gsemr\n")
cat("========================================\n\n")

# ==========================================================================
# LDSC Benchmark
# ==========================================================================
cat("--- Benchmarking LDSC ---\n")

ldsc_bench <- bench::mark(
  R_GenomicSEM = {
    GenomicSEM::ldsc(
      traits = traits,
      sample.prev = c(NA, NA, NA),
      population.prev = c(NA, NA, NA),
      ld = ld, wld = ld,
      trait.names = trait_names,
      n.blocks = 200
    )
  },
  gsemr_Rust = {
    gsemr::ldsc(
      traits = traits,
      sample.prev = c(NA, NA, NA),
      population.prev = c(NA, NA, NA),
      ld = ld, wld = ld,
      trait.names = trait_names,
      n.blocks = 200L
    )
  },
  min_iterations = 3,
  check = FALSE,
  filter_gc = FALSE
)

cat("\nLDSC results:\n")
print(ldsc_bench[, c("expression", "min", "median", "mem_alloc", "n_itr")])

# ==========================================================================
# Common Factor SEM Benchmark
# ==========================================================================
cat("\n--- Benchmarking Common Factor SEM ---\n")

# Pre-compute LDSC results for SEM input
r_covstruc <- GenomicSEM::ldsc(
  traits = traits,
  sample.prev = c(NA, NA, NA),
  population.prev = c(NA, NA, NA),
  ld = ld, wld = ld,
  trait.names = trait_names,
  n.blocks = 200
)
rust_covstruc <- gsemr::ldsc(
  traits = traits,
  sample.prev = c(NA, NA, NA),
  population.prev = c(NA, NA, NA),
  ld = ld, wld = ld,
  trait.names = trait_names,
  n.blocks = 200L
)

sem_bench <- bench::mark(
  R_GenomicSEM = {
    GenomicSEM::commonfactor(r_covstruc, estimation = "DWLS")
  },
  gsemr_Rust = {
    gsemr::commonfactor(rust_covstruc, estimation = "DWLS")
  },
  min_iterations = 5,
  check = FALSE,
  filter_gc = FALSE
)

cat("\nSEM results:\n")
print(sem_bench[, c("expression", "min", "median", "mem_alloc", "n_itr")])

# ==========================================================================
# User Model Benchmark
# ==========================================================================
cat("\n--- Benchmarking User Model ---\n")

model_str <- "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3"

usermodel_bench <- bench::mark(
  R_GenomicSEM = {
    GenomicSEM::usermodel(r_covstruc, estimation = "DWLS", model = model_str)
  },
  gsemr_Rust = {
    gsemr::usermodel(rust_covstruc, estimation = "DWLS", model = model_str)
  },
  min_iterations = 5,
  check = FALSE,
  filter_gc = FALSE
)

cat("\nUser Model results:\n")
print(usermodel_bench[, c("expression", "min", "median", "mem_alloc", "n_itr")])

# ==========================================================================
# Speedup Summary
# ==========================================================================
ldsc_speedup <- as.numeric(ldsc_bench$median[1]) / as.numeric(ldsc_bench$median[2])
sem_speedup  <- as.numeric(sem_bench$median[1]) / as.numeric(sem_bench$median[2])
um_speedup   <- as.numeric(usermodel_bench$median[1]) / as.numeric(usermodel_bench$median[2])

cat("\n========================================\n")
cat(sprintf("  LDSC         speedup (median): %.1fx\n", ldsc_speedup))
cat(sprintf("  commonfactor speedup (median): %.1fx\n", sem_speedup))
cat(sprintf("  usermodel    speedup (median): %.1fx\n", um_speedup))
cat("========================================\n")

# ==========================================================================
# Correctness Comparison
# ==========================================================================
cat("\n--- Correctness Check ---\n")

S_r <- as.matrix(r_covstruc$S)
V_r <- as.matrix(r_covstruc$V)
I_r <- as.matrix(r_covstruc$I)

S_rust <- as.matrix(rust_covstruc$S)
V_rust <- as.matrix(rust_covstruc$V)
I_rust <- as.matrix(rust_covstruc$I)

s_diff <- max(abs(S_r - S_rust))
v_diff <- max(abs(V_r - V_rust))
i_diff <- max(abs(I_r - I_rust))

cat(sprintf("  S matrix max diff: %.2e\n", s_diff))
cat(sprintf("  V matrix max diff: %.2e\n", v_diff))
cat(sprintf("  I matrix max diff: %.2e\n", i_diff))

# ==========================================================================
# Write results.md
# ==========================================================================
fmt_time <- function(x) sprintf("%.2f", as.numeric(x))

lines <- c(
  "# Benchmark Results: gsemr (Rust) vs R GenomicSEM",
  "",
  sprintf("Data: 3 simulated traits, N=50,000 | bench::mark"),
  "",
  "## Timing (median)",
  "",
  "| Step | R GenomicSEM | gsemr (Rust) | Speedup |",
  "|------|-------------|-------------|---------|",
  sprintf("| LDSC | %ss | %ss | **%.1fx** |",
    fmt_time(ldsc_bench$median[1]), fmt_time(ldsc_bench$median[2]), ldsc_speedup),
  sprintf("| commonfactor | %ss | %ss | **%.1fx** |",
    fmt_time(sem_bench$median[1]), fmt_time(sem_bench$median[2]), sem_speedup),
  sprintf("| usermodel | %ss | %ss | **%.1fx** |",
    fmt_time(usermodel_bench$median[1]), fmt_time(usermodel_bench$median[2]), um_speedup),
  "",
  "## Numerical Accuracy (max element-wise absolute diff)",
  "",
  "| Matrix | Max Diff |",
  "|--------|---------|",
  sprintf("| S (genetic cov) | %.2e |", s_diff),
  sprintf("| V (sampling cov) | %.2e |", v_diff),
  sprintf("| I (intercepts) | %.2e |", i_diff),
  "",
  "## S Matrix (R)",
  "```",
  capture.output(print(round(S_r, 4))),
  "```",
  "",
  "## S Matrix (gsemr / Rust)",
  "```",
  capture.output(print(round(S_rust, 4))),
  "```"
)

writeLines(lines, "results.md")
cat("\nResults written to bench/results.md\n")
