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

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

library(bench)
library(jsonlite)
library(GenomicSEM)

if (!requireNamespace("gsemr", quietly = TRUE)) {
  stop("gsemr not installed. Run: R CMD INSTALL ../bindings/r/")
}
library(gsemr)

# -- Configuration --
traits <- paste0("data/iter1GWAS", 1:3, ".sumstats.gz")
ld     <- "data/eur_w_ld_chr/"
trait_names <- c("V1", "V2", "V3")
model_str  <- "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3"

cat("\n========================================\n")
cat("  GenomicSEM Benchmark: R vs gsemr\n")
cat("========================================\n\n")

all_results <- list()

# ==========================================================================
# 1. LDSC
# ==========================================================================
cat("--- [1/6] LDSC ---\n")
all_results$ldsc <- bench::mark(
  R = GenomicSEM::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                       trait.names = trait_names, n.blocks = 200),
  gsemr = gsemr::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                       trait.names = trait_names, n.blocks = 200L),
  min_iterations = 3, check = FALSE, filter_gc = FALSE
)

# Pre-compute LDSC for downstream benchmarks
r_cov <- GenomicSEM::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                           trait.names = trait_names, n.blocks = 200)
rust_cov <- gsemr::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                          trait.names = trait_names, n.blocks = 200L)

# ==========================================================================
# 2. commonfactor
# ==========================================================================
cat("--- [2/6] commonfactor ---\n")
all_results$commonfactor <- bench::mark(
  R = GenomicSEM::commonfactor(r_cov, estimation = "DWLS"),
  gsemr = gsemr::commonfactor(rust_cov, estimation = "DWLS"),
  min_iterations = 5, check = FALSE, filter_gc = FALSE
)

# ==========================================================================
# 3. usermodel
# ==========================================================================
cat("--- [3/6] usermodel ---\n")
all_results$usermodel <- bench::mark(
  R = GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str),
  gsemr = gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str),
  min_iterations = 5, check = FALSE, filter_gc = FALSE
)

# ==========================================================================
# 4. paLDSC (parallel analysis)
# ==========================================================================
cat("--- [4/6] paLDSC ---\n")
# R GenomicSEM's paLDSC has `function(S = S, V = V, ...)` which causes
# "promise already under evaluation" errors. This is a known R bug in the
# original package. We benchmark gsemr only.
tryCatch({
  t_rust <- system.time(rust_pa <- gsemr::paLDSC(rust_cov, r = 500))["elapsed"]
  cat(sprintf("  gsemr: %.3fs (R GenomicSEM skipped: recursive default arg bug)\n", t_rust))
  cat(sprintf("  n_factors: %d\n", rust_pa$n_factors))
}, error = function(e) {
  cat("  SKIPPED:", conditionMessage(e), "\n")
})

# ==========================================================================
# 5. write.model
# ==========================================================================
cat("--- [5/6] write.model ---\n")
loadings <- matrix(c(0.7, 0.6, 0.5), ncol = 1)
rownames(loadings) <- trait_names

tryCatch({
  all_results$write.model <- bench::mark(
    R = GenomicSEM::write.model(loadings, r_cov$S, cutoff = 0.3),
    gsemr = gsemr::write.model(loadings, rust_cov$S, cutoff = 0.3),
    min_iterations = 50, check = FALSE, filter_gc = FALSE
  )
}, error = function(e) {
  cat("  SKIPPED:", conditionMessage(e), "\n")
})

# ==========================================================================
# 6. rgmodel
# ==========================================================================
cat("--- [6/6] rgmodel ---\n")
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
  "3 simulated traits, ~1.29M SNPs, N=50,000",
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
