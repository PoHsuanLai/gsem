#!/usr/bin/env Rscript
# Compare R and Rust LDSC outputs and generate results.md

library(jsonlite)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

cat("=== Comparing R vs Rust outputs ===\n\n")

# Load R outputs
r_ldsc <- fromJSON(readLines("out_r/ldsc_result.json"))
r_timings <- fromJSON(readLines("out_r/timings.json"))

# Load Rust outputs (faer serde format: {nrows, ncols, data: [col-major]})
rust_raw <- fromJSON(readLines("out_rust/ldsc.json"))

faer_to_matrix <- function(obj) {
  nr <- obj$nrows
  nc <- obj$ncols
  matrix(obj$data, nrow = nr, ncol = nc, byrow = FALSE)
}

S_r <- as.matrix(data.frame(r_ldsc$S))
S_rust <- faer_to_matrix(rust_raw$s)

V_r <- as.matrix(data.frame(r_ldsc$V))
V_rust <- faer_to_matrix(rust_raw$v)

I_r <- as.matrix(data.frame(r_ldsc$I))
I_rust <- faer_to_matrix(rust_raw$i_mat)

# Compute max absolute differences
s_diff <- max(abs(S_r - S_rust))
v_diff <- max(abs(V_r - V_rust))
i_diff <- max(abs(I_r - I_rust))

cat(sprintf("S matrix max diff: %.6e\n", s_diff))
cat(sprintf("V matrix max diff: %.6e\n", v_diff))
cat(sprintf("I matrix max diff: %.6e\n", i_diff))

cat("\nS (R):\n")
print(round(S_r, 4))
cat("\nS (Rust):\n")
print(round(S_rust, 4))
cat("\nI (R):\n")
print(round(I_r, 4))
cat("\nI (Rust):\n")
print(round(I_rust, 4))

# Load timings
rust_timings <- fromJSON(readLines("out_rust/timings.json"))

# Generate results.md
lines <- c(
  "# Benchmark Results: Rust vs R GenomicSEM",
  "",
  sprintf("Data: 3 simulated traits, ~1.29M SNPs, N=50,000"),
  "",
  "## Timing",
  "",
  "| Step | R (s) | Rust (s) | Speedup |",
  "|------|-------|----------|---------|",
  sprintf("| LDSC | %.2f | %.2f | **%.1fx** |",
    r_timings$ldsc, rust_timings$ldsc,
    r_timings$ldsc / rust_timings$ldsc),
  sprintf("| SEM  | %.2f | %.2f | **%.1fx** |",
    r_timings$sem, rust_timings$sem,
    r_timings$sem / rust_timings$sem),
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
  "## S Matrix (Rust)",
  "```",
  capture.output(print(round(S_rust, 4))),
  "```"
)

writeLines(lines, "results.md")
cat("\nResults written to bench/results.md\n")
