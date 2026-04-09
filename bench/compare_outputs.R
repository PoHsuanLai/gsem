#!/usr/bin/env Rscript
# Standalone correctness validation: R GenomicSEM vs gsemr (Rust).
# Runs both pipelines once and compares S, V, I matrices.
#
# Usage: cd bench && Rscript compare_outputs.R

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

library(jsonlite)
library(GenomicSEM)

if (!requireNamespace("gsemr", quietly = TRUE)) {
  stop("gsemr package not installed. Install with: R CMD INSTALL ../crates/gsem-r/")
}
library(gsemr)

traits <- paste0("data/iter1GWAS", 1:3, ".sumstats.gz")
ld     <- "data/eur_w_ld_chr/"
trait_names <- c("V1", "V2", "V3")

cat("=== Running R GenomicSEM LDSC ===\n")
r_result <- ldsc(
  traits = traits,
  sample.prev = c(NA, NA, NA),
  population.prev = c(NA, NA, NA),
  ld = ld, wld = ld,
  trait.names = trait_names,
  n.blocks = 200
)

cat("\n=== Running gsemr (Rust) LDSC ===\n")
rust_json <- ldsc_r(traits, c(NA, NA, NA), c(NA, NA, NA), ld, ld, 200L)
rust_result <- fromJSON(rust_json)

# Extract matrices
S_r <- as.matrix(r_result$S)
V_r <- as.matrix(r_result$V)
I_r <- as.matrix(r_result$I)

S_rust <- as.matrix(rust_result$s)
V_rust <- as.matrix(rust_result$v)
I_rust <- as.matrix(rust_result$i_mat)

# Compare
s_diff <- max(abs(S_r - S_rust))
v_diff <- max(abs(V_r - V_rust))
i_diff <- max(abs(I_r - I_rust))

cat("\n=== Comparison ===\n")
cat(sprintf("S matrix max diff: %.6e\n", s_diff))
cat(sprintf("V matrix max diff: %.6e\n", v_diff))
cat(sprintf("I matrix max diff: %.6e\n", i_diff))

cat("\nS (R):\n"); print(round(S_r, 4))
cat("\nS (Rust):\n"); print(round(S_rust, 4))
cat("\nI (R):\n"); print(round(I_r, 4))
cat("\nI (Rust):\n"); print(round(I_rust, 4))

# Thresholds
pass <- s_diff < 1e-4 && v_diff < 1e-3 && i_diff < 1e-4
cat(sprintf("\nOverall: %s\n", if (pass) "PASS" else "FAIL"))
if (!pass) quit(status = 1)
