#!/usr/bin/env Rscript
# Verify the upstream commonfactorGWAS patch fixes the crash on the
# pathological SNPs (rs3863622, rs4411087, rs4422948) identified earlier,
# and that the patched R output agrees numerically with gsemr.
#
# Run from bench/.

setwd("/Users/pohsuanlai/Documents/GenomicSEM-rs/bench")
suppressPackageStartupMessages({
  library(GenomicSEM)
  library(gsemr)
})

source("patch_genomicsem.R")

trait_names <- c("ANX", "OCD", "PTSD")

# --- 0. Load cached cov and sumstats (built by debug_usergwas_divergence.R) ---
if (!file.exists("cache/cov.rds") || !file.exists("cache/sumstats.rds")) {
  stop("Missing cache/cov.rds or cache/sumstats.rds; ",
       "run debug_usergwas_divergence.R first")
}
cov_r <- readRDS("cache/cov.rds")
ss    <- readRDS("cache/sumstats.rds")

# Small input: first 60 SNPs of the 1000-SNP subset (enough to include
# rs3863622 at row 51) plus the two later-row crashers rs4411087 (row 96)
# and rs4422948 (row 99) appended explicitly. Total ~62 SNPs, which is
# tiny for a bench but fully exercises the fix.
snps <- rbind(
  head(ss, 60),
  ss[ss$SNP %in% c("rs4411087", "rs4422948"), ]
)
snps <- snps[!duplicated(snps$SNP), ]
cat(sprintf("input: %d SNPs (includes %s)\n",
            nrow(snps),
            paste(intersect(c("rs3863622","rs4411087","rs4422948"), snps$SNP),
                  collapse = ", ")))

# --- 1. Baseline: confirm the crash still reproduces BEFORE patching ---------
cat("\n==============================================================\n")
cat("  [1/4] baseline: unpatched R on head(ss, 1000)\n")
cat("==============================================================\n")
baseline <- tryCatch(
  GenomicSEM::commonfactorGWAS(cov_r, SNPs = snps, parallel = FALSE),
  error = function(e) e
)
if (inherits(baseline, "error")) {
  cat("  BASELINE: crashed as expected ->",
      conditionMessage(baseline), "\n")
} else {
  cat("  BASELINE: surprise — unpatched R succeeded with nrow =",
      nrow(baseline), "\n")
}

# --- 2. Install patch and rerun ----------------------------------------------
cat("\n==============================================================\n")
cat("  [2/4] applying patch\n")
cat("==============================================================\n")
ok <- patch_genomicsem()
if (!ok) stop("patch_genomicsem failed; see warnings above")

cat("\n==============================================================\n")
cat("  [3/4] patched R on head(ss, 1000)\n")
cat("==============================================================\n")
t0 <- Sys.time()
patched <- tryCatch(
  GenomicSEM::commonfactorGWAS(cov_r, SNPs = snps, parallel = FALSE),
  error = function(e) e
)
t1 <- Sys.time()
if (inherits(patched, "error")) {
  stop("PATCHED R still crashed: ", conditionMessage(patched))
}
cat(sprintf("  PATCHED: ok, elapsed %.2fs, nrow = %d, cols = %s\n",
            as.numeric(t1 - t0, units = "secs"),
            nrow(patched), paste(names(patched), collapse = ", ")))
cat("  any se_c NA:", sum(is.na(patched$se_c)), "\n")
cat("  any Q==Not Computed:",
    sum(as.character(patched$Q) == "Not Computed", na.rm = TRUE), "\n")

targets <- c("rs3863622", "rs4411087", "rs4422948")
cat("\n  rows for previously-pathological SNPs:\n")
for (t in targets) {
  r <- patched[patched$SNP == t, ]
  if (nrow(r) == 0) { cat("    ", t, ": not in output\n"); next }
  cat(sprintf("    %-12s est=%+.4e  se_c=%s  Q=%s  warning=%s\n",
              t, r$est,
              ifelse(is.na(r$se_c), "NA    ", sprintf("%.3e", r$se_c)),
              substr(as.character(r$Q), 1, 12),
              substr(as.character(r$warning), 1, 50)))
}

# --- 4. Compare patched R vs gsemr -------------------------------------------
cat("\n==============================================================\n")
cat("  [4/4] gsemr::commonfactorGWAS on the same input\n")
cat("==============================================================\n")
t0 <- Sys.time()
rust_res <- gsemr::commonfactorGWAS(cov_r, SNPs = snps, parallel = FALSE)
t1 <- Sys.time()
cat(sprintf("  Rust: elapsed %.2fs, nrow = %d\n",
            as.numeric(t1 - t0, units = "secs"), nrow(rust_res)))

# Join on SNP. R output has cols SNP, est, se, se_c, Q, Z_Estimate, Pval_Estimate, ...
r_cmp <- data.frame(
  SNP = patched$SNP, est_r = patched$est, se_r = patched$se_c,
  stringsAsFactors = FALSE
)
rust_cmp <- data.frame(
  SNP = rust_res$SNP, est_rust = rust_res$est, se_rust = rust_res$se,
  stringsAsFactors = FALSE
)
cmp <- merge(r_cmp, rust_cmp, by = "SNP")
cmp$est_diff <- cmp$est_rust - cmp$est_r
cmp$est_ratio <- cmp$est_rust / cmp$est_r
cmp$se_ratio  <- cmp$se_rust  / cmp$se_r

# Strip rows where either side is NA (shouldn't be many after the patch)
valid <- is.finite(cmp$est_r) & is.finite(cmp$est_rust) &
         is.finite(cmp$se_r)  & is.finite(cmp$se_rust)
cat(sprintf("\n  joined rows: %d (%d fully finite)\n", nrow(cmp), sum(valid)))

finite_cmp <- cmp[valid, ]
cat(sprintf("  cor(est_r, est_rust) = %.6f\n",
            cor(finite_cmp$est_r, finite_cmp$est_rust)))
cat(sprintf("  cor(se_r,  se_rust)  = %.6f\n",
            cor(finite_cmp$se_r,  finite_cmp$se_rust)))
cat(sprintf("  max |est diff|       = %.3e\n", max(abs(finite_cmp$est_diff))))
cat(sprintf("  median |est diff|    = %.3e\n", median(abs(finite_cmp$est_diff))))
cat(sprintf("  est ratio (rust/R)   median = %.4f  range = [%.4f, %.4f]\n",
            median(finite_cmp$est_ratio),
            min(finite_cmp$est_ratio), max(finite_cmp$est_ratio)))
cat(sprintf("  se  ratio (rust/R)   median = %.4f  range = [%.4f, %.4f]\n",
            median(finite_cmp$se_ratio),
            min(finite_cmp$se_ratio), max(finite_cmp$se_ratio)))

cat("\n  head of worst est disagreements:\n")
ord <- order(-abs(finite_cmp$est_diff))
print(head(finite_cmp[ord, c("SNP","est_r","est_rust","est_diff","se_r","se_rust")], 10))

cat("\n  pathological SNPs R vs Rust:\n")
for (t in targets) {
  r <- finite_cmp[finite_cmp$SNP == t, ]
  if (nrow(r) == 0) {
    cat("    ", t, ": dropped by finite filter (R has NA?)\n")
    next
  }
  cat(sprintf("    %-12s  R: est=%+.4e se=%+.4e  |  Rust: est=%+.4e se=%+.4e\n",
              t, r$est_r, r$se_r, r$est_rust, r$se_rust))
}
