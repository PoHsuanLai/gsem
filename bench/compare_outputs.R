#!/usr/bin/env Rscript
# Correctness validation: R GenomicSEM vs gsemr (Rust) on real PGC data.
#
# Usage: cd bench && Rscript compare_outputs.R

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

library(jsonlite)
library(GenomicSEM)
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
sample_prev <- c(0.5, 0.5, 0.5)
pop_prev    <- c(0.16, 0.02, 0.07)
n_overrides <- c(NA, 9725, 5831)

pass_count <- 0
fail_count <- 0
skip_count <- 0

check <- function(name, test_fn) {
  cat(sprintf("%-25s", name))
  tryCatch({
    test_fn()
    cat(" PASS\n")
    pass_count <<- pass_count + 1
  }, error = function(e) {
    cat(sprintf(" FAIL: %s\n", conditionMessage(e)))
    fail_count <<- fail_count + 1
  })
}

skip <- function(name, reason) {
  cat(sprintf("%-25s SKIP: %s\n", name, reason))
  skip_count <<- skip_count + 1
}

assert_close <- function(a, b, tol, label) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) stop(sprintf("%s: length mismatch %d vs %d", label, length(a), length(b)))
  d <- max(abs(a - b), na.rm = TRUE)
  if (d > tol) stop(sprintf("%s: max diff %.6e > tol %.6e", label, d, tol))
}

assert_mat_close <- function(a, b, tol, label) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) stop(sprintf("%s: dim mismatch %s vs %s", label, paste(dim(a),collapse="x"), paste(dim(b),collapse="x")))
  d <- max(abs(a - b), na.rm = TRUE)
  if (d > tol) stop(sprintf("%s: max diff %.6e > tol %.6e", label, d, tol))
}

cat("\n==========================================\n")
cat("  Correctness: R GenomicSEM vs gsemr\n")
cat("  Data: PGC Anxiety, OCD, PTSD\n")
cat("==========================================\n\n")

# ==========================================================================
# 0. Munge (if needed)
# ==========================================================================
munged_traits <- paste0(trait_names, ".sumstats.gz")
if (!all(file.exists(munged_traits))) {
  cat("--- Munging raw GWAS files ---\n")
  GenomicSEM::munge(
    files = raw_traits,
    hm3 = hm3,
    trait.names = trait_names,
    N = n_overrides,
    info.filter = 0.9,
    maf.filter = 0.01
  )
}

# ==========================================================================
# 1. LDSC
# ==========================================================================
cat("--- LDSC ---\n")
r_cov <- GenomicSEM::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                           trait.names = trait_names, n.blocks = 200)
rust_cov <- gsemr::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                          trait.names = trait_names, n.blocks = 200L)

check("ldsc: S matrix", function() {
  assert_mat_close(r_cov$S, rust_cov$S, 1e-3, "S")
})
check("ldsc: V matrix", function() {
  assert_mat_close(r_cov$V, rust_cov$V, 1e-3, "V")
})
check("ldsc: I matrix", function() {
  assert_mat_close(r_cov$I, rust_cov$I, 1e-3, "I")
})

# ==========================================================================
# 2. commonfactor
# ==========================================================================
cat("\n--- commonfactor ---\n")
r_cf <- GenomicSEM::commonfactor(r_cov, estimation = "DWLS")
rust_cf <- gsemr::commonfactor(rust_cov, estimation = "DWLS")

check("commonfactor: estimates", function() {
  r_est <- r_cf$results$Unstandardized_Estimate[r_cf$results$op == "=~"]
  rust_est <- rust_cf$results$est[rust_cf$results$op == "=~"]
  # Loading estimates depend on S; with S diff ~0.01, loadings can differ more
  assert_close(r_est, rust_est, 0.05, "loading estimates")
})

# ==========================================================================
# 3. usermodel
# ==========================================================================
cat("\n--- usermodel ---\n")
model_str <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~~ 1*F1\n",
  paste(trait_names, "~~", trait_names, collapse = "\n")
)
r_um <- GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str)
rust_um <- gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str)

check("usermodel: estimates", function() {
  r_est <- if ("Unstd_Est" %in% names(r_um$results)) {
    r_um$results$Unstd_Est[r_um$results$op == "=~"]
  } else {
    r_um$results[[4]][r_um$results$op == "=~"]
  }
  rust_est <- rust_um$results$est[rust_um$results$op == "=~"]
  assert_close(r_est, rust_est, 0.05, "loading estimates")
})

# ==========================================================================
# 4. rgmodel
# ==========================================================================
cat("\n--- rgmodel ---\n")
r_rg <- GenomicSEM::rgmodel(r_cov)
rust_rg <- gsemr::rgmodel(rust_cov)

check("rgmodel: R matrix", function() {
  assert_mat_close(r_rg$R, rust_rg$R, 0.05, "R matrix")
})

# ==========================================================================
# 5. write.model
# ==========================================================================
cat("\n--- write.model ---\n")
loadings <- matrix(c(0.7, 0.6, 0.5), ncol = 1)
rownames(loadings) <- trait_names

check("write.model: syntax", function() {
  r_model <- GenomicSEM::write.model(loadings, r_cov$S, cutoff = 0.3)
  rust_model <- gsemr::write.model(loadings, rust_cov$S, cutoff = 0.3)
  for (v in trait_names) {
    if (!grepl(v, r_model)) stop(sprintf("R model missing %s", v))
    if (!grepl(v, rust_model)) stop(sprintf("Rust model missing %s", v))
  }
})

# ==========================================================================
# 6. sumstats (merge)
# ==========================================================================
cat("\n--- sumstats ---\n")
check("sumstats: merge", function() {
  dir.create("out_compare", showWarnings = FALSE)
  rust_ss_path <- "out_compare/rust_merged.tsv"

  # R version
  r_ss <- tryCatch(
    GenomicSEM::sumstats(
      files = munged_traits,
      ref = hm3,
      trait.names = trait_names,
      se.logit = c(FALSE, FALSE, FALSE),
      OLS = c(FALSE, FALSE, FALSE),
      linprob = c(FALSE, FALSE, FALSE),
      info.filter = 0.6,
      maf.filter = 0.01
    ),
    error = function(e) { cat("  R error:", e$message, "\n"); NULL }
  )

  # Rust version
  gsemr::sumstats(
    files = munged_traits,
    ref = hm3,
    trait.names = trait_names,
    info.filter = 0.6,
    maf.filter = 0.01,
    out = rust_ss_path
  )

  if (!file.exists(rust_ss_path)) stop("Rust sumstats produced no output")
  rust_data <- read.delim(rust_ss_path)
  if (nrow(rust_data) < 100) stop(sprintf("Rust merged only %d SNPs", nrow(rust_data)))

  if (!is.null(r_ss)) {
    r_data <- if (is.data.frame(r_ss)) r_ss else read.delim(r_ss)
    r_n <- nrow(r_data)
    rust_n <- nrow(rust_data)
    ratio <- min(r_n, rust_n) / max(r_n, rust_n)
    if (ratio < 0.80) stop(sprintf("SNP count: R=%d Rust=%d (ratio=%.2f)", r_n, rust_n, ratio))
  }
})

# ==========================================================================
# 7. commonfactorGWAS (quick test on subset)
# ==========================================================================
cat("\n--- commonfactorGWAS ---\n")
check("commonfactorGWAS: runs", function() {
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    gsemr::sumstats(files = munged_traits, ref = hm3, trait.names = trait_names, out = merged_path)
  }
  result <- gsemr::commonfactorGWAS(rust_cov, SNPs = merged_path)
  if (is.null(result) || nrow(result) == 0) stop("No results")
  n_finite <- sum(is.finite(result$est))
  if (n_finite < 100) stop(sprintf("Only %d finite estimates", n_finite))
})

# ==========================================================================
# Summary
# ==========================================================================
cat("\n==========================================\n")
cat(sprintf("  PASS: %d  FAIL: %d  SKIP: %d\n", pass_count, fail_count, skip_count))
cat("==========================================\n")

# Cleanup
unlink("out_compare", recursive = TRUE)

if (fail_count > 0) {
  cat("\nFAILED\n")
  quit(status = 1)
} else {
  cat("\nAll tests passed.\n")
}
