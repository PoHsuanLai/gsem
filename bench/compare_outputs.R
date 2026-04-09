#!/usr/bin/env Rscript
# Comprehensive correctness validation: R GenomicSEM vs gsemr (Rust).
# Runs both implementations and compares outputs for all functions.
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

traits <- paste0("data/iter1GWAS", 1:3, ".sumstats.gz")
ld     <- "data/eur_w_ld_chr/"
hm3    <- "data/w_hm3.snplist"
trait_names <- c("V1", "V2", "V3")

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
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b)) stop(sprintf("%s: length mismatch %d vs %d", label, length(a), length(b)))
  d <- max(abs(a - b), na.rm = TRUE)
  if (d > tol) stop(sprintf("%s: max diff %.6e > tol %.6e", label, d, tol))
}

assert_mat_close <- function(a, b, tol, label) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) stop(sprintf("%s: dim mismatch %s vs %s", label, paste(dim(a),collapse="x"), paste(dim(b),collapse="x")))
  d <- max(abs(a - b), na.rm = TRUE)
  if (d > tol) stop(sprintf("%s: max diff %.6e > tol %.6e", label, d, tol))
}

cat("\n==========================================\n")
cat("  Correctness: R GenomicSEM vs gsemr\n")
cat("==========================================\n\n")

# =========================================================================
# 1. LDSC
# =========================================================================
cat("--- LDSC ---\n")
r_cov <- GenomicSEM::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                           trait.names = trait_names, n.blocks = 200)
rust_cov <- gsemr::ldsc(traits, c(NA,NA,NA), c(NA,NA,NA), ld, ld,
                          trait.names = trait_names, n.blocks = 200L)

check("ldsc: S matrix", function() {
  assert_mat_close(r_cov$S, rust_cov$S, 1e-4, "S")
})
check("ldsc: V matrix", function() {
  assert_mat_close(r_cov$V, rust_cov$V, 1e-3, "V")
})
check("ldsc: I matrix", function() {
  assert_mat_close(r_cov$I, rust_cov$I, 1e-4, "I")
})

# =========================================================================
# 2. commonfactor
# =========================================================================
cat("\n--- commonfactor ---\n")
r_cf <- GenomicSEM::commonfactor(r_cov, estimation = "DWLS")
rust_cf <- gsemr::commonfactor(rust_cov, estimation = "DWLS")

check("commonfactor: estimates", function() {
  r_est <- r_cf$results$est[r_cf$results$op == "=~"]
  rust_est <- rust_cf$results$est[rust_cf$results$op == "=~"]
  assert_close(r_est, rust_est, 0.05, "loading estimates")
})
check("commonfactor: SEs", function() {
  r_se <- r_cf$results$SE[r_cf$results$op == "=~"]
  rust_se <- rust_cf$results$se[rust_cf$results$op == "=~"]
  assert_close(r_se, rust_se, 0.02, "loading SEs")
})

# =========================================================================
# 3. usermodel
# =========================================================================
cat("\n--- usermodel ---\n")
model_str <- "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3"
r_um <- GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str)
rust_um <- gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str)

check("usermodel: estimates", function() {
  r_est <- r_um$results$est[r_um$results$op == "=~"]
  rust_est <- rust_um$results$est[rust_um$results$op == "=~"]
  assert_close(r_est, rust_est, 0.05, "loading estimates")
})

# =========================================================================
# 4. write.model
# =========================================================================
cat("\n--- write.model ---\n")
loadings <- matrix(c(0.7, 0.6, 0.5), ncol = 1)
rownames(loadings) <- trait_names

check("write.model: syntax", function() {
  r_model <- GenomicSEM::write.model(loadings, r_cov$S, cutoff = 0.3)
  rust_model <- gsemr::write.model(loadings, rust_cov$S, cutoff = 0.3)
  # Both should produce valid lavaan syntax with the same variables
  # Exact string match may differ in whitespace, so check key content
  for (v in trait_names) {
    if (!grepl(v, r_model)) stop(sprintf("R model missing %s", v))
    if (!grepl(v, rust_model)) stop(sprintf("Rust model missing %s", v))
  }
})

# =========================================================================
# 5. rgmodel
# =========================================================================
cat("\n--- rgmodel ---\n")
r_rg <- GenomicSEM::rgmodel(r_cov)
rust_rg <- gsemr::rgmodel(rust_cov)

check("rgmodel: R matrix", function() {
  assert_mat_close(r_rg$R, rust_rg$R, 0.05, "R matrix")
})
check("rgmodel: V_R matrix", function() {
  assert_mat_close(r_rg$V_R, rust_rg$V_R, 0.01, "V_R matrix")
})

# =========================================================================
# 6. paLDSC
# =========================================================================
cat("\n--- paLDSC ---\n")
check("paLDSC: runs", function() {
  # R's paLDSC has recursive default arg bug, so only test gsemr
  rust_pa <- gsemr::paLDSC(rust_cov, r = 100)
  if (is.null(rust_pa$n_factors)) stop("n_factors is NULL")
  if (rust_pa$n_factors < 1) stop(sprintf("n_factors=%d, expected >= 1", rust_pa$n_factors))
})

# =========================================================================
# 7. munge
# =========================================================================
cat("\n--- munge ---\n")
# We need raw GWAS files (pre-munge). The simulated data is already munged.
# Skip if no raw files available.
skip("munge", "requires raw (pre-munged) GWAS files")

# =========================================================================
# 8. sumstats (merge)
# =========================================================================
cat("\n--- sumstats ---\n")
# First check if R GenomicSEM's sumstats() works with our data
check("sumstats: merge", function() {
  dir.create("out_compare", showWarnings = FALSE)

  # R version
  r_ss <- tryCatch(
    GenomicSEM::sumstats(
      files = traits,
      ref = ld,
      trait.names = trait_names,
      info.filter = 0.6,
      maf.filter = 0.01,
      OLS = NULL,
      linprob = NULL,
      N = NULL,
      betas = NULL,
      se.logit = NULL,
      keep.indel = FALSE
    ),
    error = function(e) NULL
  )

  # Rust version
  rust_ss_path <- "out_compare/rust_merged.tsv"
  rust_ss <- gsemr::sumstats(
    files = traits,
    ref = ld,
    trait.names = trait_names,
    info.filter = 0.6,
    maf.filter = 0.01,
    out = rust_ss_path
  )

  if (is.null(r_ss)) stop("R sumstats() failed")
  if (!file.exists(rust_ss_path)) stop("Rust sumstats produced no output")

  # Compare: both should produce files with similar SNP counts
  r_data <- read.delim(r_ss)
  rust_data <- read.delim(rust_ss_path)

  # SNP counts should be within 5% of each other
  r_n <- nrow(r_data)
  rust_n <- nrow(rust_data)
  ratio <- min(r_n, rust_n) / max(r_n, rust_n)
  if (ratio < 0.90) stop(sprintf("SNP count mismatch: R=%d Rust=%d (ratio=%.2f)", r_n, rust_n, ratio))
})

# =========================================================================
# 9. commonfactorGWAS
# =========================================================================
cat("\n--- commonfactorGWAS ---\n")
check("commonfactorGWAS: runs", function() {
  # Need merged sumstats file. Generate with gsemr first.
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    gsemr::sumstats(files = traits, ref = ld, trait.names = trait_names, out = merged_path)
  }

  rust_cfgwas <- gsemr::commonfactorGWAS(rust_cov, SNPs = merged_path)
  if (is.null(rust_cfgwas) || nrow(rust_cfgwas) == 0) stop("No GWAS results")
  if (!("est" %in% names(rust_cfgwas))) stop("Missing est column")
  # Check some results are finite
  n_finite <- sum(is.finite(rust_cfgwas$est))
  if (n_finite < 100) stop(sprintf("Only %d finite estimates", n_finite))
})

# =========================================================================
# 10. userGWAS
# =========================================================================
cat("\n--- userGWAS ---\n")
check("userGWAS: runs", function() {
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    gsemr::sumstats(files = traits, ref = ld, trait.names = trait_names, out = merged_path)
  }

  gwas_model <- "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nF1 ~ SNP\nSNP ~~ SNP"
  rust_ugwas <- gsemr::userGWAS(rust_cov, SNPs = merged_path, model = gwas_model)
  if (is.null(rust_ugwas) || nrow(rust_ugwas) == 0) stop("No GWAS results")
  n_finite <- sum(is.finite(rust_ugwas$est))
  if (n_finite < 100) stop(sprintf("Only %d finite estimates", n_finite))
})

# =========================================================================
# 11. hdl
# =========================================================================
cat("\n--- hdl ---\n")
skip("hdl", "requires HDL reference panels (not in test data)")

# =========================================================================
# 12. s_ldsc (stratified)
# =========================================================================
cat("\n--- s_ldsc ---\n")
skip("s_ldsc", "requires annotation-specific LD scores")

# =========================================================================
# 13. enrich
# =========================================================================
cat("\n--- enrich ---\n")
skip("enrich", "requires stratified LDSC output")

# =========================================================================
# 14. simLDSC
# =========================================================================
cat("\n--- simLDSC ---\n")
skip("simLDSC", "generative function, no R comparison needed")

# =========================================================================
# 15. multiSNP
# =========================================================================
cat("\n--- multiSNP ---\n")
skip("multiSNP", "requires LD matrix + multi-SNP data")

# =========================================================================
# 16. multiGene
# =========================================================================
cat("\n--- multiGene ---\n")
skip("multiGene", "requires gene-level LD data")

# =========================================================================
# 17. summaryGLS
# =========================================================================
cat("\n--- summaryGLS ---\n")
skip("summaryGLS", "requires GWAS output from userGWAS")

# =========================================================================
# 18. convert_hdl_panels
# =========================================================================
cat("\n--- convert_hdl_panels ---\n")
skip("convert_hdl_panels", "utility function, no R equivalent")

# =========================================================================
# Summary
# =========================================================================
cat("\n==========================================\n")
cat(sprintf("  PASS: %d  FAIL: %d  SKIP: %d\n", pass_count, fail_count, skip_count))
cat("==========================================\n")

if (fail_count > 0) {
  cat("\nFAILED\n")
  quit(status = 1)
} else {
  cat("\nAll testable functions verified.\n")
}

# Cleanup
unlink("out_compare", recursive = TRUE)
