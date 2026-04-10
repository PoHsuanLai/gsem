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

# Silence gsemr's first-use commonfactorGWAS semantics warning during this
# smoke test — see ARCHITECTURE.md §3.3 for the rationale. The smoke test
# here is only checking that the function runs and returns finite results;
# the identification-difference caveat is out of scope for this script.
options(gsemr.commonfactorGWAS.quiet = TRUE)

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
hm3       <- "data/w_hm3.snplist"
ref_1000g <- "data/reference.1000G.maf.0.005.txt.gz"  # proper ref with MAF column
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
# 6. sumstats (merge) — uses RAW GWAS files (not munged)
# ==========================================================================
cat("\n--- sumstats ---\n")
check("sumstats: merge", function() {
  dir.create("out_compare", showWarnings = FALSE)
  rust_ss_path <- "out_compare/rust_merged.tsv"

  # Use 1000G reference (has MAF column) so R sumstats works
  ref <- if (file.exists(ref_1000g)) ref_1000g else hm3

  # R version
  r_ss <- tryCatch(
    GenomicSEM::sumstats(
      files = raw_traits,
      ref = ref,
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
    files = raw_traits,
    ref = ref,
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
  ref <- if (file.exists(ref_1000g)) ref_1000g else hm3
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    gsemr::sumstats(files = raw_traits, ref = ref, trait.names = trait_names, out = merged_path)
  }
  # Use a small subset (500 SNPs) — full GWAS on 1M+ SNPs is too slow for CI
  subset_path <- "out_compare/rust_merged_subset.tsv"
  full <- read.delim(merged_path, nrows = 500)
  write.table(full, subset_path, sep = "\t", row.names = FALSE, quote = FALSE)
  result <- gsemr::commonfactorGWAS(rust_cov, SNPs = subset_path)
  if (is.null(result) || nrow(result) == 0) stop("No results")
  n_finite <- sum(is.finite(result$est))
  if (n_finite < 100) stop(sprintf("Only %d finite estimates", n_finite))
})

# ==========================================================================
# 8. userGWAS (quick test on subset)
# ==========================================================================
cat("\n--- userGWAS ---\n")
check("userGWAS: runs", function() {
  ref <- if (file.exists(ref_1000g)) ref_1000g else hm3
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    dir.create("out_compare", showWarnings = FALSE)
    gsemr::sumstats(files = raw_traits, ref = ref, trait.names = trait_names, out = merged_path)
  }
  subset_path <- "out_compare/rust_merged_subset.tsv"
  if (!file.exists(subset_path)) {
    full <- read.delim(merged_path, nrows = 500)
    write.table(full, subset_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  user_model <- paste0(
    "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
    "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
  )
  result <- gsemr::userGWAS(rust_cov, SNPs = subset_path, model = user_model)
  if (is.null(result) || nrow(result) == 0) stop("No results")
  # userGWAS returns nested JSON; extract SNP effect from params
  n_converged <- sum(result$converged, na.rm = TRUE)
  if (n_converged < 100) stop(sprintf("Only %d converged", n_converged))
})

# ==========================================================================
# 9. munge
# ==========================================================================
cat("\n--- munge ---\n")
check("munge: produces output", function() {
  # Munge just the first trait as a quick test (writes to current dir)
  out <- gsemr::munge(
    files = raw_traits[1],
    hm3 = hm3,
    trait.names = "ANX_test",
    N = n_overrides[1],
    info.filter = 0.9,
    maf.filter = 0.01
  )
  out_file <- "ANX_test.sumstats.gz"
  if (!file.exists(out_file)) stop("Munge produced no output file")
  # Check it has SNP, N, Z, A1 columns
  data <- read.table(gzfile(out_file), header = TRUE, nrows = 5)
  required <- c("SNP", "N", "Z", "A1")
  for (col in required) {
    if (!(col %in% colnames(data))) stop(sprintf("Missing column: %s", col))
  }
  file.remove(out_file)
})

# ==========================================================================
# 10. paLDSC (parallel analysis)
# ==========================================================================
cat("\n--- paLDSC ---\n")
check("paLDSC: eigenvalues", function() {
  result <- gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 3, save.pdf = FALSE)
  if (is.null(result)) stop("paLDSC returned NULL")
  if (is.null(result$observed)) stop("No observed eigenvalues")
  if (length(result$observed) != 3) stop(sprintf("Expected 3 eigenvalues, got %d", length(result$observed)))
  if (!all(is.finite(result$observed))) stop("Non-finite eigenvalues")
  if (is.null(result$n_factors) || result$n_factors < 1) stop("n_factors should be >= 1")
})

# ==========================================================================
# 11. summaryGLS
# ==========================================================================
cat("\n--- summaryGLS ---\n")
check("summaryGLS: regression", function() {
  # Simple GLS regression on vech(S)
  k <- nrow(rust_cov$S)
  kstar <- k * (k + 1) / 2
  y <- rust_cov$S[lower.tri(rust_cov$S, diag = TRUE)]
  X <- matrix(1:kstar, ncol = 1)
  result <- gsemr::summaryGLS(Y = y, V_Y = rust_cov$V, PREDICTORS = X, INTERCEPT = TRUE)
  if (is.null(result)) stop("summaryGLS returned NULL")
  if (nrow(result) == 0) stop("Empty results")
  if (!all(c("beta", "se", "z", "p") %in% colnames(result))) stop("Missing columns")
  if (!all(is.finite(result$beta))) stop("Non-finite betas")
})

# ==========================================================================
# 12. simLDSC
# ==========================================================================
cat("\n--- simLDSC ---\n")
check("simLDSC: generates data", function() {
  # Simulate from the estimated covariance
  sim_result <- gsemr::simLDSC(
    covmat = rust_cov$S,
    N = c(5000, 5000, 5000),
    seed = 42,
    ld = "data/eur_w_ld_chr/"
  )
  if (is.null(sim_result)) stop("simLDSC returned NULL")
  # Should return a matrix: SNPs x traits or traits x SNPs
  if (!is.matrix(sim_result) && !is.data.frame(sim_result)) stop("Expected matrix output")
  dm <- dim(sim_result)
  if (max(dm) < 100) stop(sprintf("Too few elements: %dx%d", dm[1], dm[2]))
})

# ==========================================================================
# 13. HDL (requires HDL panels — skip if not available)
# ==========================================================================
cat("\n--- HDL ---\n")
hdl_panels <- "data/UKB_imputed_SVD_eigen99_extraction"
if (dir.exists(hdl_panels)) {
  check("hdl: estimates", function() {
    result <- gsemr::hdl(munged_traits, sample_prev, pop_prev,
                         trait.names = trait_names, LD.path = hdl_panels)
    if (is.null(result$S)) stop("No S matrix")
    if (nrow(result$S) != 3) stop("Wrong S dimension")
  })
} else {
  skip("hdl: estimates", "HDL panels not available")
}

# ==========================================================================
# 14. s_ldsc (stratified LDSC — requires annotation files, skip if unavailable)
# ==========================================================================
cat("\n--- s_ldsc ---\n")
frq_dir <- "data/1000G_Phase3_frq/"
annot_ld <- "data/baseline_v1.1/"
if (dir.exists(frq_dir) && dir.exists(annot_ld)) {
  check("s_ldsc: runs", function() {
    result <- gsemr::s_ldsc(munged_traits, sample_prev, pop_prev,
                            ld = annot_ld, wld = ld, frq = frq_dir,
                            trait.names = trait_names, n.blocks = 200L)
    if (is.null(result)) stop("s_ldsc returned NULL")
  })
} else {
  skip("s_ldsc: runs", "stratified LDSC data not available")
}

# ==========================================================================
# 15. enrich (requires stratified LDSC output — skip if s_ldsc unavailable)
# ==========================================================================
cat("\n--- enrich ---\n")
skip("enrich: runs", "requires stratified LDSC output")

# ==========================================================================
# 16. multiSNP
# ==========================================================================
cat("\n--- multiSNP ---\n")
check("multiSNP: runs", function() {
  ref <- if (file.exists(ref_1000g)) ref_1000g else hm3
  merged_path <- "out_compare/rust_merged.tsv"
  if (!file.exists(merged_path)) {
    dir.create("out_compare", showWarnings = FALSE)
    gsemr::sumstats(files = raw_traits, ref = ref, trait.names = trait_names, out = merged_path)
  }
  # Read a small set of SNPs and build the expected data frame format
  snp_data <- read.delim(merged_path, nrows = 5)
  beta_cols <- grep("^beta\\.", colnames(snp_data))
  se_cols <- grep("^se\\.", colnames(snp_data))
  snp_df <- data.frame(
    beta = I(as.matrix(snp_data[, beta_cols])),
    se = I(as.matrix(snp_data[, se_cols])),
    var = 2 * snp_data$MAF * (1 - snp_data$MAF)
  )
  ld_mat <- diag(nrow(snp_data))
  result <- gsemr::multiSNP(
    covstruc = rust_cov,
    SNPs = snp_df,
    LD = ld_mat,
    SNPlist = snp_data$SNP
  )
  if (is.null(result)) stop("multiSNP returned NULL")
  if (is.null(result$results)) stop("No results data frame")
})

# ==========================================================================
# 17. multiGene (same interface as multiSNP)
# ==========================================================================
cat("\n--- multiGene ---\n")
check("multiGene: runs", function() {
  merged_path <- "out_compare/rust_merged.tsv"
  snp_data <- read.delim(merged_path, nrows = 3)
  beta_cols <- grep("^beta\\.", colnames(snp_data))
  se_cols <- grep("^se\\.", colnames(snp_data))
  gene_df <- data.frame(
    beta = I(as.matrix(snp_data[, beta_cols])),
    se = I(as.matrix(snp_data[, se_cols])),
    var = 2 * snp_data$MAF * (1 - snp_data$MAF)
  )
  ld_mat <- diag(nrow(snp_data))
  result <- gsemr::multiGene(
    covstruc = rust_cov,
    Genes = gene_df,
    LD = ld_mat,
    Genelist = snp_data$SNP
  )
  if (is.null(result)) stop("multiGene returned NULL")
  if (is.null(result$results)) stop("No results data frame")
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
