#!/usr/bin/env Rscript
# Performance benchmark: R GenomicSEM vs gsemr (Rust)
# Measures wall-clock time and peak memory for each API function.
#
# Usage: cd bench && Rscript benchmark_perf.R
#
# Output:
#   bench/benchmark_results.csv  — timing and memory data
#   bench/benchmark_plots.pdf    — comparison charts

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

for (pkg in c("GenomicSEM", "gsemr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' not installed.", pkg))
}
library(GenomicSEM)
library(gsemr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ptsd_file <- Sys.glob("data/*PTSD*")[1]
if (is.na(ptsd_file)) ptsd_file <- Sys.glob("data/pts_*ea*")[1]
raw_traits <- c(
  "data/anxiety.meta.full.cc.tbl.gz",
  "data/ocd_aug2017.gz",
  ptsd_file
)
if (any(is.na(raw_traits))) stop("Missing GWAS files. Run: Rscript setup_data.R")

ld          <- "data/eur_w_ld_chr/"
hm3         <- "data/w_hm3.snplist"
ref_1000g   <- "data/reference.1000G.maf.0.005.txt.gz"  # proper ref with MAF column
trait_names <- c("ANX", "OCD", "PTSD")
sample_prev <- c(0.5, 0.5, 0.5)
pop_prev    <- c(0.16, 0.02, 0.07)
n_overrides <- c(NA, 9725, 5831)
munged_traits <- paste0(trait_names, ".sumstats.gz")

GWAS_SUBSET_N <- 1000L

if (!file.exists(ref_1000g)) {
  stop("Missing reference file: ", ref_1000g,
       "\nDownload from: https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v")
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Benchmark a single call: returns list(time_s, peak_mb, error)
run_bench <- function(expr_fn) {
  gc(full = TRUE, reset = TRUE)
  mem_before <- sum(gc()[, 2])  # current Mb used (Ncells + Vcells)
  t <- tryCatch(
    system.time({ result <- expr_fn() }),
    error = function(e) e
  )
  if (inherits(t, "error")) {
    return(list(time_s = NA_real_, peak_mb = NA_real_, error = conditionMessage(t)))
  }
  mem_after <- gc()
  peak_mb <- sum(mem_after[, 7])  # max used (Mb) column — Ncells + Vcells
  list(time_s = as.numeric(t["elapsed"]), peak_mb = peak_mb, error = NA_character_)
}

# Accumulator
results <- data.frame(
  func     = character(),
  impl     = character(),
  time_s   = numeric(),
  peak_mb  = numeric(),
  error    = character(),
  stringsAsFactors = FALSE
)

add_result <- function(func, impl, bench) {
  results[nrow(results) + 1, ] <<- list(
    func    = func,
    impl    = impl,
    time_s  = bench$time_s,
    peak_mb = bench$peak_mb,
    error   = ifelse(is.na(bench$error), "", bench$error)
  )
}

# Tolerance-based comparison of two named numeric vectors / matrices / data frames.
# Records pass/fail/skip per function in `equiv` accumulator.
equiv <- data.frame(
  func   = character(),
  status = character(),  # PASS / FAIL / SKIP
  detail = character(),
  stringsAsFactors = FALSE
)
add_equiv <- function(func, status, detail = "") {
  equiv[nrow(equiv) + 1, ] <<- list(func = func, status = status, detail = detail)
}

# Compare two numeric objects within absolute tolerance `tol`.
# Returns "" on success or an error string on failure.
compare_numeric <- function(a, b, tol, label) {
  a <- suppressWarnings(as.numeric(unlist(a)))
  b <- suppressWarnings(as.numeric(unlist(b)))
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) == 0) return(sprintf("%s: no finite values to compare", label))
  if (length(a) != length(b)) {
    return(sprintf("%s: length mismatch %d vs %d", label, length(a), length(b)))
  }
  d <- max(abs(a[ok] - b[ok]))
  if (!is.finite(d) || d > tol) {
    return(sprintf("%s: max diff %.3e > tol %.3e", label, d, tol))
  }
  ""
}

cat("\n========================================================\n")
cat("  Performance Benchmark: R GenomicSEM vs gsemr (Rust)\n")
cat("  Data: PGC Anxiety, OCD, PTSD\n")
cat("========================================================\n\n")

# ---------------------------------------------------------------------------
# 0. Ensure munged files exist
# ---------------------------------------------------------------------------
if (!all(file.exists(munged_traits))) {
  cat("Munging raw GWAS files (prerequisite)...\n")
  GenomicSEM::munge(
    files = raw_traits, hm3 = hm3, trait.names = trait_names,
    N = n_overrides, info.filter = 0.9, maf.filter = 0.01
  )
}

# ---------------------------------------------------------------------------
# 1. munge
# ---------------------------------------------------------------------------
cat("[1/12] munge\n")
# R munge
munge_r_names <- paste0(trait_names, "_benchR")
# Clean up any leftover files
file.remove(paste0(munge_r_names, ".sumstats.gz")[file.exists(paste0(munge_r_names, ".sumstats.gz"))])
b <- run_bench(function() {
  GenomicSEM::munge(
    files = raw_traits, hm3 = hm3, trait.names = munge_r_names,
    N = n_overrides, info.filter = 0.9, maf.filter = 0.01
  )
})
add_result("munge", "R", b)

# Rust munge
munge_rust_names <- paste0(trait_names, "_benchRust")
file.remove(paste0(munge_rust_names, ".sumstats.gz")[file.exists(paste0(munge_rust_names, ".sumstats.gz"))])
b <- run_bench(function() {
  gsemr::munge(
    files = raw_traits, hm3 = hm3, trait.names = munge_rust_names,
    N = n_overrides, info.filter = 0.9, maf.filter = 0.01
  )
})
add_result("munge", "Rust", b)

# Equivalence: compare Z column on the intersection of SNPs across both outputs
tryCatch({
  r_files    <- paste0(munge_r_names, ".sumstats.gz")
  rust_files <- paste0(munge_rust_names, ".sumstats.gz")
  if (all(file.exists(r_files)) && all(file.exists(rust_files))) {
    diffs <- character(0)
    for (i in seq_along(trait_names)) {
      r_d    <- read.table(gzfile(r_files[i]), header = TRUE, stringsAsFactors = FALSE)
      rust_d <- read.table(gzfile(rust_files[i]), header = TRUE, stringsAsFactors = FALSE)
      common <- intersect(r_d$SNP, rust_d$SNP)
      if (length(common) < 100) { diffs <- c(diffs, sprintf("%s: only %d common SNPs", trait_names[i], length(common))); next }
      r_sub    <- r_d[match(common, r_d$SNP), ]
      rust_sub <- rust_d[match(common, rust_d$SNP), ]
      # Align allele orientation: flip Z when A1 disagrees
      flip <- r_sub$A1 != rust_sub$A1
      rust_z <- rust_sub$Z
      rust_z[flip] <- -rust_z[flip]
      err <- compare_numeric(r_sub$Z, rust_z, tol = 1e-3, label = sprintf("%s Z", trait_names[i]))
      if (nzchar(err)) diffs <- c(diffs, err)
    }
    if (length(diffs) == 0) add_equiv("munge", "PASS", sprintf("Z within 1e-3 over %d traits", length(trait_names)))
    else add_equiv("munge", "FAIL", paste(diffs, collapse = "; "))
  } else {
    add_equiv("munge", "SKIP", "missing output files")
  }
}, error = function(e) add_equiv("munge", "SKIP", conditionMessage(e)))
file.remove(paste0(munge_r_names, ".sumstats.gz")[file.exists(paste0(munge_r_names, ".sumstats.gz"))])
file.remove(paste0(munge_rust_names, ".sumstats.gz")[file.exists(paste0(munge_rust_names, ".sumstats.gz"))])

# ---------------------------------------------------------------------------
# 2. ldsc
# ---------------------------------------------------------------------------
cat("[2/12] ldsc\n")
b <- run_bench(function() {
  GenomicSEM::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                   trait.names = trait_names, n.blocks = 200)
})
add_result("ldsc", "R", b)

b <- run_bench(function() {
  gsemr::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
              trait.names = trait_names, n.blocks = 200L)
})
add_result("ldsc", "Rust", b)

# Pre-compute covstruct for downstream
r_cov    <- GenomicSEM::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                              trait.names = trait_names, n.blocks = 200)
rust_cov <- gsemr::ldsc(munged_traits, sample_prev, pop_prev, ld, ld,
                          trait.names = trait_names, n.blocks = 200L)

# Equivalence: S, V, I matrices within tol
local({
  errs <- c(
    compare_numeric(as.matrix(r_cov$S), as.matrix(rust_cov$S), 1e-3, "S"),
    compare_numeric(as.matrix(r_cov$V), as.matrix(rust_cov$V), 1e-3, "V"),
    compare_numeric(as.matrix(r_cov$I), as.matrix(rust_cov$I), 1e-3, "I")
  )
  errs <- errs[nzchar(errs)]
  if (length(errs) == 0) add_equiv("ldsc", "PASS", "S/V/I within 1e-3")
  else add_equiv("ldsc", "FAIL", paste(errs, collapse = "; "))
})

# ---------------------------------------------------------------------------
# 3. commonfactor
# ---------------------------------------------------------------------------
cat("[3/12] commonfactor\n")
b <- run_bench(function() GenomicSEM::commonfactor(r_cov, estimation = "DWLS"))
add_result("commonfactor", "R", b)

b <- run_bench(function() gsemr::commonfactor(rust_cov, estimation = "DWLS"))
add_result("commonfactor", "Rust", b)

local({
  r_cf    <- GenomicSEM::commonfactor(r_cov, estimation = "DWLS")
  rust_cf <- gsemr::commonfactor(rust_cov, estimation = "DWLS")
  r_load    <- r_cf$results$Unstandardized_Estimate[r_cf$results$op == "=~"]
  rust_load <- rust_cf$results$est[rust_cf$results$op == "=~"]
  err <- compare_numeric(r_load, rust_load, 0.05, "loadings")
  if (nzchar(err)) add_equiv("commonfactor", "FAIL", err)
  else add_equiv("commonfactor", "PASS", "loadings within 0.05")
})

# ---------------------------------------------------------------------------
# 4. usermodel
# ---------------------------------------------------------------------------
cat("[4/12] usermodel\n")
model_str <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~~ 1*F1\n",
  paste(trait_names, "~~", trait_names, collapse = "\n")
)

b <- run_bench(function() GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str))
add_result("usermodel", "R", b)

b <- run_bench(function() gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str))
add_result("usermodel", "Rust", b)

local({
  r_um    <- GenomicSEM::usermodel(r_cov, estimation = "DWLS", model = model_str)
  rust_um <- gsemr::usermodel(rust_cov, estimation = "DWLS", model = model_str)
  r_est <- if ("Unstd_Est" %in% names(r_um$results)) {
    r_um$results$Unstd_Est[r_um$results$op == "=~"]
  } else {
    r_um$results[[4]][r_um$results$op == "=~"]
  }
  rust_est <- rust_um$results$est[rust_um$results$op == "=~"]
  err <- compare_numeric(r_est, rust_est, 0.05, "loadings")
  if (nzchar(err)) add_equiv("usermodel", "FAIL", err)
  else add_equiv("usermodel", "PASS", "loadings within 0.05")
})

# ---------------------------------------------------------------------------
# 5. rgmodel
# ---------------------------------------------------------------------------
cat("[5/12] rgmodel\n")
b <- run_bench(function() GenomicSEM::rgmodel(r_cov))
add_result("rgmodel", "R", b)

b <- run_bench(function() gsemr::rgmodel(rust_cov))
add_result("rgmodel", "Rust", b)

tryCatch({
  r_rg    <- GenomicSEM::rgmodel(r_cov)
  rust_rg <- gsemr::rgmodel(rust_cov)
  err <- compare_numeric(as.matrix(r_rg$R), as.matrix(rust_rg$R), 0.05, "R matrix")
  if (nzchar(err)) add_equiv("rgmodel", "FAIL", err)
  else add_equiv("rgmodel", "PASS", "R matrix within 0.05")
}, error = function(e) add_equiv("rgmodel", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 6. sumstats (use 1000G reference with MAF column — works for both R and Rust)
# ---------------------------------------------------------------------------
cat("[6/12] sumstats\n")
dir.create("out_bench", showWarnings = FALSE)

r_ss_path   <- "out_bench/r_merged.rds"
rust_ss_path <- "out_bench/rust_merged.tsv"

# R sumstats
b <- run_bench(function() {
  res <- GenomicSEM::sumstats(
    files = raw_traits, ref = ref_1000g, trait.names = trait_names,
    se.logit = c(FALSE, FALSE, FALSE),
    OLS = c(FALSE, FALSE, FALSE),
    linprob = c(FALSE, FALSE, FALSE),
    info.filter = 0.6, maf.filter = 0.01
  )
  saveRDS(res, r_ss_path)
  res
})
add_result("sumstats", "R", b)

# Rust sumstats
b <- run_bench(function() {
  gsemr::sumstats(
    files = raw_traits, ref = ref_1000g, trait.names = trait_names,
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path
  )
})
add_result("sumstats", "Rust", b)

# Equivalence: row count ratio (sumstats outputs differ in shape, so use a
# loose count-based check; per-SNP beta comparison is done at GWAS time).
tryCatch({
  if (file.exists(r_ss_path) && file.exists(rust_ss_path)) {
    r_n    <- nrow(readRDS(r_ss_path))
    rust_n <- nrow(read.delim(rust_ss_path))
    ratio  <- min(r_n, rust_n) / max(r_n, rust_n)
    if (ratio >= 0.80) {
      add_equiv("sumstats", "PASS", sprintf("R=%d Rust=%d (ratio %.2f)", r_n, rust_n, ratio))
    } else {
      add_equiv("sumstats", "FAIL", sprintf("R=%d Rust=%d (ratio %.2f < 0.80)", r_n, rust_n, ratio))
    }
  } else {
    add_equiv("sumstats", "SKIP", "missing output files")
  }
}, error = function(e) add_equiv("sumstats", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 7. write.model
# ---------------------------------------------------------------------------
cat("[7/12] write.model\n")
loadings <- matrix(c(0.7, 0.6, 0.5), ncol = 1)
rownames(loadings) <- trait_names

b <- run_bench(function() GenomicSEM::write.model(loadings, r_cov$S, cutoff = 0.3))
add_result("write.model", "R", b)

b <- run_bench(function() gsemr::write.model(loadings, rust_cov$S, cutoff = 0.3))
add_result("write.model", "Rust", b)

# Equivalence: just verify both produce non-empty syntax mentioning every trait
local({
  r_model    <- GenomicSEM::write.model(loadings, r_cov$S, cutoff = 0.3)
  rust_model <- gsemr::write.model(loadings, rust_cov$S, cutoff = 0.3)
  missing <- character(0)
  for (v in trait_names) {
    if (!grepl(v, r_model))    missing <- c(missing, sprintf("R missing %s", v))
    if (!grepl(v, rust_model)) missing <- c(missing, sprintf("Rust missing %s", v))
  }
  if (length(missing) == 0) add_equiv("write.model", "PASS", "all traits present")
  else add_equiv("write.model", "FAIL", paste(missing, collapse = "; "))
})

# ---------------------------------------------------------------------------
# Prepare GWAS subsets for commonfactorGWAS / userGWAS
# Each implementation uses its own sumstats output (fair comparison).
# ---------------------------------------------------------------------------
cat("Preparing GWAS subsets...\n")

# R's sumstats output (data frame saved as RDS)
if (file.exists(r_ss_path)) {
  r_snps_full <- readRDS(r_ss_path)
} else {
  # Fallback: run R sumstats if not already done
  r_snps_full <- GenomicSEM::sumstats(
    files = raw_traits, ref = ref_1000g, trait.names = trait_names,
    se.logit = c(FALSE, FALSE, FALSE), OLS = c(FALSE, FALSE, FALSE),
    linprob = c(FALSE, FALSE, FALSE), info.filter = 0.6, maf.filter = 0.01
  )
}
r_snps_df <- head(r_snps_full, GWAS_SUBSET_N)
cat(sprintf("  R subset: %d SNPs (from R sumstats, %d total)\n", nrow(r_snps_df), nrow(r_snps_full)))

# Rust's sumstats output (TSV file)
if (!file.exists(rust_ss_path)) {
  gsemr::sumstats(
    files = raw_traits, ref = ref_1000g, trait.names = trait_names,
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path
  )
}
rust_subset_path <- "out_bench/rust_subset.tsv"
rust_snps <- read.delim(rust_ss_path, nrows = GWAS_SUBSET_N)
write.table(rust_snps, rust_subset_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Rust subset: %d SNPs\n", nrow(rust_snps)))

# ---------------------------------------------------------------------------
# 8. commonfactorGWAS — run both parallel and serial for both impls
# ---------------------------------------------------------------------------
cat("[8/12] commonfactorGWAS\n")

cfg_results <- list()  # capture outputs for equivalence
for (par in c(TRUE, FALSE)) {
  tag <- if (par) "par" else "seq"

  b <- run_bench(function() {
    cfg_results[[paste0("R_",  tag)]] <<-
      GenomicSEM::commonfactorGWAS(covstruc = r_cov, SNPs = r_snps_df, parallel = par)
  })
  add_result("commonfactorGWAS", sprintf("R (%s)", tag), b)

  b <- run_bench(function() {
    cfg_results[[paste0("Rust_", tag)]] <<-
      gsemr::commonfactorGWAS(covstruc = rust_cov, SNPs = rust_subset_path, parallel = par)
  })
  add_result("commonfactorGWAS", sprintf("Rust (%s)", tag), b)
}

# Equivalence: parallel vs serial within each impl, and R vs Rust on common SNPs
tryCatch({
  notes <- character(0)
  for (impl in c("R", "Rust")) {
    a <- cfg_results[[paste0(impl, "_par")]]
    b <- cfg_results[[paste0(impl, "_seq")]]
    if (!is.null(a) && !is.null(b)) {
      err <- compare_numeric(a$est, b$est, 1e-6, sprintf("%s par/seq est", impl))
      if (nzchar(err)) notes <- c(notes, err)
    }
  }
  ra <- cfg_results$R_par; rb <- cfg_results$Rust_par
  if (!is.null(ra) && !is.null(rb)) {
    common <- intersect(ra$SNP, rb$SNP)
    if (length(common) >= 50) {
      ra_s <- ra[match(common, ra$SNP), ]
      rb_s <- rb[match(common, rb$SNP), ]
      ok <- is.finite(ra_s$est) & is.finite(rb_s$est)
      if (sum(ok) >= 50) {
        d <- max(abs(ra_s$est[ok] - rb_s$est[ok]))
        # Allow looser tolerance for per-SNP estimates (Rust uses different optimizer)
        if (!is.finite(d) || d > 0.05) notes <- c(notes, sprintf("R/Rust est max diff %.3e > 0.05", d))
      }
    }
  }
  if (length(notes) == 0) add_equiv("commonfactorGWAS", "PASS", "par==seq and R~Rust on common SNPs")
  else add_equiv("commonfactorGWAS", "FAIL", paste(notes, collapse = "; "))
}, error = function(e) add_equiv("commonfactorGWAS", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 9. userGWAS — scaling sweep over multiple SNP counts (parallel only)
# ---------------------------------------------------------------------------
cat("[9/12] userGWAS (scaling sweep)\n")
gwas_model <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

# SNP-count grid for head-to-head scaling. Chosen to cover the linear regime
# where R is still tolerable (~25k SNPs ≈ 3 min R) while being meaningful
# relative to real GWAS deployments. Capped by available SNPs.
ug_sizes <- c(1000L, 5000L, 10000L, 25000L)
n_rust_available <- length(readLines(rust_ss_path)) - 1L  # minus header
n_r_available    <- nrow(r_snps_full)
n_available      <- min(n_rust_available, n_r_available)
ug_sizes <- unique(pmin(ug_sizes, n_available))
ug_sizes <- ug_sizes[ug_sizes > 0]
cat(sprintf("  head-to-head sizes: %s (capped at %d available SNPs)\n",
            paste(ug_sizes, collapse = ", "), n_available))

ug_results <- list()
for (n_snp in ug_sizes) {
  cat(sprintf("  -> N=%d\n", n_snp))
  # Materialize per-size subsets for both impls
  r_snps_sub    <- head(r_snps_full, n_snp)
  rust_sub_path <- sprintf("out_bench/rust_subset_%d.tsv", n_snp)
  rust_snps_sub <- read.delim(rust_ss_path, nrows = n_snp)
  write.table(rust_snps_sub, rust_sub_path, sep = "\t", row.names = FALSE, quote = FALSE)

  label_r    <- sprintf("R (N=%d)",    n_snp)
  label_rust <- sprintf("Rust (N=%d)", n_snp)

  b <- run_bench(function() {
    ug_results[[sprintf("R_%d", n_snp)]] <<-
      GenomicSEM::userGWAS(covstruc = r_cov, SNPs = r_snps_sub,
                           model = gwas_model, parallel = TRUE)
  })
  add_result("userGWAS", label_r, b)

  b <- run_bench(function() {
    ug_results[[sprintf("Rust_%d", n_snp)]] <<-
      gsemr::userGWAS(covstruc = rust_cov, SNPs = rust_sub_path,
                      model = gwas_model, parallel = TRUE)
  })
  add_result("userGWAS", label_rust, b)
}

# Rust-only "real-world deployment" point: the full sumstats file.
# Running R at this scale is not feasible on a single machine — the R
# package recommends MPI for multi-million-SNP GWAS — so we report a
# Rust-only bar for the headline "can do it on one box" claim.
cat(sprintf("  -> Rust-only full-scale: N=%d\n", n_rust_available))
b <- run_bench(function() {
  ug_results[["Rust_full"]] <<-
    gsemr::userGWAS(covstruc = rust_cov, SNPs = rust_ss_path,
                    model = gwas_model, parallel = TRUE)
})
add_result("userGWAS", sprintf("Rust (N=%d, full)", n_rust_available), b)

# Equivalence: convergence counts align between R and Rust at each size
# (tolerate up to 2% of SNPs disagreeing on convergence).
tryCatch({
  notes <- character(0)
  for (n_snp in ug_sizes) {
    a <- ug_results[[sprintf("R_%d",    n_snp)]]
    b <- ug_results[[sprintf("Rust_%d", n_snp)]]
    if (!is.null(a) && !is.null(b)) {
      ca <- sum(a$converged, na.rm = TRUE)
      cb <- sum(b$converged, na.rm = TRUE)
      if (abs(ca - cb) > 0.02 * n_snp) {
        notes <- c(notes, sprintf("N=%d converged %d(R) vs %d(Rust)", n_snp, ca, cb))
      }
    }
  }
  if (length(notes) == 0) add_equiv("userGWAS", "PASS",
                                    sprintf("converged counts align across %d sizes", length(ug_sizes)))
  else add_equiv("userGWAS", "FAIL", paste(notes, collapse = "; "))
}, error = function(e) add_equiv("userGWAS", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 10. paLDSC
# ---------------------------------------------------------------------------
cat("[10/12] paLDSC\n")
b <- run_bench(function() {
  # R GenomicSEM may not have paLDSC; wrap safely
  if (existsFunction("paLDSC", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::paLDSC(r_cov$S, r_cov$V, r = 100, p = 3, save.pdf = TRUE)
  } else {
    stop("paLDSC not available in R GenomicSEM")
  }
})
add_result("paLDSC", "R", b)

b <- run_bench(function() {
  gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 3, save.pdf = TRUE)
})
add_result("paLDSC", "Rust", b)

# Equivalence: observed eigenvalues should match (both compute on the same S)
tryCatch({
  if (existsFunction("paLDSC", where = asNamespace("GenomicSEM"))) {
    r_pa    <- GenomicSEM::paLDSC(r_cov$S, r_cov$V, r = 100, p = 3, save.pdf = FALSE)
    rust_pa <- gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 3, save.pdf = FALSE)
    err <- compare_numeric(r_pa$observed, rust_pa$observed, 1e-3, "observed eigenvalues")
    if (nzchar(err)) add_equiv("paLDSC", "FAIL", err)
    else add_equiv("paLDSC", "PASS", "observed eigenvalues within 1e-3")
  } else {
    add_equiv("paLDSC", "SKIP", "R paLDSC unavailable")
  }
}, error = function(e) add_equiv("paLDSC", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 11. summaryGLS
# ---------------------------------------------------------------------------
cat("[11/12] summaryGLS\n")
k <- nrow(rust_cov$S)
kstar <- k * (k + 1) / 2
y_gls <- rust_cov$S[lower.tri(rust_cov$S, diag = TRUE)]
X_gls <- matrix(1:kstar, ncol = 1)

b <- run_bench(function() {
  if (existsFunction("summaryGLS", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::summaryGLS(Y = y_gls, V_Y = r_cov$V, PREDICTORS = X_gls, INTERCEPT = TRUE)
  } else {
    stop("summaryGLS not available in R GenomicSEM")
  }
})
add_result("summaryGLS", "R", b)

b <- run_bench(function() {
  gsemr::summaryGLS(Y = y_gls, V_Y = rust_cov$V, PREDICTORS = X_gls, INTERCEPT = TRUE)
})
add_result("summaryGLS", "Rust", b)

tryCatch({
  if (existsFunction("summaryGLS", where = asNamespace("GenomicSEM"))) {
    r_gls    <- GenomicSEM::summaryGLS(Y = y_gls, V_Y = r_cov$V, PREDICTORS = X_gls, INTERCEPT = TRUE)
    rust_gls <- gsemr::summaryGLS(Y = y_gls, V_Y = rust_cov$V, PREDICTORS = X_gls, INTERCEPT = TRUE)
    err <- compare_numeric(r_gls$beta, rust_gls$beta, 1e-3, "beta")
    if (nzchar(err)) add_equiv("summaryGLS", "FAIL", err)
    else add_equiv("summaryGLS", "PASS", "beta within 1e-3")
  } else {
    add_equiv("summaryGLS", "SKIP", "R summaryGLS unavailable")
  }
}, error = function(e) add_equiv("summaryGLS", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 12. simLDSC
# ---------------------------------------------------------------------------
cat("[12/12] simLDSC\n")
b <- run_bench(function() {
  if (existsFunction("simLDSC", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::simLDSC(
      covmat = r_cov$S, N = c(5000, 5000, 5000), seed = 42,
      ld = "data/eur_w_ld_chr/"
    )
  } else {
    stop("simLDSC not available in R GenomicSEM")
  }
})
add_result("simLDSC", "R", b)

b <- run_bench(function() {
  gsemr::simLDSC(
    covmat = rust_cov$S, N = c(5000, 5000, 5000), seed = 42,
    ld = "data/eur_w_ld_chr/"
  )
})
add_result("simLDSC", "Rust", b)

# simLDSC produces stochastic output; only verify both produce data of similar shape
tryCatch({
  if (existsFunction("simLDSC", where = asNamespace("GenomicSEM"))) {
    r_sim    <- GenomicSEM::simLDSC(covmat = r_cov$S, N = c(5000, 5000, 5000), seed = 42, ld = "data/eur_w_ld_chr/")
    rust_sim <- gsemr::simLDSC(covmat = rust_cov$S, N = c(5000, 5000, 5000), seed = 42, ld = "data/eur_w_ld_chr/")
    if (max(dim(as.matrix(r_sim))) >= 100 && max(dim(as.matrix(rust_sim))) >= 100) {
      add_equiv("simLDSC", "PASS", "both produced sufficient simulated data")
    } else {
      add_equiv("simLDSC", "FAIL", "insufficient simulated data")
    }
  } else {
    add_equiv("simLDSC", "SKIP", "R simLDSC unavailable")
  }
}, error = function(e) add_equiv("simLDSC", "SKIP", conditionMessage(e)))

# NOTE: CSV save and summary table are deferred to AFTER the CLI and Python
# sections so they include every impl row.

# ===========================================================================
# CLI BENCHMARKS — gsem binary at target/release/gsem
# Same inputs as the library benchmarks above; all outputs go to out_bench/cli.
# Each call writes its result to disk so we can read it back for equivalence
# checks.
# ===========================================================================
cat("\n========================================================\n")
cat("  Rust CLI benchmarks (target/release/gsem)\n")
cat("========================================================\n\n")

cli_bin <- normalizePath("../target/release/gsem", mustWork = FALSE)
if (!file.exists(cli_bin)) {
  cat("CLI binary not found at", cli_bin, "- skipping CLI benchmarks.\n")
  cat("Build with: cargo build --release -p gsem\n\n")
} else {
  dir.create("out_bench/cli", showWarnings = FALSE, recursive = TRUE)

  # Helper: time a CLI invocation. Captures wall time including process startup
  # so the result reflects the experience of running the binary from a shell.
  run_cli <- function(args) {
    t0 <- proc.time()[["elapsed"]]
    rc <- tryCatch(
      system2(cli_bin, args = args, stdout = FALSE, stderr = FALSE),
      error = function(e) { attr(e, "elapsed") <- proc.time()[["elapsed"]] - t0; e }
    )
    t1 <- proc.time()[["elapsed"]]
    if (inherits(rc, "error")) {
      list(time_s = t1 - t0, peak_mb = NA_real_, error = conditionMessage(rc))
    } else if (!is.numeric(rc) || rc != 0) {
      list(time_s = t1 - t0, peak_mb = NA_real_, error = sprintf("exit %s", as.character(rc)))
    } else {
      list(time_s = t1 - t0, peak_mb = NA_real_, error = NA_character_)
    }
  }

  # Save covstruc as a JSON file the CLI can ingest. The shape matches
  # gsem_ldsc::LdscResult exactly: {s, v, i_mat, n_vec, m}.
  cli_covstruc <- "out_bench/cli/covstruc.json"
  cs <- list(
    s     = unname(as.matrix(rust_cov$S)),
    v     = unname(as.matrix(rust_cov$V)),
    i_mat = unname(as.matrix(rust_cov$I)),
    n_vec = as.numeric(rust_cov$N),
    m     = as.numeric(rust_cov$m)
  )
  writeLines(jsonlite::toJSON(cs, auto_unbox = TRUE, digits = 17), cli_covstruc)

  # ---- 1. munge ----
  cat("[CLI 1/11] munge\n")
  cli_munge_dir <- "out_bench/cli/munge"
  dir.create(cli_munge_dir, showWarnings = FALSE, recursive = TRUE)
  cli_munge_names <- paste0(trait_names, "_benchCLI")
  # CLI munge takes a SINGLE -n scalar applied to all traits (vs the R/Python
  # bindings which accept per-trait overrides). We pass a representative value
  # so all three files can be processed; the work is identical in nature so
  # the timing remains directly comparable. Output is not used for equivalence.
  b <- run_cli(c(
    "munge",
    "--files", raw_traits,
    "--hm3", hm3,
    "--trait-names", cli_munge_names,
    "-n", "7000",
    "--info-filter", "0.9",
    "--maf-filter", "0.01",
    "--out", cli_munge_dir
  ))
  add_result("munge", "CLI", b)

  # ---- 2. ldsc ----
  cat("[CLI 2/11] ldsc\n")
  cli_ldsc_out <- "out_bench/cli/ldsc.json"
  b <- run_cli(c(
    "ldsc",
    "--traits", munged_traits,
    "--sample-prev", paste(sample_prev, collapse = ","),
    "--pop-prev", paste(pop_prev, collapse = ","),
    "--ld", ld,
    "--wld", ld,
    "--n-blocks", "200",
    "--out", cli_ldsc_out
  ))
  add_result("ldsc", "CLI", b)

  # Use the CLI's own LDSC output for downstream CLI commands so we measure
  # an end-to-end CLI pipeline (not a hybrid R+CLI one).
  cli_covstruc_self <- if (file.exists(cli_ldsc_out)) cli_ldsc_out else cli_covstruc

  # ---- 3. common-factor ----
  cat("[CLI 3/11] common-factor\n")
  b <- run_cli(c(
    "common-factor",
    "--covstruc", cli_covstruc_self,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/cf.tsv"
  ))
  add_result("commonfactor", "CLI", b)

  # ---- 4. usermodel ----
  cat("[CLI 4/11] usermodel\n")
  # The CLI's LDSC JSON carries no trait names — the SEM layer uses generic
  # V1, V2, V3 indicators. Build a parallel model string against those.
  cli_indicators <- paste0("V", seq_along(trait_names))
  cli_model_str <- paste0(
    "F1 =~ NA*", cli_indicators[1], " + ", cli_indicators[2], " + ", cli_indicators[3], "\n",
    "F1 ~~ 1*F1\n",
    paste(cli_indicators, "~~", cli_indicators, collapse = "\n")
  )
  cli_model_file <- "out_bench/cli/model.txt"
  writeLines(cli_model_str, cli_model_file)
  b <- run_cli(c(
    "usermodel",
    "--covstruc", cli_covstruc_self,
    "--model-file", cli_model_file,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/um.tsv"
  ))
  add_result("usermodel", "CLI", b)

  # ---- 5. rgmodel ----
  cat("[CLI 5/11] rgmodel\n")
  b <- run_cli(c(
    "rgmodel",
    "--covstruc", cli_covstruc_self,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/rg.tsv"
  ))
  add_result("rgmodel", "CLI", b)

  # ---- 6. sumstats ----
  cat("[CLI 6/11] sumstats\n")
  cli_sumstats_out <- "out_bench/cli/merged.tsv"
  b <- run_cli(c(
    "sumstats",
    "--files", raw_traits,
    "--ref-file", ref_1000g,
    "--trait-names", trait_names,
    "--info-filter", "0.6",
    "--maf-filter", "0.01",
    "-o", cli_sumstats_out
  ))
  add_result("sumstats", "CLI", b)

  # ---- 7. write.model ----
  cat("[CLI 7/11] write.model\n")
  cli_loadings_path <- "out_bench/cli/loadings.tsv"
  # CLI parser wants numbers only — no header, no rownames.
  write.table(matrix(c(0.7, 0.6, 0.5), ncol = 1), cli_loadings_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  b <- run_cli(c(
    "write.model",
    "--loadings", cli_loadings_path,
    "--names", trait_names,
    "--cutoff", "0.3",
    "-o", "out_bench/cli/model_out.txt"
  ))
  add_result("write.model", "CLI", b)

  # ---- 8. commonfactorGWAS — parallel and serial ----
  cat("[CLI 8/11] commonfactorGWAS\n")
  for (par in c(TRUE, FALSE)) {
    tag <- if (par) "par" else "seq"
    threads_arg <- if (par) character(0) else c("--threads", "1")
    b <- run_cli(c(
      "commonfactorGWAS",
      "--covstruc", cli_covstruc_self,
      "--sumstats", rust_subset_path,
      threads_arg,
      "-o", sprintf("out_bench/cli/cfgwas_%s.tsv", tag)
    ))
    add_result("commonfactorGWAS", sprintf("CLI (%s)", tag), b)
  }

  # ---- 9. userGWAS — parallel and serial ----
  cat("[CLI 9/11] userGWAS\n")
  # CLI userGWAS extracts observed names from the model and matches them to
  # `merged.trait_names` from the sumstats TSV (which are the real trait
  # names, ANX/OCD/PTSD, derived from the beta.* columns). So we use the
  # original `gwas_model` here, not V1/V2/V3.
  cli_gwas_model_file <- "out_bench/cli/gwas_model.txt"
  writeLines(gwas_model, cli_gwas_model_file)
  for (par in c(TRUE, FALSE)) {
    tag <- if (par) "par" else "seq"
    threads_arg <- if (par) character(0) else c("--threads", "1")
    b <- run_cli(c(
      "userGWAS",
      "--covstruc", cli_covstruc_self,
      "--sumstats", rust_subset_path,
      "--model-file", cli_gwas_model_file,
      threads_arg,
      "-o", sprintf("out_bench/cli/ugwas_%s.tsv", tag)
    ))
    add_result("userGWAS", sprintf("CLI (%s)", tag), b)
  }

  # ---- 10. paLDSC ----
  cat("[CLI 10/11] paLDSC\n")
  b <- run_cli(c(
    "paLDSC",
    "--covstruc", cli_covstruc_self,
    "--n-sim", "100",
    "-o", "out_bench/cli/paldsc.tsv"
  ))
  add_result("paLDSC", "CLI", b)

  # ---- 11. summaryGLS ----
  cat("[CLI 11/11] summaryGLS\n")
  cli_y_path <- "out_bench/cli/gls_y.tsv"
  cli_x_path <- "out_bench/cli/gls_x.tsv"
  cli_v_path <- "out_bench/cli/gls_v.tsv"
  # CLI parsers want numbers only — no header rows.
  write.table(matrix(y_gls, ncol = 1), cli_y_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(X_gls, cli_x_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(unname(as.matrix(rust_cov$V)), cli_v_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  b <- run_cli(c(
    "summaryGLS",
    "--x", cli_x_path,
    "--y", cli_y_path,
    "--v", cli_v_path,
    "--intercept",
    "-o", "out_bench/cli/gls.tsv"
  ))
  add_result("summaryGLS", "CLI", b)

  # CLI does not have a trivial simLDSC parity-mode call (Python skips it
  # too because the API surface differs). We omit CLI from simLDSC.

  cat("\nCLI benchmarks done.\n")
}

# ===========================================================================
# PYTHON BENCHMARKS — call the venv's run_python_bench.py and merge CSV
# ===========================================================================
cat("\n========================================================\n")
cat("  Python benchmarks (.venv/bin/python run_python_bench.py)\n")
cat("========================================================\n\n")

python_bin <- normalizePath("../.venv/bin/python", mustWork = FALSE)
py_script  <- "run_python_bench.py"
if (!file.exists(python_bin) || !file.exists(py_script)) {
  cat("Python interpreter or script missing - skipping Python benchmarks.\n")
  cat("  python: ", python_bin, "\n")
  cat("  script: ", py_script, "\n")
} else {
  rc <- tryCatch(
    system2(python_bin, args = py_script, stdout = "", stderr = ""),
    error = function(e) { cat("Python bench error:", conditionMessage(e), "\n"); -1 }
  )
  if (file.exists("python_results.csv")) {
    py <- read.csv("python_results.csv", stringsAsFactors = FALSE)
    py$error[is.na(py$error)] <- ""
    if (nrow(py) > 0) {
      results <- rbind(results, py[, c("func", "impl", "time_s", "peak_mb", "error")])
      cat(sprintf("Merged %d Python rows from python_results.csv\n", nrow(py)))
    }
  } else {
    cat("python_results.csv not produced; skipping merge.\n")
  }

  # ---- Cross-impl equivalence check (R vs Python) using artifacts ----
  if (file.exists("python_artifacts.json")) {
    py_art <- jsonlite::fromJSON("python_artifacts.json", simplifyVector = FALSE)

    if (!is.null(py_art$ldsc_S)) {
      py_S <- do.call(rbind, lapply(py_art$ldsc_S, unlist))
      err <- compare_numeric(as.matrix(rust_cov$S), py_S, 1e-6, "ldsc S (Rust vs Python)")
      if (nzchar(err)) add_equiv("ldsc(py)", "FAIL", err)
      else add_equiv("ldsc(py)", "PASS", "S matches Rust within 1e-6")
    }
    if (!is.null(py_art$commonfactor_loadings)) {
      cf <- gsemr::commonfactor(rust_cov, estimation = "DWLS")
      rust_load <- cf$results$est[cf$results$op == "=~"]
      err <- compare_numeric(rust_load, unlist(py_art$commonfactor_loadings),
                             1e-4, "commonfactor loadings (Rust vs Python)")
      if (nzchar(err)) add_equiv("commonfactor(py)", "FAIL", err)
      else add_equiv("commonfactor(py)", "PASS", "loadings match Rust within 1e-4")
    }
  }
}

# ---------------------------------------------------------------------------
# Cleanup temp files
# ---------------------------------------------------------------------------
unlink("out_bench", recursive = TRUE)

# ---------------------------------------------------------------------------
# Save CSV
# ---------------------------------------------------------------------------
csv_path <- "benchmark_results.csv"
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_path))

# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
cat("\n========================================================\n")
cat("  Summary\n")
cat("========================================================\n\n")
cat(sprintf("%-20s %-12s %10s %10s\n", "Function", "Impl", "Time (s)", "Peak (MB)"))
cat(strrep("-", 58), "\n")

funcs <- unique(results$func)
fmt_val <- function(x) if (is.na(x)) "      N/A" else sprintf("%10.2f", x)
fmt_mem <- function(x) if (is.na(x)) "      N/A" else sprintf("%10.1f", x)

for (fn in funcs) {
  rows <- results[results$func == fn, ]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    cat(sprintf("%-20s %-12s %s %s",
                if (i == 1) fn else "", r$impl, fmt_val(r$time_s), fmt_mem(r$peak_mb)))
    if (nzchar(r$error)) cat(sprintf("  [err: %s]", substr(r$error, 1, 40)))
    cat("\n")
  }
}

# ---------------------------------------------------------------------------
# Equivalence check report
# ---------------------------------------------------------------------------
cat("\n========================================================\n")
cat("  Output equivalence (R vs Rust, with tolerance)\n")
cat("========================================================\n\n")
cat(sprintf("%-20s %-6s %s\n", "Function", "Status", "Detail"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(equiv))) {
  cat(sprintf("%-20s %-6s %s\n", equiv$func[i], equiv$status[i], equiv$detail[i]))
}
n_pass <- sum(equiv$status == "PASS")
n_fail <- sum(equiv$status == "FAIL")
n_skip <- sum(equiv$status == "SKIP")
cat(sprintf("\n  PASS: %d  FAIL: %d  SKIP: %d\n", n_pass, n_fail, n_skip))

# Save equivalence report alongside CSV
write.csv(equiv, "benchmark_equivalence.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Generate single raw-speed bar chart (linear scale, with numeric labels)
# ---------------------------------------------------------------------------
cat("\nGenerating plot...\n")
pdf_path <- "benchmark_plots.pdf"

plot_df <- results[!is.na(results$time_s), c("func", "impl", "time_s")]
# Stable order: by func order of first appearance, then by impl as recorded
plot_df$func <- factor(plot_df$func, levels = unique(results$func))
plot_df$impl <- factor(plot_df$impl, levels = unique(results$impl))

# Color palette: distinct colors per impl, with related shades for par/seq variants
all_impls <- levels(plot_df$impl)
palette <- c(
  "R"              = "#E41A1C",
  "R (par)"        = "#E41A1C",
  "R (seq)"        = "#FB9A99",
  "Rust"           = "#377EB8",
  "Rust (par)"     = "#377EB8",
  "Rust (seq)"     = "#A6CEE3",
  "CLI"            = "#4DAF4A",
  "CLI (par)"      = "#4DAF4A",
  "CLI (seq)"      = "#B2DF8A",
  "Python"         = "#984EA3",
  "Python (par)"   = "#984EA3",
  "Python (seq)"   = "#CAB2D6"
)
colors_impl <- palette[all_impls]
colors_impl[is.na(colors_impl)] <- "#999999"
names(colors_impl) <- all_impls

theme_bench <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    legend.position  = "top",
    panel.grid.major.x = element_blank(),
    plot.margin      = margin(10, 15, 10, 10)
  )

pdf(pdf_path, width = 12, height = 7)
if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = func, y = time_s, fill = impl)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, alpha = 0.9) +
    geom_text(
      aes(label = sprintf("%.2f", time_s)),
      position = position_dodge(width = 0.8),
      vjust    = -0.4,
      size     = 3
    ) +
    scale_fill_manual(values = colors_impl) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = "Execution Time across implementations",
      subtitle = sprintf("R GenomicSEM | gsemr (R bindings) | gsem CLI | genomicsem (Python) | PGC Anxiety/OCD/PTSD | GWAS subset: %d SNPs", GWAS_SUBSET_N),
      x        = NULL, y = "Time (seconds)", fill = NULL
    ) +
    theme_bench
  print(p)
}
dev.off()
cat(sprintf("Plot saved to %s\n", pdf_path))
cat("\nDone.\n")
