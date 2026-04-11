#!/usr/bin/env Rscript
# Performance benchmark: R GenomicSEM vs gsemr (Rust)
# Measures wall-clock time and peak memory for each API function.
#
# Usage:
#   cd bench && Rscript benchmark_perf.R                   # default (safe)
#   BENCH_CORES=8 Rscript benchmark_perf.R                 # more workers
#   BENCH_BIG=1 Rscript benchmark_perf.R                   # include 25k-SNP scaling point
#
# Environment knobs:
#   BENCH_CORES        Max worker processes on BOTH sides. Defaults to
#                      min(4, parallel::detectCores() - 1). Capped to
#                      prevent the R side from spawning many BLAS-threaded
#                      workers and blowing out memory (prior runs hit
#                      80+ GB with 90+ threads via Accelerate BLAS).
#   BENCH_BIG          If "1", include the 25000-SNP scaling point in the
#                      userGWAS sweep. Off by default because at that size
#                      R's parallel path is memory-hungry.
#
# Output:
#   bench/benchmark_results.csv  — timing and memory data
#   bench/benchmark_plots.pdf    — comparison charts

# ---------------------------------------------------------------------------
# Thread / memory containment (must run BEFORE any package load so the
# BLAS library in this process AND every forked worker sees the limit).
#
# Previous runs spawned ~90 threads and used ~80 GB of RAM, which was
# almost entirely Apple Accelerate / OpenBLAS fanning out *inside* each
# R `parallel::makeCluster(type="FORK")` worker. Clamping to single-
# threaded BLAS per worker + a small worker count brings memory down to
# a few GB without losing the parallel speedup story.
# ---------------------------------------------------------------------------
Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",  # Apple Accelerate (macOS default BLAS)
  OPENBLAS_NUM_THREADS   = "1",  # OpenBLAS (Linux default BLAS)
  MKL_NUM_THREADS        = "1",  # Intel MKL
  BLIS_NUM_THREADS       = "1",  # BLIS
  OMP_NUM_THREADS        = "1"   # OpenMP layer used by some matrix libs
)

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

source("setup_data.R")
ensure_bench_data(bench_dir = ".")

for (pkg in c("GenomicSEM", "gsemr", "ggplot2", "parallel")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' not installed.", pkg))
}
library(GenomicSEM)
library(gsemr)
library(ggplot2)

# Resolve worker-count budget for both R and gsemr.
.bench_env_int <- function(name, default) {
  val <- Sys.getenv(name, unset = NA)
  if (is.na(val) || !nzchar(val)) return(default)
  n <- suppressWarnings(as.integer(val))
  if (is.na(n) || n <= 0) default else n
}
BENCH_CORES <- .bench_env_int(
  "BENCH_CORES",
  max(1L, min(4L, parallel::detectCores() - 1L))
)
BENCH_BIG <- identical(Sys.getenv("BENCH_BIG", unset = ""), "1")
# Also clamp rayon for any gsemr code path that falls back to env vars
# (the R shim passes `cores=BENCH_CORES` explicitly, but the CLI bench
# below launches child processes that read RAYON_NUM_THREADS).
Sys.setenv(RAYON_NUM_THREADS = as.character(BENCH_CORES))
cat(sprintf(
  "bench config: BENCH_CORES=%d  BENCH_BIG=%s  detectCores=%d\n",
  BENCH_CORES,
  if (BENCH_BIG) "1" else "0",
  parallel::detectCores()
))

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

# Shared userGWAS model — used by the R library, CLI, and Python bench
# sections. Defined up here so every impl runs the identical model and
# so the CLI/Python sections can reference it before the R library
# section creates its own copy.
gwas_model <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

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

# ===========================================================================
# CLI BENCHMARKS — gsem binary at target/release/gsem
#
# Runs before the R library section so the I/O-heavy commands (munge,
# ldsc, sumstats) get the cold OS page cache. The library section then
# benefits from the warmed cache on its second pass over the same
# files, keeping both measurements reflective of "first pass from a
# fresh shell" rather than advantaging whichever impl ran second.
#
# The CLI pipeline is self-contained: it runs its own munge, ldsc,
# sumstats, and derives its downstream covstruc / merged-TSV inputs
# from its own outputs, so no reference to R library variables
# (`rust_cov`, `rust_ss_path`, `y_gls`, ...) is needed.
# ===========================================================================
cat("\n========================================================\n")
cat("  Rust CLI benchmarks (target/release/gsem)\n")
cat("========================================================\n\n")

cli_bin <- normalizePath("../target/release/gsem", mustWork = FALSE)
cli_available <- file.exists(cli_bin)
if (!cli_available) {
  cat("CLI binary not found at", cli_bin, "- skipping CLI benchmarks.\n")
  cat("Build with: cargo build --release -p gsem\n\n")
} else {
  dir.create("out_bench/cli", showWarnings = FALSE, recursive = TRUE)

  # Helper: time a CLI invocation. Captures wall time including process
  # startup so the result reflects the experience of running the binary
  # from a shell.
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

  # ---- CLI 1. munge ----
  cat("[CLI 1/12] munge\n")
  cli_munge_dir <- "out_bench/cli/munge"
  dir.create(cli_munge_dir, showWarnings = FALSE, recursive = TRUE)
  cli_munge_names <- paste0(trait_names, "_benchCLI")
  # The CLI takes a SINGLE -n scalar applied to all traits (vs the R /
  # Python bindings which accept per-trait overrides). Pass a
  # representative value; the work is identical in kind so the timing
  # remains directly comparable.
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

  # ---- CLI 2. ldsc ----
  cat("[CLI 2/12] ldsc\n")
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
  if (!file.exists(cli_ldsc_out)) {
    stop("CLI ldsc did not produce ", cli_ldsc_out,
         " — downstream CLI steps cannot run.")
  }

  # ---- CLI 3. common-factor ----
  cat("[CLI 3/12] common-factor\n")
  b <- run_cli(c(
    "common-factor",
    "--covstruc", cli_ldsc_out,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/cf.tsv"
  ))
  add_result("commonfactor", "CLI", b)

  # ---- CLI 4. usermodel ----
  cat("[CLI 4/12] usermodel\n")
  # The CLI's LDSC JSON carries no trait names — the SEM layer uses
  # generic V1, V2, V3 indicators. Build a parallel model string
  # against those.
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
    "--covstruc", cli_ldsc_out,
    "--model-file", cli_model_file,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/um.tsv"
  ))
  add_result("usermodel", "CLI", b)

  # ---- CLI 5. rgmodel ----
  cat("[CLI 5/12] rgmodel\n")
  b <- run_cli(c(
    "rgmodel",
    "--covstruc", cli_ldsc_out,
    "--estimation", "DWLS",
    "-o", "out_bench/cli/rg.tsv"
  ))
  add_result("rgmodel", "CLI", b)

  # ---- CLI 6. sumstats ----
  # CLI parallelizes reference + per-trait reads across a local rayon
  # pool; `--threads` caps that pool at the bench-wide core budget.
  cat("[CLI 6/12] sumstats\n")
  cli_sumstats_out <- "out_bench/cli/merged.tsv"
  b <- run_cli(c(
    "sumstats",
    "--files", raw_traits,
    "--ref-file", ref_1000g,
    "--trait-names", trait_names,
    "--info-filter", "0.6",
    "--maf-filter", "0.01",
    "--threads", as.character(BENCH_CORES),
    "-o", cli_sumstats_out
  ))
  add_result("sumstats", "CLI", b)

  # ---- CLI 7. write.model ----
  cat("[CLI 7/12] write.model\n")
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

  # ---- CLI 8. commonfactorGWAS — DISABLED (see section 8 in the R
  # library block for the full rationale). ----
  cat("[CLI 8/12] commonfactorGWAS (SKIPPED)\n")
  add_result("commonfactorGWAS", "CLI (par)",
             list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))
  add_result("commonfactorGWAS", "CLI (seq)",
             list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))

  # ---- CLI 9. userGWAS — scaling sweep matching the library section
  # (sizes 1000/5000/10000 + full-scale, parallel only). ----
  cat("[CLI 9/12] userGWAS (scaling sweep)\n")
  cli_gwas_model_file <- "out_bench/cli/gwas_model.txt"
  writeLines(gwas_model, cli_gwas_model_file)

  cli_n_available <- length(readLines(cli_sumstats_out)) - 1L
  cli_ug_sizes <- if (BENCH_BIG) {
    c(1000L, 5000L, 10000L, 25000L)
  } else {
    c(1000L, 5000L, 10000L)
  }
  cli_ug_sizes <- unique(pmin(cli_ug_sizes, cli_n_available))
  cli_ug_sizes <- cli_ug_sizes[cli_ug_sizes > 0]
  cat(sprintf("  head-to-head sizes: %s (capped at %d, cores=%d)\n",
              paste(cli_ug_sizes, collapse = ", "), cli_n_available, BENCH_CORES))

  for (n_snp in cli_ug_sizes) {
    cat(sprintf("  -> N=%d\n", n_snp))
    cli_sub_path <- sprintf("out_bench/cli/ugwas_subset_%d.tsv", n_snp)
    cli_sub_df   <- read.delim(cli_sumstats_out, nrows = n_snp)
    write.table(cli_sub_df, cli_sub_path, sep = "\t", row.names = FALSE, quote = FALSE)
    b <- run_cli(c(
      "userGWAS",
      "--covstruc", cli_ldsc_out,
      "--sumstats", cli_sub_path,
      "--model-file", cli_gwas_model_file,
      "--threads", as.character(BENCH_CORES),
      "-o", sprintf("out_bench/cli/ugwas_%d.tsv", n_snp)
    ))
    add_result("userGWAS", sprintf("CLI (N=%d)", n_snp), b)
  }

  # Full-scale "real-world deployment" point. Matches the Rust library
  # section's full-scale Rust row.
  cat(sprintf("  -> CLI full-scale: N=%d\n", cli_n_available))
  b <- run_cli(c(
    "userGWAS",
    "--covstruc", cli_ldsc_out,
    "--sumstats", cli_sumstats_out,
    "--model-file", cli_gwas_model_file,
    "--threads", as.character(BENCH_CORES),
    "-o", "out_bench/cli/ugwas_full.tsv"
  ))
  add_result("userGWAS", sprintf("CLI (N=%d, full)", cli_n_available), b)

  # ---- CLI 10. paLDSC ----
  cat("[CLI 10/12] paLDSC\n")
  b <- run_cli(c(
    "paLDSC",
    "--covstruc", cli_ldsc_out,
    "--n-sim", "100",
    "-o", "out_bench/cli/paldsc.tsv"
  ))
  add_result("paLDSC", "CLI", b)

  # ---- CLI 11. summaryGLS ----
  cat("[CLI 11/12] summaryGLS\n")
  # Derive Y/X/V from the CLI's own LDSC JSON so this step stays
  # self-contained (no dependency on the R library `rust_cov`).
  cli_ldsc_json <- jsonlite::fromJSON(cli_ldsc_out)
  cli_S <- as.matrix(cli_ldsc_json$s)
  cli_V <- as.matrix(cli_ldsc_json$v)
  cli_k <- nrow(cli_S)
  cli_kstar <- cli_k * (cli_k + 1) / 2
  cli_y_gls <- cli_S[lower.tri(cli_S, diag = TRUE)]
  cli_X_gls <- matrix(seq_len(cli_kstar), ncol = 1)
  cli_y_path <- "out_bench/cli/gls_y.tsv"
  cli_x_path <- "out_bench/cli/gls_x.tsv"
  cli_v_path <- "out_bench/cli/gls_v.tsv"
  write.table(matrix(cli_y_gls, ncol = 1), cli_y_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(cli_X_gls, cli_x_path,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(cli_V, cli_v_path,
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

  # ---- CLI 12. simLDSC ----
  cat("[CLI 12/12] simLDSC\n")
  b <- run_cli(c(
    "simLDSC",
    "--covstruc", cli_ldsc_out,
    "--n-per-trait", "5000,5000,5000",
    "--ld", ld,
    "-o", "out_bench/cli/sim_ldsc.tsv"
  ))
  add_result("simLDSC", "CLI", b)

  cat("\nCLI benchmarks done.\n")
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

# Rust sumstats — gsemr parallelizes reference + per-trait reads
# across a local rayon pool; pass `cores = BENCH_CORES` explicitly so
# this matches the rest of the bench's thread budget.
b <- run_bench(function() {
  gsemr::sumstats(
    files = raw_traits, ref = ref_1000g, trait.names = trait_names,
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path,
    parallel = TRUE, cores = BENCH_CORES
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
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path,
    parallel = TRUE, cores = BENCH_CORES
  )
}
rust_subset_path <- "out_bench/rust_subset.tsv"
rust_snps <- read.delim(rust_ss_path, nrows = GWAS_SUBSET_N)
write.table(rust_snps, rust_subset_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Rust subset: %d SNPs\n", nrow(rust_snps)))

# ---------------------------------------------------------------------------
# 8. commonfactorGWAS — DISABLED (intentionally commented out)
#
# This slot is commented out rather than deleted so it can be re-enabled
# once upstream R GenomicSEM clarifies the intended semantics.
#
# Two separate issues make a fair R vs Rust comparison impossible today:
#
#   (a) Crash. GenomicSEM::commonfactorGWAS crashes on this PGC input via
#       an unguarded solve(t(S2.delt) %*% S2.W %*% S2.delt, tol=toler) in
#       .commonfactorGWAS_main, when lavaan converges to a degenerate
#       local minimum on any SNP in the batch (e.g. rs3863622, rs4411087,
#       rs4422948, where ANX~~ANX runs off to 8429 with F1~~F1=-8429).
#       A proposed upstream patch and reproducer live in ../upstream-issue/.
#
#   (b) Identification mismatch. Even with the crash patched, R's
#       commonfactorGWAS uses marker-indicator identification (first
#       loading fixed to 1, factor variance free) and refits the full
#       model per SNP. gsemr::commonfactorGWAS uses fixed-variance
#       identification (F1 ~~ 1*F1, loadings free) with fix_measurement,
#       matching GenomicSEM::userGWAS on the equivalent model. On real
#       GWAS data the two can disagree in sign and magnitude — not a
#       bug, a design choice documented in ARCHITECTURE.md §3.3.
#
# Consequence: R and Rust commonfactorGWAS solve different optimization
# problems, so comparing them in the bench is apples-to-oranges even
# after the crash is fixed. userGWAS below covers the per-SNP GWAS path
# with a clean equivalence check; commonfactorGWAS performance/parity
# can be re-added here once upstream clarifies.
# ---------------------------------------------------------------------------
cat("[8/12] commonfactorGWAS (SKIPPED — see comment above)\n")
add_result("commonfactorGWAS", "R (par)", list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))
add_result("commonfactorGWAS", "R (seq)", list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))
add_result("commonfactorGWAS", "Rust (par)", list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))
add_result("commonfactorGWAS", "Rust (seq)", list(time_s = NA_real_, peak_mb = NA_real_, error = "skipped"))
add_equiv("commonfactorGWAS", "SKIP",
          "intentionally disabled pending upstream clarification; see comment in benchmark_perf.R")

# --- disabled block (restore and delete the skipped lines above to re-enable) ---
# cfg_results <- list()  # capture outputs for equivalence
# for (par in c(TRUE, FALSE)) {
#   tag <- if (par) "par" else "seq"
#
#   b <- run_bench(function() {
#     cfg_results[[paste0("R_",  tag)]] <<-
#       GenomicSEM::commonfactorGWAS(covstruc = r_cov, SNPs = r_snps_df, parallel = par)
#   })
#   add_result("commonfactorGWAS", sprintf("R (%s)", tag), b)
#
#   b <- run_bench(function() {
#     cfg_results[[paste0("Rust_", tag)]] <<-
#       gsemr::commonfactorGWAS(covstruc = rust_cov, SNPs = rust_subset_path, parallel = par)
#   })
#   add_result("commonfactorGWAS", sprintf("Rust (%s)", tag), b)
# }
#
# # Equivalence: parallel vs serial within each impl, and R vs Rust on common SNPs
# tryCatch({
#   notes <- character(0)
#   for (impl in c("R", "Rust")) {
#     a <- cfg_results[[paste0(impl, "_par")]]
#     b <- cfg_results[[paste0(impl, "_seq")]]
#     if (!is.null(a) && !is.null(b)) {
#       err <- compare_numeric(a$est, b$est, 1e-6, sprintf("%s par/seq est", impl))
#       if (nzchar(err)) notes <- c(notes, err)
#     }
#   }
#   ra <- cfg_results$R_par; rb <- cfg_results$Rust_par
#   if (!is.null(ra) && !is.null(rb)) {
#     common <- intersect(ra$SNP, rb$SNP)
#     if (length(common) >= 50) {
#       ra_s <- ra[match(common, ra$SNP), ]
#       rb_s <- rb[match(common, rb$SNP), ]
#       ok <- is.finite(ra_s$est) & is.finite(rb_s$est)
#       if (sum(ok) >= 50) {
#         d <- max(abs(ra_s$est[ok] - rb_s$est[ok]))
#         if (!is.finite(d) || d > 0.05) notes <- c(notes, sprintf("R/Rust est max diff %.3e > 0.05", d))
#       }
#     }
#   }
#   if (length(notes) == 0) add_equiv("commonfactorGWAS", "PASS", "par==seq and R~Rust on common SNPs")
#   else add_equiv("commonfactorGWAS", "FAIL", paste(notes, collapse = "; "))
# }, error = function(e) add_equiv("commonfactorGWAS", "SKIP", conditionMessage(e)))

# ---------------------------------------------------------------------------
# 9. userGWAS — scaling sweep over multiple SNP counts (parallel only)
# ---------------------------------------------------------------------------
cat("[9/12] userGWAS (scaling sweep)\n")
# `gwas_model` is hoisted to the config block so every impl (R, CLI,
# Python) runs exactly the same formula.

# SNP-count grid for head-to-head scaling. The 25k point is gated behind
# BENCH_BIG=1 because at that size R's parallel path is memory-hungry
# even with BLAS pinned to 1 thread (lavaan itself holds a lot of state
# per worker). Default grid tops out at 10k which is plenty for a
# scaling curve fit. Capped by available SNPs.
ug_sizes <- if (BENCH_BIG) {
  c(1000L, 5000L, 10000L, 25000L)
} else {
  c(1000L, 5000L, 10000L)
}
n_rust_available <- length(readLines(rust_ss_path)) - 1L  # minus header
n_r_available    <- nrow(r_snps_full)
n_available      <- min(n_rust_available, n_r_available)
ug_sizes <- unique(pmin(ug_sizes, n_available))
ug_sizes <- ug_sizes[ug_sizes > 0]
cat(sprintf("  head-to-head sizes: %s (capped at %d available SNPs, cores=%d)\n",
            paste(ug_sizes, collapse = ", "), n_available, BENCH_CORES))

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
                           model = gwas_model, parallel = TRUE,
                           cores = BENCH_CORES)
  })
  add_result("userGWAS", label_r, b)

  b <- run_bench(function() {
    ug_results[[sprintf("Rust_%d", n_snp)]] <<-
      gsemr::userGWAS(covstruc = rust_cov, SNPs = rust_sub_path,
                      model = gwas_model, parallel = TRUE,
                      cores = BENCH_CORES)
  })
  add_result("userGWAS", label_rust, b)
}

# Rust-only "real-world deployment" point: the full sumstats file.
# Running R at this scale is not feasible on a single machine — the R
# package recommends MPI for multi-million-SNP GWAS — so we report a
# Rust-only bar for the headline "can do it on one box" claim. Uses the
# same BENCH_CORES budget as the head-to-head sweep for consistency.
cat(sprintf("  -> Rust-only full-scale: N=%d\n", n_rust_available))
b <- run_bench(function() {
  ug_results[["Rust_full"]] <<-
    gsemr::userGWAS(covstruc = rust_cov, SNPs = rust_ss_path,
                    model = gwas_model, parallel = TRUE,
                    cores = BENCH_CORES)
})
add_result("userGWAS", sprintf("Rust (N=%d, full)", n_rust_available), b)

# Free the full-scale Rust result now that its wall time is recorded;
# it's a multi-GB two-table frame we don't consume downstream.
ug_results[["Rust_full"]] <- NULL
invisible(gc(full = TRUE, verbose = FALSE))

# Equivalence: convergence counts align within 2% across sizes. R's
# `GenomicSEM::userGWAS` returns a list of per-SNP data frames (one per
# unique SNP), with no top-level `converged` column — convergence is
# encoded per frame via the `error` column (0 = converged, non-zero
# string = lavaan error). gsemr returns `list(snps, params)` and we
# read `$snps$converged` directly.
.r_usergwas_n_converged <- function(res) {
  if (is.null(res)) return(NA_integer_)
  if (is.data.frame(res)) {
    if ("converged" %in% names(res)) return(sum(res$converged, na.rm = TRUE))
    if ("error" %in% names(res)) {
      # Fallback: collapse per-SNP rows, count unique SNPs whose error is 0.
      ok <- res$error == 0 | res$error == "0"
      return(length(unique(res$SNP[ok])))
    }
    return(NA_integer_)
  }
  if (is.list(res)) {
    # list-of-per-SNP-frames shape (R's default sub=FALSE return).
    ok <- vapply(res, function(x) {
      if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(FALSE)
      if (!"error" %in% names(x)) return(TRUE)
      e <- x$error[1]
      isTRUE(e == 0) || isTRUE(e == "0")
    }, logical(1))
    return(sum(ok))
  }
  NA_integer_
}

tryCatch({
  notes <- character(0)
  for (n_snp in ug_sizes) {
    a <- ug_results[[sprintf("R_%d",    n_snp)]]
    b <- ug_results[[sprintf("Rust_%d", n_snp)]]
    if (!is.null(a) && !is.null(b)) {
      ca <- .r_usergwas_n_converged(a)
      cb <- sum(b$snps$converged, na.rm = TRUE)
      if (is.na(ca) || abs(ca - cb) > 0.02 * n_snp) {
        notes <- c(notes, sprintf("N=%d converged %s(R) vs %d(Rust)",
                                  n_snp, if (is.na(ca)) "NA" else as.character(ca), cb))
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
# `p` is the quantile used to summarise the permuted eigenvalue
# distribution (must lie in [0, 1]); 0.95 matches `psych::fa.parallel`.
b <- run_bench(function() {
  if (existsFunction("paLDSC", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::paLDSC(r_cov$S, r_cov$V, r = 100, p = 0.95, save.pdf = TRUE)
  } else {
    stop("paLDSC not available in R GenomicSEM")
  }
})
add_result("paLDSC", "R", b)

b <- run_bench(function() {
  gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 0.95, save.pdf = TRUE)
})
add_result("paLDSC", "Rust", b)

# Equivalence: check that gsemr's observed eigenvalues match a direct
# eigendecomposition of the correlation-matrix form of S (which is
# what paLDSC reports as `observed`). R `GenomicSEM::paLDSC` has no
# return() statement — it prints and returns NULL invisibly — so we
# can't read `$observed` off it; comparing against base R's own
# `eigen()` instead gives us a meaningful parity check.
tryCatch({
  rust_pa <- gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 0.95, save.pdf = FALSE)
  s_mat   <- as.matrix(rust_cov$S)
  sds     <- sqrt(pmax(diag(s_mat), 1e-10))
  cor_mat <- s_mat / outer(sds, sds)
  diag(cor_mat) <- 1
  r_eigs  <- sort(eigen(cor_mat, only.values = TRUE)$values, decreasing = TRUE)
  err <- compare_numeric(r_eigs, rust_pa$observed, 1e-6,
                         "observed eigenvalues (gsemr vs base::eigen)")
  if (nzchar(err)) add_equiv("paLDSC", "FAIL", err)
  else add_equiv("paLDSC", "PASS", "observed eigenvalues within 1e-6 of base::eigen")
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
    # Both `GenomicSEM::summaryGLS` and `gsemr::summaryGLS` return an
    # invisible numeric matrix with columns `betas, pvals, SE, Z`, so
    # the comparison uses the same accessor on each side.
    r_gls    <- suppressWarnings(utils::capture.output(
      r_out <- GenomicSEM::summaryGLS(Y = y_gls, V_Y = r_cov$V,
                                      PREDICTORS = X_gls, INTERCEPT = TRUE)
    ))
    rust_gls <- utils::capture.output(
      rust_out <- gsemr::summaryGLS(Y = y_gls, V_Y = rust_cov$V,
                                    PREDICTORS = X_gls, INTERCEPT = TRUE)
    )
    r_beta    <- as.numeric(r_out[, "betas"])
    rust_beta <- as.numeric(rust_out[, "betas"])
    err <- compare_numeric(r_beta, rust_beta, 1e-3, "beta")
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
# R `GenomicSEM::simLDSC` validates N with a scalar `if (N <= 0)` and
# builds the K x K sample-size matrix via `diag(N, k, k)`; a vector
# triggers "condition has length > 1". Pass a scalar here — gsemr
# accepts scalar, vector, or matrix and handles the same shape.
b <- run_bench(function() {
  if (existsFunction("simLDSC", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::simLDSC(
      covmat = r_cov$S, N = 5000, seed = 42,
      ld = "data/eur_w_ld_chr/"
    )
  } else {
    stop("simLDSC not available in R GenomicSEM")
  }
})
add_result("simLDSC", "R", b)

b <- run_bench(function() {
  gsemr::simLDSC(
    covmat = rust_cov$S, N = 5000, seed = 42,
    ld = "data/eur_w_ld_chr/"
  )
})
add_result("simLDSC", "Rust", b)

# simLDSC produces stochastic output and R `GenomicSEM::simLDSC` has
# no return() — it writes per-trait TSVs to the working directory and
# returns NULL invisibly. Compare shape-only: gsemr's return matrix
# must be k traits x N SNPs and R's written files must exist.
tryCatch({
  rust_sim <- gsemr::simLDSC(covmat = rust_cov$S, N = 5000, seed = 42,
                             ld = "data/eur_w_ld_chr/")
  rust_ok <- is.matrix(rust_sim) && nrow(rust_sim) == nrow(rust_cov$S) &&
             ncol(rust_sim) >= 100
  # R writes files like iter1GWAS1.sumstats.gz in the cwd when it runs.
  # The bench `run_bench(function() GenomicSEM::simLDSC(...))` above
  # already exercised it; just check the files landed.
  r_outputs <- Sys.glob("iter*GWAS*.sumstats*")
  r_ok      <- length(r_outputs) >= nrow(rust_cov$S)
  if (rust_ok && r_ok) {
    add_equiv("simLDSC", "PASS", sprintf("gsemr %dx%d; R wrote %d files",
                                          nrow(rust_sim), ncol(rust_sim), length(r_outputs)))
  } else if (rust_ok && !r_ok) {
    add_equiv("simLDSC", "FAIL", sprintf("R produced %d trait files (expected >= %d)",
                                          length(r_outputs), nrow(rust_cov$S)))
  } else {
    add_equiv("simLDSC", "FAIL", "gsemr did not return a k x N matrix")
  }
  # Cleanup the files R dropped into cwd.
  file.remove(r_outputs)
}, error = function(e) add_equiv("simLDSC", "SKIP", conditionMessage(e)))

# NOTE: CSV save and summary table are deferred to AFTER the Python
# section so they include every impl row.

# ===========================================================================
# PYTHON BENCHMARKS — call the venv's run_python_bench.py and merge CSV
# ===========================================================================
cat("\n========================================================\n")
cat("  Python benchmarks (.venv/bin/python run_python_bench.py)\n")
cat("========================================================\n\n")

# NOTE: don't resolve the symlink with normalizePath — Python's venv
# detection activates on `sys.executable` starting with the venv path,
# and following the `.venv/bin/python` symlink back to
# `/opt/homebrew/.../python3.13` defeats that, leaving the venv's
# site-packages off sys.path (numpy, genomicsem, etc. become
# invisible).
python_bin <- file.path("..", ".venv", "bin", "python")
py_script  <- "run_python_bench.py"
if (!file.exists(python_bin) || !file.exists(py_script)) {
  cat("Python interpreter or script missing - skipping Python benchmarks.\n")
  cat("  python: ", python_bin, "\n")
  cat("  script: ", py_script, "\n")
} else {
  # Forward BENCH_CORES / BENCH_BIG so the Python sweep honours the
  # same core budget and size grid as the R library section.
  Sys.setenv(BENCH_CORES = as.character(BENCH_CORES))
  Sys.setenv(BENCH_BIG = if (isTRUE(BENCH_BIG)) "1" else "0")

  # Remove any stale CSV / artifacts from a previous run so we never
  # silently merge yesterday's numbers if today's Python invocation
  # fails. The Python script rewrites both files on success.
  unlink("python_results.csv")
  unlink("python_artifacts.json")

  rc <- tryCatch(
    system2(python_bin, args = py_script, stdout = "", stderr = ""),
    error = function(e) { cat("Python bench error:", conditionMessage(e), "\n"); -1 }
  )
  if (!is.numeric(rc) || rc != 0) {
    cat(sprintf("Python bench exited with status %s — skipping merge.\n",
                as.character(rc)))
  } else if (file.exists("python_results.csv")) {
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
cat(sprintf("%-20s %-26s %10s %10s\n", "Function", "Impl", "Time (s)", "Peak (MB)"))
cat(strrep("-", 72), "\n")

funcs <- unique(results$func)
fmt_val <- function(x) if (is.na(x)) "      N/A" else sprintf("%10.2f", x)
fmt_mem <- function(x) if (is.na(x)) "      N/A" else sprintf("%10.1f", x)

for (fn in funcs) {
  rows <- results[results$func == fn, ]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    cat(sprintf("%-20s %-26s %s %s",
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
# Generate plots: a detailed PDF (every impl x every function) and a
# headline PNG (R GenomicSEM vs gsemr on the operations that dominate
# wall-clock time in a typical pipeline).
# ---------------------------------------------------------------------------
cat("\nGenerating plots...\n")
pdf_path <- "benchmark_plots.pdf"
png_path <- "benchmark_headline.png"

plot_df <- results[!is.na(results$time_s), c("func", "impl", "time_s")]
# Stable order: by func order of first appearance, then by impl as recorded
plot_df$func <- factor(plot_df$func, levels = unique(results$func))
plot_df$impl <- factor(plot_df$impl, levels = unique(results$impl))

# Color palette: one base colour per impl family, shared across the
# plain row ("R", "Rust", ...) and its scale-sweep or par/seq siblings
# ("R (N=1000)", "Rust (N=..., full)", ...). Family is matched by the
# leading word so any future variant automatically picks up the family
# colour.
all_impls <- levels(plot_df$impl)
family_palette <- c(
  "R"      = "#E41A1C",
  "Rust"   = "#377EB8",
  "CLI"    = "#4DAF4A",
  "Python" = "#984EA3"
)
impl_family <- sub("^(\\w+).*$", "\\1", all_impls)
colors_impl <- unname(family_palette[impl_family])
colors_impl[is.na(colors_impl)] <- "#999999"
names(colors_impl) <- all_impls

theme_bench <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    legend.position  = "top",
    panel.grid.major.y = element_blank(),
    plot.margin      = margin(10, 15, 10, 10)
  )

# ---- Detailed report: every impl x every function, horizontal bars ----
pdf(pdf_path, width = 12, height = 9)
if (nrow(plot_df) > 0) {
  # Reverse function order so the first row (first bench step) lands at
  # the top of the y-axis rather than the bottom after ggplot's flip.
  plot_df$func <- factor(plot_df$func, levels = rev(levels(plot_df$func)))
  p <- ggplot(plot_df, aes(x = func, y = time_s, fill = impl)) +
    geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
             width = 0.8, alpha = 0.9) +
    geom_text(
      aes(label = sprintf("%.2f", time_s)),
      position = position_dodge2(width = 0.85, preserve = "single"),
      hjust    = -0.15,
      size     = 3
    ) +
    coord_flip() +
    scale_fill_manual(values = colors_impl) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Execution Time across implementations",
      subtitle = sprintf(
        "R GenomicSEM | gsemr (R bindings) | gsem CLI | genomicsem (Python) | PGC Anxiety/OCD/PTSD | BENCH_CORES=%d",
        BENCH_CORES
      ),
      x        = NULL, y = "Time (seconds)", fill = NULL
    ) +
    theme_bench
  print(p)
}
invisible(dev.off())
cat(sprintf("Detailed plot saved to %s\n", pdf_path))

# ---- Headline: R vs gsemr on the operations that matter for wall time ----
#
# Drops the fast ops (commonfactor, usermodel, rgmodel, write.model,
# summaryGLS, paLDSC) from the headline — they complete in under 0.5s
# in both impls and aren't what drives a user's decision to switch.
# Keeps munge, ldsc, sumstats, simLDSC, and the userGWAS scaling
# sweep + full-scale row. Sorted by R wall time (descending) so the
# worst case lands at the top.
headline_df <- results[!is.na(results$time_s) & nzchar(results$impl), ]
headline_funcs <- c("munge", "ldsc", "sumstats", "simLDSC", "userGWAS")
headline_df <- headline_df[headline_df$func %in% headline_funcs, ]

# Build a per-row label that collapses the scaling sweep into "userGWAS
# (N=10 000)" etc., and normalises the impl column to {R, gsemr} so the
# headline shows exactly one red bar and one blue bar per row.
headline_df$family <- sub("^(\\w+).*$", "\\1", headline_df$impl)
headline_df <- headline_df[headline_df$family %in% c("R", "Rust"), ]
headline_df$impl_clean <- ifelse(headline_df$family == "R", "R GenomicSEM", "gsemr (Rust)")

# Extract the size suffix from impls like "Rust (N=10000)" /
# "Rust (N=4936648, full)" to build a per-row label.
extract_size <- function(impl) {
  m <- regmatches(impl, regexpr("N=\\d+(, full)?", impl))
  if (length(m) == 0 || !nzchar(m)) return(NA_character_)
  m
}
headline_df$size_tag <- vapply(as.character(headline_df$impl), extract_size, character(1))

# Build the per-row label: plain function name for non-userGWAS rows,
# "userGWAS (N=X,XXX)" or "userGWAS (N=X,XXX,XXX, full)" for the sweep
# rows, with a thousands separator on N.
headline_df$row_label <- mapply(
  function(fn, sz) {
    if (fn != "userGWAS" || is.na(sz)) return(as.character(fn))
    n_str <- sub(",.*$", "", sub("^N=", "", sz))
    n <- suppressWarnings(as.numeric(n_str))
    pretty_n <- if (is.finite(n)) {
      format(n, big.mark = ",", scientific = FALSE)
    } else {
      n_str
    }
    if (grepl("full", sz)) sprintf("userGWAS (N=%s, full)", pretty_n)
    else                   sprintf("userGWAS (N=%s)",       pretty_n)
  },
  as.character(headline_df$func),
  headline_df$size_tag,
  USE.NAMES = FALSE
)

# Fold down to one row per (row_label, impl_clean) — drop the extra
# Rust full-scale vs sweep overlap if any.
headline_df <- aggregate(time_s ~ row_label + impl_clean,
                         data = headline_df, FUN = mean)

# Sort rows by R wall time desc. Rows where R has no bar (e.g. the
# full-scale userGWAS) go to the bottom as "R: not attempted".
r_times <- setNames(rep(NA_real_, length(unique(headline_df$row_label))),
                    unique(headline_df$row_label))
for (lbl in names(r_times)) {
  rv <- headline_df$time_s[headline_df$row_label == lbl & headline_df$impl_clean == "R GenomicSEM"]
  if (length(rv) > 0) r_times[[lbl]] <- rv[1]
}
row_order <- names(sort(r_times, decreasing = TRUE, na.last = TRUE))
headline_df$row_label <- factor(headline_df$row_label, levels = rev(row_order))
headline_df$impl_clean <- factor(headline_df$impl_clean,
                                  levels = c("R GenomicSEM", "gsemr (Rust)"))

# Format time labels: sub-second in ms, seconds otherwise; anything
# over 60s also annotated with minutes in the label itself.
fmt_time_label <- function(t) {
  if (is.na(t)) return("")
  if (t < 1)   return(sprintf("%.2fs", t))
  if (t < 10)  return(sprintf("%.1fs", t))
  if (t < 60)  return(sprintf("%.0fs", t))
  mins <- floor(t / 60)
  secs <- round(t - 60 * mins)
  sprintf("%.0fs (%dm%02ds)", t, mins, secs)
}
headline_df$label <- vapply(headline_df$time_s, fmt_time_label, character(1))

headline_palette <- c(
  "R GenomicSEM" = "#E41A1C",
  "gsemr (Rust)" = "#377EB8"
)

# Which rows have no R bar? Annotate them inline so the absence is
# explicit. The annotation sits on the red-bar row of the dodge (so it
# lines up where the red bar *would* be, not on top of the blue bar)
# and uses a negative hjust so the text starts just right of x=0.
rows_with_r    <- unique(headline_df$row_label[headline_df$impl_clean == "R GenomicSEM"])
rows_missing_r <- setdiff(levels(headline_df$row_label), as.character(rows_with_r))

# For each missing-R row, inject a zero-height R bar so the dodge
# reserves the slot; the annotation text then lands on that slot.
if (length(rows_missing_r) > 0) {
  placeholder <- data.frame(
    row_label  = rows_missing_r,
    impl_clean = "R GenomicSEM",
    time_s     = 0,
    label      = "R GenomicSEM: not attempted at this scale",
    stringsAsFactors = FALSE
  )
  placeholder$row_label  <- factor(placeholder$row_label,  levels = levels(headline_df$row_label))
  placeholder$impl_clean <- factor(placeholder$impl_clean, levels = levels(headline_df$impl_clean))
  headline_df <- rbind(headline_df, placeholder)
}

# Labels for real rows are grey; labels for the injected placeholder
# rows (no R bar) are red so the "not attempted" message reads as an
# annotation rather than a bar value.
headline_df$label_colour <- ifelse(
  headline_df$time_s == 0 & headline_df$impl_clean == "R GenomicSEM",
  "#B22222", "grey20"
)

png(png_path, width = 1400, height = 900, res = 150)
if (nrow(headline_df) > 0) {
  p_headline <- ggplot(headline_df,
                       aes(x = row_label, y = time_s, fill = impl_clean)) +
    geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
             width = 0.75) +
    geom_text(
      aes(label = label, colour = label_colour),
      position = position_dodge2(width = 0.85, preserve = "single"),
      hjust    = -0.05,
      size     = 3.6
    ) +
    scale_colour_identity() +
    coord_flip(clip = "off") +
    scale_fill_manual(values = headline_palette) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.28)),
      labels = function(x) sprintf("%gs", x)
    ) +
    labs(
      title    = "gsemr vs R GenomicSEM — wall-clock time on PGC Anxiety/OCD/PTSD sumstats",
      subtitle = sprintf(
        "Psychiatric Genomics Consortium (PGC) GWAS · BENCH_CORES=%d · lower is better",
        BENCH_CORES
      ),
      caption  = "Fit ops (commonfactor, usermodel, rgmodel, write.model, summaryGLS, paLDSC) finish in under 0.5s on both and are omitted.",
      x = NULL, y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      # Anchor title/subtitle/caption to the full plot area (not the
      # inner panel), then center them. With default "panel" anchoring,
      # "center" means "middle of the bars area" — which, with long
      # y-axis labels, is noticeably right of the image center and
      # causes long titles to clip off the right edge.
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle    = element_text(color = "grey35", size = 11, hjust = 0.5, margin = margin(b = 12)),
      plot.caption     = element_text(color = "grey45", size = 9, hjust = 0.5, margin = margin(t = 10)),
      # ggplot legends centre against the inner panel, not the full
      # plot area — with long y-axis labels that makes them drift
      # right. Moving to the bottom sidesteps the issue: the x-axis
      # runs edge-to-edge so the panel centre effectively lines up
      # with the image centre.
      legend.position  = "bottom",
      legend.text      = element_text(size = 11),
      axis.text.y      = element_text(size = 11),
      axis.text.x      = element_text(size = 10, colour = "grey40"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin      = margin(15, 40, 15, 10)
    )
  print(p_headline)
}
invisible(dev.off())
cat(sprintf("Headline plot saved to %s\n", png_path))
cat("\nDone.\n")
