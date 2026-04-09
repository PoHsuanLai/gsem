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
trait_names <- c("ANX", "OCD", "PTSD")
sample_prev <- c(0.5, 0.5, 0.5)
pop_prev    <- c(0.16, 0.02, 0.07)
n_overrides <- c(NA, 9725, 5831)
munged_traits <- paste0(trait_names, ".sumstats.gz")

GWAS_SUBSET_N <- 1000L

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Benchmark a single call: returns list(time_s, peak_mb, error)
run_bench <- function(expr_fn) {
  gc(reset = TRUE)
  t <- tryCatch(
    system.time({ result <- expr_fn() }),
    error = function(e) e
  )
  if (inherits(t, "error")) {
    return(list(time_s = NA_real_, peak_mb = NA_real_, error = conditionMessage(t)))
  }
  mem <- gc()
  peak_mb <- mem[2, 6]  # max used Mb (Vcells)
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
file.remove(paste0(munge_r_names, ".sumstats.gz")[file.exists(paste0(munge_r_names, ".sumstats.gz"))])

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

# ---------------------------------------------------------------------------
# 3. commonfactor
# ---------------------------------------------------------------------------
cat("[3/12] commonfactor\n")
b <- run_bench(function() GenomicSEM::commonfactor(r_cov, estimation = "DWLS"))
add_result("commonfactor", "R", b)

b <- run_bench(function() gsemr::commonfactor(rust_cov, estimation = "DWLS"))
add_result("commonfactor", "Rust", b)

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

# ---------------------------------------------------------------------------
# 5. rgmodel
# ---------------------------------------------------------------------------
cat("[5/12] rgmodel\n")
b <- run_bench(function() GenomicSEM::rgmodel(r_cov))
add_result("rgmodel", "R", b)

b <- run_bench(function() gsemr::rgmodel(rust_cov))
add_result("rgmodel", "Rust", b)

# ---------------------------------------------------------------------------
# 6. sumstats
# ---------------------------------------------------------------------------
cat("[6/12] sumstats\n")
dir.create("out_bench", showWarnings = FALSE)

# R sumstats — may error on w_hm3.snplist (no MAF column)
b <- run_bench(function() {
  GenomicSEM::sumstats(
    files = raw_traits, ref = hm3, trait.names = trait_names,
    se.logit = c(FALSE, FALSE, FALSE),
    OLS = c(FALSE, FALSE, FALSE),
    linprob = c(FALSE, FALSE, FALSE),
    info.filter = 0.6, maf.filter = 0.01
  )
})
add_result("sumstats", "R", b)

# Rust sumstats
rust_ss_path <- "out_bench/rust_merged.tsv"
b <- run_bench(function() {
  gsemr::sumstats(
    files = raw_traits, ref = hm3, trait.names = trait_names,
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path
  )
})
add_result("sumstats", "Rust", b)

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

# ---------------------------------------------------------------------------
# Prepare GWAS subset for commonfactorGWAS / userGWAS
# ---------------------------------------------------------------------------
cat("Preparing GWAS subset...\n")
if (!file.exists(rust_ss_path)) {
  gsemr::sumstats(
    files = raw_traits, ref = hm3, trait.names = trait_names,
    info.filter = 0.6, maf.filter = 0.01, out = rust_ss_path
  )
}
subset_path <- "out_bench/merged_subset.tsv"
full_snps <- read.delim(rust_ss_path, nrows = GWAS_SUBSET_N)
write.table(full_snps, subset_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Using %d SNPs for GWAS benchmarks.\n", nrow(full_snps)))

# ---------------------------------------------------------------------------
# 8. commonfactorGWAS
# ---------------------------------------------------------------------------
cat("[8/12] commonfactorGWAS\n")
b <- run_bench(function() {
  GenomicSEM::commonfactorGWAS(covstruc = r_cov, SNPs = subset_path)
})
add_result("commonfactorGWAS", "R", b)

b <- run_bench(function() {
  gsemr::commonfactorGWAS(covstruc = rust_cov, SNPs = subset_path)
})
add_result("commonfactorGWAS", "Rust", b)

# ---------------------------------------------------------------------------
# 9. userGWAS
# ---------------------------------------------------------------------------
cat("[9/12] userGWAS\n")
gwas_model <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

b <- run_bench(function() {
  GenomicSEM::userGWAS(covstruc = r_cov, SNPs = subset_path, model = gwas_model)
})
add_result("userGWAS", "R", b)

b <- run_bench(function() {
  gsemr::userGWAS(covstruc = rust_cov, SNPs = subset_path, model = gwas_model)
})
add_result("userGWAS", "Rust", b)

# ---------------------------------------------------------------------------
# 10. paLDSC
# ---------------------------------------------------------------------------
cat("[10/12] paLDSC\n")
b <- run_bench(function() {
  # R GenomicSEM may not have paLDSC; wrap safely
  if (existsFunction("paLDSC", where = asNamespace("GenomicSEM"))) {
    GenomicSEM::paLDSC(r_cov$S, r_cov$V, r = 100, p = 3, save.pdf = FALSE)
  } else {
    stop("paLDSC not available in R GenomicSEM")
  }
})
add_result("paLDSC", "R", b)

b <- run_bench(function() {
  gsemr::paLDSC(rust_cov$S, rust_cov$V, r = 100, p = 3, save.pdf = FALSE)
})
add_result("paLDSC", "Rust", b)

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
cat(sprintf("%-20s %10s %10s %10s %10s %10s\n",
            "Function", "R (s)", "Rust (s)", "Speedup", "R (MB)", "Rust (MB)"))
cat(strrep("-", 72), "\n")

funcs <- unique(results$func)
for (fn in funcs) {
  r_row    <- results[results$func == fn & results$impl == "R", ]
  rust_row <- results[results$func == fn & results$impl == "Rust", ]
  r_time   <- if (nrow(r_row) && !is.na(r_row$time_s[1])) r_row$time_s[1] else NA
  rust_time <- if (nrow(rust_row) && !is.na(rust_row$time_s[1])) rust_row$time_s[1] else NA
  speedup  <- if (!is.na(r_time) && !is.na(rust_time) && rust_time > 0) r_time / rust_time else NA
  r_mem    <- if (nrow(r_row) && !is.na(r_row$peak_mb[1])) r_row$peak_mb[1] else NA
  rust_mem <- if (nrow(rust_row) && !is.na(rust_row$peak_mb[1])) rust_row$peak_mb[1] else NA

  fmt_val <- function(x) if (is.na(x)) "   N/A" else sprintf("%10.2f", x)
  fmt_spd <- function(x) if (is.na(x)) "   N/A" else sprintf("%9.1fx", x)

  # Print errors if any
  r_err <- if (nrow(r_row) && nzchar(r_row$error[1])) r_row$error[1] else ""
  rust_err <- if (nrow(rust_row) && nzchar(rust_row$error[1])) rust_row$error[1] else ""

  cat(sprintf("%-20s %s %s %s %s %s",
              fn, fmt_val(r_time), fmt_val(rust_time), fmt_spd(speedup),
              fmt_val(r_mem), fmt_val(rust_mem)))
  if (nzchar(r_err)) cat(sprintf("  [R err: %s]", substr(r_err, 1, 40)))
  if (nzchar(rust_err)) cat(sprintf("  [Rust err: %s]", substr(rust_err, 1, 40)))
  cat("\n")
}

# ---------------------------------------------------------------------------
# Generate plots
# ---------------------------------------------------------------------------
cat("\nGenerating plots...\n")
pdf_path <- "benchmark_plots.pdf"

# Prepare data for plotting: only functions where both R and Rust succeeded
plot_data <- do.call(rbind, lapply(funcs, function(fn) {
  r_row    <- results[results$func == fn & results$impl == "R", ]
  rust_row <- results[results$func == fn & results$impl == "Rust", ]
  r_time   <- if (nrow(r_row) && !is.na(r_row$time_s[1])) r_row$time_s[1] else NA
  rust_time <- if (nrow(rust_row) && !is.na(rust_row$time_s[1])) rust_row$time_s[1] else NA
  r_mem    <- if (nrow(r_row) && !is.na(r_row$peak_mb[1])) r_row$peak_mb[1] else NA
  rust_mem <- if (nrow(rust_row) && !is.na(rust_row$peak_mb[1])) rust_row$peak_mb[1] else NA
  data.frame(
    func      = fn,
    r_time    = r_time,
    rust_time = rust_time,
    r_mem     = r_mem,
    rust_mem  = rust_mem,
    stringsAsFactors = FALSE
  )
}))

# Filter to functions with both times available
time_data <- plot_data[!is.na(plot_data$r_time) & !is.na(plot_data$rust_time), ]
time_data$speedup <- time_data$r_time / time_data$rust_time

# Long format for time comparison
time_long <- rbind(
  data.frame(func = time_data$func, impl = "R GenomicSEM", time_s = time_data$r_time,
             stringsAsFactors = FALSE),
  data.frame(func = time_data$func, impl = "gsemr (Rust)", time_s = time_data$rust_time,
             stringsAsFactors = FALSE)
)
# Order functions by speedup descending
func_order <- time_data$func[order(time_data$speedup, decreasing = TRUE)]
time_long$func <- factor(time_long$func, levels = func_order)
time_data$func_f <- factor(time_data$func, levels = func_order)

# Theme
theme_bench <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

colors_impl <- c("R GenomicSEM" = "#E41A1C", "gsemr (Rust)" = "#377EB8")

pdf(pdf_path, width = 11, height = 8)

# --- Plot 1: Speed comparison (log scale) ---
if (nrow(time_long) > 0) {
  p1 <- ggplot(time_long, aes(x = func, y = time_s, fill = impl)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
    scale_y_log10(labels = function(x) ifelse(x >= 1, sprintf("%.1fs", x), sprintf("%.0fms", x * 1000))) +
    scale_fill_manual(values = colors_impl) +
    labs(
      title = "Execution Time: R GenomicSEM vs gsemr (Rust)",
      subtitle = sprintf("PGC Anxiety/OCD/PTSD data | GWAS subset: %d SNPs | Log scale", GWAS_SUBSET_N),
      x = NULL, y = "Time (log scale)", fill = NULL
    ) +
    theme_bench
  print(p1)
}

# --- Plot 2: Speedup factor ---
if (nrow(time_data) > 0) {
  p2 <- ggplot(time_data, aes(x = func_f, y = speedup)) +
    geom_col(fill = "#377EB8", alpha = 0.85, width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = sprintf("%.1fx", speedup)),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Speedup Factor: gsemr (Rust) over R GenomicSEM",
      subtitle = "Values > 1 mean Rust is faster",
      x = NULL, y = "Speedup (R time / Rust time)"
    ) +
    theme_bench
  print(p2)
}

# --- Plot 3: Memory comparison ---
mem_data <- plot_data[!is.na(plot_data$r_mem) & !is.na(plot_data$rust_mem), ]
if (nrow(mem_data) > 0) {
  mem_long <- rbind(
    data.frame(func = mem_data$func, impl = "R GenomicSEM", peak_mb = mem_data$r_mem,
               stringsAsFactors = FALSE),
    data.frame(func = mem_data$func, impl = "gsemr (Rust)", peak_mb = mem_data$rust_mem,
               stringsAsFactors = FALSE)
  )
  mem_long$func <- factor(mem_long$func, levels = func_order[func_order %in% mem_data$func])

  p3 <- ggplot(mem_long, aes(x = func, y = peak_mb, fill = impl)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
    scale_fill_manual(values = colors_impl) +
    labs(
      title = "Peak Memory Usage: R GenomicSEM vs gsemr (Rust)",
      subtitle = "Measured via gc() max used Mb (Vcells)",
      x = NULL, y = "Peak Memory (MB)", fill = NULL
    ) +
    theme_bench
  print(p3)
}

dev.off()
cat(sprintf("Plots saved to %s\n", pdf_path))
cat("\nDone.\n")
