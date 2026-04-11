#!/usr/bin/env Rscript
# Reality-check the numbers in upstream-issue/ISSUE-parameterization.md.
#
# Goal: on the same 1000-SNP batch, the same covstruc, and the same
# conceptual single-factor model, measure how much
# GenomicSEM::commonfactorGWAS() and GenomicSEM::userGWAS() actually
# disagree. Produces the summary() / cor() / max-diff numbers we'd
# paste into the upstream issue, using only upstream R GenomicSEM (the
# comparison is R-vs-R, gsemr is not involved).
#
# Requires the upstream crash patch to be hot-applied, otherwise
# commonfactorGWAS crashes partway through the batch.
#
# Run from bench/.

setwd("/Users/pohsuanlai/Documents/GenomicSEM-rs/bench")
suppressPackageStartupMessages({
  library(GenomicSEM)
})
source("patch_genomicsem.R")

trait_names <- c("ANX", "OCD", "PTSD")
sample_prev <- c(0.5, 0.5, 0.5)
pop_prev    <- c(0.16, 0.02, 0.07)
ld          <- "data/eur_w_ld_chr/"
ref_1000g   <- "data/reference.1000G.maf.0.005.txt.gz"
munged      <- paste0(trait_names, ".sumstats.gz")
raw_traits  <- c("data/anxiety.meta.full.cc.tbl.gz",
                 "data/ocd_aug2017.gz",
                 "data/SORTED_PTSD_EA9_ALL_study_specific_PCs1.txt")

cache_dir <- "cache"
dir.create(cache_dir, showWarnings = FALSE)
cov_path <- file.path(cache_dir, "cov.rds")
ss_path  <- file.path(cache_dir, "sumstats.rds")

# --- 0. covstruc (cached) -----------------------------------------------------
if (file.exists(cov_path)) {
  cat("[cache] loading covstruc from", cov_path, "\n")
  covstruc <- readRDS(cov_path)
} else {
  if (!all(file.exists(munged))) {
    stop("Missing munged files: ", paste(munged[!file.exists(munged)], collapse=", "),
         " — run benchmark_perf.R first to generate them")
  }
  cat("[compute] ldsc on 3 traits...\n")
  t0 <- Sys.time()
  covstruc <- GenomicSEM::ldsc(
    munged, sample_prev, pop_prev, ld, ld,
    trait.names = trait_names, n.blocks = 200
  )
  cat(sprintf("[compute] ldsc took %.1fs\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  saveRDS(covstruc, cov_path)
}

# --- 1. sumstats (cached) — 3-trait merge -------------------------------------
if (file.exists(ss_path)) {
  cat("[cache] loading sumstats from", ss_path, "\n")
  ss_all <- readRDS(ss_path)
} else {
  # Use R's own sumstats() so there's no gsemr involvement in the inputs.
  # sumstats() wants the RAW GWAS files (with P and SE columns), not the
  # munged Z/N files that ldsc() consumes.
  cat("[compute] sumstats merge on 3 traits (this takes a few minutes)...\n")
  t0 <- Sys.time()
  ss_all <- GenomicSEM::sumstats(
    files = raw_traits, ref = ref_1000g,
    trait.names = trait_names,
    se.logit    = c(FALSE, FALSE, FALSE),
    OLS         = c(FALSE, FALSE, FALSE),
    linprob     = c(FALSE, FALSE, FALSE),
    info.filter = 0.6, maf.filter = 0.01
  )
  cat(sprintf("[compute] sumstats took %.1fs; %d rows\n",
              as.numeric(Sys.time() - t0, units = "secs"), nrow(ss_all)))
  saveRDS(ss_all, ss_path)
}

# --- 2. Apply the crash patch so commonfactorGWAS doesn't die mid-batch -------
cat("\n[patch] applying upstream crash fix...\n")
ok <- patch_genomicsem()
if (!ok) stop("patch_genomicsem failed")

# --- 3. Pick a 1000-SNP batch -------------------------------------------------
# Use a deterministic slice so the numbers are reproducible.
batch <- head(ss_all, 1000)
cat(sprintf("\n[input] batch = head(ss, 1000), %d SNPs\n", nrow(batch)))

# --- 4. Run commonfactorGWAS --------------------------------------------------
cat("\n==============================================================\n")
cat("  commonfactorGWAS (patched R)\n")
cat("==============================================================\n")
t0 <- Sys.time()
cf_out <- GenomicSEM::commonfactorGWAS(covstruc, SNPs = batch, parallel = FALSE)
cat(sprintf("elapsed: %.1fs, nrow = %d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(cf_out)))
cat("colnames:", paste(names(cf_out), collapse=", "), "\n")

# --- 5. Run userGWAS with the explicit single-factor model -------------------
# Same single-factor model the bench uses for userGWAS — fixed factor
# variance, all three loadings free. This is the fixed-variance
# parameterization that commonfactorGWAS does NOT use internally.
single_factor_model <- paste0(
  "F1 =~ NA*", trait_names[1], " + ", trait_names[2], " + ", trait_names[3], "\n",
  "F1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP"
)

cat("\n==============================================================\n")
cat("  userGWAS (fixed-variance parameterization)\n")
cat("==============================================================\n")
cat("model:\n", single_factor_model, "\n")
t0 <- Sys.time()
ug_out <- GenomicSEM::userGWAS(
  covstruc, SNPs = batch, model = single_factor_model,
  parallel = FALSE
)
cat(sprintf("elapsed: %.1fs, nrow(list) = %d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            length(ug_out)))

# --- 6. Extract F1 ~ SNP effect from userGWAS nested list ---------------------
ug_rows <- do.call(rbind, lapply(seq_along(ug_out), function(i) {
  x <- ug_out[[i]]
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
    return(data.frame(SNP = batch$SNP[i], est = NA_real_, se = NA_real_,
                      stringsAsFactors = FALSE))
  }
  r <- x[x$lhs == "F1" & x$op == "~" & x$rhs == "SNP", ]
  if (nrow(r) == 0) {
    return(data.frame(SNP = batch$SNP[i], est = NA_real_, se = NA_real_,
                      stringsAsFactors = FALSE))
  }
  # upstream R userGWAS columns: lhs, op, rhs, est, SE, Z_Estimate, ...
  data.frame(
    SNP = as.character(r$SNP[1]),
    est = as.numeric(r$est[1]),
    se  = as.numeric(r$SE[1]),
    stringsAsFactors = FALSE
  )
}))

# --- 7. Join and compare ------------------------------------------------------
cat("\n==============================================================\n")
cat("  commonfactorGWAS vs userGWAS on the same batch\n")
cat("==============================================================\n")

# commonfactorGWAS columns: SNP, est, se, se_c, Z_Estimate, Pval_Estimate, ...
# We compare `est` (the raw estimate) and `se` vs the userGWAS-extracted
# Unstand_Est / Unstand_SE for the F1 ~ SNP row.

cmp <- merge(
  data.frame(SNP = cf_out$SNP, cf_est = cf_out$est, cf_se = cf_out$se),
  ug_rows,
  by = "SNP", all = FALSE
)
names(cmp) <- c("SNP", "cf_est", "cf_se", "ug_est", "ug_se")

valid <- is.finite(cmp$cf_est) & is.finite(cmp$ug_est) &
         is.finite(cmp$cf_se)  & is.finite(cmp$ug_se)
cat(sprintf("\njoined rows: %d  (of which finite on both sides: %d)\n",
            nrow(cmp), sum(valid)))

v <- cmp[valid, ]
v$est_diff <- v$ug_est - v$cf_est
v$se_diff  <- v$ug_se  - v$cf_se

cat("\n--- est: cf vs ug ---\n")
cat("  summary(ug_est - cf_est):\n")
print(summary(v$est_diff))
cat(sprintf("\n  cor(cf_est, ug_est) = %.4f\n", cor(v$cf_est, v$ug_est)))
cat(sprintf("  max |diff|          = %.4e\n", max(abs(v$est_diff))))
cat(sprintf("  median |diff|       = %.4e\n", median(abs(v$est_diff))))
cat(sprintf("  mean |diff|         = %.4e\n", mean(abs(v$est_diff))))
cat(sprintf("  n SNPs where sign differs: %d / %d\n",
            sum(sign(v$cf_est) != sign(v$ug_est)), nrow(v)))

cat("\n--- se: cf vs ug ---\n")
cat("  summary(ug_se - cf_se):\n")
print(summary(v$se_diff))
cat(sprintf("  cor(cf_se, ug_se) = %.4f\n", cor(v$cf_se, v$ug_se)))

cat("\n--- 10 worst |est| disagreements ---\n")
ord <- order(-abs(v$est_diff))
print(head(v[ord, c("SNP","cf_est","ug_est","cf_se","ug_se","est_diff")], 10))

# --- 8. Save a compact JSON-ish summary the issue can quote -------------------
results_path <- "diag_cfgwas_vs_usergwas_summary.txt"
sink(results_path)
cat("Reality-check summary for upstream-issue/ISSUE-parameterization.md\n")
cat("==================================================================\n\n")
cat("Input: head(ss, 1000) from the 3-trait PGC merge (ANX/OCD/PTSD)\n")
cat("covstruc: GenomicSEM::ldsc with n.blocks=200, sample.prev=0.5 each,\n")
cat("  pop.prev=c(0.16, 0.02, 0.07)\n")
cat("R GenomicSEM: upstream master (commit 123ec96) + crash-fix hot-patch\n\n")
cat("Model (for userGWAS):\n")
cat("  F1 =~ NA*V1 + V2 + V3\n")
cat("  F1 ~ SNP\n")
cat("  F1 ~~ 1*F1\n")
cat("  SNP ~~ SNP\n\n")
cat(sprintf("SNPs joined on both sides: %d of %d\n", sum(valid), nrow(cmp)))
cat("\n--- ug_est - cf_est ---\n")
print(summary(v$est_diff))
cat(sprintf("\ncor(cf_est, ug_est) = %.4f\n", cor(v$cf_est, v$ug_est)))
cat(sprintf("max |diff|          = %.4e\n", max(abs(v$est_diff))))
cat(sprintf("median |diff|       = %.4e\n", median(abs(v$est_diff))))
cat(sprintf("mean |diff|         = %.4e\n", mean(abs(v$est_diff))))
cat(sprintf("sign disagreements  = %d / %d\n",
            sum(sign(v$cf_est) != sign(v$ug_est)), nrow(v)))
cat("\n--- ug_se - cf_se ---\n")
print(summary(v$se_diff))
cat(sprintf("cor(cf_se, ug_se) = %.4f\n", cor(v$cf_se, v$ug_se)))
cat("\n--- top 10 worst |est| disagreements ---\n")
print(head(v[ord, c("SNP","cf_est","ug_est","cf_se","ug_se","est_diff")], 10),
      row.names = FALSE)
sink()
cat(sprintf("\nwrote %s\n", results_path))
