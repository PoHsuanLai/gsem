#!/usr/bin/env Rscript
# Debug LDSC: print intermediate values at each step for trait 1 (h2).
library(GenomicSEM)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))[1]
if (!is.na(script_dir) && nzchar(script_dir)) setwd(script_dir)

# Read just trait 1
traits <- "data/iter1GWAS1.sumstats.gz"
ld_path <- "data/eur_w_ld_chr/"

# Manually replicate ldsc() internals for trait 1
suppressMessages(library(data.table))

# Read LD scores
ld_files <- paste0(ld_path, 1:22, ".l2.ldscore.gz")
ld_all <- do.call(rbind, lapply(ld_files, fread))
# Read weight LD
wld_all <- ld_all  # same dir

# Read M counts
m_vals <- sapply(1:22, function(i) {
  as.numeric(strsplit(readLines(paste0(ld_path, i, ".l2.M_5_50")), "\\s+")[[1]])
})
M.tot <- sum(unlist(m_vals))

cat(sprintf("M.tot = %f\n", M.tot))
cat(sprintf("n_ld_snps = %d\n", nrow(ld_all)))

# Read trait 1
trait1 <- fread(traits)
cat(sprintf("n_trait1_snps = %d\n", nrow(trait1)))

# Merge
merged <- merge(trait1, ld_all, by = "SNP")
cat(sprintf("n_merged = %d\n", nrow(merged)))

# Use L2 from LD file (may be L2.y after merge if trait also has L2)
if ("L2.y" %in% names(merged)) {
  merged$L2 <- merged$L2.y
} else if (!"L2" %in% names(merged)) {
  stop("No L2 column found")
}
merged$chi1 <- merged$Z^2
merged$N_trait <- merged$N

# Chi-sq filter
chisq_max <- max(0.001 * max(merged$N_trait), 80)
cat(sprintf("chisq_max = %f\n", chisq_max))
n_before <- nrow(merged)
merged <- merged[merged$chi1 < chisq_max, ]
cat(sprintf("n_after_chisq_filter = %d (removed %d)\n", nrow(merged), n_before - nrow(merged)))

n.snps <- nrow(merged)
N.bar <- mean(merged$N_trait)
cat(sprintf("N.bar = %f\n", N.bar))
cat(sprintf("mean(chi2) = %f\n", mean(merged$chi1)))
cat(sprintf("mean(L2) = %f\n", mean(merged$L2)))
cat(sprintf("mean(L2*N) = %f\n", mean(merged$L2 * merged$N_trait)))

# Weights
tot.agg <- (M.tot * (mean(merged$chi1) - 1)) / mean(merged$L2 * merged$N_trait)
tot.agg <- max(tot.agg, 0)
tot.agg <- min(tot.agg, 1)
cat(sprintf("tot.agg = %.10f\n", tot.agg))

merged$ld <- pmax(merged$L2, 1)
merged$w.ld <- pmax(merged$L2, 1)  # using same as LD for weights
merged$c <- tot.agg * merged$N_trait / M.tot
merged$het.w <- 1 / (2 * (1 + (merged$c * merged$ld))^2)
merged$oc.w <- 1 / merged$w.ld
merged$w <- merged$het.w * merged$oc.w
merged$initial.w <- sqrt(merged$w)
merged$weights <- merged$initial.w / sum(merged$initial.w)

cat(sprintf("mean(weights) = %.15e\n", mean(merged$weights)))
cat(sprintf("sum(weights) = %.15e\n", sum(merged$weights)))
cat(sprintf("weights[1:5] = %s\n", paste(sprintf("%.10e", merged$weights[1:5]), collapse=", ")))

# Weighted regression
weighted.LD <- cbind(merged$L2, 1) * merged$weights
weighted.chi <- merged$chi1 * merged$weights

xtx <- t(weighted.LD) %*% weighted.LD
xty <- t(weighted.LD) %*% weighted.chi
reg <- solve(xtx, xty)

cat(sprintf("xtx = [[%e, %e], [%e, %e]]\n", xtx[1,1], xtx[1,2], xtx[2,1], xtx[2,2]))
cat(sprintf("xty = [%e, %e]\n", xty[1], xty[2]))
cat(sprintf("reg (slope, intercept) = [%e, %e]\n", reg[1], reg[2]))

intercept <- reg[2]
coefs <- reg[1] / N.bar
h2 <- coefs * M.tot

cat(sprintf("\n=== FINAL ===\n"))
cat(sprintf("slope (raw) = %.10f\n", reg[1]))
cat(sprintf("intercept   = %.10f\n", intercept))
cat(sprintf("h2 = slope / N.bar * M = %.10f / %.1f * %.1f = %.10f\n",
            reg[1], N.bar, M.tot, h2))
