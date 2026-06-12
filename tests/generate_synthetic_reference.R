#!/usr/bin/env Rscript
#
# Generate synthetic-data reference fixtures for the LDSC half of the
# pipeline (munge, sumstats, ldsc) by running R GenomicSEM on small,
# fully reproducible synthetic inputs.
#
# Unlike generate_reference.R (matrix/SEM ops via lavaan), this script
# exercises the GWAS-summary-statistics I/O + LD-score-regression path.
# All inputs are simulated from the standard LDSC generative model so the
# committed input files are small (~hundreds of KB) and CI-runnable, while
# still flowing through the exact same code paths as real PGC data.
#
# Requires: GenomicSEM, jsonlite
# Usage: cd tests && Rscript generate_synthetic_reference.R
#
# Output:
#   tests/fixtures/synth/            committed synthetic INPUT files
#   tests/fixtures/{ldsc,munge,sumstats}_synth.json   R reference OUTPUTS

suppressMessages({
  library(jsonlite)
  library(GenomicSEM)
})

set.seed(20240611)

outdir   <- "fixtures"
synthdir <- file.path(outdir, "synth")
lddir    <- file.path(synthdir, "eur_w_ld_chr")
dir.create(lddir, showWarnings = FALSE, recursive = TRUE)

mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))
write_fixture <- function(data, name) {
  path <- file.path(outdir, paste0(name, ".json"))
  writeLines(toJSON(data, auto_unbox = TRUE, digits = 15), path)
  cat("Wrote", path, "\n")
}

# ---------------------------------------------------------------------------
# 1. Simulate SNPs + LD scores across a couple of chromosomes
# ---------------------------------------------------------------------------
chr           <- 2
snps_per_chr  <- 900
n_traits      <- 3
trait_names   <- c("T1", "T2", "T3")

bases <- c("A", "C", "G", "T")
make_snp_table <- function() {
  rows <- list()
  for (c in seq_len(chr)) {
    bp <- sort(sample(1e6:9e6, snps_per_chr))
    for (j in seq_len(snps_per_chr)) {
      a1 <- sample(bases, 1)
      a2 <- sample(setdiff(bases, a1), 1)
      rows[[length(rows) + 1]] <- data.frame(
        CHR = c,
        SNP = sprintf("rs%d", length(rows) + 1 + 1000000),
        BP  = bp[j],
        A1  = a1,
        A2  = a2,
        L2  = round(runif(1, 1, 60), 4),   # LD score
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
snp_tab <- make_snp_table()
M_snps  <- nrow(snp_tab)
M_5_50  <- M_snps   # normalization count (all SNPs "common" in this sim)

# Per-chromosome LD score files (CHR SNP BP L2) + M files
for (c in seq_len(chr)) {
  sub <- snp_tab[snp_tab$CHR == c, c("CHR", "SNP", "BP", "L2")]
  gz <- gzfile(file.path(lddir, paste0(c, ".l2.ldscore.gz")), "w")
  write.table(sub, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
  m_c <- sum(snp_tab$CHR == c)
  writeLines(as.character(m_c), file.path(lddir, paste0(c, ".l2.M_5_50")))
  writeLines(as.character(m_c), file.path(lddir, paste0(c, ".l2.M")))
}

# HapMap3-style SNP list (SNP A1 A2) and reference panel (SNP A1 A2 MAF)
hm3 <- snp_tab[, c("SNP", "A1", "A2")]
write.table(hm3, file.path(synthdir, "w_hm3.snplist"),
            sep = "\t", quote = FALSE, row.names = FALSE)

maf <- round(runif(M_snps, 0.05, 0.5), 4)
ref <- data.frame(SNP = snp_tab$SNP, A1 = snp_tab$A1, A2 = snp_tab$A2, MAF = maf,
                  stringsAsFactors = FALSE)
gzr <- gzfile(file.path(synthdir, "reference.txt.gz"), "w")
write.table(ref, gzr, sep = "\t", quote = FALSE, row.names = FALSE)
close(gzr)

# ---------------------------------------------------------------------------
# 2. Simulate GWAS z-scores under the LDSC model, parameterized directly by
#    the per-trait LD-score regression slope so the synthetic data has a
#    realistic mean chi-square (~1 + slope*mean(L2) ~ 1.3) rather than the
#    absurd inflation a tiny synthetic M would otherwise produce:
#       diag(Sigma_i)    = 1 + slope_j * L2_i
#       offdiag(Sigma_i) = rg_jk * sqrt(slope_j slope_k) * L2_i
#    LDSC recovers genetic covariance S_jk = sqrt(Nj Nk)/M * slope_jk.
# ---------------------------------------------------------------------------
slope <- c(0.012, 0.010, 0.014)   # ~mean(L2)~30 -> mean chi-square ~1.3-1.4
Nj    <- c(50000, 40000, 60000)
rg    <- matrix(c(1.0, 0.40, 0.30,
                  0.40, 1.0, 0.50,
                  0.30, 0.50, 1.0), 3, 3)

Z <- matrix(0, M_snps, n_traits)
for (i in seq_len(M_snps)) {
  L2 <- snp_tab$L2[i]
  Sig <- matrix(0, n_traits, n_traits)
  for (a in seq_len(n_traits)) for (b in seq_len(n_traits)) {
    if (a == b) {
      Sig[a, b] <- 1 + slope[a] * L2
    } else {
      Sig[a, b] <- rg[a, b] * sqrt(slope[a] * slope[b]) * L2
    }
  }
  # symmetric sqrt for the MVN draw
  e <- eigen(Sig, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 1e-8))) %*% t(e$vectors)
  Z[i, ] <- A %*% rnorm(n_traits)
}

# ---------------------------------------------------------------------------
# 3. Write munged sumstats (SNP A1 A2 Z N) for the ldsc test, and raw GWAS
#    (SNP A1 A2 BETA P N MAF) for the munge test.
# ---------------------------------------------------------------------------
munged_paths <- character(n_traits)
for (j in seq_len(n_traits)) {
  ss <- data.frame(SNP = snp_tab$SNP, A1 = snp_tab$A1, A2 = snp_tab$A2,
                   Z = round(Z[, j], 6), N = Nj[j], stringsAsFactors = FALSE)
  p <- file.path(synthdir, paste0(trait_names[j], ".sumstats.gz"))
  gz <- gzfile(p, "w")
  write.table(ss, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
  munged_paths[j] <- p
}

# Raw GWAS for the munge test (trait 1 only is enough to validate munge).
# beta carries the sign of Z; P is the two-sided p-value -> munge recovers
# Z = sign(beta) * |qnorm(P/2)| == original Z.
raw_paths <- character(n_traits)
for (j in seq_len(n_traits)) {
  beta <- Z[, j] / sqrt(Nj[j])          # plausible effect-size scale
  se   <- rep(1 / sqrt(Nj[j]), M_snps)  # so beta/se == Z exactly
  pval <- 2 * pnorm(-abs(Z[, j]))
  pval[pval < 1e-300] <- 1e-300
  # Full precision: R derives Z from P, gsem from BETA/SE (different columns),
  # so rounding either would cap the achievable agreement artificially.
  raw <- data.frame(SNP = snp_tab$SNP, A1 = snp_tab$A1, A2 = snp_tab$A2,
                    BETA = beta, SE = se, P = pval,
                    N = Nj[j], MAF = maf, stringsAsFactors = FALSE)
  p <- file.path(synthdir, paste0(trait_names[j], ".raw.gz"))
  gz <- gzfile(p, "w")
  write.table(raw, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
  raw_paths[j] <- p
}

# ---------------------------------------------------------------------------
# 4. Run R GenomicSEM munge -> reference output
# ---------------------------------------------------------------------------
cat("=== munge ===\n")
owd <- getwd()
munge_tmp <- file.path(tempdir(), "munge_out")
dir.create(munge_tmp, showWarnings = FALSE)
# GenomicSEM::munge writes <trait>.sumstats.gz into the working dir
raw1_abs <- normalizePath(raw_paths[1])
hm3_abs  <- normalizePath(file.path(synthdir, "w_hm3.snplist"))
setwd(munge_tmp)
suppressMessages(
  GenomicSEM::munge(
    files       = raw1_abs,
    hm3         = hm3_abs,
    trait.names = "MUNGE1",
    N           = Nj[1],
    info.filter = 0.0,
    maf.filter  = 0.01
  )
)
setwd(owd)
munged_r <- suppressWarnings(
  read.table(file.path(munge_tmp, "MUNGE1.sumstats.gz"), header = TRUE,
             stringsAsFactors = FALSE))
# Sort by SNP for a deterministic comparison
munged_r <- munged_r[order(munged_r$SNP), ]
write_fixture(list(
  raw_file   = "synth/T1.raw.gz",
  hm3_file   = "synth/w_hm3.snplist",
  n          = Nj[1],
  snp        = munged_r$SNP,
  a1         = toupper(munged_r$A1),
  a2         = toupper(munged_r$A2),
  z          = as.numeric(munged_r$Z),
  n_out      = as.numeric(munged_r$N)
), "munge_synth")

# ---------------------------------------------------------------------------
# 5. Run R GenomicSEM ldsc -> S, V, I reference
# ---------------------------------------------------------------------------
cat("=== ldsc ===\n")
ld <- lddir
r_ldsc <- suppressMessages(GenomicSEM::ldsc(
  traits          = munged_paths,
  sample.prev     = rep(NA, n_traits),
  population.prev = rep(NA, n_traits),
  ld              = ld,
  wld             = ld,
  trait.names     = trait_names,
  chr             = chr,
  n.blocks        = 20,
  stand           = FALSE
))
write_fixture(list(
  munged_files = paste0("synth/", trait_names, ".sumstats.gz"),
  ld_dir       = "synth/eur_w_ld_chr",
  chr          = chr,
  n_blocks     = 20,
  trait_names  = trait_names,
  m_total      = M_5_50,
  s            = mat_to_list(r_ldsc$S),
  v            = mat_to_list(r_ldsc$V),
  i            = mat_to_list(r_ldsc$I)
), "ldsc_synth")

# ---------------------------------------------------------------------------
# 6. Run R GenomicSEM sumstats -> merged per-SNP betas/SEs reference
# ---------------------------------------------------------------------------
cat("=== sumstats (all standardization modes) ===\n")
# Run R sumstats once per standardization mode on the SAME raw files. The
# synthetic effects are continuous, so linprob/se.logit/none are not
# biologically meaningful here — but R and gsem apply identical formulas, so
# this still validates formula-for-formula equivalence of every mode.
run_mode <- function(ols, linprob, se_logit) {
  r <- suppressMessages(GenomicSEM::sumstats(
    files       = raw_paths,
    ref         = file.path(synthdir, "reference.txt.gz"),
    trait.names = trait_names,
    se.logit    = rep(se_logit, n_traits),
    OLS         = rep(ols, n_traits),
    linprob     = rep(linprob, n_traits),
    N           = Nj,
    betas       = NULL,
    info.filter = 0.0,
    maf.filter  = 0.01
  ))
  r <- r[order(r$SNP), ]
  beta_cols <- grep("^beta\\.", colnames(r), value = TRUE)
  se_cols   <- grep("^se\\.",   colnames(r), value = TRUE)
  list(r = r, beta = beta_cols, se = se_cols)
}

modes <- list(
  ols     = run_mode(TRUE,  FALSE, FALSE),
  linprob = run_mode(FALSE, TRUE,  FALSE),
  se_logit= run_mode(FALSE, FALSE, TRUE),
  none    = run_mode(FALSE, FALSE, FALSE)
)
# All modes operate on the same merged SNP set / alleles.
base <- modes$ols$r
mode_fixture <- function(m) list(
  beta = mat_to_list(as.matrix(m$r[, m$beta])),
  se   = mat_to_list(as.matrix(m$r[, m$se]))
)
write_fixture(list(
  raw_files    = paste0("synth/", trait_names, ".raw.gz"),
  ref_file     = "synth/reference.txt.gz",
  trait_names  = trait_names,
  n            = Nj,
  snp          = base$SNP,
  a1           = toupper(base$A1),
  a2           = toupper(base$A2),
  # Back-compat: top-level beta/se are the OLS mode.
  beta         = mat_to_list(as.matrix(modes$ols$r[, modes$ols$beta])),
  se           = mat_to_list(as.matrix(modes$ols$r[, modes$ols$se])),
  modes        = list(
    ols      = mode_fixture(modes$ols),
    linprob  = mode_fixture(modes$linprob),
    se_logit = mode_fixture(modes$se_logit),
    none     = mode_fixture(modes$none)
  )
), "sumstats_synth")

# GenomicSEM munge/sumstats write *.log files into the working dir; remove
# them so re-running the generator leaves a clean tree.
unlink(list.files(".", pattern = "_(munge|sumstats)\\.log$", full.names = TRUE))
unlink(list.files(synthdir, pattern = "\\.log$", full.names = TRUE))

cat("\n=== Synthetic reference fixtures generated ===\n")
