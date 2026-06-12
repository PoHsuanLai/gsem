#!/usr/bin/env Rscript
#
# Reference fixture for stratified LDSC (s_ldsc) with OVERLAPPING annotations.
# Validates the overlap-aware partitioned heritability: per-annotation
# S = overlap.matrix %*% (tau * M). Both R and gsem run the identical
# partitioned regression on the same annotation LD scores, so equivalence
# holds by construction — the data need only be well-conditioned and the
# annotations overlapping (base covering all SNPs + two binary categories).
#
# Requires: GenomicSEM, jsonlite
# Usage: cd tests && Rscript generate_sldsc_reference.R

suppressMessages({ library(jsonlite); library(GenomicSEM) })
set.seed(303)

outdir   <- "fixtures"
synthdir <- file.path(outdir, "sldsc")
annotdir <- file.path(synthdir, "annot")     # ld (annotation ldscores + M + annot.gz + frq)
wlddir   <- file.path(synthdir, "weights")   # wld
dir.create(annotdir, showWarnings = FALSE, recursive = TRUE)
dir.create(wlddir,   showWarnings = FALSE, recursive = TRUE)

mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))
write_fixture <- function(data, name) {
  writeLines(toJSON(data, auto_unbox = TRUE, digits = 15),
             file.path(outdir, paste0(name, ".json")))
  cat("Wrote", file.path(outdir, paste0(name, ".json")), "\n")
}

chr          <- 2
snps_per_chr <- 700
n_traits     <- 2
trait_names  <- c("T1", "T2")
annot_names  <- c("base", "A", "B")        # base = all SNPs; A,B overlapping binary
n_annot      <- length(annot_names)
bases <- c("A", "C", "G", "T")

# ---- SNP table with annotation membership + per-annotation LD scores --------
rows <- list()
gid <- 0
for (c in seq_len(chr)) {
  bp <- sort(sample(1e6:9e6, snps_per_chr))
  for (j in seq_len(snps_per_chr)) {
    gid <- gid + 1
    a1 <- sample(bases, 1); a2 <- sample(setdiff(bases, a1), 1)
    inA <- as.integer(runif(1) < 0.45)
    inB <- as.integer(runif(1) < 0.35)
    # Per-annotation LD scores: base is continuous; A/B LD scores are nonzero
    # mostly where the SNP is in that annotation, plus a little leakage so the
    # design matrix is well-conditioned (not perfectly collinear with membership).
    baseL2 <- round(runif(1, 10, 50), 4)
    AL2 <- round(inA * baseL2 * runif(1, 0.3, 0.6) + runif(1, 0, 1.5), 4)
    BL2 <- round(inB * baseL2 * runif(1, 0.3, 0.6) + runif(1, 0, 1.5), 4)
    rows[[gid]] <- data.frame(
      CHR = c, SNP = sprintf("rs%d", 2000000 + gid), BP = bp[j],
      A1 = a1, A2 = a2, base = 1L, A = inA, B = inB,
      baseL2 = baseL2, AL2 = AL2, BL2 = BL2,
      wLD = round(runif(1, 1, 40), 4), MAF = round(runif(1, 0.06, 0.5), 4),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
M_snps <- nrow(tab)

# ---- Write per-chromosome files in the formats s_ldsc expects ---------------
for (c in seq_len(chr)) {
  sub <- tab[tab$CHR == c, ]
  # annotation LD scores: CHR SNP BP <annot>L2 ...
  gz <- gzfile(file.path(annotdir, paste0(c, ".l2.ldscore.gz")), "w")
  write.table(sub[, c("CHR","SNP","BP","baseL2","AL2","BL2")], gz,
              sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
  # M_5_50: per-annotation count of (MAF 5-50%) SNPs in this chr
  m_c <- colSums(sub[, c("base","A","B")])
  writeLines(paste(m_c, collapse = "\t"), file.path(annotdir, paste0(c, ".l2.M_5_50")))
  # annot.gz: CHR BP SNP CM <annot membership>
  ag <- gzfile(file.path(annotdir, paste0(c, ".annot.gz")), "w")
  annot_df <- data.frame(CHR = sub$CHR, BP = sub$BP, SNP = sub$SNP, CM = 0,
                         base = sub$base, A = sub$A, B = sub$B)
  write.table(annot_df, ag, sep = "\t", quote = FALSE, row.names = FALSE)
  close(ag)
  # frq (PLINK): CHR SNP A1 A2 MAF NCHROBS
  frq_df <- data.frame(CHR = sub$CHR, SNP = sub$SNP, A1 = sub$A1, A2 = sub$A2,
                       MAF = sub$MAF, NCHROBS = 1000L)
  write.table(frq_df, file.path(annotdir, paste0(c, ".frq")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # weights (wld): CHR SNP BP L2
  wz <- gzfile(file.path(wlddir, paste0(c, ".l2.ldscore.gz")), "w")
  wdf <- data.frame(CHR = sub$CHR, SNP = sub$SNP, BP = sub$BP, L2 = sub$wLD)
  write.table(wdf, wz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(wz)
  writeLines("0", file.path(wlddir, paste0(c, ".l2.M_5_50")))
}

# ---- Simulate GWAS z-scores so the partitioned regression is well-posed -----
# E[z_j^2] = 1 + N_j * (annotLD %*% tau_j); cross via shared latent.
Nj  <- c(50000, 45000)
tau <- matrix(c(8e-7, 1.2e-6, 9e-7,    # trait 1 per-annotation coefficients
                7e-7, 1.0e-6, 1.1e-6), # trait 2
              nrow = n_annot)
annotLD <- as.matrix(tab[, c("baseL2","AL2","BL2")])
rg_e <- 0.4  # cross-trait correlation of the per-SNP signal
Z <- matrix(0, M_snps, n_traits)
for (i in seq_len(M_snps)) {
  v1 <- 1 + Nj[1] * sum(annotLD[i, ] * tau[, 1])
  v2 <- 1 + Nj[2] * sum(annotLD[i, ] * tau[, 2])
  cv <- rg_e * sqrt((v1 - 1) * (v2 - 1))
  Sig <- matrix(c(v1, cv, cv, v2), 2, 2)
  e <- eigen(Sig, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 1e-8))) %*% t(e$vectors)
  Z[i, ] <- A %*% rnorm(2)
}
munged_paths <- character(n_traits)
for (j in seq_len(n_traits)) {
  ss <- data.frame(SNP = tab$SNP, A1 = tab$A1, A2 = tab$A2,
                   Z = round(Z[, j], 6), N = Nj[j])
  p <- file.path(synthdir, paste0(trait_names[j], ".sumstats.gz"))
  gz <- gzfile(p, "w"); write.table(ss, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  munged_paths[j] <- p
}

# ---- Run R s_ldsc -----------------------------------------------------------
cat("=== s_ldsc ===\n")
res <- GenomicSEM::s_ldsc(
  traits          = munged_paths,
  sample.prev     = rep(NA, n_traits),
  population.prev = rep(NA, n_traits),
  ld              = paste0(annotdir, "/"),
  wld             = paste0(wlddir, "/"),
  frq             = paste0(annotdir, "/"),
  trait.names     = trait_names,
  n.blocks        = 20,
  exclude_cont    = FALSE
)
# res: S_Tau / S_List style — capture per-annotation S and V.
str(res, max.level = 1)
saveRDS(res, file.path(synthdir, "r_sldsc.rds"))
cat("\n=== done ===\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
