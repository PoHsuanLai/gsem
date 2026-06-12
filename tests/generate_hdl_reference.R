#!/usr/bin/env Rscript
#
# Reference fixture for hdl (R GenomicSEM's multivariable High-Definition
# Likelihood). The real UKB LD reference panels are large and externally hosted
# (blocked here), so this synthesises a small, self-contained panel in the exact
# .rda / .bim format R's hdl() expects, runs R hdl() on simulated sumstats, and
# dumps both the panel (per-piece eigenvalues `lam`, eigenvectors `V`, LD scores
# `LDsc`) and the resulting S / I / V so gsem can reproduce them.
#
# hdl() hardcodes a `for(chr in 1:22)` loop, so we lay out 22 single-piece
# "chromosomes" of a few SNPs each. The HDL likelihood is evaluated in the
# eigenspace (bstar = V' bhat, lam = eigenvalues), which gsem now mirrors.
#
# Requires: GenomicSEM, jsonlite, MASS. Usage: cd tests && Rscript generate_hdl_reference.R

suppressMessages({ library(jsonlite); library(GenomicSEM); library(MASS) })
set.seed(7)
outdir <- "fixtures"
paneldir <- file.path(tempdir(), "hdl_panel")
dir.create(paneldir, showWarnings = FALSE, recursive = TRUE)
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

n_chr <- 22
snps_per_piece <- 4
alleles <- c("A", "C", "G", "T")

# True per-SNP standardised effects for 2 correlated traits.
S_true <- matrix(c(0.05, 0.02, 0.02, 0.06), 2, 2)   # genetic covariance
L_true <- t(chol(S_true))
N1 <- 60000; N2 <- 80000

nsnps.list.imputed <- vector("list", n_chr)
all_snps <- c(); all_a1 <- c(); all_a2 <- c()
panel_pieces <- list()
# accumulate true per-SNP betas to build Z later
beta_true <- matrix(0, 0, 2)

snp_counter <- 0
for (chr in 1:n_chr) {
  m <- snps_per_piece
  ids <- paste0("rs", chr, "_", 1:m)
  a1 <- sample(alleles, m, replace = TRUE)
  a2 <- sapply(a1, function(x) sample(setdiff(alleles, x), 1))

  # Build a PSD block LD correlation matrix (unit diagonal).
  B <- matrix(rnorm(m * m, 0, 0.3), m, m)
  R <- cov2cor(B %*% t(B) + diag(m))
  e <- eigen(R, symmetric = TRUE)
  lam <- e$values
  V <- e$vectors
  LDsc <- rowSums(R^2)

  save(LDsc, lam, V, file = file.path(paneldir, paste0("chr", chr, ".1.rda")))
  bim <- data.frame(chr = chr, id = ids, non = 0,
                    pos = (1:m) * 1000, A1 = a1, A2 = a2)
  write.table(bim, file = file.path(paneldir, paste0("chr", chr, ".1.bim")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  nsnps.list.imputed[[chr]] <- m
  all_snps <- c(all_snps, ids); all_a1 <- c(all_a1, a1); all_a2 <- c(all_a2, a2)
  # true per-SNP betas ~ N(0, S_true) so that LDSC structure ~ S_true
  bt <- t(L_true %*% matrix(rnorm(2 * m), 2, m))
  beta_true <- rbind(beta_true, bt)
  panel_pieces[[chr]] <- list(
    snps = ids, a1 = a1, a2 = a2,
    ldsc = as.numeric(LDsc), lam = as.numeric(lam), v = mat_to_list(V)
  )
  snp_counter <- snp_counter + m
}

snps.list.imputed.vector <- all_snps
nsnps.list <- nsnps.list.imputed
save(nsnps.list.imputed, nsnps.list,
     file = file.path(paneldir, "UKB_snp_counter_demo.rda"))
save(snps.list.imputed.vector,
     file = file.path(paneldir, "UKB_snp_list_demo.rda"))

# Simulate per-trait Z from true betas: Z = beta*sqrt(N) + noise.
M_tot <- nrow(beta_true)
write_sumstats <- function(beta_col, N, fname) {
  z <- beta_col * sqrt(N) + rnorm(M_tot, 0, 1)
  df <- data.frame(SNP = all_snps, A1 = all_a1, A2 = all_a2, N = N, Z = z)
  f <- file.path(paneldir, fname)
  write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  list(snp = all_snps, a1 = all_a1, a2 = all_a2, n = rep(N, M_tot), z = z)
}
t1 <- write_sumstats(beta_true[, 1], N1, "trait1.sumstats")
t2 <- write_sumstats(beta_true[, 2], N2, "trait2.sumstats")

res <- hdl(
  traits = c(file.path(paneldir, "trait1.sumstats"),
             file.path(paneldir, "trait2.sumstats")),
  trait.names = c("T1", "T2"),
  LD.path = paneldir, Nref = 1000, method = "piecewise"
)

cat("\nHDL S:\n"); print(res$S)
cat("HDL I:\n"); print(res$I)

out <- list(
  n_ref = 1000,
  panel = panel_pieces,
  trait1 = list(snp = t1$snp, a1 = t1$a1, a2 = t1$a2, n = t1$n, z = t1$z),
  trait2 = list(snp = t2$snp, a1 = t2$a1, a2 = t2$a2, n = t2$n, z = t2$z),
  s = mat_to_list(res$S),
  i = mat_to_list(res$I),
  v = mat_to_list(res$V)
)
writeLines(toJSON(out, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "hdl_synth.json"))
cat("\nWrote fixtures/hdl_synth.json\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
