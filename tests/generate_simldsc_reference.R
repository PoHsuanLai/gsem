#!/usr/bin/env Rscript
#
# Reference fixture for simLDSC's deterministic per-SNP Z covariance (Sigma).
#
# simLDSC writes random GWAS sumstats; the random draw (MASS::mvrnorm on R's RNG
# stream) cannot be reproduced bit-for-bit cross-platform. The scientifically
# meaningful, deterministic core is the per-SNP variance-covariance matrix of Z
# statistics, Sigma(ld), which simLDSC builds before drawing. This script runs
# R GenomicSEM simLDSC's OWN construction block (lines ~109-170 of simLDSC.R,
# copied verbatim) on synthetic inputs and dumps the resulting Sigma for several
# LD-score values, so gsem must reproduce R's exact varZ/covZ algebra:
#   varZ[i]  = (N_i * S_ii / M) * ld + int_i
#   covZ[ij] = (sqrt(N_i N_j) * S_ij / M) * ld + rPheno_ij * N_ij / sqrt(N_i N_j)
#
# Requires: jsonlite. Usage: cd tests && Rscript generate_simldsc_reference.R

suppressMessages({ library(jsonlite) })
outdir <- "fixtures"
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

# ── synthetic model: k = 3 traits ───────────────────────────────────────────
covMatrix <- matrix(c(0.100, 0.050, 0.025,
                      0.050, 0.080, 0.030,
                      0.025, 0.030, 0.090), 3, 3)
phenos <- nrow(covMatrix)
rownames(covMatrix) <- colnames(covMatrix) <- paste("Pheno_", 1:phenos, sep = "")
covMatrix <- as.data.frame(covMatrix)

# N matrix: diagonal sample sizes, off-diagonal sample overlaps.
N <- diag(c(75000, 300000, 120000), phenos, phenos)
N_overlap <- 0.50
N[lower.tri(N)] <- N[upper.tri(N)] <- 0
for (a in 1:phenos) for (b in 1:phenos) if (a != b)
  N[a, b] <- sqrt(diag(N)[a] * diag(N)[b]) * N_overlap
colnames(N) <- rownames(N) <- paste("Pheno_", 1:phenos, sep = "")

int <- c(1.02, 1.01, 1.03)               # LDSC intercepts (per trait)
rG <- cov2cor(as.matrix(covMatrix))
rPheno <- matrix(0.40, phenos, phenos); diag(rPheno) <- 1
colnames(rPheno) <- rownames(rPheno) <- colnames(rG)
M <- 1.5e6

# ── simLDSC construction block (verbatim algebra from simLDSC.R) ─────────────
sigma_for_ld <- function(ld) {
  varZ <- list()
  for (i in 1:phenos) {
    varZ[[i]] <- (diag(N)[i] * covMatrix[i, i] / M) * ld + int[i]
    names(varZ)[i] <- paste(colnames(covMatrix)[i], ",", colnames(covMatrix)[i], sep = "")
  }
  covZ <- list()
  covs <- t(unique(combn(colnames(rG), 2)))
  for (i in 1:nrow(covs)) {
    P1 <- covs[i, 1]; P2 <- covs[i, 2]
    covZ[[i]] <- (sqrt(diag(N)[P1] * diag(N)[P2]) * covMatrix[P1, P2] / M) * ld +
      rPheno[P1, P2] * N[P1, P2] / sqrt(diag(N)[P1] * diag(N)[P2])
    names(covZ)[i] <- paste0(P1, ",", P2)
  }
  SigmaNames <- matrix(NA, phenos, phenos)
  for (a in 1:phenos) for (c in 1:phenos)
    SigmaNames[a, c] <- paste(rownames(covMatrix)[a], ",", colnames(covMatrix)[c], sep = "")
  makeSymm <- function(m) { m[upper.tri(m)] <- t(m)[upper.tri(m)]; m }
  SigmaNames <- makeSymm(t(SigmaNames))
  varcovarZ <- c(varZ, covZ)
  Sigma <- matrix(NA, phenos, phenos)
  for (a in 1:phenos) for (c in 1:phenos)
    Sigma[a, c] <- varcovarZ[[SigmaNames[a, c]]]
  Sigma
}

ld_values <- c(1.0, 5.0, 12.5, 30.0, 80.0)
sigmas <- lapply(ld_values, function(l) mat_to_list(sigma_for_ld(l)))
cat("Sigma at ld=5:\n"); print(round(sigma_for_ld(5.0), 6))

out <- list(
  s = mat_to_list(as.matrix(covMatrix)),
  n_diag = as.numeric(diag(N)),
  n_overlap = N_overlap,
  intercepts = int,
  r_pheno = mat_to_list(rPheno),
  m = M,
  ld_values = ld_values,
  sigma = sigmas
)
writeLines(toJSON(out, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "simldsc_synth.json"))
cat("\nWrote fixtures/simldsc_synth.json\n")
