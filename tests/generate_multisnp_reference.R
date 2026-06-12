#!/usr/bin/env Rscript
#
# Reference fixture for multiSNP (R GenomicSEM's joint multi-SNP S/V builder).
#
# multiSNP(covstruc, SNPs, LD, ...) expands an LDSC covariance structure
# (k traits) with f SNPs into a (k+f) "S_Full" observed-covariance matrix and
# its sampling covariance "V_Full", accounting for:
#   - SNP variances (2*MAF*(1-MAF)) on the SNP diagonal,
#   - SNP-SNP covariances from the LD correlation matrix,
#   - SNP-trait covariances (varSNP * beta),
#   - the trait-trait LDSC block (S_LD),
# and a fully-populated V_Full with cross-trait/cross-SNP sampling covariances.
#
# Everything here is deterministic given the inputs (no RNG), so gsem must
# reproduce S_Full and V_Full to numerical precision.
#
# Requires: GenomicSEM, jsonlite. Usage: cd tests && Rscript generate_multisnp_reference.R

suppressMessages({ library(jsonlite); library(GenomicSEM) })
set.seed(2024)
outdir <- "fixtures"
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

# ── covstruc: k = 2 traits, well-conditioned S_LD / V_LD / I_LD ──────────────
k <- 2
traits <- c("T1", "T2")
S_LD <- matrix(c(0.40, 0.18,
                 0.18, 0.50), 2, 2)
colnames(S_LD) <- rownames(S_LD) <- traits
kstar_ld <- k * (k + 1) / 2                     # 3
V_LD <- diag(c(4e-4, 2.5e-4, 5e-4))             # sampling var of vech(S_LD)
# small off-diagonal sampling covariances to exercise the full V machinery
V_LD[1, 2] <- V_LD[2, 1] <- 0.4e-4
V_LD[1, 3] <- V_LD[3, 1] <- 0.3e-4
V_LD[2, 3] <- V_LD[3, 2] <- 0.35e-4
I_LD <- matrix(c(1.02, 0.05,
                 0.05, 1.03), 2, 2)             # LDSC intercept (off-diag = bivar)

covstruc <- list(V_LD = V_LD, S_LD = S_LD, I_LD = I_LD)

# ── SNPs data frame: f = 3 SNPs ──────────────────────────────────────────────
f <- 3
SNPs <- data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  CHR = c(1L, 1L, 1L),
  BP  = c(1000L, 2000L, 3000L),
  MAF = c(0.25, 0.40, 0.15),
  A1  = c("A", "A", "A"),
  A2  = c("G", "C", "T"),
  beta.T1 = c(0.030, -0.020, 0.015),
  se.T1   = c(0.006,  0.005, 0.008),
  beta.T2 = c(0.018,  0.025, -0.010),
  se.T2   = c(0.007,  0.006, 0.009),
  stringsAsFactors = FALSE
)

# ── LD correlation matrix, rownames "<rs>_<A1>" as multiSNP expects ──────────
# NOTE: equal off-diagonal LD (0.1) is used deliberately. R GenomicSEM's
# multiSNP has a bug in the cross-SNP cross-trait sampling-covariance block: it
# weights every such cell by a CONSTANT LD value (LD2[(f^2-f)/2], the last
# lower-triangle entry) instead of the actual SNP-pair LD. With equal
# off-diagonal LD, that buggy constant coincides with the correct per-pair LD,
# so R's full V_Full equals the semantically-correct matrix that gsem builds.
# (gsem uses the correct per-pair LD; the divergence under unequal LD is
# documented in CHANGELOG and unit-tested separately.)
LD <- matrix(c(1.00, 0.10, 0.10,
               0.10, 1.00, 0.10,
               0.10, 0.10, 1.00), 3, 3)
rownames(LD) <- c("rs1_A", "rs2_A", "rs3_A")
colnames(LD) <- rownames(LD)

res <- multiSNP(covstruc, SNPs = SNPs, LD = LD, SNPSE = FALSE)

S_Full <- as.matrix(res$S_Full)
V_Full <- as.matrix(res$V_Full)
cat("S_Full dim:", dim(S_Full), "  V_Full dim:", dim(V_Full), "\n")
print(round(S_Full, 5))

out <- list(
  traits = traits, snp_names = SNPs$SNP,
  maf = SNPs$MAF, a1 = SNPs$A1,
  beta = mat_to_list(as.matrix(SNPs[, c("beta.T1", "beta.T2")])),
  se   = mat_to_list(as.matrix(SNPs[, c("se.T1",   "se.T2")])),
  s_ld = mat_to_list(S_LD), v_ld = mat_to_list(V_LD), i_ld = mat_to_list(I_LD),
  ld   = mat_to_list(LD),
  s_full = mat_to_list(S_Full),
  v_full = mat_to_list(V_Full)
)
writeLines(toJSON(out, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "multisnp_synth.json"))
cat("\nWrote fixtures/multisnp_synth.json\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
