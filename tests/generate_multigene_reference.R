#!/usr/bin/env Rscript
#
# Reference fixture for multiGene (R GenomicSEM's joint multi-gene S/V builder).
#
# IMPORTANT: stock R GenomicSEM's multiGene() is BROKEN for k >= 2 traits: at
# the cross-trait within-gene step it assigns into `V_SNP[y,x]`, a variable that
# does not exist in multiGene (the surrounding matrix is named `V_Gene`), so the
# function aborts with "object 'V_SNP' not found". We therefore source a
# MINIMALLY-PATCHED copy (the single `V_SNP` -> `V_Gene` typo fixed) to obtain
# R's *intended* output. gsem implements this corrected algorithm.
#
# multiGene is otherwise algorithmically identical to multiSNP, with gene
# heritabilities (Genes$HSQ) as the "variances" and a much smaller fixed SE
# floor (1e-8). As with multiSNP, R's cross-gene cross-trait block has a
# constant-LD indexing bug, so the fixture uses equal off-diagonal LD.
#
# Requires: GenomicSEM (for nearPD/deps) + dplyr, jsonlite.
# Usage: cd tests && Rscript generate_multigene_reference.R

suppressMessages({
  library(jsonlite); library(GenomicSEM); library(dplyr); library(Matrix)
})
source("/tmp/multiGene_patched.R")   # patched multiGene (V_SNP -> V_Gene)
outdir <- "fixtures"
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

k <- 2
traits <- c("T1", "T2")
S_LD <- matrix(c(0.40, 0.18, 0.18, 0.50), 2, 2)
colnames(S_LD) <- rownames(S_LD) <- traits
V_LD <- diag(c(4e-4, 2.5e-4, 5e-4))
V_LD[1, 2] <- V_LD[2, 1] <- 0.4e-4
V_LD[1, 3] <- V_LD[3, 1] <- 0.3e-4
V_LD[2, 3] <- V_LD[3, 2] <- 0.35e-4
I_LD <- matrix(c(1.02, 0.05, 0.05, 1.03), 2, 2)
covstruc <- list(V_LD = V_LD, S_LD = S_LD, I_LD = I_LD)

# Genes: f = 3, heritabilities in HSQ.
Genes <- data.frame(
  Gene = c("G1", "G2", "G3"),
  ID   = c("G1", "G2", "G3"),
  CHR  = c(1L, 1L, 1L),
  HSQ  = c(0.020, 0.030, 0.015),
  beta.T1 = c(0.30, -0.20, 0.15), se.T1 = c(0.06, 0.05, 0.08),
  beta.T2 = c(0.18,  0.25, -0.10), se.T2 = c(0.07, 0.06, 0.09),
  stringsAsFactors = FALSE
)
# Equal off-diagonal LD neutralises R's cross-gene cross-trait constant-LD bug.
LD <- matrix(0.10, 3, 3); diag(LD) <- 1
rownames(LD) <- colnames(LD) <- Genes$Gene

res <- multiGene(covstruc, Genes = Genes, LD = LD)
S_Full <- as.matrix(res$S_Full)
V_Full <- as.matrix(res$V_Full)
cat("S_Full dim:", dim(S_Full), "  V_Full dim:", dim(V_Full), "\n")
print(round(S_Full, 5))

out <- list(
  traits = traits, gene_names = Genes$Gene,
  hsq = Genes$HSQ,
  beta = mat_to_list(as.matrix(Genes[, c("beta.T1", "beta.T2")])),
  se   = mat_to_list(as.matrix(Genes[, c("se.T1",   "se.T2")])),
  s_ld = mat_to_list(S_LD), v_ld = mat_to_list(V_LD), i_ld = mat_to_list(I_LD),
  ld   = mat_to_list(LD),
  s_full = mat_to_list(S_Full),
  v_full = mat_to_list(V_Full)
)
writeLines(toJSON(out, auto_unbox = TRUE, digits = 15),
           file.path(outdir, "multigene_synth.json"))
cat("\nWrote fixtures/multigene_synth.json\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
