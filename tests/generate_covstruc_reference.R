#!/usr/bin/env Rscript
#
# Reference fixtures for the covstruc-derived functions: summaryGLS and
# paLDSC. These take an LDSC covstruc (S genetic covariance + V sampling
# covariance) — or, for GLS, a design + outcome + weight — and produce
# deterministic output we can check against gsem.
#
# Requires: GenomicSEM, jsonlite
# Usage: cd tests && Rscript generate_covstruc_reference.R

suppressMessages({
  library(jsonlite)
  library(GenomicSEM)
})

set.seed(7)
outdir <- "fixtures"
mat_to_list <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))
write_fixture <- function(data, name) {
  path <- file.path(outdir, paste0(name, ".json"))
  writeLines(toJSON(data, auto_unbox = TRUE, digits = 15), path)
  cat("Wrote", path, "\n")
}

# ---------------------------------------------------------------------------
# A clean, well-conditioned 4-trait genetic covariance S (PD by construction)
# and a kstar x kstar sampling covariance V (diagonal, PD).
# ---------------------------------------------------------------------------
k <- 4
lambda <- c(0.8, 0.6, 0.7, 0.5)          # single-factor loadings
S <- outer(lambda, lambda)               # common part
diag(S) <- diag(S) + c(0.30, 0.40, 0.25, 0.45)  # add uniqueness -> PD
colnames(S) <- rownames(S) <- paste0("V", 1:k)
kstar <- k * (k + 1) / 2                  # 10
V <- diag(kstar) * 0.002

# ---------------------------------------------------------------------------
# summaryGLS: beta = (X'Ω⁻¹X)⁻¹ X'Ω⁻¹ y, deterministic GLS.
# ---------------------------------------------------------------------------
cat("=== summaryGLS ===\n")
ny <- 6
set.seed(11)
y_gls <- c(0.42, 0.31, 0.55, 0.19, 0.48, 0.27)
# A PD weight matrix Ω (6x6): diagonal + small off-diagonal.
Omega <- diag(ny) * 0.01
Omega[1, 2] <- Omega[2, 1] <- 0.002
Omega[3, 4] <- Omega[4, 3] <- 0.0015
Omega[5, 6] <- Omega[6, 5] <- 0.0025
predictors <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), ncol = 1)  # one predictor

gls_out <- GenomicSEM::summaryGLS(
  Y = y_gls, V_Y = Omega, PREDICTORS = predictors, INTERCEPT = TRUE
)
# gls_out columns: betas, pvals, SE, Z
write_fixture(list(
  y       = as.numeric(y_gls),
  v       = mat_to_list(Omega),
  # Design matrix gsem expects (intercept + predictor), matching R INTERCEPT=T.
  x       = mat_to_list(cbind(1, predictors)),
  betas   = as.numeric(gls_out[, "betas"]),
  se      = as.numeric(gls_out[, "SE"]),
  z       = as.numeric(gls_out[, "Z"]),
  pvals   = as.numeric(gls_out[, "pvals"])
), "summary_gls")

# ---------------------------------------------------------------------------
# paLDSC: observed eigenvalues of cor(S) are deterministic; the simulated
# threshold and suggested nfactors use the LDSC sampling distribution (RNG).
# We commit the observed eigenvalues (exact) plus nfactors (large r -> stable).
# ---------------------------------------------------------------------------
cat("=== paLDSC (observed eigenvalues) ===\n")
# paLDSC's RNG-independent output is the observed eigenvalue spectrum of the
# correlation matrix of S (Horn's parallel analysis observed component). The
# simulated thresholds / suggested nfactors depend on RNG and on plotting
# packages, so we pin only the deterministic observed eigenvalues here.
obs_eig <- sort(eigen(cov2cor(S), symmetric = TRUE, only.values = TRUE)$values,
                decreasing = TRUE)
write_fixture(list(
  s              = mat_to_list(S),
  v              = mat_to_list(V),
  observed_eig   = as.numeric(obs_eig)
), "paldsc")

# ---------------------------------------------------------------------------
# rgmodel: saturated genetic-correlation model. Returns R (the genetic
# correlation matrix) and V_R (its sampling covariance) — matching gsem's
# run_rgmodel RgModelResult.{r, v_r}. R GenomicSEM's rgmodel reads a full
# covstruc list(V, S, I, N, m).
# ---------------------------------------------------------------------------
cat("=== rgmodel ===\n")
I_mat <- diag(k)
N_vec <- rep(50000, k)
m_snps <- 1000
covstruc <- list(V = V, S = S, I = I_mat, N = N_vec, m = m_snps)
rg <- GenomicSEM::rgmodel(covstruc)
write_fixture(list(
  s   = mat_to_list(S),
  v   = mat_to_list(V),
  r   = mat_to_list(as.matrix(rg$R)),
  v_r = mat_to_list(as.matrix(rg$V_R))
), "rgmodel")

cat("\n=== covstruc reference fixtures generated ===\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
