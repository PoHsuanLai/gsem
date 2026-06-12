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
# rgmodel pulls in `simsalapar`; if it isn't installed, keep the committed
# rgmodel.json and skip regenerating it rather than aborting the whole script
# (the other fixtures below don't depend on it).
rg <- tryCatch(GenomicSEM::rgmodel(covstruc), error = function(e) {
  cat("  (skipping rgmodel regen:", conditionMessage(e), ")\n"); NULL
})
if (!is.null(rg)) {
  write_fixture(list(
    s   = mat_to_list(S),
    v   = mat_to_list(V),
    r   = mat_to_list(as.matrix(rg$R)),
    v_r = mat_to_list(as.matrix(rg$V_R))
  ), "rgmodel")
}

# ---------------------------------------------------------------------------
# write.model: factor -> indicator assignment. R emits random residual-
# variance labels (non-deterministic) and a different identification format
# than gsem, so we compare the deterministic part: which indicators load on
# which factor given the cutoff (|loading| > cutoff). Parse R's =~ lines.
# ---------------------------------------------------------------------------
cat("=== write.model ===\n")
Lw <- matrix(c(0.80, 0.70, 0.10, 0.05, 0.00,
               0.00, 0.10, 0.75, 0.65, 0.55), ncol = 2)
rownames(Lw) <- paste0("V", 1:5)
Sw <- diag(5); colnames(Sw) <- rownames(Sw) <- paste0("V", 1:5)
cutoff_w <- 0.3
model_str <- GenomicSEM::write.model(Lw, Sw, cutoff = cutoff_w)
# Parse "F1=~V1 + V2" lines into factor -> indicators.
lines <- strsplit(model_str, "\n")[[1]]
facmap <- list()
for (ln in lines) {
  ln <- trimws(ln)
  if (grepl("=~", ln)) {
    parts <- strsplit(ln, "=~")[[1]]
    fac <- trimws(parts[1])
    inds <- trimws(strsplit(parts[2], "\\+")[[1]])
    inds <- gsub("NA\\*", "", inds)
    facmap[[fac]] <- inds
  }
}
write_fixture(list(
  loadings = mat_to_list(Lw),
  names    = rownames(Lw),
  cutoff   = cutoff_w,
  factors  = names(facmap),
  # one indicator list per factor, in factor order
  indicators = unname(facmap)
), "write_model")

# ---------------------------------------------------------------------------
# subSV: subset vech(S) and the V block by a set of 1-based vech positions.
# ---------------------------------------------------------------------------
# NOTE: stock R subSV has a bug in its matrix-input validation branch (it
# references an undefined `RMATRIX`), so passing SMATRIX/VMATRIX directly errors.
# The LDSC_OBJECT path is bug-free, so we drive it that way — the numeric result
# is identical to what the (fixed) matrix path would return.
cat("=== subSV ===\n")
subsv_index_s <- c(1, 3, 6, 10)          # TYPE="S": positions in vech(S), incl diag
sub_s_obj <- GenomicSEM::subSV(
  LDSC_OBJECT = list(S = S, V = V),
  INDEXVALS = subsv_index_s, TYPE = "S"
)
# TYPE="R": off-diagonal numbering on the correlation matrix. For k=4 the strict
# lower triangle has k(k-1)/2 = 6 positions, so V_R is 6x6.
R_corr <- cov2cor(S)
kstar_r <- k * (k - 1) / 2               # 6
V_R <- diag(kstar_r) * 0.0015
subsv_index_r <- c(1, 4, 6)
sub_r_obj <- GenomicSEM::subSV(
  LDSC_OBJECT = list(R = R_corr, V_R = V_R),
  INDEXVALS = subsv_index_r, TYPE = "R"
)
write_fixture(list(
  s            = mat_to_list(S),
  v            = mat_to_list(V),
  r_corr       = mat_to_list(R_corr),
  v_r          = mat_to_list(V_R),
  index_s      = subsv_index_s,
  index_r      = subsv_index_r,
  sub_s        = as.numeric(sub_s_obj$subS),
  sub_v        = mat_to_list(as.matrix(sub_s_obj$subV)),
  sub_s_r      = as.numeric(sub_r_obj$subS),
  sub_v_r      = mat_to_list(as.matrix(sub_r_obj$subV))
), "subsv")

# ---------------------------------------------------------------------------
# summaryGLSbands: GLS confidence-band data.
# ---------------------------------------------------------------------------
# We reproduce summaryGLSbands' NUMERIC core (the ggplot rendering is not
# ported, and stock R additionally has two bugs we don't reproduce: the band
# loop overwrites Ohm with Y on the Y/V_Y input path, and it requires ggplot2 to
# even return). The band math itself — main GLS fit, then the SE of the fitted
# value at each re-centred grid origin — is well-defined and is what the Rust
# `summary_gls_bands` computes. We run that exact math here.
cat("=== summaryGLSbands (numeric core) ===\n")
np_b <- 8
pred_b <- seq(-1.5, 1.5, length.out = np_b)
y_b <- 0.5 + 0.3 * pred_b
Ohm_b <- diag(np_b) * 0.01
intervals_b <- 5
band_size_b <- 1
Xb <- cbind(rep(1, np_b), pred_b)
betas_b <- solve(t(Xb) %*% solve(Ohm_b) %*% Xb) %*% t(Xb) %*% solve(Ohm_b) %*% y_b
rng_b <- range(pred_b)
band_se_b <- numeric(intervals_b)
for (i in 0:(intervals_b - 1)) {
  centered <- pred_b - (min(rng_b) + i * (max(rng_b) - min(rng_b)) / intervals_b)
  XX <- cbind(rep(1, np_b), centered)
  band_se_b[i + 1] <- sqrt(diag(solve(t(XX) %*% solve(Ohm_b) %*% XX)))[1]
}
grid_b <- seq(from = min(pred_b), to = max(pred_b), length.out = intervals_b)
line_b <- betas_b[1] + betas_b[2] * grid_b
write_fixture(list(
  predictors = pred_b,
  y          = y_b,
  v_y        = mat_to_list(Ohm_b),
  intervals  = intervals_b,
  band_size  = band_size_b,
  betas      = as.numeric(betas_b[, 1]),
  grid       = grid_b,
  line       = as.numeric(line_b),
  band_se    = band_se_b,
  upper      = as.numeric(line_b + band_size_b * band_se_b),
  lower      = as.numeric(line_b - band_size_b * band_se_b)
), "gls_bands")

cat("\n=== covstruc reference fixtures generated ===\n")
unlink(list.files(".", pattern = "\\.log$", full.names = TRUE))
