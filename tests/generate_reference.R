#!/usr/bin/env Rscript
#
# Generate reference fixtures for GenomicSEM-rs validation tests.
#
# Requires: Matrix, lavaan, jsonlite
# Optional: GenomicSEM (for full LDSC pipeline tests)
#
# Usage: Rscript generate_reference.R
#
# Output: tests/fixtures/*.json

library(jsonlite)
library(Matrix)
library(lavaan)

outdir <- "fixtures"
dir.create(outdir, showWarnings = FALSE)

mat_to_list <- function(m) {
  lapply(1:nrow(m), function(i) as.numeric(m[i,]))
}

write_fixture <- function(data, name) {
  path <- file.path(outdir, paste0(name, ".json"))
  writeLines(toJSON(data, auto_unbox = TRUE, digits = 15), path)
  cat("Wrote", path, "\n")
}

# ============================================================
# Test Case 1: nearPD
# ============================================================
cat("=== nearPD ===\n")
input_npd <- matrix(c(
  1.0,  0.9,  0.9,
  0.9,  1.0, -0.1,
  0.9, -0.1,  1.0
), 3, 3, byrow = TRUE)

result_npd <- as.matrix(nearPD(input_npd, corr = FALSE, keepDiag = TRUE)$mat)

write_fixture(list(
  input = mat_to_list(input_npd),
  expected = mat_to_list(result_npd)
), "near_pd")

# ============================================================
# Test Case 2: vech (half-vectorization)
# ============================================================
cat("=== vech ===\n")
m3 <- matrix(c(
  1.0, 0.5, 0.3,
  0.5, 2.0, 0.4,
  0.3, 0.4, 3.0
), 3, 3, byrow = TRUE)

# lavaan's vech: column-major lower triangle
v3 <- lav_matrix_vech(m3)
# reverse
m3_rev <- lav_matrix_vech_reverse(v3)

write_fixture(list(
  input = mat_to_list(m3),
  vech = as.numeric(v3),
  reverse = mat_to_list(m3_rev)
), "vech_3x3")

m4 <- matrix(c(
  1.0, 0.2, 0.3, 0.1,
  0.2, 2.0, 0.4, 0.5,
  0.3, 0.4, 3.0, 0.6,
  0.1, 0.5, 0.6, 4.0
), 4, 4, byrow = TRUE)

v4 <- lav_matrix_vech(m4)
m4_rev <- lav_matrix_vech_reverse(v4)

write_fixture(list(
  input = mat_to_list(m4),
  vech = as.numeric(v4),
  reverse = mat_to_list(m4_rev)
), "vech_4x4")

# ============================================================
# Test Case 3: cov2cor
# ============================================================
cat("=== cov2cor ===\n")
cov_mat <- matrix(c(
  4.0, 1.2, 0.8,
  1.2, 9.0, 2.1,
  0.8, 2.1, 16.0
), 3, 3, byrow = TRUE)

cor_mat <- cov2cor(cov_mat)

write_fixture(list(
  input = mat_to_list(cov_mat),
  expected = mat_to_list(cor_mat)
), "cov_to_cor")

# ============================================================
# Test Case 4: V_SNP construction (.get_V_SNP equivalent)
# ============================================================
cat("=== V_SNP ===\n")

# Inputs: k=2 traits, single SNP
k <- 2
SE_SNP <- matrix(c(0.05, 0.08), nrow = 1) # 1 SNP, k traits
I_LD <- matrix(c(1.05, 0.03, 0.03, 1.08), 2, 2) # intercept matrix
varSNP <- c(0.25) # 2*MAF*(1-MAF) for MAF=0.15 approx
coords <- which(!is.na(I_LD), arr.ind = TRUE)

# Source the utility functions
source("../GenomicSEM/R/utils.R")

V_SNP_standard <- .get_V_SNP(SE_SNP, I_LD, varSNP, "standard", coords, k, 1)
V_SNP_conserv <- .get_V_SNP(SE_SNP, I_LD, varSNP, "conserv", coords, k, 1)
V_SNP_none <- .get_V_SNP(SE_SNP, I_LD, varSNP, "none", coords, k, 1)

v_snp_input <- list(
  se_snp = as.numeric(SE_SNP[1,]),
  i_ld = mat_to_list(I_LD),
  var_snp = varSNP[1],
  k = k
)

write_fixture(c(v_snp_input, list(expected = mat_to_list(V_SNP_standard))), "v_snp_standard")
write_fixture(c(v_snp_input, list(expected = mat_to_list(V_SNP_conserv))), "v_snp_conservative")
write_fixture(c(v_snp_input, list(expected = mat_to_list(V_SNP_none))), "v_snp_none")

# ============================================================
# Test Case 5: S_Full construction (.get_S_Full equivalent)
# ============================================================
cat("=== S_Full ===\n")

S_LD <- matrix(c(0.30, 0.08, 0.08, 0.25), 2, 2)
colnames(S_LD) <- rownames(S_LD) <- c("V1", "V2")
beta_SNP <- matrix(c(0.04, -0.02), nrow = 1)

S_Full <- .get_S_Full(k, S_LD, varSNP, beta_SNP, TWAS = FALSE, i = 1)

write_fixture(list(
  s_ld = mat_to_list(S_LD),
  beta_snp = as.numeric(beta_SNP[1,]),
  var_snp = varSNP[1],
  expected = mat_to_list(unname(S_Full))
), "s_full")

# ============================================================
# Test Case 6: V_Full construction (.get_V_full equivalent)
# ============================================================
cat("=== V_Full ===\n")

kstar_ld <- k * (k + 1) / 2
V_LD <- diag(kstar_ld) * 0.01  # simple diagonal V_LD
V_LD[1, 2] <- 0.002; V_LD[2, 1] <- 0.002
V_LD[1, 3] <- 0.001; V_LD[3, 1] <- 0.001
V_LD[2, 3] <- 0.003; V_LD[3, 2] <- 0.003

varSNPSE2 <- (0.0005)^2 # standard default

V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP_standard)

write_fixture(list(
  v_ld = mat_to_list(V_LD),
  v_snp = mat_to_list(V_SNP_standard),
  var_snp_se2 = varSNPSE2,
  k = k,
  expected = mat_to_list(V_Full)
), "v_full")

# ============================================================
# Test Case 7: Z_pre computation (.get_Z_pre equivalent)
# ============================================================
cat("=== Z_pre ===\n")

z_standard <- .get_Z_pre(1, beta_SNP, SE_SNP, I_LD, "standard")
z_conserv <- .get_Z_pre(1, beta_SNP, SE_SNP, I_LD, "conserv")
z_none <- .get_Z_pre(1, beta_SNP, SE_SNP, I_LD, "none")

write_fixture(list(
  beta = as.numeric(beta_SNP[1,]),
  se = as.numeric(SE_SNP[1,]),
  i_ld = mat_to_list(I_LD),
  standard = as.numeric(z_standard),
  conservative = as.numeric(z_conserv),
  none = as.numeric(z_none)
), "z_pre")

# ============================================================
# Test Case 8: SEM fitting with known S/V
# ============================================================
cat("=== SEM fitting ===\n")

# 3-trait genetic covariance from a 1-factor model
# True: F1 =~ 0.7*V1 + 0.6*V2 + 0.5*V3, F1~~1*F1
# S = Lambda * Psi * Lambda' + Theta
S_sem <- matrix(c(
  0.60, 0.42, 0.35,
  0.42, 0.50, 0.30,
  0.35, 0.30, 0.40
), 3, 3, byrow = TRUE)
colnames(S_sem) <- rownames(S_sem) <- c("V1", "V2", "V3")

# Sampling covariance (diagonal for simplicity)
kstar_sem <- 3 * 4 / 2  # = 6
V_sem <- diag(kstar_sem) * 0.001

# Fit with lavaan DWLS
model_str <- "F1 =~ NA*V1 + V2 + V3
F1 ~~ 1*F1
V1 ~~ V1
V2 ~~ V2
V3 ~~ V3"

# Weight matrix: inverse of diagonal of V (DWLS)
V_diag_mat <- diag(as.numeric(diag(V_sem)), nrow = kstar_sem)
W_sem <- solve(V_diag_mat)

fit_sem <- tryCatch({
  sem(model_str,
      sample.cov = S_sem,
      estimator = "DWLS",
      WLS.V = W_sem,
      sample.nobs = 200,
      se = "standard",
      optim.dx.tol = 0.01)
}, error = function(e) {
  cat("SEM fitting with DWLS failed:", e$message, "\n")
  cat("Trying ML instead...\n")
  tryCatch({
    sem(model_str,
        sample.cov = S_sem,
        estimator = "ML",
        sample.nobs = 200)
  }, error = function(e2) {
    cat("ML also failed:", e2$message, "\n")
    NULL
  })
})

if (!is.null(fit_sem)) {
  pe <- parameterEstimates(fit_sem)
  free_params <- pe[pe$op %in% c("=~", "~~") & pe$se > 0, ]

  fit_measures <- fitMeasures(fit_sem, c("chisq", "df", "pvalue", "cfi", "aic", "srmr"))

  # Get sandwich SEs (matching GenomicSEM's commonfactor.R approach)
  Delta_raw <- lavInspect(fit_sem, "delta")
  # In newer lavaan, delta may be a matrix directly rather than a list
  if (is.list(Delta_raw)) { Delta <- Delta_raw[[1]] } else { Delta <- as.matrix(Delta_raw) }
  W_lav <- lavInspect(fit_sem, "WLS.V")
  bread <- solve(t(Delta) %*% W_lav %*% Delta)
  lettuce <- W_lav %*% Delta

  # Reorder V
  order_idx <- .rearrange(k = 3, fit = fit_sem, names = c("V1", "V2", "V3"))
  V_reorder <- V_sem[order_idx, order_idx]

  Ohtt <- bread %*% t(lettuce) %*% V_reorder %*% lettuce %*% bread
  sandwich_se <- sqrt(diag(Ohtt))

  write_fixture(list(
    s = mat_to_list(S_sem),
    v = mat_to_list(V_sem),
    v_diag = as.numeric(diag(V_sem)),
    model = model_str,
    estimates = lapply(1:nrow(free_params), function(i) {
      list(
        lhs = free_params$lhs[i],
        op = free_params$op[i],
        rhs = free_params$rhs[i],
        est = free_params$est[i],
        se = free_params$se[i]
      )
    }),
    sandwich_se = as.numeric(sandwich_se),
    fit_indices = list(
      chisq = as.numeric(fit_measures["chisq"]),
      df = as.numeric(fit_measures["df"]),
      p_chisq = as.numeric(fit_measures["pvalue"]),
      cfi = as.numeric(fit_measures["cfi"]),
      aic = as.numeric(fit_measures["aic"]),
      srmr = as.numeric(fit_measures["srmr"])
    ),
    reorder_perm = as.numeric(order_idx)
  ), "sem_1factor")

  # Also save reordered V for reorder test
  write_fixture(list(
    v = mat_to_list(V_sem),
    user_order = c("V1", "V2", "V3"),
    model_order = rownames(lavInspect(fit_sem, "theta")),
    perm = as.numeric(order_idx),
    v_reordered = mat_to_list(V_reorder)
  ), "reorder")
} else {
  cat("Skipping SEM fixtures (lavaan fitting failed)\n")
}

# ============================================================
# Test Case 10: 2-factor SEM with multiple fixed rows
# Tests that parameter indexing is correct when fixed rows
# appear in various positions in the partable.
# ============================================================
cat("=== 2-factor SEM ===\n")

S_2f <- matrix(c(
  0.60, 0.42, 0.10, 0.05,
  0.42, 0.50, 0.08, 0.04,
  0.10, 0.08, 0.45, 0.30,
  0.05, 0.04, 0.30, 0.55
), 4, 4, byrow = TRUE)
colnames(S_2f) <- rownames(S_2f) <- c("V1", "V2", "V3", "V4")

kstar_2f <- 4 * 5 / 2  # = 10
V_2f <- diag(kstar_2f) * 0.001

# 2 factors, cross-loading fixed to 0, factor variances fixed to 1
# partable will have: loadings (free), fixed factor variances, fixed zero
# cross-loadings, free residual variances, free factor covariance
model_2f <- "
F1 =~ NA*V1 + V2
F2 =~ NA*V3 + V4
F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~~ F2
V1 ~~ V1
V2 ~~ V2
V3 ~~ V3
V4 ~~ V4
"

V_diag_2f <- diag(as.numeric(diag(V_2f)), nrow = kstar_2f)
W_2f <- solve(V_diag_2f)

fit_2f <- tryCatch({
  sem(model_2f,
      sample.cov = S_2f,
      estimator = "DWLS",
      WLS.V = W_2f,
      sample.nobs = 200,
      se = "standard",
      optim.dx.tol = 0.01)
}, error = function(e) {
  cat("2-factor SEM failed:", e$message, "\n")
  NULL
})

if (!is.null(fit_2f)) {
  pe_2f <- parameterEstimates(fit_2f)
  free_2f <- pe_2f[pe_2f$op %in% c("=~", "~~") & pe_2f$se > 0, ]

  write_fixture(list(
    s = mat_to_list(S_2f),
    v = mat_to_list(V_2f),
    v_diag = as.numeric(diag(V_2f)),
    model = model_2f,
    estimates = lapply(1:nrow(free_2f), function(i) {
      list(
        lhs = free_2f$lhs[i],
        op = free_2f$op[i],
        rhs = free_2f$rhs[i],
        est = free_2f$est[i],
        se = free_2f$se[i]
      )
    })
  ), "sem_2factor")
} else {
  cat("Skipping 2-factor fixture\n")
}

cat("\n=== All fixtures generated ===\n")
cat("Files in", outdir, ":\n")
cat(paste(" ", list.files(outdir)), sep = "\n")
