# Parity: gsemr::summaryGLS() must reproduce R GenomicSEM's GLS coefficients,
# SEs, Z and p (committed as summary_gls.json). Pure matrix in/out — exercises
# the binding's numeric-matrix marshalling with no file I/O.

test_that("summaryGLS() reproduces R-reference GLS estimates", {
  skip_if_no_fixtures()
  fx <- load_fixture("summary_gls")

  x <- as.matrix(fx$x) # predictors incl. the intercept column
  y <- as.numeric(fx$y)
  v <- as.matrix(fx$v)

  # x already carries the intercept column, so INTERCEPT = FALSE to avoid
  # double-adding it (matches how the fixture was generated).
  res <- gsemr::summaryGLS(Y = y, V_Y = v, PREDICTORS = x, INTERCEPT = FALSE)

  # The wrapper returns a numeric matrix with columns betas, pvals, SE, Z.
  betas <- res[, "betas"]
  se <- res[, "SE"]
  z <- res[, "Z"]
  pvals <- res[, "pvals"]

  expect_equal(unname(betas), unname(fx$betas), tolerance = 1e-8)
  expect_equal(unname(se), unname(fx$se), tolerance = 1e-8)
  expect_equal(unname(z), unname(fx$z), tolerance = 1e-8)
  expect_equal(unname(pvals), unname(fx$pvals), tolerance = 1e-8)
})
