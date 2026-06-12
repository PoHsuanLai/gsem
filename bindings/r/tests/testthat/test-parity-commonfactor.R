# Parity: gsemr::commonfactor() must reproduce R GenomicSEM's parameter
# estimates and fit indices (committed as commonfactor.json). Exercises the
# binding's covstruc list -> Rust -> results data.frame round-trip, the
# trickiest marshalling path.

test_that("commonfactor() reproduces R-reference parameter estimates", {
  skip_if_no_fixtures()
  fx <- load_fixture("commonfactor")

  # A real covstruc from ldsc() carries S, V and I; commonfactor only consumes
  # S and V, but the binding requires the I slot to be present, so supply an
  # identity placeholder (unused by the common-factor fit).
  k <- nrow(as.matrix(fx$s))
  covstruc <- list(
    S = as.matrix(fx$s),
    V = as.matrix(fx$v),
    I = diag(k)
  )
  nm <- paste0("V", seq_len(k))
  dimnames(covstruc$S) <- list(nm, nm)
  dimnames(covstruc$I) <- list(nm, nm)

  cf <- gsemr::commonfactor(covstruc, estimation = "DWLS")

  # Match each expected parameter by (lhs, op, rhs). The common factor's
  # orientation is unidentified (F1 vs -F1), so loadings are compared on |est|.
  ref <- fx$parameters
  got <- cf$results
  expect_true(nrow(got) >= nrow(ref), "fewer parameters than reference")

  for (i in seq_len(nrow(ref))) {
    row <- ref[i, ]
    match_idx <- which(
      got$lhs == row$lhs & got$op == row$op & got$rhs == row$rhs
    )
    expect_length(match_idx, 1)
    expect_equal(
      abs(got$est[match_idx]), abs(row$est),
      tolerance = 1e-4,
      info = sprintf("%s %s %s", row$lhs, row$op, row$rhs)
    )
  }

  # Fit indices: a saturated 1-factor / 3-indicator model has df = 0, chisq ~ 0.
  expect_equal(unname(cf$modelfit["chisq"]), as.numeric(fx$chisq), tolerance = 1e-4)
  expect_equal(unname(cf$modelfit["SRMR"]), as.numeric(fx$srmr), tolerance = 1e-4)
})
