# Parity: gsemr::ldsc() must reproduce the R-reference S/V/I that R GenomicSEM
# produced for the same synthetic inputs (committed as ldsc_synth.json). This
# is the binding-layer counterpart of crates/gsem/tests/r_validation_ldsc.rs.

test_that("ldsc() reproduces the R-reference S/V/I matrices", {
  skip_if_no_fixtures()
  fx <- load_fixture("ldsc_synth")

  # The fixture records the exact call R made: munged files, LD dir, chr,
  # n.blocks and trait names — replay it through the Rust-backed wrapper.
  munged <- fixture_path(fx$munged_files)
  ld <- fixture_path(fx$ld_dir)

  res <- gsemr::ldsc(
    traits = munged,
    sample.prev = rep(NA, length(munged)),
    population.prev = rep(NA, length(munged)),
    ld = ld,
    wld = ld,
    trait.names = fx$trait_names,
    chr = fx$chr,
    n.blocks = fx$n_blocks
  )

  # S (genetic covariance) and I (intercepts) to tight relative tolerance;
  # V (sampling covariance) is noisier from jackknife so a touch looser.
  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-4)
  expect_equal(unname(as.matrix(res$I)), unname(fx$i), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$V)), unname(fx$v), tolerance = 1e-4)
})

test_that("ldsc(stand = TRUE) reproduces R's S_Stand / V_Stand", {
  skip_if_no_fixtures()
  fx <- load_fixture("ldsc_stand")
  munged <- fixture_path(fx$munged_files)
  ld <- fixture_path(fx$ld_dir)

  res <- gsemr::ldsc(
    traits = munged,
    sample.prev = rep(NA, length(munged)),
    population.prev = rep(NA, length(munged)),
    ld = ld, wld = ld,
    trait.names = fx$trait_names,
    chr = fx$chr, n.blocks = fx$n_blocks,
    stand = TRUE
  )

  # Unstandardized pieces still match, plus the correlation-scale outputs that
  # only appear under stand = TRUE (binding ldsc_result_to_list_stand path).
  #
  # gsemr's S_Stand == cov2cor(gsemr's S) exactly (verified), as does R's. The
  # only gap is the ~1e-7 run-to-run noise in the S off-diagonals, which the
  # correlation scaling of these tiny (O(1e-4)) covariances amplifies to ~2e-4
  # in S_Stand and ~3e-3 in V_Stand. Tolerances are sized to that amplified
  # noise floor, not loosened arbitrarily.
  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-4)
  expect_false(is.null(res$S_Stand), info = "stand=TRUE must return S_Stand")
  expect_equal(unname(as.matrix(res$S_Stand)), unname(fx$s_stand), tolerance = 5e-4)
  expect_equal(unname(as.matrix(res$V_Stand)), unname(fx$v_stand), tolerance = 5e-3)
  # S_Stand is a correlation matrix: unit diagonal.
  expect_equal(unname(diag(as.matrix(res$S_Stand))), rep(1, length(fx$trait_names)), tolerance = 1e-6)
})

test_that("ldsc(select = 'ODD') restricts to odd chromosomes", {
  skip_if_no_fixtures()
  fx <- load_fixture("ldsc_select_odd")
  munged <- fixture_path(fx$munged_files)
  ld <- fixture_path(fx$ld_dir)

  res <- gsemr::ldsc(
    traits = munged,
    sample.prev = rep(NA, length(munged)),
    population.prev = rep(NA, length(munged)),
    ld = ld, wld = ld,
    trait.names = fx$trait_names,
    chr = fx$chr, n.blocks = fx$n_blocks,
    select = "ODD"
  )

  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-4)
  expect_equal(unname(as.matrix(res$I)), unname(fx$i), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$V)), unname(fx$v), tolerance = 1e-4)
})
