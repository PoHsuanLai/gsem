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
    n.blocks = fx$n_blocks,
    # Single-threaded: the parallel jackknife reduction order is
    # non-deterministic and perturbs the S off-diagonals at ~1e-7, which is a
    # large *relative* shift on these O(1e-4) entries. parallel=FALSE makes the
    # parity check reproducible run-to-run.
    parallel = FALSE
  )

  # Binding-layer parity tolerances. The S off-diagonals (and V) carry ~1e-7
  # absolute run-to-run noise from the block jackknife / parallel reduction
  # ordering; since these S entries are O(1e-4), that is up to ~5e-4 *relative*.
  # 1e-3 sits comfortably above that noise floor. (Bit-level S/V/I validation
  # lives in the core crates/gsem/tests/r_validation_ldsc.rs single-threaded
  # tests; here we only check the binding round-trips R's values.)
  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$I)), unname(fx$i), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$V)), unname(fx$v), tolerance = 1e-3)
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
    stand = TRUE, parallel = FALSE
  )

  # Unstandardized pieces still match, plus the correlation-scale outputs that
  # only appear under stand = TRUE (binding ldsc_result_to_list_stand path).
  #
  # gsemr's S_Stand == cov2cor(gsemr's S) exactly (verified), as does R's. The
  # only gap is the ~1e-7 run-to-run noise in the S off-diagonals, which the
  # correlation scaling of these tiny (O(1e-4)) covariances amplifies to ~2e-4
  # in S_Stand and ~3e-3 in V_Stand. Tolerances are sized to that amplified
  # noise floor, not loosened arbitrarily.
  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-3)
  expect_false(is.null(res$S_Stand), info = "stand=TRUE must return S_Stand")
  expect_equal(unname(as.matrix(res$S_Stand)), unname(fx$s_stand), tolerance = 5e-4)
  # V_Stand is the correlation-scaled sampling covariance — doubly noisy (V
  # jackknife noise carried through the standardization). Some entries are
  # near zero, so `tolerance` (which testthat scores *relative*ly) is the wrong
  # gauge; check the max ABSOLUTE deviation instead, against the observed
  # ~3e-3 noise floor.
  expect_lt(max(abs(unname(as.matrix(res$V_Stand)) - unname(fx$v_stand))), 5e-3)
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
    select = "ODD", parallel = FALSE
  )

  # Same jackknife/parallel noise floor as the base test (see comment there).
  expect_equal(unname(as.matrix(res$S)), unname(fx$s), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$I)), unname(fx$i), tolerance = 1e-3)
  expect_equal(unname(as.matrix(res$V)), unname(fx$v), tolerance = 1e-3)
})
