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
