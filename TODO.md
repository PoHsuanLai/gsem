# Post-0.1.0 TODO

## Coverage gaps

- [ ] **`convert_hdl_panels`: add Python binding.** Currently R-only. Requires an R `.rda` parser; consider `rdata-rs` or a thin Rust helper that shells out to R's conversion logic.
- [ ] **`convert_hdl_panels`: add CLI subcommand.** Blocked on the Python binding / Rust helper above.
- [ ] **`multiGene`: add CLI subcommand.** The engine exists in `crates/gsem/src/gwas/multi_snp.rs` and is exposed via R (`gsemr::multiGene`) and Python (`genomicsem.multi_gene`).  Wiring a clap subcommand is ~30 lines in `crates/gsem/src/main.rs` — clone `MultiSnpArgs` to `MultiGeneArgs` and route through the same dispatcher.

## Known issues in v0.1.0 (to fix in 0.1.1)

- [ ] **PyPI 0.1.0 project page missing description.** The `0.1.0` wheels were published before `readme = "README.md"` was added to `pyproject.toml`, so `pypi.org/project/genomicsem/` shows "The author of this package has not provided a project description". Fixed in master for 0.1.1 onwards; the 0.1.0 page can't be updated because PyPI doesn't allow re-uploading the same version.

## Also noted while fixing Windows R CMD check

Non-blocking warnings that still surface in `R CMD check --as-cran`
after the v0.1.0 fix pass — they're reported but don't fail the build
(we use `error_on = "error"`). Worth cleaning up for a possible CRAN
submission:

- [ ] **`_exit` / `abort` / `exit` symbols in `libgsemr.a`** — Rust's stdlib pulls these in for panic handling. CRAN's "portable packages" rule flags them but they're unavoidable without dropping to `no_std`. The standard workaround (used across Rust-backed R packages) is to ship as a non-CRAN source tarball, which is what we do today. Revisit if CRAN becomes a target.
- [ ] **`unlockBinding(".gsemr_cfgwas_warned", ns)` in `commonfactorGWAS.R`** — flagged as "possibly unsafe call". Rewrite the first-use warning flag to use an environment (`.gsemr_env$cfgwas_warned`) instead of mutating the namespace, which R CMD check accepts silently.

## Other follow-ups noted during 0.1.0 release audit

- [ ] Consider extending the NumPy-wrapper modules (`_pa_ldsc.py`, `_enrich.py`, `_sim.py`) with `Examples` sections in the same style as the PyO3 `#[pyfunction]` docstrings so `help(genomicsem.pa_ldsc_extended)` etc. are consistent.

(Add items here as they come up — keep this file short and cross off completed work.)
