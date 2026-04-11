# Post-0.1.0 TODO

## Coverage gaps

- [ ] **`convert_hdl_panels`: add Python binding.** Currently R-only. Requires an R `.rda` parser; consider `rdata-rs` or a thin Rust helper that shells out to R's conversion logic.
- [ ] **`convert_hdl_panels`: add CLI subcommand.** Blocked on the Python binding / Rust helper above.
- [ ] **`multiGene`: add CLI subcommand.** The engine exists in `crates/gsem/src/gwas/multi_snp.rs` and is exposed via R (`gsemr::multiGene`) and Python (`genomicsem.multi_gene`).  Wiring a clap subcommand is ~30 lines in `crates/gsem/src/main.rs` — clone `MultiSnpArgs` to `MultiGeneArgs` and route through the same dispatcher.

## Other follow-ups noted during 0.1.0 release audit

- [ ] `cargo fmt --all --check` has 6 hunks in 2 files (`sumstats_reader.rs`, `sumstats.rs`) on master — independent of the example-adding pass.  Run `cargo fmt --all` and commit.
- [ ] Consider extending the NumPy-wrapper modules (`_pa_ldsc.py`, `_enrich.py`, `_sim.py`) with `Examples` sections in the same style as the PyO3 `#[pyfunction]` docstrings so `help(genomicsem.pa_ldsc_extended)` etc. are consistent.

(Add items here as they come up — keep this file short and cross off completed work.)
