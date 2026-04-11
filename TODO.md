# Post-0.1.0 TODO

## Coverage gaps

- [ ] **`convert_hdl_panels`: add Python binding.** Currently R-only. Requires an R `.rda` parser; consider `rdata-rs` or a thin Rust helper that shells out to R's conversion logic.
- [ ] **`convert_hdl_panels`: add CLI subcommand.** Blocked on the Python binding / Rust helper above.
- [ ] **`multiGene`: add CLI subcommand.** The engine exists in `crates/gsem/src/gwas/multi_snp.rs` and is exposed via R (`gsemr::multiGene`) and Python (`genomicsem.multi_gene`).  Wiring a clap subcommand is ~30 lines in `crates/gsem/src/main.rs` — clone `MultiSnpArgs` to `MultiGeneArgs` and route through the same dispatcher.

## Known issues in v0.1.0 (to fix in 0.1.1)

- [ ] **PyPI 0.1.0 project page missing description.** The `0.1.0` wheels were published before `readme = "README.md"` was added to `pyproject.toml`, so `pypi.org/project/genomicsem/` shows "The author of this package has not provided a project description". Fixed in master for 0.1.1 onwards; the 0.1.0 page can't be updated because PyPI doesn't allow re-uploading the same version.
- [ ] **R on Windows**: first-pass `configure.win` + `src/Makevars.win` were added in master and the manual publish.yml re-trigger uploaded a Windows-compatible tarball. If further Windows issues surface, the likely culprits are (a) Rtools `sh.exe` path detection, (b) MinGW native-static-libs drift between rustc versions (run `RUSTFLAGS=--print=native-static-libs cargo build --target=x86_64-pc-windows-gnu` in a Windows Rtools shell to re-check the link flags), or (c) long-path issues with `target/` inside R's temp install dir.

## Other follow-ups noted during 0.1.0 release audit

- [ ] Consider extending the NumPy-wrapper modules (`_pa_ldsc.py`, `_enrich.py`, `_sim.py`) with `Examples` sections in the same style as the PyO3 `#[pyfunction]` docstrings so `help(genomicsem.pa_ldsc_extended)` etc. are consistent.

(Add items here as they come up — keep this file short and cross off completed work.)
