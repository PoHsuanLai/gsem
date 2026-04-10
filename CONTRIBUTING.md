# Contributing to gsem

## Getting Started

```sh
git clone https://github.com/PoHsuanLai/gsem.git
cd gsem
cargo build
cargo test --workspace
```

## Development Workflow

1. Create a branch from `master`
2. Make changes
3. Run `cargo fmt --all` and `cargo clippy --workspace`
4. Run `cargo test --workspace` (all tests must pass)
5. Open a pull request

## Project Structure

Workspace crates (under `crates/`):

| Crate | Purpose |
|-------|---------|
| `gsem-matrix` | Matrix utilities (nearPD, vech, smoothing) |
| `gsem-ldsc` | LD Score Regression engine |
| `gsem-sem` | SEM fitting engine (DWLS/ML) |
| `gsem` | Main binary, I/O, munge, GWAS pipeline |

Language bindings (standalone crates, not part of the workspace):

| Path | Crate | Purpose |
|------|-------|---------|
| `bindings/r` | `gsemr` | R package (extendr) |
| `bindings/python` | `genomicsem` | Python package (PyO3 / maturin) |

To build the bindings locally:

```sh
# R
Rscript -e 'install.packages("bindings/r", repos=NULL, type="source")'

# Python
cd bindings/python && maturin develop --release
```

Before hacking on the estimator, sandwich-SE path, or per-SNP GWAS
loops, read [`ARCHITECTURE.md`](./ARCHITECTURE.md) — it documents the
design choices and numerical-parity constraints that test changes in
those areas need to respect.

## Code Style

- `cargo fmt` enforces formatting
- `cargo clippy` enforces lints
- CI runs both with `-Dwarnings` (warnings are errors)
- Prefer `anyhow::Result` for application code, `thiserror` for library errors
- Use `log::info!` / `log::warn!` for diagnostics, `eprintln!` for CLI user output

## Testing

- Unit tests live next to the code in `#[cfg(test)] mod tests { ... }`
- `crates/gsem/tests/r_validation.rs` pins numerical parity against R
  GenomicSEM at the primitive level (`vech`, `nearPD`, `V_SNP`, `S_Full`,
  `V_Full`, SEM fit, commonfactor/userGWAS per-SNP). Extend it whenever
  you touch a numerical path.
- `bench/` contains end-to-end parity harnesses against R GenomicSEM on
  real PGC data. Useful for catching regressions the unit tests miss.
- Use `approx` crate for floating-point comparisons in tests
