# Contributing to gsem

## Getting Started

```sh
git clone https://github.com/PoHsuanLai/gsem.git
cd gsem
cargo build
cargo test --workspace
```

## Development Workflow

1. Create a branch from `main`
2. Make changes
3. Run `cargo fmt --all` and `cargo clippy --workspace`
4. Run `cargo test --workspace` (all tests must pass)
5. Open a pull request

## Project Structure

| Crate | Purpose |
|-------|---------|
| `gsem-matrix` | Matrix utilities (nearPD, vech, smoothing) |
| `gsem-ldsc` | LD Score Regression engine |
| `gsem-sem` | SEM fitting engine (DWLS/ML) |
| `gsem` | Main binary, I/O, munge, GWAS pipeline |
| `gsem-r` | R bindings (feature-gated) |
| `gsem-py` | Python bindings (feature-gated) |

## Code Style

- `cargo fmt` enforces formatting
- `cargo clippy` enforces lints
- CI runs both with `-Dwarnings` (warnings are errors)
- Prefer `anyhow::Result` for application code, `thiserror` for library errors
- Use `log::info!` / `log::warn!` for diagnostics, `eprintln!` for CLI user output

## Testing

- Unit tests live next to the code in `#[cfg(test)] mod tests { ... }`
- Test against known R GenomicSEM outputs when possible
- Use `approx` crate for floating-point comparisons in tests
