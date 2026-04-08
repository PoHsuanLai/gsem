# gsem

[![CI](https://github.com/PoHsuanLai/gsem/actions/workflows/ci.yml/badge.svg)](https://github.com/PoHsuanLai/gsem/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Rust](https://img.shields.io/badge/rust-2024_edition-orange.svg)](https://www.rust-lang.org)

A Rust rewrite of [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) for multivariate genomic structural equation modeling on GWAS summary statistics.

## Why Rust?

- **Performance**: Per-SNP GWAS loop parallelized with rayon (vs R's foreach/doParallel)
- **No system deps**: Pure-Rust linear algebra via [faer](https://github.com/sarah-quinones/faer-rs) (no LAPACK/cmake)
- **Portable**: Single static binary, plus R and Python bindings
- **Correctness**: 91 unit tests, algorithms validated against R GenomicSEM output

## Architecture

```
gsem/
  crates/
    gsem-matrix/   Matrix utilities: nearPD, vech, PSD smoothing
    gsem-ldsc/     LD Score Regression with block jackknife
    gsem-sem/      Minimal SEM engine (DWLS/ML, L-BFGS, sandwich SE)
    gsem/          Main binary + I/O, munge, GWAS pipeline
    gsem-r/        R bindings via extendr
    gsem-py/       Python bindings via PyO3
```

Each core crate (`gsem-matrix`, `gsem-ldsc`, `gsem-sem`) is independently useful outside of gsem.

## Installation

### CLI

```sh
cargo install --path crates/gsem
```

### As a Rust library

```toml
[dependencies]
gsem-ldsc = { git = "https://github.com/PoHsuanLai/gsem" }
gsem-sem  = { git = "https://github.com/PoHsuanLai/gsem" }
```

### Python

```sh
cd crates/gsem-py
maturin develop --release --features pyo3
```

```python
import genomicsem as gsem

result = gsem.ldsc(
    traits=["trait1.sumstats.gz", "trait2.sumstats.gz"],
    sample_prev=[None, 0.5],
    pop_prev=[None, 0.01],
    ld="eur_w_ld_chr/",
)
result.S  # np.ndarray
result.V  # np.ndarray
```

### R

Requires the extendr toolchain. Build with:

```sh
cd crates/gsem-r
cargo build --release --features extendr
```

Then wrap with an R package (see `crates/gsem-r/` for API details).

## CLI Usage

### 1. Munge GWAS summary statistics

```sh
genomicsem munge \
  --files trait1.txt trait2.txt \
  --hm3 w_hm3.snplist \
  --trait-names Trait1 Trait2 \
  --out ./munged/
```

### 2. Run multivariate LDSC

```sh
genomicsem ldsc \
  --traits Trait1.sumstats.gz Trait2.sumstats.gz \
  --sample-prev NA,0.5 \
  --pop-prev NA,0.01 \
  --ld eur_w_ld_chr/ \
  --out ldsc_result.json
```

### 3. Fit SEM model

```sh
genomicsem sem \
  --covstruc ldsc_result.json \
  --model "F1 =~ NA*V1 + V2 + V3; F1 ~~ 1*F1" \
  --estimation DWLS \
  --out sem_result.tsv
```

### 4. Run multivariate GWAS

```sh
genomicsem gwas \
  --covstruc ldsc_result.json \
  --sumstats merged_sumstats.tsv \
  --model "F1 =~ NA*V1 + V2; F1 ~ SNP; F1 ~~ 1*F1; SNP ~~ SNP" \
  --estimation DWLS \
  --gc standard \
  --threads 8 \
  --out gwas_result.tsv
```

## Pipeline

```
Raw GWAS files ──► munge ──► .sumstats.gz ──► ldsc ──► S, V, I (JSON)
                                                           │
                                              ┌────────────┼────────────┐
                                              ▼            ▼            ▼
                                            sem          gwas      enrichment
                                         (usermodel)  (per-SNP)   (stratified)
```

## Crate Details

### gsem-matrix

| Feature | Description |
|---------|-------------|
| `nearest_pd()` | Higham (2002) alternating projections with Dykstra correction |
| `vech()` / `vech_reverse()` | Half-vectorization (column-major lower triangle) |
| `smooth_if_needed()` | PSD check via Cholesky, auto-smooth with nearPD |
| `cov_to_cor()` / `cor_to_cov()` | Covariance-correlation conversions |

### gsem-ldsc

| Feature | Description |
|---------|-------------|
| `ldsc()` | Full multivariate LDSC: weights, block WLS, jackknife V, liability scale |
| `estimate_h2()` | Single-trait heritability via LD regression |
| `estimate_gcov()` | Pairwise genetic covariance |
| Serialization | `LdscResult` supports JSON via faer's serde feature |

### gsem-sem

| Feature | Description |
|---------|-------------|
| Syntax parser | Lavaan-compatible: `=~`, `~~`, `~`, `:=`, labels, constraints |
| DWLS / ML | L-BFGS optimizer with Armijo line search and projected bounds |
| Sandwich SE | Corrected standard errors via `(Delta' W Delta)^{-1} Delta' W V W Delta (Delta' W Delta)^{-1}` |
| Fit indices | Chi-square, CFI, AIC, SRMR |
| V reorder | Port of `.rearrange()` for model variable ordering |

## Input Formats

### Munged sumstats (`.sumstats.gz`)

Tab-delimited, gzipped: `SNP  N  Z  A1  A2`

### Merged sumstats (for GWAS)

Tab-delimited: `SNP  A1  A2  MAF  beta.Trait1  se.Trait1  beta.Trait2  se.Trait2  ...`

### LD scores

Standard LDSC format: `{chr}.l2.ldscore.gz` files + `{chr}.l2.M_5_50` files in a directory.

## Differences from R GenomicSEM

| Aspect | R GenomicSEM | gsem |
|--------|-------------|------|
| SEM engine | lavaan | Built-in L-BFGS |
| Linear algebra | R's LAPACK bindings | faer |
| GWAS parallelism | foreach + doParallel | rayon |
| Serialization | R `.rds` files | JSON |
| Distribution | R package | CLI + Rust lib + R/Python bindings |

## Development

```sh
# Build
cargo build --release

# Test
cargo test --workspace

# Lint
cargo clippy --workspace

# Run with logging
RUST_LOG=info genomicsem ldsc --traits ...

# Build Python wheel
cd crates/gsem-py && maturin build --release --features pyo3

# Build with R bindings
cargo build --release --features gsem-r/extendr
```

## Citation

If you use this software, please cite the original GenomicSEM paper:

> Grotzinger, A.D., Rhemtulla, M., de Vlaming, R. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513-520 (2019). https://doi.org/10.1038/s41588-019-0360-8

## License

GPL-3.0 -- see [LICENSE](LICENSE).
