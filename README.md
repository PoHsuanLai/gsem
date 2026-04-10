# gsem

[![CI](https://github.com/PoHsuanLai/gsem/actions/workflows/ci.yml/badge.svg)](https://github.com/PoHsuanLai/gsem/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/gsem.svg)](https://crates.io/crates/gsem)
[![PyPI](https://img.shields.io/pypi/v/genomicsem.svg)](https://pypi.org/project/genomicsem/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Rust implementation of [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) — multivariate LD Score Regression and Structural Equation Modeling on GWAS summary statistics.

Available as an **R package** (`gsemr`), **Python package** (`genomicsem`), **Rust crates**, and a **CLI**.

## Performance

Benchmarked against R GenomicSEM on 3 simulated traits (~1.29M SNPs, N=50,000):

| Function | R | gsemr (Rust) | Speedup |
|----------|---|-------------|---------|
| ldsc | 8.8s | 3.2s | 2.8x |
| commonfactor | 115ms | 11ms | 10x |
| usermodel | 91ms | 11ms | 8x |
| rgmodel | 127ms | 1ms | 124x |

S matrix max diff: 1.7e-5. V matrix max diff: 8.2e-8.

See [`ARCHITECTURE.md §2`](./ARCHITECTURE.md#2-why-its-faster) for the
full performance breakdown.

## Documentation

- **[`API_COMPAT.md`](./API_COMPAT.md)** — which R GenomicSEM parameters
  are implemented in the Rust port, and which are accepted as stubs.
- **[`ARCHITECTURE.md`](./ARCHITECTURE.md)** — why the Rust port runs
  faster than R GenomicSEM, and where its algorithms and numerical
  outputs diverge from R.
- **[`CONTRIBUTING.md`](./CONTRIBUTING.md)** — dev setup and
  contribution guide.

## Install

### R (gsemr)

```r
# Requires Rust toolchain (https://rustup.rs)
remotes::install_github("PoHsuanLai/gsem", subdir = "bindings/r")
```

### Python (genomicsem)

```sh
# Requires Rust toolchain
pip install genomicsem
# or from source:
cd bindings/python && maturin develop --release
```

### CLI

```sh
cargo install --path crates/gsem
```

### Rust crates

```toml
[dependencies]
gsem-ldsc = { git = "https://github.com/PoHsuanLai/gsem" }
gsem-sem  = { git = "https://github.com/PoHsuanLai/gsem" }
```

## Usage

### R

`gsemr` is a drop-in replacement for GenomicSEM. Same function names, same arguments.

```r
library(gsemr)

# LDSC
covstruc <- ldsc(
  traits = c("trait1.sumstats.gz", "trait2.sumstats.gz", "trait3.sumstats.gz"),
  sample.prev = c(NA, NA, NA),
  population.prev = c(NA, NA, NA),
  ld = "eur_w_ld_chr/",
  wld = "eur_w_ld_chr/",
  trait.names = c("V1", "V2", "V3")
)

# Common factor
cf <- commonfactor(covstruc, estimation = "DWLS")

# User-specified model
um <- usermodel(covstruc, model = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3")

# Munge
munge(files = c("gwas1.txt", "gwas2.txt"), hm3 = "w_hm3.snplist", trait.names = c("T1", "T2"))

# Merge sumstats
sumstats(files = c("T1.sumstats.gz", "T2.sumstats.gz"), ref = "eur_w_ld_chr/", trait.names = c("V1", "V2"))

# GWAS
commonfactorGWAS(covstruc, SNPs = "merged_sumstats.tsv")
userGWAS(covstruc, SNPs = "merged_sumstats.tsv", model = "F1 =~ NA*V1 + V2\nF1 ~ SNP")

# Other
paLDSC(covstruc, r = 500)
rgmodel(covstruc)
write.model(loadings_matrix, covstruc$S)
hdl(traits, sample.prev, population.prev, LD.path = "hdl_panels/")
```

All 18 R GenomicSEM functions are implemented:
`ldsc`, `commonfactor`, `usermodel`, `munge`, `sumstats`, `commonfactorGWAS`, `userGWAS`, `paLDSC`, `write.model`, `rgmodel`, `hdl`, `s_ldsc`, `enrich`, `simLDSC`, `multiSNP`, `multiGene`, `summaryGLS`, `convert_hdl_panels`

### Python

```python
import genomicsem as gsem

result = gsem.ldsc(
    traits=["trait1.sumstats.gz", "trait2.sumstats.gz"],
    sample_prev=[None, None],
    pop_prev=[None, None],
    ld="eur_w_ld_chr/",
)
result.S   # np.ndarray
result.V   # np.ndarray
result.I   # np.ndarray

cf = gsem.commonfactor(result.to_json())
um = gsem.usermodel(result.to_json(), model="F1 =~ NA*V1 + V2\nF1 ~~ 1*F1")
```

### CLI

```sh
# Munge
gsem munge --files trait1.txt trait2.txt --hm3 w_hm3.snplist --trait-names T1 T2

# LDSC
gsem ldsc --traits T1.sumstats.gz T2.sumstats.gz --ld eur_w_ld_chr/ --out ldsc.json

# SEM
gsem usermodel --covstruc ldsc.json --model "F1 =~ NA*V1 + V2; F1 ~~ 1*F1" --out sem.tsv

# GWAS
gsem userGWAS --covstruc ldsc.json --sumstats merged.tsv --model "F1 =~ NA*V1 + V2; F1 ~ SNP"
```

## Architecture

```
crates/
  gsem-matrix/     nearPD, vech, PSD smoothing
  gsem-ldsc/       LDSC, HDL, stratified LDSC, block jackknife
  gsem-sem/        SEM engine (DWLS/ML, L-BFGS, sandwich SE, fit indices)
  gsem/            CLI, I/O, munge, sumstats merge, GWAS pipeline
bindings/
  r/               R package (gsemr) via extendr
  python/          Python package (genomicsem) via PyO3 + maturin
```

## Development

```sh
cargo build --release
cargo test --workspace      # 117 tests including R validation
cargo bench                 # criterion micro-benchmarks
cd bench && Rscript benchmark.R  # R vs gsemr comparison
```

## Citation

If you use this software, please cite the original GenomicSEM paper:

> Grotzinger, A.D., Rhemtulla, M., de Vlaming, R. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513–520 (2019). https://doi.org/10.1038/s41588-019-0360-8

## License

GPL-3.0 — see [LICENSE](LICENSE).
