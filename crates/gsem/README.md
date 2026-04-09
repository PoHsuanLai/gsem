# gsem

Rust implementation of [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) — multivariate LD Score Regression and Structural Equation Modeling on GWAS summary statistics.

## Install

```sh
cargo install gsem
```

## CLI

```sh
# Munge GWAS summary statistics
gsem munge --files trait1.txt trait2.txt --hm3 w_hm3.snplist --trait-names T1 T2

# Run multivariate LDSC
gsem ldsc --traits T1.sumstats.gz T2.sumstats.gz --ld eur_w_ld_chr/ --out ldsc.json

# Fit SEM model
gsem usermodel --covstruc ldsc.json --model "F1 =~ NA*V1 + V2; F1 ~~ 1*F1" --out sem.tsv

# Common factor
gsem commonfactor --covstruc ldsc.json --out cf.tsv

# Multivariate GWAS
gsem userGWAS --covstruc ldsc.json --sumstats merged.tsv --model "F1 =~ NA*V1 + V2; F1 ~ SNP"
gsem commonfactorGWAS --covstruc ldsc.json --sumstats merged.tsv

# Merge sumstats
gsem sumstats --files T1.sumstats.gz T2.sumstats.gz --ref eur_w_ld_chr/ --out merged.tsv
```

## As a library

```toml
[dependencies]
gsem-ldsc = "0.1"
gsem-sem  = "0.1"
```

```rust
use gsem_ldsc::{ldsc, LdscConfig, TraitSumstats};
use gsem_sem::{syntax, model::Model, estimator};
```

## Workspace crates

| Crate | Description |
|-------|-------------|
| `gsem-matrix` | nearPD, vech, PSD smoothing |
| `gsem-ldsc` | LDSC, HDL, stratified LDSC, block jackknife |
| `gsem-sem` | SEM engine (DWLS/ML, L-BFGS, sandwich SE) |
| `gsem` | CLI, I/O, munge, sumstats merge, GWAS pipeline |

## Citation

> Grotzinger, A.D. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513–520 (2019).

## License

GPL-3.0
