# genomicsem

Rust-accelerated [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) for Python. Multivariate LD Score Regression and Structural Equation Modeling on GWAS summary statistics.

> See [`/API_COMPAT.md`](../../API_COMPAT.md) for parameter compatibility
> with R GenomicSEM and [`/ARCHITECTURE.md`](../../ARCHITECTURE.md) for
> algorithmic and performance notes.

## Install

```sh
pip install genomicsem
```

From source (requires Rust toolchain):

```sh
pip install maturin
cd bindings/python && maturin develop --release
```

## Usage

```python
import genomicsem as gsem

# LDSC
result = gsem.ldsc(
    traits=["trait1.sumstats.gz", "trait2.sumstats.gz"],
    sample_prev=[None, None],
    pop_prev=[None, None],
    ld="eur_w_ld_chr/",
)
result.S       # np.ndarray — genetic covariance matrix
result.V       # np.ndarray — sampling covariance matrix
result.I       # np.ndarray — intercept matrix
result.n       # list — sample sizes
result.m_total # float — number of SNPs

# Common factor model
cf = gsem.commonfactor(result.to_json())

# User-specified SEM
um = gsem.usermodel(
    result.to_json(),
    model="F1 =~ NA*V1 + V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2",
    estimation="DWLS",
)

# Munge
gsem.munge(
    files=["gwas1.txt", "gwas2.txt"],
    hm3="w_hm3.snplist",
    trait_names=["T1", "T2"],
)

# Merge sumstats
gsem.sumstats(
    files=["T1.sumstats.gz", "T2.sumstats.gz"],
    ref_dir="eur_w_ld_chr/",
    trait_names=["V1", "V2"],
)

# GWAS
gsem.commonfactor_gwas(result.to_json(), "merged_sumstats.tsv")
gsem.user_gwas(result.to_json(), "merged_sumstats.tsv", model="F1 =~ NA*V1 + V2\nF1 ~ SNP")

# Parallel analysis
gsem.parallel_analysis(result.to_json(), n_sim=500)

# Auto-generate model syntax
gsem.write_model(loadings=[[0.7, 0.0], [0.6, 0.0], [0.0, 0.8]], names=["V1", "V2", "V3"])

# Genetic correlation
gsem.rgmodel(result.to_json())
```

## Functions

`ldsc`, `commonfactor`, `usermodel`, `munge`, `sumstats`, `commonfactor_gwas`, `user_gwas`, `parallel_analysis`, `write_model`, `rgmodel`

## Citation

> Grotzinger, A.D. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513–520 (2019).

## License

GPL-3.0
