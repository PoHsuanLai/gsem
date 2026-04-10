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
result.s       # np.ndarray — genetic covariance matrix
result.v       # np.ndarray — sampling covariance matrix
result.i_mat   # np.ndarray — intercept matrix
result.n       # list — sample sizes
result.m_total # float — number of SNPs

# Every downstream function accepts the `result` object directly — no
# need to call `result.to_json()`. You can also pass a plain dict with
# `s`, `v`, `i_mat`, `n_vec`, `m` fields, or a dict using the uppercase
# `S`, `V`, `I`, `N`, `m` keys that match the R GenomicSEM convention.

# Common factor model
cf = gsem.commonfactor(result)
# cf is a dict: {parameters: {lhs, op, rhs, est, se, z, p}, chisq, df, ...}

# User-specified SEM
um = gsem.usermodel(
    result,
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

# GWAS — results come back as a columnar dict. Feed it straight into pandas:
#   import pandas as pd; df = pd.DataFrame({k: v for k, v in res.items() if k != "params"})
res = gsem.commonfactor_gwas(result, "merged_sumstats.tsv")
res = gsem.user_gwas(result, "merged_sumstats.tsv", model="F1 =~ NA*V1 + V2\nF1 ~ SNP")

# Parallel analysis
gsem.parallel_analysis(result.s, result.v, r=500)

# Auto-generate model syntax
gsem.write_model(loadings=[[0.7, 0.0], [0.6, 0.0], [0.0, 0.8]], names=["V1", "V2", "V3"])

# Genetic correlation
gsem.rgmodel(result, model="")
```

`result.to_json()` / `gsem.LdscResult.from_json(s)` are still available as
optional compatibility helpers when you need to checkpoint an LDSC run to
disk, but they are no longer on the hot path between functions.

## Functions

`ldsc`, `commonfactor`, `usermodel`, `munge`, `sumstats`, `commonfactor_gwas`, `user_gwas`, `parallel_analysis`, `write_model`, `rgmodel`

## Citation

> Grotzinger, A.D. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513–520 (2019).

## License

GPL-3.0
