# gsemr

Rust-accelerated drop-in replacement for [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM). Same functions, same arguments, 2–100x faster.

## Install

```r
# Requires Rust toolchain (https://rustup.rs)
remotes::install_github("PoHsuanLai/gsem", subdir = "bindings/r")
```

## Usage

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

## Functions

All 18 R GenomicSEM functions are implemented:
`ldsc`, `commonfactor`, `usermodel`, `munge`, `sumstats`, `commonfactorGWAS`, `userGWAS`, `paLDSC`, `write.model`, `rgmodel`, `hdl`, `s_ldsc`, `enrich`, `simLDSC`, `multiSNP`, `multiGene`, `summaryGLS`, `convert_hdl_panels`

## Performance

Benchmarked against R GenomicSEM (3 traits, ~1.29M SNPs, N=50,000):

| Function | R GenomicSEM | gsemr | Speedup |
|----------|-------------|-------|---------|
| ldsc | 8.8s | 3.2s | 2.8x |
| commonfactor | 115ms | 11ms | 10x |
| usermodel | 91ms | 11ms | 8x |
| rgmodel | 127ms | 1ms | 124x |

Numerically identical results (S matrix max diff: 1.7e-5).

## Citation

> Grotzinger, A.D. et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. *Nat Genet* 51, 513–520 (2019).

## License

GPL-3.0
