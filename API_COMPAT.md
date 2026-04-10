# API Compatibility with R GenomicSEM

This document tracks which R GenomicSEM parameters are accepted by the
Rust port, and which are accepted as stubs.

It applies to all three frontends:

- `gsemr` — R package (`library(gsemr)`)
- `genomicsem` — Python package (`import genomicsem`)
- `gsem` — CLI (`gsem userGWAS --...`)

Parameter names match R GenomicSEM exactly in the R binding. The Python
binding uses snake_case equivalents (`fix.measurement` → `fix_measurement`,
`sample.prev` → `sample_prev`, etc.). The CLI uses `--kebab-case` flags
with the same meanings.

For the **algorithmic and numerical differences** between the Rust port
and R GenomicSEM — why `commonfactorGWAS` uses a different parameterization,
why SEs differ, why optimizer behavior differs on Heywood cases — see
[`ARCHITECTURE.md`](./ARCHITECTURE.md).

## Behavioral notes (summary)

For the bulk of the pipeline (LDSC, S/V/I matrices, single-factor and
2-factor SEM, `usermodel`, `userGWAS` per-SNP point estimates) the Rust
port produces results identical to R GenomicSEM within numerical
tolerance (~1e-5 for S, ~1e-8 for V, ~1e-4 for SEM point estimates).

The known places where outputs differ from R:

- **`commonfactorGWAS` per-SNP signs and magnitudes** do not match R's
  `commonfactorGWAS` by default — they match R's `userGWAS`. See
  [`ARCHITECTURE.md §3.3`](./ARCHITECTURE.md#33-commonfactorgwas-parameterization).
- **Standard errors** are sandwich (robust) SEs, not lavaan's
  information-matrix SEs. See
  [`ARCHITECTURE.md §3.2`](./ARCHITECTURE.md#32-standard-errors-sandwich-rust-vs-information-matrix-r--lavaan).
- **SEM optimizer** is L-BFGS, not nlminb; converges to the same minimum
  on well-conditioned problems but can diverge on Heywood cases. See
  [`ARCHITECTURE.md §3.1`](./ARCHITECTURE.md#31-sem-optimizer-l-bfgs-rust-vs-nlminb-r--lavaan).
- **Heywood cases** are allowed by default in both packages; see
  [`ARCHITECTURE.md §3.5`](./ARCHITECTURE.md#35-heywood-cases-negative-residual-variances).

## Fully implemented parameters

All core parameters for every function work identically to R GenomicSEM,
including:

- **ldsc**: traits, sample.prev, population.prev, ld, wld, trait.names, n.blocks, chr, stand, select, chisq.max, sep_weights, ldsc.log, parallel, cores
- **commonfactor**: covstruc, estimation
- **usermodel**: covstruc, estimation, model, std.lv, fix_resid, imp_cov, Q_Factor, toler, CFIcalc
- **munge**: files, hm3, trait.names, N, info.filter, maf.filter, column.names, overwrite, log.name
- **sumstats**: files, ref, trait.names, se.logit, OLS, linprob, N, info.filter, maf.filter, keep.indel, out, ambig, betas, direct.filter
- **userGWAS**: covstruc, SNPs, model, estimation, GC, sub, SNPSE, smooth_check, std.lv, fix_measurement, Q_SNP, printwarn, TWAS, parallel, cores
- **commonfactorGWAS**: covstruc, SNPs, estimation, GC, SNPSE, smooth_check, TWAS, identification, parallel, cores
- **paLDSC**: covstruc, r, p, diag, save.pdf, fa, fm, nfactors, parallel, cores
- **write.model**: Loadings, S_LD, cutoff, fix_resid, bifactor, mustload, common
- **rgmodel**: LDSCoutput, model, std.lv, estimation, sub
- **hdl**: traits, sample.prev, population.prev, LD.path, Nref, method
- **s_ldsc**: traits, sample.prev, population.prev, ld, wld, frq, trait.names, n.blocks, exclude_cont, ldsc.log
- **enrich**: s_baseline, s_annot, v_annot, model, params, fix, std.lv, toler, fixparam, tau, rm_flank
- **simLDSC**: covmat, N, ld, rPheno, int, N_overlap
- **multiSNP**: covstruc, model, beta, se, var_snp, ld_matrix, snp_names, SNPSE
- **multiGene**: covstruc, model, beta, se, var_gene, ld_matrix, gene_names, GeneSE, Genelist
- **summaryGLS**: covstruc, results
- **cores/parallel**: per-call thread budget for `ldsc`, `userGWAS`, `commonfactorGWAS`, and `paLDSC`. Each call builds its own local rayon pool, so concurrent calls do not share thread state and `parallel=FALSE` is fully scoped. Currently a no-op for `munge`, `sumstats`, and `simLDSC` (their underlying implementations are serial).

## Not yet implemented

These parameters are accepted but have no effect. An informational
message is printed when they are used.

### MPI

- `userGWAS(MPI=TRUE)` and `commonfactorGWAS(MPI=TRUE)` — MPI
  distributed computing. Not applicable to the Rust backend, which uses
  shared-memory parallelism via rayon. For distributed workloads, split
  the sumstats file and run independent `gsem` CLI processes on each
  chunk, then concatenate the TSV outputs.
