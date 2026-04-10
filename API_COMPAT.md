# gsemr API Compatibility with R GenomicSEM

gsemr accepts all the same function signatures as R GenomicSEM. Most parameters are fully implemented. A small number of advanced parameters are accepted for compatibility but not yet functional — these print an informational message when used.

## Behavioral differences

For the bulk of the pipeline (LDSC, S/V/I matrices, single-factor and 2-factor SEM, `usermodel`, `userGWAS` per-SNP point estimates) gsemr produces results identical to R GenomicSEM within numerical tolerance (~1e-5 for S, ~1e-8 for V, ~1e-4 for SEM point estimates).

The known differences are listed below.

### `commonfactorGWAS` parameterization

R `GenomicSEM::commonfactorGWAS` runs lavaan internally with **marker-indicator** identification (first loading fixed to 1, factor variance free), and on degenerate `sample.nobs=2` per-SNP surfaces it appears to free `F1~~F1` jointly with `F1~SNP` per SNP. R's own `commonfactorGWAS` and `userGWAS` outputs disagree in sign for ~80% of SNPs on the same data — the two parameterizations are not numerically equivalent on this kind of surface.

gsemr `commonfactorGWAS` defaults to `Identification::FixedVariance` (factor variance fixed to 1, all loadings free), which matches R's `userGWAS` per-SNP estimates exactly. An opt-in `Identification::MarkerIndicator` mode is available for users who want R-style anchoring; it produces output that is mathematically consistent with the `FixedVariance` path under the loading-rescaling identity but does **not** numerically match R's `commonfactorGWAS` per-SNP estimates because the per-SNP free-parameter set differs.

**Recommendation:** for parity with R `userGWAS`, use either gsemr `userGWAS` or gsemr `commonfactorGWAS` (default). For parity with R `commonfactorGWAS` per-SNP signs and magnitudes, no exact replacement is currently provided.

### Standard errors are sandwich-corrected

gsemr always reports sandwich (robust) standard errors for SEM and per-SNP GWAS parameters. R GenomicSEM passes through lavaan's default `se = "standard"` (information-matrix SE). On well-conditioned problems the two agree closely; on degenerate `sample.nobs=2` per-SNP surfaces they routinely differ by tens of percent. The point estimates and z-statistics have the same interpretation in both packages, but a per-SNP `SE` value will not numerically match R.

### Optimizer

gsemr uses L-BFGS for SEM and per-SNP GWAS fits; R/lavaan uses nlminb (trust-region quasi-Newton). On well-identified problems they reach the same minimum to many decimals. On non-convex objectives (e.g., a Heywood case where a residual variance is negative), they can pick different local minima. The single-factor SEM and `usermodel` test fixtures cover both situations and pass within tolerance.

### Heywood cases (negative residual variances)

Both packages allow negative residual variances by default — neither bounds `theta` away from zero. gsemr's `commonfactor` and `userGWAS` pipelines will land on the same Heywood solutions R does, but the per-SNP fits on top of a Heywood baseline are sensitive to starting point; see the `commonfactorGWAS` note above.

## Fully implemented parameters

All core parameters for every function work identically to R GenomicSEM, including:

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

These parameters are accepted but have no effect. An informational message is printed when they are used.

### MPI
- `userGWAS(MPI=TRUE)` and `commonfactorGWAS(MPI=TRUE)` — MPI distributed computing. Not applicable to the Rust backend which uses shared-memory parallelism via rayon.

