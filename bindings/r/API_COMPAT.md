# gsemr API Compatibility with R GenomicSEM

gsemr accepts all the same function signatures as R GenomicSEM. Most parameters are fully implemented. A small number of advanced parameters are accepted for compatibility but not yet functional — these print an informational message when used.

## Behavioral differences

No known behavioral differences. gsemr produces identical results to R GenomicSEM within numerical tolerance (~1e-5 for S, ~1e-8 for V).

## Fully implemented parameters

All core parameters for every function work identically to R GenomicSEM, including:

- **ldsc**: traits, sample.prev, population.prev, ld, wld, trait.names, n.blocks, chr, stand, select, chisq.max, sep_weights, ldsc.log
- **commonfactor**: covstruc, estimation
- **usermodel**: covstruc, estimation, model, std.lv, fix_resid, imp_cov, Q_Factor, toler, CFIcalc
- **munge**: files, hm3, trait.names, N, info.filter, maf.filter, column.names, overwrite, log.name
- **sumstats**: files, ref, trait.names, se.logit, OLS, linprob, N, info.filter, maf.filter, keep.indel, out, ambig, betas, direct.filter
- **userGWAS**: covstruc, SNPs, model, estimation, GC, sub, SNPSE, smooth_check, std.lv, fix_measurement, Q_SNP, printwarn, TWAS
- **commonfactorGWAS**: covstruc, SNPs, estimation, GC, SNPSE, smooth_check, TWAS
- **paLDSC**: covstruc, r, p, diag, save.pdf, fa, fm, nfactors
- **write.model**: Loadings, S_LD, cutoff, fix_resid, bifactor, mustload, common
- **rgmodel**: LDSCoutput, model, std.lv, estimation, sub
- **hdl**: traits, sample.prev, population.prev, LD.path, Nref, method
- **s_ldsc**: traits, sample.prev, population.prev, ld, wld, frq, trait.names, n.blocks, exclude_cont, ldsc.log
- **enrich**: s_baseline, s_annot, v_annot, model, params, fix, std.lv, toler, fixparam, tau, rm_flank
- **simLDSC**: covmat, N, ld, rPheno, int, N_overlap
- **multiSNP**: covstruc, model, beta, se, var_snp, ld_matrix, snp_names, SNPSE
- **multiGene**: covstruc, model, beta, se, var_gene, ld_matrix, gene_names, GeneSE, Genelist
- **summaryGLS**: covstruc, results
- **cores/parallel**: controls rayon thread count (parallel=FALSE → single-threaded)

## Not yet implemented

These parameters are accepted but have no effect. An informational message is printed when they are used.

### MPI
- `userGWAS(MPI=TRUE)` and `commonfactorGWAS(MPI=TRUE)` — MPI distributed computing. Not applicable to the Rust backend which uses shared-memory parallelism via rayon.

