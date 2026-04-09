# gsemr API Compatibility with R GenomicSEM

gsemr accepts all the same function signatures as R GenomicSEM. Most parameters are fully implemented. Some advanced parameters are accepted for compatibility but not yet functional — these print an informational message when used.

## Fully implemented parameters

All core parameters for every function work identically to R GenomicSEM, including:

- **ldsc**: traits, sample.prev, population.prev, ld, wld, trait.names, n.blocks, chr, stand, select, chisq.max, sep_weights, ldsc.log
- **commonfactor**: covstruc, estimation
- **usermodel**: covstruc, estimation, model, std.lv, fix_resid
- **munge**: files, hm3, trait.names, N, info.filter, maf.filter, column.names, overwrite, log.name
- **sumstats**: files, ref, trait.names, se.logit, OLS, linprob, N, info.filter, maf.filter, keep.indel, out
- **userGWAS**: covstruc, SNPs, model, estimation, GC, sub, SNPSE, smooth_check, std.lv, fix_measurement, Q_SNP, printwarn
- **commonfactorGWAS**: covstruc, SNPs, estimation, GC, SNPSE, smooth_check
- **paLDSC**: covstruc, r
- **write.model**: Loadings, S_LD, cutoff, fix_resid, bifactor
- **rgmodel**: LDSCoutput, model, std.lv, estimation, sub
- **hdl**: traits, sample.prev, population.prev, LD.path, Nref, method
- **s_ldsc**: traits, sample.prev, population.prev, ld, wld, frq, trait.names, n.blocks, exclude_cont, ldsc.log
- **enrich**: s_baseline, s_annot, v_annot
- **simLDSC**: covmat, N, ld
- **multiSNP**: covstruc, model, beta, se, var_snp, ld_matrix, snp_names
- **multiGene**: covstruc, model, beta, se, var_gene, ld_matrix, gene_names
- **summaryGLS**: covstruc, results
- **cores/parallel**: controls rayon thread count (parallel=FALSE → single-threaded)

## Not yet implemented

These parameters are accepted but have no effect. An informational message is printed when they are used.

### TWAS mode
- `userGWAS(TWAS=TRUE)` and `commonfactorGWAS(TWAS=TRUE)` — TWAS gene-level analysis mode. The `snp_label` infrastructure exists but the full TWAS pipeline (gene expression weights, panel integration) is not yet ported.

### MPI
- `userGWAS(MPI=TRUE)` and `commonfactorGWAS(MPI=TRUE)` — MPI distributed computing. Not applicable to the Rust backend which uses shared-memory parallelism via rayon.

### Enrichment model options
- `enrich(params, fix, std.lv, rm_flank, tau, toler, fixparam)` — The Rust enrichment implementation uses a simplified test. The full lavaan-based enrichment model with custom parameter constraints is not yet ported.

### Parallel analysis variants
- `paLDSC(p, diag, fa, fm, nfactors)` — Custom percentile thresholds, diagonal-only mode, factor analysis mode, and factor extraction methods. The Rust implementation uses the default 95th percentile with full V matrix.
- `paLDSC(save.pdf)` — PDF plot generation. Not supported in the Rust backend.

### Simulation parameters
- `simLDSC(rPheno, int, N_overlap)` — Phenotypic correlation, intercept inflation, and sample overlap parameters. The Rust simulation generates Z-statistics directly from the genetic covariance matrix.

### Model options
- `usermodel(imp_cov)` — Return model-implied covariance matrix. Not yet implemented.
- `usermodel(Q_Factor)` — Q factor heterogeneity statistic. Not yet implemented.
- `usermodel(CFIcalc)` — CFI is always computed in gsemr; this flag has no effect.

### Sumstats options
- `sumstats(betas)` — Custom beta column name overrides. gsemr uses automatic column detection.
- `sumstats(ambig)` — Keep strand-ambiguous SNPs. gsemr always removes them.
- `sumstats(direct.filter)` — Direct allele frequency filter. Not yet implemented.

### write.model options
- `write.model(mustload)` — Require all variables to load on at least one factor. Not yet implemented.
- `write.model(common)` — Include a common factor. Not yet implemented.

### Multi-SNP/Gene options
- `multiSNP(SNPSE)` — Custom SNP standard error. Not yet implemented.
- `multiGene(GeneSE, Genelist)` — Custom gene SE and gene subsetting. Not yet implemented.

### Convergence tolerance
- `toler` in usermodel, userGWAS, commonfactorGWAS — The L-BFGS optimizer uses its own convergence criteria. Custom tolerance is not yet exposed.
