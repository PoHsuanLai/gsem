use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use faer::Mat;

use gsem::io::gwas_reader;
use gsem::munge;

#[derive(Parser)]
#[command(name = "genomicsem", about = "Genomic Structural Equation Modeling")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// QC and munge raw GWAS summary statistics
    Munge {
        /// Input GWAS summary statistics files
        #[arg(short, long, num_args = 1..)]
        files: Vec<PathBuf>,

        /// HapMap3 reference SNP list
        #[arg(long)]
        hm3: PathBuf,

        /// Trait names
        #[arg(long, num_args = 1..)]
        trait_names: Option<Vec<String>>,

        /// INFO score filter threshold
        #[arg(long, default_value = "0.9")]
        info_filter: f64,

        /// MAF filter threshold
        #[arg(long, default_value = "0.01")]
        maf_filter: f64,

        /// Override sample size
        #[arg(short, long)]
        n: Option<f64>,

        /// Column name overrides (KEY=VALUE pairs, e.g., SNP=RSID,P=PVAL)
        #[arg(long)]
        column_names: Option<String>,

        /// Output directory
        #[arg(short, long, default_value = ".")]
        out: PathBuf,
    },

    /// Run multivariate LD Score Regression
    Ldsc {
        /// Munged summary statistics files
        #[arg(short, long, num_args = 1..)]
        traits: Vec<PathBuf>,

        /// Sample prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        sample_prev: Option<String>,

        /// Population prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        pop_prev: Option<String>,

        /// LD score directory
        #[arg(long)]
        ld: PathBuf,

        /// Weight LD score directory
        #[arg(long)]
        wld: Option<PathBuf>,

        /// Number of jackknife blocks
        #[arg(long, default_value = "200")]
        n_blocks: usize,

        /// Maximum chi-square for SNP filtering (default: auto = 0.001*max_N or 80)
        #[arg(long)]
        chisq_max: Option<f64>,

        /// Number of chromosomes (default: 22)
        #[arg(long, default_value = "22")]
        chr: usize,

        /// Chromosome selection: "all", "odd", "even", or comma-separated list (e.g., "1,3,5")
        #[arg(long, default_value = "all")]
        select: String,

        /// Return standardized genetic correlation matrix
        #[arg(long)]
        stand: bool,

        /// Output file (JSON)
        #[arg(short, long, default_value = "ldsc_result.json")]
        out: PathBuf,
    },

    /// Fit structural equation model
    Sem {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Model specification (lavaan syntax)
        #[arg(long)]
        model: Option<String>,

        /// Model specification file
        #[arg(long)]
        model_file: Option<PathBuf>,

        /// Estimation method (DWLS or ML)
        #[arg(long, default_value = "DWLS")]
        estimation: String,

        /// Standardize latent variables (fix first loading to free, set variance to 1)
        #[arg(long)]
        std_lv: bool,

        /// Auto-add lower bounds on residual variances and retry if model fails to converge
        #[arg(long)]
        fix_resid: bool,

        /// Output file
        #[arg(short, long, default_value = "sem_result.tsv")]
        out: PathBuf,
    },

    /// Merge GWAS summary statistics for multivariate GWAS
    Sumstats {
        /// Input GWAS summary statistics files
        #[arg(short, long, num_args = 1..)]
        files: Vec<PathBuf>,

        /// LD score reference directory (for SNP filtering)
        #[arg(long, name = "ref")]
        ref_dir: PathBuf,

        /// Trait names
        #[arg(long, num_args = 1..)]
        trait_names: Option<Vec<String>>,

        /// INFO score filter threshold
        #[arg(long, default_value = "0.6")]
        info_filter: f64,

        /// MAF filter threshold
        #[arg(long, default_value = "0.01")]
        maf_filter: f64,

        /// Keep indels (multi-character alleles)
        #[arg(long)]
        keep_indel: bool,

        /// Output file
        #[arg(short, long, default_value = "merged_sumstats.tsv")]
        out: PathBuf,
    },

    /// Fit common factor model (auto-generated 1-factor CFA)
    CommonFactor {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Estimation method (DWLS or ML)
        #[arg(long, default_value = "DWLS")]
        estimation: String,

        /// Output file
        #[arg(short, long, default_value = "commonfactor_result.tsv")]
        out: PathBuf,
    },

    /// Run multivariate GWAS
    Gwas {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Merged summary statistics file
        #[arg(long)]
        sumstats: PathBuf,

        /// Model specification
        #[arg(long)]
        model: Option<String>,

        /// Model specification file
        #[arg(long)]
        model_file: Option<PathBuf>,

        /// Estimation method
        #[arg(long, default_value = "DWLS")]
        estimation: String,

        /// Genomic control mode (conservative, standard, none)
        #[arg(long, default_value = "standard")]
        gc: String,

        /// Number of threads
        #[arg(long)]
        threads: Option<usize>,

        /// Standardize latent variables
        #[arg(long)]
        std_lv: bool,

        /// Filter output to specific parameters (comma-separated, e.g., "F1~SNP,F2~SNP")
        #[arg(long)]
        sub: Option<String>,

        /// Log warnings when covariance matrix requires smoothing
        #[arg(long)]
        smooth_check: bool,

        /// Compute Q_SNP heterogeneity statistic per SNP
        #[arg(long)]
        q_snp: bool,

        /// Fix measurement model parameters from baseline fit (estimate only SNP paths)
        #[arg(long)]
        fix_measurement: bool,

        /// Output file
        #[arg(short, long, default_value = "gwas_result.tsv")]
        out: PathBuf,
    },

    /// Run common factor GWAS (auto-generated 1-factor model per SNP)
    CommonfactorGwas {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Merged summary statistics file
        #[arg(long)]
        sumstats: PathBuf,

        /// Genomic control mode (conservative, standard, none)
        #[arg(long, default_value = "standard")]
        gc: String,

        /// Number of threads
        #[arg(long)]
        threads: Option<usize>,

        /// Output file
        #[arg(short, long, default_value = "commonfactor_gwas_result.tsv")]
        out: PathBuf,
    },

    /// Auto-generate model syntax from factor loadings
    WriteModel {
        /// TSV file of factor loadings (rows=phenotypes, cols=factors)
        #[arg(long)]
        loadings: PathBuf,

        /// Phenotype / variable names
        #[arg(long, num_args = 1..)]
        names: Vec<String>,

        /// Minimum absolute loading to include an indicator
        #[arg(long, default_value = "0.3")]
        cutoff: f64,

        /// Fix residual variances with positivity constraints
        #[arg(long)]
        fix_resid: bool,

        /// Generate bifactor model
        #[arg(long)]
        bifactor: bool,

        /// Output file
        #[arg(short, long, default_value = "model.txt")]
        out: PathBuf,
    },

    /// Parallel analysis to determine number of factors
    ParallelAnalysis {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Number of Monte Carlo simulations
        #[arg(long, default_value = "500")]
        n_sim: usize,

        /// Output file
        #[arg(short, long, default_value = "pa_result.tsv")]
        out: PathBuf,
    },

    /// Enrichment analysis using stratified LDSC results
    Enrich {
        /// JSON file with enrichment input data
        #[arg(long)]
        input: PathBuf,

        /// Output file
        #[arg(short, long, default_value = "enrich_result.tsv")]
        out: PathBuf,
    },

    /// Run stratified (partitioned) LD Score Regression
    SLdsc {
        /// Munged summary statistics files
        #[arg(short, long, num_args = 1..)]
        traits: Vec<PathBuf>,

        /// Sample prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        sample_prev: Option<String>,

        /// Population prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        pop_prev: Option<String>,

        /// Annotation LD score directory
        #[arg(long)]
        ld: PathBuf,

        /// Weight LD score directory
        #[arg(long)]
        wld: Option<PathBuf>,

        /// Number of jackknife blocks
        #[arg(long, default_value = "200")]
        n_blocks: usize,

        /// Number of chromosomes (default: 22)
        #[arg(long, default_value = "22")]
        chr: usize,

        /// Chromosome selection: "all", "odd", "even", or comma-separated list
        #[arg(long, default_value = "all")]
        select: String,

        /// Output file (JSON)
        #[arg(short, long, default_value = "s_ldsc_result.json")]
        out: PathBuf,
    },

    /// Simulate GWAS summary statistics
    Simulate {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Comma-separated per-trait sample sizes
        #[arg(long)]
        n_per_trait: String,

        /// LD score directory
        #[arg(long)]
        ld: PathBuf,

        /// Output file
        #[arg(short, long, default_value = "simulated_sumstats.tsv")]
        out: PathBuf,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Munge {
            files,
            hm3,
            trait_names,
            info_filter,
            maf_filter,
            n,
            column_names,
            out,
        } => run_munge(
            &files,
            &hm3,
            trait_names.as_deref(),
            info_filter,
            maf_filter,
            n,
            column_names,
            &out,
        ),
        Commands::Ldsc {
            traits,
            sample_prev,
            pop_prev,
            ld,
            wld,
            n_blocks,
            chisq_max,
            chr,
            select,
            stand,
            out,
        } => run_ldsc(
            &traits,
            sample_prev,
            pop_prev,
            &ld,
            wld,
            n_blocks,
            chisq_max,
            chr,
            &select,
            stand,
            &out,
        ),
        Commands::Sem {
            covstruc,
            model,
            model_file,
            estimation,
            std_lv,
            fix_resid,
            out,
        } => run_sem(&covstruc, model, model_file, &estimation, std_lv, fix_resid, &out),
        Commands::Sumstats {
            files,
            ref_dir,
            trait_names,
            info_filter,
            maf_filter,
            keep_indel,
            out,
        } => run_sumstats(&files, &ref_dir, trait_names, info_filter, maf_filter, keep_indel, &out),
        Commands::CommonFactor {
            covstruc,
            estimation,
            out,
        } => run_commonfactor_cmd(&covstruc, &estimation, &out),
        Commands::Gwas {
            covstruc,
            sumstats,
            model,
            model_file,
            estimation,
            gc,
            threads,
            std_lv,
            sub,
            smooth_check,
            q_snp,
            fix_measurement,
            out,
        } => run_gwas(
            &covstruc,
            &sumstats,
            model,
            model_file,
            &estimation,
            &gc,
            threads,
            std_lv,
            sub,
            smooth_check,
            q_snp,
            fix_measurement,
            &out,
        ),
        Commands::CommonfactorGwas {
            covstruc,
            sumstats,
            gc,
            threads,
            out,
        } => run_commonfactor_gwas(&covstruc, &sumstats, &gc, threads, &out),
        Commands::WriteModel {
            loadings,
            names,
            cutoff,
            fix_resid,
            bifactor,
            out,
        } => run_write_model(&loadings, &names, cutoff, fix_resid, bifactor, &out),
        Commands::ParallelAnalysis {
            covstruc,
            n_sim,
            out,
        } => run_parallel_analysis(&covstruc, n_sim, &out),
        Commands::SLdsc {
            traits,
            sample_prev,
            pop_prev,
            ld,
            wld,
            n_blocks,
            chr,
            select,
            out,
        } => run_s_ldsc(
            &traits,
            sample_prev,
            pop_prev,
            &ld,
            wld,
            n_blocks,
            chr,
            &select,
            &out,
        ),
        Commands::Enrich { input, out } => run_enrich(&input, &out),
        Commands::Simulate {
            covstruc,
            n_per_trait,
            ld,
            out,
        } => run_simulate(&covstruc, &n_per_trait, &ld, &out),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_munge(
    files: &[PathBuf],
    hm3: &Path,
    trait_names: Option<&[String]>,
    info_filter: f64,
    maf_filter: f64,
    n_override: Option<f64>,
    column_names: Option<String>,
    out_dir: &Path,
) -> Result<()> {
    // Read reference
    eprintln!("Reading reference panel: {}", hm3.display());
    let reference = munge::read_reference(hm3).context("failed to read HapMap3 reference")?;
    eprintln!("Loaded {} reference SNPs", reference.len());

    // Parse column name overrides from comma-separated KEY=VALUE pairs
    let column_overrides = column_names.map(|s| {
        s.split(',')
            .filter_map(|pair| {
                let mut parts = pair.splitn(2, '=');
                let key = parts.next()?.trim().to_string();
                let value = parts.next()?.trim().to_string();
                if key.is_empty() || value.is_empty() {
                    None
                } else {
                    Some((key, value))
                }
            })
            .collect::<HashMap<String, String>>()
    });

    let config = munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override,
        column_overrides,
    };

    for (i, file) in files.iter().enumerate() {
        let trait_name = trait_names
            .and_then(|names| names.get(i))
            .map(|s| s.as_str())
            .or_else(|| file.file_stem().and_then(|s| s.to_str()))
            .unwrap_or("trait");

        let out_path = out_dir.join(format!("{trait_name}.sumstats.gz"));
        eprintln!("Munging: {} -> {}", file.display(), out_path.display());

        munge::munge_and_write(file, &reference, &config, &out_path)
            .with_context(|| format!("failed to munge {}", file.display()))?;
    }

    eprintln!("Done.");
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_sumstats(
    files: &[PathBuf],
    ref_dir: &Path,
    trait_names: Option<Vec<String>>,
    info_filter: f64,
    maf_filter: f64,
    keep_indel: bool,
    out: &Path,
) -> Result<()> {
    let names: Vec<String> = trait_names.unwrap_or_else(|| {
        files
            .iter()
            .map(|f| {
                f.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("trait")
                    .to_string()
            })
            .collect()
    });

    let k = files.len();
    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        n_overrides: vec![None; k],
        se_logit: vec![false; k],
        ols: vec![false; k],
        linprob: vec![false; k],
        keep_indel,
    };

    let file_refs: Vec<&Path> = files.iter().map(|p| p.as_path()).collect();
    eprintln!("Merging {} GWAS files...", k);
    gsem::sumstats::merge_sumstats(&file_refs, ref_dir, &names, &config, out)?;
    eprintln!("Done.");
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_ldsc(
    traits: &[PathBuf],
    sample_prev: Option<String>,
    pop_prev: Option<String>,
    ld: &Path,
    wld: Option<PathBuf>,
    n_blocks: usize,
    chisq_max: Option<f64>,
    chr: usize,
    select: &str,
    stand: bool,
    out: &Path,
) -> Result<()> {
    eprintln!("Reading {} trait files...", traits.len());

    // Read summary stats
    let mut trait_data = Vec::new();
    for path in traits {
        let records = gwas_reader::read_sumstats(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        let n_snps = records.len();
        eprintln!("  {}: {} SNPs", path.display(), n_snps);

        trait_data.push(gsem_ldsc::TraitSumstats {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

    // Read LD scores
    let wld_dir = wld.as_deref().unwrap_or(ld);
    let chromosomes = parse_chromosome_selection(select, chr);
    let ld_data = gsem::io::ld_reader::read_ld_scores(ld, wld_dir, &chromosomes)
        .context("failed to read LD scores")?;

    eprintln!(
        "Loaded {} LD score SNPs, M={}",
        ld_data.records.len(),
        ld_data.total_m
    );

    // Parse prevalences
    let k = traits.len();
    let sp = parse_prevalences(&sample_prev, k);
    let pp = parse_prevalences(&pop_prev, k);

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max,
    };

    eprintln!("Running LDSC...");
    let result = gsem_ldsc::ldsc(
        &trait_data,
        &sp,
        &pp,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
    )?;

    // Write JSON output
    let json = result.to_json_string()?;
    std::fs::write(out, &json).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("LDSC complete. Results written to {}", out.display());

    // Print S matrix
    let k = result.s.nrows();
    eprintln!("\nGenetic covariance matrix (S):");
    for i in 0..k {
        let row: Vec<String> = (0..k)
            .map(|j| format!("{:8.4}", result.s[(i, j)]))
            .collect();
        eprintln!("  {}", row.join(" "));
    }

    eprintln!("\nIntercepts (I):");
    for i in 0..k {
        let row: Vec<String> = (0..k)
            .map(|j| format!("{:8.4}", result.i_mat[(i, j)]))
            .collect();
        eprintln!("  {}", row.join(" "));
    }

    if stand {
        eprintln!("\nStandardized genetic correlation matrix (S_Stand):");
        for i in 0..k {
            let row: Vec<String> = (0..k)
                .map(|j| {
                    let d_i = result.s[(i, i)].sqrt();
                    let d_j = result.s[(j, j)].sqrt();
                    if d_i > 0.0 && d_j > 0.0 {
                        format!("{:8.4}", result.s[(i, j)] / (d_i * d_j))
                    } else {
                        format!("{:8.4}", 0.0)
                    }
                })
                .collect();
            eprintln!("  {}", row.join(" "));
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_s_ldsc(
    traits: &[PathBuf],
    sample_prev: Option<String>,
    pop_prev: Option<String>,
    ld: &Path,
    wld: Option<PathBuf>,
    n_blocks: usize,
    chr: usize,
    select: &str,
    out: &Path,
) -> Result<()> {
    eprintln!("Reading {} trait files...", traits.len());

    // Read summary stats
    let mut trait_data = Vec::new();
    for path in traits {
        let records = gwas_reader::read_sumstats(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        let n_snps = records.len();
        eprintln!("  {}: {} SNPs", path.display(), n_snps);

        trait_data.push(gsem_ldsc::TraitSumstats {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

    // Read annotation LD scores
    let wld_dir = wld.as_deref().unwrap_or(ld);
    let chromosomes = parse_chromosome_selection(select, chr);
    let annot_data = gsem_ldsc::annot_reader::read_annot_ld_scores(ld, wld_dir, &chromosomes)
        .context("failed to read annotation LD scores")?;

    eprintln!(
        "Loaded {} LD score SNPs, {} annotations, M_total={}",
        annot_data.snps.len(),
        annot_data.annotation_names.len(),
        annot_data.m_annot.iter().sum::<f64>()
    );

    // Parse prevalences
    let k = traits.len();
    let sp = parse_prevalences(&sample_prev, k);
    let pp = parse_prevalences(&pop_prev, k);

    let config = gsem_ldsc::stratified::StratifiedLdscConfig { n_blocks };

    eprintln!("Running stratified LDSC...");
    let result = gsem_ldsc::stratified::s_ldsc(
        &trait_data,
        &sp,
        &pp,
        &annot_data.annot_ld,
        &annot_data.w_ld,
        &annot_data.snps,
        &annot_data.annotation_names,
        &annot_data.m_annot,
        &config,
    )?;

    // Write JSON output
    let json = result.to_json_string()?;
    std::fs::write(out, &json).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!(
        "Stratified LDSC complete. Results written to {}",
        out.display()
    );

    // Print summary
    let n_annot = result.annotations.len();
    eprintln!("\nAnnotations ({n_annot}):");
    for (a, name) in result.annotations.iter().enumerate() {
        let h2_diag: Vec<String> = (0..result.s_annot[a].nrows())
            .map(|i| format!("{:.6}", result.s_annot[a][(i, i)]))
            .collect();
        eprintln!(
            "  {name}: M={:.0}, prop={:.4}, h2_diag=[{}]",
            result.m_annot[a],
            result.prop[a],
            h2_diag.join(", ")
        );
    }

    eprintln!("\nIntercepts (I):");
    let ki = result.i_mat.nrows();
    for i in 0..ki {
        let row: Vec<String> = (0..ki)
            .map(|j| format!("{:8.4}", result.i_mat[(i, j)]))
            .collect();
        eprintln!("  {}", row.join(" "));
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_sem(
    covstruc: &Path,
    model: Option<String>,
    model_file: Option<PathBuf>,
    estimation: &str,
    std_lv: bool,
    fix_resid: bool,
    out: &Path,
) -> Result<()> {
    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    // Get model string
    let model_str = if let Some(m) = model {
        m
    } else if let Some(f) = model_file {
        std::fs::read_to_string(&f).with_context(|| format!("failed to read {}", f.display()))?
    } else {
        anyhow::bail!("must provide --model or --model-file");
    };

    eprintln!("Fitting SEM model (estimation={estimation}, std_lv={std_lv})...");

    // Parse and fit
    let mut pt =
        gsem_sem::syntax::parse_model(&model_str, std_lv).map_err(|e| anyhow::anyhow!("{e}"))?;
    let k = ldsc_result.s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    let mut model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut model, &ldsc_result.s, 1000)
    } else {
        gsem_sem::estimator::fit_dwls(&mut model, &ldsc_result.s, &v_diag, 1000)
    };

    // If model failed to converge and fix_resid is set, add lower bounds on
    // residual variances (op == "~~" && lhs == rhs) and retry.
    let fit = if !fit.converged && fix_resid {
        eprintln!("Model did not converge; retrying with residual variance lower bounds (fix_resid)...");
        for row in &mut pt.rows {
            if row.op == gsem_sem::syntax::Op::Covariance
                && row.lhs == row.rhs
                && row.free > 0
                && row.lower_bound.is_none()
            {
                row.lower_bound = Some(0.0001);
            }
        }
        // Rebuild model with updated lower bounds
        let mut model = gsem_sem::model::Model::from_partable(&pt, &obs_names);
        if estimation.to_uppercase() == "ML" {
            gsem_sem::estimator::fit_ml(&mut model, &ldsc_result.s, 1000)
        } else {
            gsem_sem::estimator::fit_dwls(&mut model, &ldsc_result.s, &v_diag, 1000)
        }
    } else {
        fit
    };

    eprintln!("Converged: {}", fit.converged);
    eprintln!("Objective: {:.6}", fit.objective);

    // Write results
    let mut output = String::from("lhs\top\trhs\test\n");
    for (i, row) in pt.rows.iter().enumerate() {
        if row.free > 0 {
            let est = fit.params.get(i).copied().unwrap_or(0.0);
            output.push_str(&format!(
                "{}\t{}\t{}\t{:.6}\n",
                row.lhs, row.op, row.rhs, est
            ));
        }
    }
    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Results written to {}", out.display());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_gwas(
    covstruc: &Path,
    sumstats: &Path,
    model: Option<String>,
    model_file: Option<PathBuf>,
    estimation: &str,
    gc: &str,
    threads: Option<usize>,
    std_lv: bool,
    sub: Option<String>,
    smooth_check: bool,
    _q_snp: bool,
    _fix_measurement: bool,
    out: &Path,
) -> Result<()> {
    // Set thread count
    if let Some(t) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .context("failed to initialize thread pool")?;
    }

    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    // Get model string
    let model_str = if let Some(m) = model {
        m
    } else if let Some(f) = model_file {
        std::fs::read_to_string(&f)
            .with_context(|| format!("failed to read model file {}", f.display()))?
    } else {
        anyhow::bail!("must provide --model or --model-file");
    };

    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().expect("infallible");
    let k = ldsc_result.s.nrows();

    // Read merged sumstats
    eprintln!("Reading merged sumstats: {}", sumstats.display());
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(sumstats)
        .with_context(|| format!("failed to read {}", sumstats.display()))?;

    let n_snps = merged.snps.len();
    eprintln!("GWAS: {n_snps} SNPs, {k} traits, estimation={estimation}, gc={gc}");

    if merged.trait_names.len() != k {
        anyhow::bail!(
            "trait count mismatch: LDSC has {k} traits, sumstats has {} (beta.* columns)",
            merged.trait_names.len()
        );
    }

    // Extract per-SNP arrays for user_gwas
    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged
        .snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    // Clamp intercept diagonals to >= 1 (matching R)
    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k {
        if i_ld[(i, i)] < 1.0 {
            i_ld[(i, i)] = 1.0;
        }
    }

    // TODO: --fix-measurement: fit baseline model without SNP, then fix non-SNP params
    // per SNP. This requires pre-fitting, identifying measurement params, and constraining
    // them in the per-SNP loop. Not yet implemented.

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: model_str,
        estimation: estimation.to_string(),
        gc: gc_mode,
        max_iter: 500,
        std_lv,
        smooth_check,
        snp_se: None,
    };

    eprintln!("Running GWAS across {n_snps} SNPs...");
    let results = gsem::gwas::user_gwas::run_user_gwas(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
    );

    // Parse --sub filter list
    let sub_filters: Option<Vec<String>> = sub.map(|s| {
        s.split(',')
            .map(|entry| entry.trim().to_string())
            .filter(|entry| !entry.is_empty())
            .collect()
    });

    // Write TSV output
    let mut output = String::from("SNP\tlhs\top\trhs\test\tse\tz\tp\tchisq\tdf\tconverged\n");
    for snp_result in &results {
        let snp_name = &merged.snps[snp_result.snp_idx].snp;
        for param in &snp_result.params {
            // Apply --sub filter: only include params matching the filter list
            if let Some(ref filters) = sub_filters {
                let param_key = format!("{}{}{}", param.lhs, param.op, param.rhs);
                if !filters.iter().any(|f| f == &param_key) {
                    continue;
                }
            }

            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\t{:.4}\t{}\t{}\n",
                snp_name,
                param.lhs,
                param.op,
                param.rhs,
                param.est,
                param.se,
                param.z_stat,
                param.p_value,
                snp_result.chisq,
                snp_result.chisq_df,
                snp_result.converged,
            ));
        }
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    let n_converged = results.iter().filter(|r| r.converged).count();
    eprintln!(
        "GWAS complete: {n_snps} SNPs, {n_converged} converged. Results: {}",
        out.display()
    );
    Ok(())
}

fn run_commonfactor_cmd(covstruc: &Path, estimation: &str, out: &Path) -> Result<()> {
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    eprintln!("Fitting common factor model (estimation={estimation})...");

    let result = gsem_sem::commonfactor::run_commonfactor(
        &ldsc_result.s,
        &ldsc_result.v,
        estimation,
    )?;

    eprintln!(
        "Converged. Chi-sq={:.4}, df={}, CFI={:.4}, SRMR={:.4}",
        result.fit.chisq, result.fit.df, result.fit.cfi, result.fit.srmr
    );

    // Write TSV output
    let mut output = String::from("lhs\top\trhs\test\tse\tz\tp\n");
    for param in &result.parameters {
        output.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\n",
            param.lhs, param.op, param.rhs, param.est, param.se, param.z, param.p
        ));
    }

    // Append fit indices
    output.push_str(&format!(
        "\n# Model fit: chisq={:.4} df={} p={:.6e} AIC={:.4} CFI={:.4} SRMR={:.4}\n",
        result.fit.chisq, result.fit.df, result.fit.p_chisq,
        result.fit.aic, result.fit.cfi, result.fit.srmr
    ));

    std::fs::write(out, &output)
        .with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Results written to {}", out.display());
    Ok(())
}

fn run_commonfactor_gwas(
    covstruc: &Path,
    sumstats: &Path,
    gc: &str,
    threads: Option<usize>,
    out: &Path,
) -> Result<()> {
    // Set thread count
    if let Some(t) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .context("failed to initialize thread pool")?;
    }

    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().expect("infallible");
    let k = ldsc_result.s.nrows();

    // Read merged sumstats
    eprintln!("Reading merged sumstats: {}", sumstats.display());
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(sumstats)
        .with_context(|| format!("failed to read {}", sumstats.display()))?;

    let n_snps = merged.snps.len();
    eprintln!("Common factor GWAS: {n_snps} SNPs, {k} traits, gc={gc}");

    if merged.trait_names.len() != k {
        anyhow::bail!(
            "trait count mismatch: LDSC has {k} traits, sumstats has {} (beta.* columns)",
            merged.trait_names.len()
        );
    }

    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged
        .snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    // Clamp intercept diagonals to >= 1 (matching R)
    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k {
        if i_ld[(i, i)] < 1.0 {
            i_ld[(i, i)] = 1.0;
        }
    }

    eprintln!("Running common factor GWAS across {n_snps} SNPs...");
    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
        gc_mode,
    );

    // Write TSV output (same format as user GWAS)
    let mut output = String::from("SNP\tlhs\top\trhs\test\tse\tz\tp\tchisq\tdf\tconverged\n");
    for snp_result in &results {
        let snp_name = &merged.snps[snp_result.snp_idx].snp;
        for param in &snp_result.params {
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\t{:.4}\t{}\t{}\n",
                snp_name,
                param.lhs,
                param.op,
                param.rhs,
                param.est,
                param.se,
                param.z_stat,
                param.p_value,
                snp_result.chisq,
                snp_result.chisq_df,
                snp_result.converged,
            ));
        }
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    let n_converged = results.iter().filter(|r| r.converged).count();
    eprintln!(
        "Common factor GWAS complete: {n_snps} SNPs, {n_converged} converged. Results: {}",
        out.display()
    );
    Ok(())
}

fn run_write_model(
    loadings_path: &Path,
    names: &[String],
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
    out: &Path,
) -> Result<()> {
    // Read loadings TSV: rows are phenotypes, columns are factors, values are floats
    let content = std::fs::read_to_string(loadings_path)
        .with_context(|| format!("failed to read {}", loadings_path.display()))?;

    let mut rows: Vec<Vec<f64>> = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let vals: Vec<f64> = trimmed
            .split(['\t', ' '])
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<f64>()
                    .with_context(|| format!("failed to parse float: {s:?}"))
            })
            .collect::<Result<Vec<_>>>()?;
        rows.push(vals);
    }

    if rows.is_empty() {
        anyhow::bail!("loadings file is empty");
    }

    let n_rows = rows.len();
    let n_cols = rows[0].len();
    let loadings = Mat::from_fn(n_rows, n_cols, |i, j| rows[i][j]);

    eprintln!(
        "Generating model from {}x{} loadings matrix (cutoff={cutoff}, fix_resid={fix_resid}, bifactor={bifactor})",
        n_rows, n_cols
    );

    let model_str =
        gsem_sem::write_model::write_model(&loadings, names, cutoff, fix_resid, bifactor);

    std::fs::write(out, &model_str)
        .with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Model written to {}", out.display());
    Ok(())
}

fn run_parallel_analysis(covstruc: &Path, n_sim: usize, out: &Path) -> Result<()> {
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    let k = ldsc_result.s.nrows();
    eprintln!("Running parallel analysis ({k} traits, {n_sim} simulations)...");

    let result =
        gsem::stats::parallel_analysis::parallel_analysis(&ldsc_result.s, &ldsc_result.v, n_sim);

    eprintln!("Suggested number of factors: {}", result.n_factors);

    // Write TSV output
    let mut output = String::from("factor\tobserved\tsimulated_95\n");
    for (i, (obs, sim)) in result
        .observed
        .iter()
        .zip(result.simulated_95.iter())
        .enumerate()
    {
        output.push_str(&format!("{}\t{:.6}\t{:.6}\n", i + 1, obs, sim));
    }
    output.push_str(&format!("\n# n_factors={}\n", result.n_factors));

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Results written to {}", out.display());
    Ok(())
}

fn run_enrich(input: &Path, out: &Path) -> Result<()> {
    // Read JSON input containing all enrichment data
    let json = std::fs::read_to_string(input)
        .with_context(|| format!("failed to read {}", input.display()))?;

    #[derive(serde::Deserialize)]
    struct EnrichInput {
        s_baseline: Vec<Vec<f64>>,
        s_annot: Vec<Vec<Vec<f64>>>,
        v_annot: Vec<Vec<Vec<f64>>>,
        annotation_names: Vec<String>,
        m_annot: Vec<f64>,
        m_total: f64,
    }

    let data: EnrichInput =
        serde_json::from_str(&json).context("failed to parse enrichment input JSON")?;

    // Convert nested vecs to faer::Mat
    let k = data.s_baseline.len();
    let s_baseline = Mat::from_fn(k, k, |i, j| data.s_baseline[i][j]);

    let s_annot: Vec<Mat<f64>> = data
        .s_annot
        .iter()
        .map(|m| {
            let rows = m.len();
            let cols = if rows > 0 { m[0].len() } else { 0 };
            Mat::from_fn(rows, cols, |i, j| m[i][j])
        })
        .collect();

    let v_annot: Vec<Mat<f64>> = data
        .v_annot
        .iter()
        .map(|m| {
            let rows = m.len();
            let cols = if rows > 0 { m[0].len() } else { 0 };
            Mat::from_fn(rows, cols, |i, j| m[i][j])
        })
        .collect();

    eprintln!(
        "Running enrichment analysis ({} annotations, {k} traits)...",
        data.annotation_names.len()
    );

    let result = gsem::stats::enrich::enrichment_test(
        &s_baseline,
        &s_annot,
        &v_annot,
        &data.annotation_names,
        &data.m_annot,
        data.m_total,
    );

    // Write TSV output
    let mut output = String::from("annotation\tenrichment\tse\tp\n");
    for i in 0..result.annotations.len() {
        output.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{:.6e}\n",
            result.annotations[i], result.enrichment[i], result.se[i], result.p[i]
        ));
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Enrichment results written to {}", out.display());
    Ok(())
}

fn run_simulate(covstruc: &Path, n_per_trait: &str, ld_dir: &Path, out: &Path) -> Result<()> {
    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    let k = ldsc_result.s.nrows();

    // Parse per-trait sample sizes
    let n_vec: Vec<f64> = n_per_trait
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<f64>()
                .with_context(|| format!("failed to parse sample size: {s:?}"))
        })
        .collect::<Result<Vec<_>>>()?;

    if n_vec.len() != k {
        anyhow::bail!(
            "expected {k} sample sizes (one per trait), got {}",
            n_vec.len()
        );
    }

    // Read LD scores (all 22 chromosomes by default)
    let chromosomes: Vec<usize> = (1..=22).collect();
    let ld_data = gsem::io::ld_reader::read_ld_scores(ld_dir, ld_dir, &chromosomes)
        .context("failed to read LD scores")?;

    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();
    let n_snps = ld_scores.len();

    eprintln!(
        "Simulating GWAS: {k} traits, {n_snps} SNPs, M={:.0}",
        ld_data.total_m
    );

    let z_all = gsem::stats::simulation::simulate_sumstats(
        &ldsc_result.s,
        &n_vec,
        &ld_scores,
        ld_data.total_m,
    );

    // Write TSV: columns are trait z-scores
    let mut output = String::new();
    // Header
    let headers: Vec<String> = (0..k).map(|i| format!("Z{}", i + 1)).collect();
    output.push_str(&headers.join("\t"));
    output.push('\n');
    // Rows (one per SNP)
    for s_idx in 0..n_snps {
        let vals: Vec<String> = (0..k).map(|t| format!("{:.6}", z_all[t][s_idx])).collect();
        output.push_str(&vals.join("\t"));
        output.push('\n');
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!(
        "Simulated {n_snps} SNPs x {k} traits. Results: {}",
        out.display()
    );
    Ok(())
}

fn parse_prevalences(s: &Option<String>, k: usize) -> Vec<Option<f64>> {
    match s {
        Some(s) => s
            .split(',')
            .map(|x| {
                let trimmed = x.trim();
                if trimmed.eq_ignore_ascii_case("na") || trimmed.is_empty() {
                    None
                } else {
                    trimmed.parse().ok()
                }
            })
            .collect(),
        None => vec![None; k],
    }
}

fn parse_chromosome_selection(select: &str, n_chr: usize) -> Vec<usize> {
    match select.to_lowercase().as_str() {
        "all" => (1..=n_chr).collect(),
        "odd" => (1..=n_chr).filter(|c| c % 2 == 1).collect(),
        "even" => (1..=n_chr).filter(|c| c % 2 == 0).collect(),
        other => {
            // Parse comma-separated list
            other
                .split(',')
                .filter_map(|s| s.trim().parse::<usize>().ok())
                .collect()
        }
    }
}
