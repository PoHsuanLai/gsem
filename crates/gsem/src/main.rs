use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use faer::Mat;

use gsem::io::gwas_reader;
use gsem::munge;

#[derive(Parser)]
#[command(name = "gsem", about = "Genomic Structural Equation Modeling")]
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
    #[command(name = "usermodel")]
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

        /// Compute Q_Factor heterogeneity test for cross-factor indicator pairs
        #[arg(long)]
        q_factor: bool,

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
    #[command(name = "userGWAS")]
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

        /// Run in TWAS mode (input uses Gene/Panel/HSQ instead of SNP/A1/A2/MAF)
        #[arg(long)]
        twas: bool,

        /// Output file
        #[arg(short, long, default_value = "gwas_result.tsv")]
        out: PathBuf,
    },

    /// Run common factor GWAS (auto-generated 1-factor model per SNP)
    #[command(name = "commonfactorGWAS")]
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
    #[command(name = "write.model")]
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
    #[command(name = "paLDSC")]
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

    /// Compute model-implied genetic correlation matrix
    Rgmodel {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Estimation method (DWLS or ML)
        #[arg(long, default_value = "DWLS")]
        estimation: String,

        /// User-specified model (lavaan syntax). If omitted, fits common factor.
        #[arg(long)]
        model: Option<String>,

        /// Standardize latent variables
        #[arg(long, default_value_t = false)]
        std_lv: bool,

        /// Output file
        #[arg(short, long, default_value = "rgmodel_result.tsv")]
        out: PathBuf,
    },

    /// Joint analysis of multiple SNPs with LD
    #[command(name = "multiSNP")]
    MultiSnp {
        /// LDSC result JSON file
        #[arg(long)]
        covstruc: PathBuf,

        /// Tab-delimited file with columns: SNP, A1, A2, MAF, beta.T1, se.T1, ...
        #[arg(long)]
        sumstats: PathBuf,

        /// LD matrix file (tab-delimited, symmetric, SNP x SNP correlations)
        #[arg(long)]
        ld_matrix: PathBuf,

        /// Model specification
        #[arg(long)]
        model: Option<String>,

        /// Model specification file
        #[arg(long)]
        model_file: Option<PathBuf>,

        /// Estimation method
        #[arg(long, default_value = "DWLS")]
        estimation: String,

        /// Output file
        #[arg(short, long, default_value = "multi_snp_result.tsv")]
        out: PathBuf,
    },

    /// Simulate GWAS summary statistics
    #[command(name = "simLDSC")]
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

    /// High-Definition Likelihood estimation of genetic covariance
    Hdl {
        /// Munged summary statistics files
        #[arg(short, long, num_args = 1..)]
        traits: Vec<PathBuf>,

        /// Sample prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        sample_prev: Option<String>,

        /// Population prevalences (comma-separated, NA for continuous)
        #[arg(long)]
        pop_prev: Option<String>,

        /// LD reference panel directory (text format)
        #[arg(long)]
        ld_path: PathBuf,

        /// Reference panel sample size (default: 335265 for UKB)
        #[arg(long, default_value = "335265")]
        n_ref: f64,

        /// Method: piecewise or jackknife
        #[arg(long, default_value = "piecewise")]
        method: String,

        /// Output file (JSON)
        #[arg(short, long, default_value = "hdl_result.json")]
        out: PathBuf,
    },

    /// Generalized Least Squares regression on genetic parameters
    #[command(name = "summaryGLS")]
    SummaryGls {
        /// TSV file with predictor matrix X (rows=observations, cols=predictors)
        #[arg(long)]
        x: PathBuf,

        /// TSV file with outcome vector Y (one value per line)
        #[arg(long)]
        y: PathBuf,

        /// TSV file with covariance matrix V (square, same rows as Y)
        #[arg(long)]
        v: PathBuf,

        /// Add intercept column (default true)
        #[arg(long, default_value = "true")]
        intercept: bool,

        /// Output file
        #[arg(short, long, default_value = "gls_result.tsv")]
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
            q_factor,
            out,
        } => run_sem(
            &covstruc,
            model,
            model_file,
            &estimation,
            std_lv,
            fix_resid,
            q_factor,
            &out,
        ),
        Commands::Sumstats {
            files,
            ref_dir,
            trait_names,
            info_filter,
            maf_filter,
            keep_indel,
            out,
        } => run_sumstats(
            &files,
            &ref_dir,
            trait_names,
            info_filter,
            maf_filter,
            keep_indel,
            &out,
        ),
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
            twas,
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
            twas,
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
        Commands::Rgmodel {
            covstruc,
            estimation,
            model,
            std_lv,
            out,
        } => run_rgmodel_cmd(&covstruc, &estimation, model.as_deref(), std_lv, &out),
        Commands::MultiSnp {
            covstruc,
            sumstats,
            ld_matrix,
            model,
            model_file,
            estimation,
            out,
        } => run_multi_snp_cmd(
            &covstruc,
            &sumstats,
            &ld_matrix,
            model,
            model_file,
            &estimation,
            &out,
        ),
        Commands::Simulate {
            covstruc,
            n_per_trait,
            ld,
            out,
        } => run_simulate(&covstruc, &n_per_trait, &ld, &out),
        Commands::Hdl {
            traits,
            sample_prev,
            pop_prev,
            ld_path,
            n_ref,
            method,
            out,
        } => run_hdl(
            &traits,
            sample_prev,
            pop_prev,
            &ld_path,
            n_ref,
            &method,
            &out,
        ),
        Commands::SummaryGls {
            x,
            y,
            v,
            intercept,
            out,
        } => run_summary_gls(&x, &y, &v, intercept, &out),
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
        keep_ambig: false,
        beta_overrides: Vec::new(),
        direct_filter: false,
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
    q_factor: bool,
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
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, None)
    } else {
        gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, None)
    };

    // If model failed to converge and fix_resid is set, add lower bounds on
    // residual variances (op == "~~" && lhs == rhs) and retry.
    let fit = if !fit.converged && fix_resid {
        eprintln!(
            "Model did not converge; retrying with residual variance lower bounds (fix_resid)..."
        );
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
        sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);
        if estimation.to_uppercase() == "ML" {
            gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, None)
        } else {
            gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, None)
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

    // Q_Factor heterogeneity test
    if q_factor {
        let sigma_hat = sem_model.implied_cov();
        let fi = gsem_sem::q_factor::factor_indicators(&pt, &obs_names);
        if fi.len() >= 2 {
            let q_results = gsem_sem::q_factor::compute_q_factor(
                &ldsc_result.s,
                &sigma_hat,
                &ldsc_result.v,
                &fi,
            );
            if !q_results.is_empty() {
                output.push_str("\n# Q_Factor heterogeneity test\n");
                output.push_str("factor1\tfactor2\tQ_chisq\tQ_df\tQ_pval\n");
                for qr in &q_results {
                    output.push_str(&format!(
                        "{}\t{}\t{:.6}\t{}\t{:.6e}\n",
                        qr.factor1, qr.factor2, qr.q_chisq, qr.q_df, qr.q_p
                    ));
                    eprintln!(
                        "Q_Factor({}, {}): chisq={:.4}, df={}, p={:.4e}",
                        qr.factor1, qr.factor2, qr.q_chisq, qr.q_df, qr.q_p
                    );
                }
            }
        } else {
            eprintln!("Q_Factor: skipped (fewer than 2 factors in model)");
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
    twas: bool,
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

    // Clamp intercept diagonals to >= 1 (matching R)
    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k {
        if i_ld[(i, i)] < 1.0 {
            i_ld[(i, i)] = 1.0;
        }
    }

    // Parse --sub filter list
    let sub_filters: Option<Vec<String>> = sub.map(|s| {
        s.split(',')
            .map(|entry| entry.trim().to_string())
            .filter(|entry| !entry.is_empty())
            .collect()
    });

    if twas {
        run_gwas_twas(
            sumstats,
            &model_str,
            estimation,
            gc,
            gc_mode,
            k,
            std_lv,
            smooth_check,
            &ldsc_result,
            &i_ld,
            &sub_filters,
            out,
        )
    } else {
        run_gwas_snp(
            sumstats,
            &model_str,
            estimation,
            gc,
            gc_mode,
            k,
            std_lv,
            smooth_check,
            &ldsc_result,
            &i_ld,
            &sub_filters,
            out,
        )
    }
}

/// Standard SNP-based GWAS path.
#[allow(clippy::too_many_arguments)]
fn run_gwas_snp(
    sumstats: &Path,
    model_str: &str,
    estimation: &str,
    gc: &str,
    gc_mode: gsem::gwas::gc_correction::GcMode,
    k: usize,
    std_lv: bool,
    smooth_check: bool,
    ldsc_result: &gsem_ldsc::LdscResult,
    i_ld: &Mat<f64>,
    sub_filters: &Option<Vec<String>>,
    out: &Path,
) -> Result<()> {
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

    // fix_measurement is handled inside run_user_gwas when config.fix_measurement is true

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: model_str.to_string(),
        estimation: estimation.to_string(),
        gc: gc_mode,
        max_iter: 500,
        std_lv,
        smooth_check,
        snp_se: None,
        snp_label: "SNP".to_string(),
        q_snp: false,
        fix_measurement: false,
    };

    eprintln!("Running GWAS across {n_snps} SNPs...");
    let results = gsem::gwas::user_gwas::run_user_gwas(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
    );

    // Write TSV output
    let mut output = String::from("SNP\tlhs\top\trhs\test\tse\tz\tp\tchisq\tdf\tconverged\n");
    for snp_result in &results {
        let snp_name = &merged.snps[snp_result.snp_idx].snp;
        for param in &snp_result.params {
            // Apply --sub filter: only include params matching the filter list
            if let Some(filters) = sub_filters {
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

/// TWAS mode: uses Gene/Panel/HSQ instead of SNP/A1/A2/MAF.
///
/// Key differences from standard GWAS:
/// - Reads TWAS-format sumstats (Gene/Panel/HSQ columns)
/// - Uses HSQ (heritability of expression) as var_snp instead of 2*MAF*(1-MAF)
/// - Replaces "SNP" with "Gene" in model syntax
/// - Outputs Gene/Panel/HSQ columns instead of SNP column
#[allow(clippy::too_many_arguments)]
fn run_gwas_twas(
    sumstats: &Path,
    model_str: &str,
    estimation: &str,
    gc: &str,
    gc_mode: gsem::gwas::gc_correction::GcMode,
    k: usize,
    std_lv: bool,
    smooth_check: bool,
    ldsc_result: &gsem_ldsc::LdscResult,
    i_ld: &Mat<f64>,
    sub_filters: &Option<Vec<String>>,
    out: &Path,
) -> Result<()> {
    // Read TWAS sumstats
    eprintln!("Reading TWAS sumstats: {}", sumstats.display());
    let twas_data = gsem::io::twas_reader::read_twas_sumstats(sumstats)
        .with_context(|| format!("failed to read TWAS sumstats {}", sumstats.display()))?;

    let n_genes = twas_data.genes.len();
    eprintln!("TWAS: {n_genes} genes, {k} traits, estimation={estimation}, gc={gc}");

    if twas_data.trait_names.len() != k {
        anyhow::bail!(
            "trait count mismatch: LDSC has {k} traits, TWAS sumstats has {} (beta.* columns)",
            twas_data.trait_names.len()
        );
    }

    // Extract per-gene arrays
    let beta_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.beta.clone()).collect();
    let se_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.se.clone()).collect();
    // In TWAS mode, var_snp = HSQ (heritability of expression)
    let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

    // Replace "SNP" with "Gene" in model syntax so that the observed variable
    // name matches what userGWAS expects in the first position
    let twas_model = model_str.replace("SNP", "Gene");

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: twas_model,
        estimation: estimation.to_string(),
        gc: gc_mode,
        max_iter: 500,
        std_lv,
        smooth_check,
        snp_se: None,
        snp_label: "Gene".to_string(),
        q_snp: false,
        fix_measurement: false,
    };

    eprintln!("Running TWAS across {n_genes} genes...");
    let results = gsem::gwas::user_gwas::run_user_gwas(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        i_ld,
        &beta_gene,
        &se_gene,
        &var_gene,
    );

    // Write TSV output with TWAS columns
    let mut output =
        String::from("Gene\tPanel\tHSQ\tlhs\top\trhs\test\tse\tz\tp\tchisq\tdf\tconverged\n");
    for gene_result in &results {
        let gene = &twas_data.genes[gene_result.snp_idx];
        for param in &gene_result.params {
            // Apply --sub filter: only include params matching the filter list
            if let Some(filters) = sub_filters {
                let param_key = format!("{}{}{}", param.lhs, param.op, param.rhs);
                if !filters.iter().any(|f| f == &param_key) {
                    continue;
                }
            }

            output.push_str(&format!(
                "{}\t{}\t{:.6}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\t{:.4}\t{}\t{}\n",
                gene.gene,
                gene.panel,
                gene.hsq,
                param.lhs,
                param.op,
                param.rhs,
                param.est,
                param.se,
                param.z_stat,
                param.p_value,
                gene_result.chisq,
                gene_result.chisq_df,
                gene_result.converged,
            ));
        }
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    let n_converged = results.iter().filter(|r| r.converged).count();
    eprintln!(
        "TWAS complete: {n_genes} genes, {n_converged} converged. Results: {}",
        out.display()
    );
    Ok(())
}

fn run_commonfactor_cmd(covstruc: &Path, estimation: &str, out: &Path) -> Result<()> {
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    eprintln!("Fitting common factor model (estimation={estimation})...");

    let result =
        gsem_sem::commonfactor::run_commonfactor(&ldsc_result.s, &ldsc_result.v, estimation)?;

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
        result.fit.chisq,
        result.fit.df,
        result.fit.p_chisq,
        result.fit.aic,
        result.fit.cfi,
        result.fit.srmr
    ));

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

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
    let cf_config = gsem::gwas::common_factor::CommonFactorGwasConfig {
        gc: gc_mode,
        ..Default::default()
    };
    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
        &cf_config,
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
        gsem_sem::write_model::write_model(&loadings, names, cutoff, fix_resid, bifactor, false, false);

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
        gsem::stats::parallel_analysis::parallel_analysis(&ldsc_result.s, &ldsc_result.v, n_sim, 0.95, false);

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

fn run_rgmodel_cmd(
    covstruc: &Path,
    estimation: &str,
    model: Option<&str>,
    std_lv: bool,
    out: &Path,
) -> Result<()> {
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    eprintln!("Fitting rgmodel (estimation={estimation})...");

    let result = gsem_sem::rgmodel::run_rgmodel_with_model(
        &ldsc_result.s, &ldsc_result.v, estimation, model, std_lv,
    )?;

    let k = result.r.nrows();
    eprintln!(
        "Converged. R is {k}x{k}. Underlying model chi-sq={:.4}, df={}, CFI={:.4}",
        result.sem_result.fit.chisq, result.sem_result.fit.df, result.sem_result.fit.cfi
    );

    // Write output: R matrix, then V_R matrix, then SEM parameters
    let mut output = String::new();
    output.push_str("# Model-implied genetic correlation matrix R\n");
    let col_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    output.push_str(&format!("\t{}\n", col_names.join("\t")));
    for i in 0..k {
        output.push_str(&format!("V{}", i + 1));
        for j in 0..k {
            output.push_str(&format!("\t{:.6}", result.r[(i, j)]));
        }
        output.push('\n');
    }

    output.push_str("\n# Sampling covariance of vech(R) - V_R\n");
    let kstar = k * (k + 1) / 2;
    for i in 0..kstar {
        for j in 0..kstar {
            if j > 0 {
                output.push('\t');
            }
            output.push_str(&format!("{:.6e}", result.v_r[(i, j)]));
        }
        output.push('\n');
    }

    output.push_str("\n# SEM parameter estimates\n");
    output.push_str("lhs\top\trhs\test\tse\tz\tp\n");
    for param in &result.sem_result.parameters {
        output.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\n",
            param.lhs, param.op, param.rhs, param.est, param.se, param.z, param.p
        ));
    }

    output.push_str(&format!(
        "\n# Model fit: chisq={:.4} df={} p={:.6e} AIC={:.4} CFI={:.4} SRMR={:.4}\n",
        result.sem_result.fit.chisq,
        result.sem_result.fit.df,
        result.sem_result.fit.p_chisq,
        result.sem_result.fit.aic,
        result.sem_result.fit.cfi,
        result.sem_result.fit.srmr
    ));

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Results written to {}", out.display());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_multi_snp_cmd(
    covstruc: &Path,
    sumstats: &Path,
    ld_matrix_path: &Path,
    model: Option<String>,
    model_file: Option<PathBuf>,
    estimation: &str,
    out: &Path,
) -> Result<()> {
    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;
    let k = ldsc_result.s.nrows();

    // Get model string
    let model_str = if let Some(m) = model {
        m
    } else if let Some(f) = model_file {
        std::fs::read_to_string(&f)
            .with_context(|| format!("failed to read model file {}", f.display()))?
    } else {
        anyhow::bail!("must provide --model or --model-file");
    };

    // Read merged sumstats
    eprintln!("Reading merged sumstats: {}", sumstats.display());
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(sumstats)
        .with_context(|| format!("failed to read {}", sumstats.display()))?;

    if merged.trait_names.len() != k {
        anyhow::bail!(
            "trait count mismatch: LDSC has {k} traits, sumstats has {} (beta.* columns)",
            merged.trait_names.len()
        );
    }

    // Read LD matrix
    eprintln!("Reading LD matrix: {}", ld_matrix_path.display());
    let (ld_mat, _ld_names) = gsem::gwas::multi_snp::read_ld_matrix(ld_matrix_path)?;
    let n_snps = ld_mat.nrows();
    eprintln!("LD matrix: {n_snps} x {n_snps}");

    if n_snps > merged.snps.len() {
        anyhow::bail!(
            "LD matrix has {n_snps} SNPs but sumstats only has {} SNPs",
            merged.snps.len()
        );
    }

    // Use the first n_snps from the sumstats
    let beta_snp: Vec<Vec<f64>> = merged.snps[..n_snps]
        .iter()
        .map(|s| s.beta.clone())
        .collect();
    let se_snp: Vec<Vec<f64>> = merged.snps[..n_snps].iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged.snps[..n_snps]
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();
    let snp_names: Vec<String> = merged.snps[..n_snps]
        .iter()
        .map(|s| s.snp.clone())
        .collect();

    let config = gsem::gwas::multi_snp::MultiSnpConfig {
        model: model_str,
        estimation: estimation.to_string(),
        max_iter: 500,
        snp_var_se: None,
    };

    eprintln!("Running multi-SNP analysis with {n_snps} SNPs, {k} traits...");
    let result = gsem::gwas::multi_snp::run_multi_snp(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &beta_snp,
        &se_snp,
        &var_snp,
        &ld_mat,
        &snp_names,
    );

    eprintln!(
        "Done. converged={}, chisq={:.4}, df={}",
        result.converged, result.chisq, result.chisq_df
    );

    // Write TSV output
    let mut output = String::from("lhs\top\trhs\test\tse\tz\tp\n");
    for param in &result.params {
        output.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\n",
            param.lhs, param.op, param.rhs, param.est, param.se, param.z_stat, param.p_value
        ));
    }
    output.push_str(&format!(
        "\n# chisq={:.4} df={} converged={}\n",
        result.chisq, result.chisq_df, result.converged
    ));

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Results written to {}", out.display());
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
        &gsem::stats::simulation::SimConfig::default(),
    );

    // Write TSV: columns are trait z-scores
    let mut output = String::new();
    // Header
    let headers: Vec<String> = (0..k).map(|i| format!("Z{}", i + 1)).collect();
    output.push_str(&headers.join("\t"));
    output.push('\n');
    // Rows (one per SNP)
    for s_idx in 0..n_snps {
        let vals: Vec<String> = z_all
            .iter()
            .map(|z_t| format!("{:.6}", z_t[s_idx]))
            .collect();
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

#[allow(clippy::too_many_arguments)]
fn run_hdl(
    traits: &[PathBuf],
    sample_prev: Option<String>,
    pop_prev: Option<String>,
    ld_path: &Path,
    n_ref: f64,
    method: &str,
    out: &Path,
) -> Result<()> {
    use gsem_ldsc::hdl::{HdlConfig, HdlMethod, HdlTraitData, LdPiece};

    eprintln!("Reading {} trait files for HDL...", traits.len());

    // Read summary stats (same pattern as LDSC)
    let mut trait_data = Vec::new();
    for path in traits {
        let records = gwas_reader::read_sumstats(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        let n_snps = records.len();
        eprintln!("  {}: {} SNPs", path.display(), n_snps);

        trait_data.push(HdlTraitData {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

    // Parse HDL method
    let hdl_method = match method.to_lowercase().as_str() {
        "jackknife" => HdlMethod::Jackknife,
        _ => HdlMethod::Piecewise,
    };

    let config = HdlConfig {
        method: hdl_method,
        n_ref,
    };

    // Load LD pieces from text format directory.
    // Expected format:
    //   {ld_path}/pieces.txt — tab-delimited: piece_idx, chr, n_snps
    //   {ld_path}/piece.{N}.snps.txt — per-piece SNP info: SNP\tA1\tA2\tLDscore
    let pieces_file = ld_path.join("pieces.txt");
    if !pieces_file.exists() {
        anyhow::bail!(
            "HDL LD reference directory not found or missing pieces.txt at {}.\n\
             HDL requires a text-format LD reference panel directory containing:\n\
             - pieces.txt: tab-delimited index file with columns: piece_idx, chr, n_snps\n\
             - piece.{{N}}.snps.txt: per-piece SNP data with columns: SNP, A1, A2, LDscore\n\n\
             The R version of GenomicSEM uses .rda files which cannot be read directly.\n\
             Convert R LD reference panels to text format first, e.g. in R:\n\
             \n  load('UKB_SVD_eigen99_extraction/LDmatrix1.RData')\n\
             \n  write.table(data.frame(SNP=snps, A1=a1, A2=a2, LDscore=ld), \
             'piece.1.snps.txt', sep='\\t', row.names=FALSE, quote=FALSE)",
            pieces_file.display()
        );
    }

    let pieces_content = std::fs::read_to_string(&pieces_file)
        .with_context(|| format!("failed to read {}", pieces_file.display()))?;

    let mut ld_pieces = Vec::new();
    for line in pieces_content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("piece") {
            continue; // skip header or comments
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let piece_idx: usize = fields[0]
            .parse()
            .with_context(|| format!("invalid piece index in pieces.txt: '{}'", fields[0]))?;

        // Read per-piece SNP file
        let snp_file = ld_path.join(format!("piece.{piece_idx}.snps.txt"));
        if !snp_file.exists() {
            eprintln!(
                "Warning: piece file {} not found, skipping",
                snp_file.display()
            );
            continue;
        }

        let snp_content = std::fs::read_to_string(&snp_file)
            .with_context(|| format!("failed to read {}", snp_file.display()))?;

        let mut snps = Vec::new();
        let mut a1 = Vec::new();
        let mut a2 = Vec::new();
        let mut ld_scores = Vec::new();

        for sline in snp_content.lines() {
            let sline = sline.trim();
            if sline.is_empty() || sline.starts_with('#') || sline.starts_with("SNP") {
                continue;
            }
            let sf: Vec<&str> = sline.split('\t').collect();
            if sf.len() < 4 {
                continue;
            }
            snps.push(sf[0].to_string());
            a1.push(sf[1].to_string());
            a2.push(sf[2].to_string());
            if let Ok(ld) = sf[3].parse::<f64>() {
                ld_scores.push(ld);
            } else {
                ld_scores.push(0.0);
            }
        }

        let m = snps.len();
        if m > 0 {
            ld_pieces.push(LdPiece {
                snps,
                a1,
                a2,
                ld_scores,
                m,
            });
        }
    }

    if ld_pieces.is_empty() {
        anyhow::bail!("no valid LD pieces loaded from {}", ld_path.display());
    }

    eprintln!(
        "Loaded {} LD pieces ({} total SNPs)",
        ld_pieces.len(),
        ld_pieces.iter().map(|p| p.m).sum::<usize>()
    );

    // Parse prevalences
    let k = traits.len();
    let sp = parse_prevalences(&sample_prev, k);
    let pp = parse_prevalences(&pop_prev, k);

    eprintln!("Running HDL...");
    let result = gsem_ldsc::hdl::hdl(&trait_data, &sp, &pp, &ld_pieces, &config)?;

    // Write JSON output
    let json = result.to_json_string()?;
    std::fs::write(out, &json).with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("HDL complete. Results written to {}", out.display());

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

    Ok(())
}

fn run_summary_gls(
    x_path: &Path,
    y_path: &Path,
    v_path: &Path,
    intercept: bool,
    out: &Path,
) -> Result<()> {
    let x_content = std::fs::read_to_string(x_path)
        .with_context(|| format!("failed to read {}", x_path.display()))?;
    let x_rows: Vec<Vec<f64>> = x_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(['\t', ' '])
                .filter(|s| !s.is_empty())
                .map(|s| s.parse::<f64>().unwrap_or(0.0))
                .collect()
        })
        .collect();

    let n = x_rows.len();
    if n == 0 {
        anyhow::bail!("empty X matrix");
    }
    let p_base = x_rows[0].len();
    let p = if intercept { p_base + 1 } else { p_base };

    let x = faer::Mat::from_fn(n, p, |i, j| {
        if intercept && j == 0 {
            1.0
        } else {
            let col = if intercept { j - 1 } else { j };
            x_rows[i][col]
        }
    });

    let y_content = std::fs::read_to_string(y_path)
        .with_context(|| format!("failed to read {}", y_path.display()))?;
    let y: Vec<f64> = y_content
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();

    if y.len() != n {
        anyhow::bail!("Y length {} != X rows {}", y.len(), n);
    }

    let v_content = std::fs::read_to_string(v_path)
        .with_context(|| format!("failed to read {}", v_path.display()))?;
    let v_rows: Vec<Vec<f64>> = v_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(['\t', ' '])
                .filter(|s| !s.is_empty())
                .map(|s| s.parse::<f64>().unwrap_or(0.0))
                .collect()
        })
        .collect();

    if v_rows.len() != n {
        anyhow::bail!("V rows {} != Y length {}", v_rows.len(), n);
    }
    let v = faer::Mat::from_fn(n, n, |i, j| v_rows[i][j]);

    eprintln!("GLS: {} observations, {} predictors", n, p);

    let result = gsem::stats::gls::summary_gls(&x, &y, &v)
        .ok_or_else(|| anyhow::anyhow!("GLS computation failed"))?;

    let mut output = String::from("predictor\tbeta\tse\tz\tp\n");
    for i in 0..result.beta.len() {
        let name = if intercept && i == 0 {
            "intercept".to_string()
        } else {
            format!("X{}", if intercept { i } else { i + 1 })
        };
        output.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{:.4}\t{:.6e}\n",
            name, result.beta[i], result.se[i], result.z[i], result.p[i]
        ));
    }

    std::fs::write(out, &output).with_context(|| format!("failed to write {}", out.display()))?;
    eprintln!("Results written to {}", out.display());
    Ok(())
}
