use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

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

        /// Trait names
        #[arg(long, num_args = 1..)]
        trait_names: Option<Vec<String>>,

        /// Number of jackknife blocks
        #[arg(long, default_value = "200")]
        n_blocks: usize,

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

        /// Output file
        #[arg(short, long, default_value = "sem_result.tsv")]
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

        /// Output file
        #[arg(short, long, default_value = "gwas_result.tsv")]
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
            out,
        } => run_munge(
            &files,
            &hm3,
            trait_names.as_deref(),
            info_filter,
            maf_filter,
            n,
            &out,
        ),
        Commands::Ldsc {
            traits,
            sample_prev,
            pop_prev,
            ld,
            wld,
            trait_names,
            n_blocks,
            out,
        } => run_ldsc(
            &traits,
            sample_prev,
            pop_prev,
            &ld,
            wld,
            trait_names,
            n_blocks,
            &out,
        ),
        Commands::Sem {
            covstruc,
            model,
            model_file,
            estimation,
            out,
        } => run_sem(&covstruc, model, model_file, &estimation, &out),
        Commands::Gwas {
            covstruc,
            sumstats,
            model,
            model_file,
            estimation,
            gc,
            threads,
            out,
        } => run_gwas(
            &covstruc,
            &sumstats,
            model,
            model_file,
            &estimation,
            &gc,
            threads,
            &out,
        ),
    }
}

fn run_munge(
    files: &[PathBuf],
    hm3: &PathBuf,
    trait_names: Option<&[String]>,
    info_filter: f64,
    maf_filter: f64,
    n_override: Option<f64>,
    out_dir: &PathBuf,
) -> Result<()> {
    // Read reference
    eprintln!("Reading reference panel: {}", hm3.display());
    let reference = munge::read_reference(hm3).context("failed to read HapMap3 reference")?;
    eprintln!("Loaded {} reference SNPs", reference.len());

    let config = munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override,
    };

    for (i, file) in files.iter().enumerate() {
        let trait_name = trait_names
            .and_then(|names| names.get(i))
            .map(|s| s.as_str())
            .unwrap_or_else(|| file.file_stem().unwrap().to_str().unwrap());

        let out_path = out_dir.join(format!("{trait_name}.sumstats.gz"));
        eprintln!("Munging: {} -> {}", file.display(), out_path.display());

        munge::munge_and_write(file, &reference, &config, &out_path)
            .with_context(|| format!("failed to munge {}", file.display()))?;
    }

    eprintln!("Done.");
    Ok(())
}

fn run_ldsc(
    traits: &[PathBuf],
    sample_prev: Option<String>,
    pop_prev: Option<String>,
    ld: &PathBuf,
    wld: Option<PathBuf>,
    _trait_names: Option<Vec<String>>,
    n_blocks: usize,
    out: &PathBuf,
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
    let wld_dir = wld.as_ref().unwrap_or(ld);
    let ld_data =
        gsem::io::ld_reader::read_ld_scores(ld, wld_dir, 22).context("failed to read LD scores")?;

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
        chisq_max: None,
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

    Ok(())
}

fn run_sem(
    covstruc: &PathBuf,
    model: Option<String>,
    model_file: Option<PathBuf>,
    estimation: &str,
    out: &PathBuf,
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

    eprintln!("Fitting SEM model (estimation={estimation})...");

    // Parse and fit
    let pt =
        gsem_sem::syntax::parse_model(&model_str, false).map_err(|e| anyhow::anyhow!("{e}"))?;
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

fn run_gwas(
    covstruc: &PathBuf,
    _sumstats: &PathBuf,
    model: Option<String>,
    model_file: Option<PathBuf>,
    estimation: &str,
    gc: &str,
    threads: Option<usize>,
    out: &PathBuf,
) -> Result<()> {
    // Set thread count
    if let Some(t) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }

    // Read LDSC result
    let json = std::fs::read_to_string(covstruc)
        .with_context(|| format!("failed to read {}", covstruc.display()))?;
    let _ldsc_result = gsem_ldsc::LdscResult::from_json_string(&json)?;

    // Get model string
    let model_str = if let Some(m) = model {
        m
    } else if let Some(f) = model_file {
        std::fs::read_to_string(&f)?
    } else {
        anyhow::bail!("must provide --model or --model-file");
    };

    let gc_mode = gsem::gwas::gc_correction::GcMode::from_str(gc);
    eprintln!(
        "GWAS: estimation={estimation}, gc={gc:?}, model length={} chars",
        model_str.len()
    );
    eprintln!("Note: full GWAS requires merged sumstats with SNP-level betas/SEs.");
    eprintln!(
        "Use the library API (gsem::gwas::user_gwas::run_user_gwas) for programmatic access."
    );

    // Placeholder output
    std::fs::write(
        out,
        "# GWAS results (run with library API for full output)\n",
    )
    .with_context(|| format!("failed to write {}", out.display()))?;

    eprintln!("Output: {}", out.display());
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
