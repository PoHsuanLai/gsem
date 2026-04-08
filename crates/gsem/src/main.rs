use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

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
        } => {
            run_munge(&files, &hm3, trait_names.as_deref(), info_filter, maf_filter, n, &out)
        }
        Commands::Ldsc { .. } => {
            eprintln!("LDSC command: use the gsem-ldsc library programmatically for now.");
            eprintln!("Full CLI integration coming soon.");
            Ok(())
        }
        Commands::Sem { .. } => {
            eprintln!("SEM command: use the gsem-sem library programmatically for now.");
            eprintln!("Full CLI integration coming soon.");
            Ok(())
        }
        Commands::Gwas { .. } => {
            eprintln!("GWAS command: use the gsem library programmatically for now.");
            eprintln!("Full CLI integration coming soon.");
            Ok(())
        }
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
    let reference = munge::read_reference(hm3)
        .context("failed to read HapMap3 reference")?;
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
