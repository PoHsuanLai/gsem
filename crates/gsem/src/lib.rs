//! gsem: multivariate genomic structural equation modeling.
//!
//! This is the main library crate providing the full GenomicSEM pipeline:
//! - [`munge`] -- QC and prepare raw GWAS summary statistics
//! - [`io`] -- Read/write GWAS files, LD scores, merged sumstats
//! - [`gwas`] -- Per-SNP multivariate GWAS with parallel processing
//! - [`stats`] -- Enrichment, parallel analysis, GLS, simulation
//!
//! For LDSC regression, see the [`gsem_ldsc`] crate.
//! For SEM fitting, see the [`gsem_sem`] crate.

pub mod gwas;
pub mod io;
pub mod munge;
pub mod stats;
