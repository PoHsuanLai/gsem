//! R bindings for GenomicSEM-rs via extendr.
//!
//! This module provides the conversion layer between R data structures and Rust.
//! When extendr is enabled, functions are exported as R-callable via `.Call()`.
//!
//! To use: uncomment `extendr-api` in Cargo.toml and add `#[extendr]` annotations.

pub mod conversions;

// When extendr is available, the actual R-exported functions go here.
// For now, we provide the conversion layer and document the intended API.
//
// Intended R API (matching original GenomicSEM):
//
// ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, ...)
//   -> list(S=matrix, V=matrix, I=matrix, N=vector, m=scalar)
//
// commonfactor(covstruc, estimation="DWLS")
//   -> list(results=data.frame, modelfit=data.frame)
//
// usermodel(covstruc, estimation="DWLS", model="", ...)
//   -> list(results=data.frame, modelfit=data.frame)
//
// userGWAS(covstruc, SNPs, estimation="DWLS", model="", ...)
//   -> data.frame (one row per SNP × parameter)
//
// munge(files, hm3, trait.names, ...)
//   -> character vector of output file paths
