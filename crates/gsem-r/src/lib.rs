//! R bindings for gsem via extendr.
//!
//! Build with `--features extendr` when R toolchain is available.
//! Without the feature, only the conversion layer is available.
//!
//! ## R API (matching original GenomicSEM)
//!
//! ```r
//! library(gsemr)
//! result <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
//! # result$S, result$V, result$I, result$N, result$m
//!
//! cf <- commonfactor(result, estimation="DWLS")
//! um <- usermodel(result, estimation="DWLS", model="F1 =~ V1 + V2")
//! ```

pub mod conversions;

#[cfg(feature = "extendr")]
mod r_exports;

#[cfg(feature = "extendr")]
pub use r_exports::*;
