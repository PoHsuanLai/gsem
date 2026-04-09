//! Minimal structural equation modeling engine.
//!
//! Standalone crate for fitting SEMs from covariance matrices, designed
//! as a pure-Rust replacement for the subset of lavaan used by GenomicSEM.
//!
//! # Key components
//!
//! - [`syntax`] -- Lavaan-compatible model parser (`=~`, `~~`, `~`, `:=`)
//! - [`model`] -- LISREL all-y parameterization (Lambda, Psi, Theta, Beta)
//! - [`estimator`] -- DWLS and ML fitting via L-BFGS optimizer
//! - [`sandwich`] -- Corrected standard errors
//! - [`fit_indices`] -- Chi-square, CFI, AIC, SRMR

pub mod commonfactor;
pub mod delta;
pub mod error;
pub mod estimator;
pub mod fit_indices;
pub mod model;
pub mod q_factor;
pub mod reorder;
pub mod rgmodel;
pub mod sandwich;
pub mod syntax;
pub mod write_model;

use faer::Mat;

/// Result of SEM fitting.
#[derive(Debug, Clone)]
pub struct SemResult {
    /// Parameter estimates
    pub parameters: Vec<ParamEstimate>,
    /// Model fit statistics
    pub fit: ModelFit,
    /// Model-implied covariance matrix
    pub implied_cov: Mat<f64>,
}

/// A single parameter estimate.
#[derive(Debug, Clone)]
pub struct ParamEstimate {
    pub lhs: String,
    pub op: String,
    pub rhs: String,
    pub est: f64,
    pub se: f64,
    pub z: f64,
    pub p: f64,
}

/// Model fit statistics.
#[derive(Debug, Clone)]
pub struct ModelFit {
    pub chisq: f64,
    pub df: usize,
    pub p_chisq: f64,
    pub aic: f64,
    pub cfi: f64,
    pub srmr: f64,
}
