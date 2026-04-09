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
pub mod jacobian;
pub mod fit_indices;
pub mod model;
pub mod q_factor;
pub mod reorder;
pub mod rgmodel;
pub mod sandwich;
pub mod syntax;
pub mod write_model;

use faer::Mat;

use crate::syntax::Op;

/// SEM estimation method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimationMethod {
    /// Diagonally Weighted Least Squares
    Dwls,
    /// Maximum Likelihood
    Ml,
}

impl EstimationMethod {
    /// Parse from a string (case-insensitive). Defaults to DWLS for unrecognized values.
    pub fn from_str_lossy(s: &str) -> Self {
        if s.eq_ignore_ascii_case("ML") {
            Self::Ml
        } else {
            Self::Dwls
        }
    }
}

impl std::fmt::Display for EstimationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Dwls => write!(f, "DWLS"),
            Self::Ml => write!(f, "ML"),
        }
    }
}

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
    /// Left-hand side variable name
    pub lhs: String,
    /// Operator type
    pub op: Op,
    /// Right-hand side variable name
    pub rhs: String,
    /// Point estimate
    pub est: f64,
    /// Standard error
    pub se: f64,
    /// Z-statistic (est / se)
    pub z: f64,
    /// P-value (two-tailed)
    pub p: f64,
}

/// Model fit statistics.
#[derive(Debug, Clone)]
pub struct ModelFit {
    /// Chi-square test statistic
    pub chisq: f64,
    /// Degrees of freedom
    pub df: usize,
    /// P-value for chi-square test
    pub p_chisq: f64,
    /// Akaike information criterion
    pub aic: f64,
    /// Comparative fit index
    pub cfi: f64,
    /// Standardized root mean square residual
    pub srmr: f64,
}
