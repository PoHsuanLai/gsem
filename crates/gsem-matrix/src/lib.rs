//! Matrix utilities for genomic structural equation modeling.
//!
//! Standalone crate providing:
//! - [`near_pd::nearest_pd`] -- Nearest positive-definite matrix (Higham 2002)
//! - [`vech`] -- Half-vectorization and reverse for symmetric matrices
//! - [`smooth`] -- PSD checking, smoothing, covariance-correlation conversions
//!
//! These are general-purpose matrix operations useful beyond GenomicSEM.

pub mod error;
pub mod near_pd;
pub mod smooth;
pub mod vech;
