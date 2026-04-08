//! Python bindings for gsem via PyO3.
//!
//! Build with `--features pyo3` and maturin:
//! ```sh
//! maturin develop --release --features pyo3
//! ```
//!
//! ## Python API
//!
//! ```python
//! import genomicsem as gsem
//!
//! # LDSC
//! result = gsem.ldsc(
//!     traits=["Trait1.sumstats.gz", "Trait2.sumstats.gz"],
//!     sample_prev=[0.5, 0.5],
//!     pop_prev=[0.01, 0.02],
//!     ld="eur_w_ld_chr/",
//!     wld="eur_w_ld_chr/",
//! )
//! result.S  # np.ndarray
//! result.V  # np.ndarray
//!
//! # SEM
//! cf = gsem.usermodel(result, model="F1 =~ V1 + V2", estimation="DWLS")
//! ```

pub mod conversions;

#[cfg(feature = "pyo3")]
mod py_exports;
