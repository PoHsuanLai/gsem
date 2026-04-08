//! Python bindings for GenomicSEM-rs via PyO3.
//!
//! When PyO3 is enabled, this crate builds a native Python extension module
//! installable via `maturin develop --release` or `pip install .`.
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
//! cf = gsem.commonfactor(result, estimation="DWLS")
//! cf.results     # pandas DataFrame
//! cf.model_fit   # dict
//!
//! # GWAS
//! gwas_df = gsem.user_gwas(result, sumstats="merged.tsv", model="...", gc="standard")
//! ```

pub mod conversions;

// PyO3 module definition (uncomment when pyo3 is enabled):
//
// #[pymodule]
// fn genomicsem(m: &Bound<'_, PyModule>) -> PyResult<()> {
//     m.add_class::<LdscResult>()?;
//     m.add_class::<SemResult>()?;
//     m.add_function(wrap_pyfunction!(ldsc, m)?)?;
//     m.add_function(wrap_pyfunction!(commonfactor, m)?)?;
//     m.add_function(wrap_pyfunction!(usermodel, m)?)?;
//     m.add_function(wrap_pyfunction!(user_gwas, m)?)?;
//     m.add_function(wrap_pyfunction!(munge, m)?)?;
//     Ok(())
// }
