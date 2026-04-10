use faer::Mat;
use gsem_ldsc::LdscResult;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;

// ---------------------------------------------------------------------------
// Matrix <-> NumPy / flat vec
// ---------------------------------------------------------------------------

/// Convert a faer Mat to a flat row-major Vec for NumPy.
///
/// Kept for the `PyLdscResult::to_json` / `from_json` compatibility path.
pub fn mat_to_flat(mat: &Mat<f64>) -> (Vec<f64>, usize, usize) {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let mut data = Vec::with_capacity(nrows * ncols);
    for i in 0..nrows {
        for j in 0..ncols {
            data.push(mat[(i, j)]);
        }
    }
    (data, nrows, ncols)
}

/// Convert JSON string back to LdscResult. Used by `PyLdscResult::from_json`.
pub fn json_to_ldsc(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}

/// Convert a faer `Mat<f64>` into an owned NumPy `PyArray2<f64>`.
pub fn mat_to_pyarray<'py>(py: Python<'py>, m: &Mat<f64>) -> Bound<'py, PyArray2<f64>> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    let mut data = Vec::with_capacity(nrows * ncols);
    for i in 0..nrows {
        for j in 0..ncols {
            data.push(m[(i, j)]);
        }
    }
    Array2::from_shape_vec((nrows, ncols), data)
        .expect("shape mismatch building PyArray2")
        .into_pyarray(py)
}

/// Copy a read-only NumPy 2D array into a faer `Mat<f64>`.
pub fn pyarray_to_mat(arr: &PyReadonlyArray2<'_, f64>) -> Mat<f64> {
    let view = arr.as_array();
    let (nrows, ncols) = (view.nrows(), view.ncols());
    Mat::from_fn(nrows, ncols, |i, j| view[[i, j]])
}

// ---------------------------------------------------------------------------
// LdscResult <-> Python objects
// ---------------------------------------------------------------------------

/// Read a field from an arbitrary Python object. Accepts either a dict
/// lookup (e.g. `d["s"]` / `d["S"]`) or an attribute access
/// (`obj.s` / `obj.S`). Returns `None` when the field is missing or
/// explicitly set to `None`/null.
fn pyany_field<'py>(
    obj: &Bound<'py, PyAny>,
    primary: &str,
    alt: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    // Try as a dict first.
    if let Ok(dict) = obj.cast::<PyDict>() {
        if let Some(v) = dict.get_item(primary)?
            && !v.is_none()
        {
            return Ok(Some(v));
        }
        if let Some(v) = dict.get_item(alt)?
            && !v.is_none()
        {
            return Ok(Some(v));
        }
        return Ok(None);
    }
    // Fall back to attribute access (covers PyLdscResult and any user class).
    if let Ok(v) = obj.getattr(primary)
        && !v.is_none()
    {
        return Ok(Some(v));
    }
    if let Ok(v) = obj.getattr(alt)
        && !v.is_none()
    {
        return Ok(Some(v));
    }
    Ok(None)
}

/// Convert any Python object (a `PyLdscResult`, a dict with `s`/`v`/`i_mat`
/// keys, a dict with `S`/`V`/`I` keys, or a user-defined class with those
/// attributes) into an owned `LdscResult`.
pub fn pyany_to_ldsc_result(obj: &Bound<'_, PyAny>) -> PyResult<LdscResult> {
    let s_obj = pyany_field(obj, "s", "S")?
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("covstruc missing 's' (or 'S')"))?;
    let v_obj = pyany_field(obj, "v", "V")?
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("covstruc missing 'v' (or 'V')"))?;
    // i_mat is exposed as `i_mat` on PyLdscResult; fall back to `I` for dict users.
    let i_obj = pyany_field(obj, "i_mat", "I")?.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("covstruc missing 'i_mat' (or 'I')")
    })?;

    let s = robj_to_mat(&s_obj).map_err(prepend_err("covstruc.s"))?;
    let v = robj_to_mat(&v_obj).map_err(prepend_err("covstruc.v"))?;
    let i_mat = robj_to_mat(&i_obj).map_err(prepend_err("covstruc.i_mat"))?;

    // n_vec: accept either `n_vec` (dict) or `n` (PyLdscResult getter) or `N` (dict).
    let n_vec: Vec<f64> = match pyany_field(obj, "n_vec", "N")? {
        Some(v) => v.extract::<Vec<f64>>().unwrap_or_default(),
        None => match pyany_field(obj, "n", "n")? {
            Some(v) => v.extract::<Vec<f64>>().unwrap_or_default(),
            None => Vec::new(),
        },
    };

    // m: accept `m`, `M`, or `m_total` (the PyLdscResult getter).
    let m: f64 = match pyany_field(obj, "m", "M")? {
        Some(v) => v.extract::<f64>().unwrap_or(0.0),
        None => match pyany_field(obj, "m_total", "m_total")? {
            Some(v) => v.extract::<f64>().unwrap_or(0.0),
            None => 0.0,
        },
    };

    Ok(LdscResult {
        s,
        v,
        i_mat,
        n_vec,
        m,
    })
}

/// Helper: convert any numpy-array-ish Python object into a faer `Mat<f64>`.
fn robj_to_mat(obj: &Bound<'_, PyAny>) -> PyResult<Mat<f64>> {
    let arr: PyReadonlyArray2<'_, f64> = obj.extract()?;
    Ok(pyarray_to_mat(&arr))
}

/// Decorate an `extract`/`downcast` error with the parameter name so the
/// caller sees `covstruc.s: ...` instead of a bare type error.
fn prepend_err(label: &'static str) -> impl FnOnce(PyErr) -> PyErr {
    move |e| pyo3::exceptions::PyValueError::new_err(format!("{label}: {e}"))
}
