"""SEM-based enrichment analysis with per-annotation model fitting."""

from __future__ import annotations

import json
from typing import Optional

import numpy as np

from genomicsem.genomicsem import enrich as enrich_rust, usermodel


def enrich_model(
    s_covstruc: dict,
    model: str = "",
    params: Optional[list[str]] = None,
    fix: str = "regressions",
    std_lv: bool = False,
    tau: bool = False,
    toler: Optional[float] = None,
    fixparam: Optional[dict[str, float]] = None,
) -> list[dict]:
    """Enrichment analysis with optional SEM model fitting per annotation.

    When model is empty, uses the fast Rust proportional enrichment test.
    When model is provided, fits the SEM per annotation and returns parameter estimates.

    Parameters
    ----------
    s_covstruc : dict
        Stratified LDSC result with keys: S_baseline, S_annot, V_annot,
        annotation_names, m_annot, m_total.
    model : str
        lavaan-style model syntax. Empty string uses basic enrichment test.
    params : list of str, optional
        Parameter names to test (e.g. ["F1~SNP"]). None = all free params.
    fix : str
        Which parameters to fix from baseline: "regressions", "loadings", or "none".
    std_lv : bool
        Standardize latent variables.
    tau : bool
        Use tau (per-SNP) parameterization instead of total S_annot.
    toler : float, optional
        Gradient tolerance for optimizer.
    fixparam : dict, optional
        Parameters to fix at specific values.

    Returns
    -------
    list of dicts with per-annotation results.
    """
    if not model:
        # Fast Rust path
        s_baseline = json.dumps(np.asarray(s_covstruc["S_baseline"]).tolist())
        s_annot = json.dumps([np.asarray(s).tolist() for s in s_covstruc["S_annot"]])
        v_annot = json.dumps([np.asarray(v).tolist() for v in s_covstruc["V_annot"]])
        annotation_names = list(s_covstruc["annotation_names"])
        m_annot = list(s_covstruc["m_annot"])
        m_total = float(s_covstruc["m_total"])

        result_json = enrich_rust(s_baseline, s_annot, v_annot, annotation_names, m_annot, m_total)
        return json.loads(result_json)

    # SEM-based enrichment
    n_annot = len(s_covstruc["S_annot"])
    annotation_names = list(s_covstruc["annotation_names"])
    m_annot = list(s_covstruc["m_annot"])

    # Tau parameterization
    s_annot_use = []
    for a in range(n_annot):
        s_a = np.asarray(s_covstruc["S_annot"][a])
        if tau and m_annot[a] > 0:
            s_a = s_a / m_annot[a]
        s_annot_use.append(s_a)

    # Fit baseline
    k = np.asarray(s_covstruc["S_baseline"]).shape[0]
    baseline_covstruc = json.dumps({
        "s": np.asarray(s_covstruc["S_baseline"]).tolist(),
        "v": np.asarray(s_covstruc["V_annot"][0]).tolist(),
        "i_mat": np.eye(k).tolist(),
        "n_vec": [0.0] * k,
        "m": float(s_covstruc["m_total"]),
    })

    baseline_json = usermodel(baseline_covstruc, model, "DWLS", False, std_lv, False, False, False)
    baseline = json.loads(baseline_json)
    if "error" in baseline:
        raise RuntimeError(f"Baseline fit failed: {baseline['error']}")
    baseline_params = baseline["parameters"]

    # Build fixed-parameter model
    model_with_fixes = model
    if fixparam:
        for pname, pval in fixparam.items():
            model_with_fixes += f"\n{pname} == {pval}"

    # Determine which params to fix from baseline
    fix_rows = []
    if fix == "regressions":
        fix_rows = [p for p in baseline_params if p["op"] == "~"]
    elif fix == "loadings":
        fix_rows = [p for p in baseline_params if p["op"] == "=~"]

    annot_model = model_with_fixes
    for p in fix_rows:
        annot_model += f"\n{p['lhs']} {p['op']} {p['est']:.10f}*{p['rhs']}"

    # Fit per annotation
    results = []
    for a in range(n_annot):
        covstruc_json = json.dumps({
            "s": s_annot_use[a].tolist(),
            "v": np.asarray(s_covstruc["V_annot"][a]).tolist(),
            "i_mat": np.eye(k).tolist(),
            "n_vec": [0.0] * k,
            "m": 1.0 if tau else m_annot[a],
        })

        try:
            result_json = usermodel(covstruc_json, annot_model, "DWLS", False, std_lv, False, False, False)
            result = json.loads(result_json)
            if "error" in result:
                results.append({"annotation": annotation_names[a]})
                continue

            annot_params = result["parameters"]
            if params:
                annot_params = [
                    p for p in annot_params
                    if f"{p['lhs']}{p['op']}{p['rhs']}" in params
                ]

            for p in annot_params:
                p["annotation"] = annotation_names[a]
            results.extend(annot_params if annot_params else [{"annotation": annotation_names[a]}])
        except Exception:
            results.append({"annotation": annotation_names[a]})

    return results
