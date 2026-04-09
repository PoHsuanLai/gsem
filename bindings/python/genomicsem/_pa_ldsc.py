"""Extended parallel analysis with optional factor analysis and PDF output."""

from __future__ import annotations

import json
from typing import Optional

import numpy as np

from genomicsem.genomicsem import parallel_analysis


def pa_ldsc_extended(
    S: np.ndarray,
    V: np.ndarray,
    r: int = 500,
    p: Optional[float] = None,
    save_pdf: bool | str = False,
    diag: bool = False,
    fa: bool = False,
    fm: Optional[str] = None,
    nfactors: Optional[int] = None,
) -> dict:
    """Parallel analysis with optional factor analysis and scree plot.

    Parameters
    ----------
    S : np.ndarray
        Genetic covariance matrix.
    V : np.ndarray
        Sampling covariance matrix.
    r : int
        Number of Monte Carlo simulations.
    p : float, optional
        Percentile threshold (default 0.95).
    save_pdf : bool or str
        Save scree plot. True uses default filename, str specifies path.
    diag : bool
        Use only diagonal of V for simulation.
    fa : bool
        Use factor analysis (requires scikit-learn or factor_analyzer).
    fm : str, optional
        Factor method (e.g. "minres", "ml", "pa").
    nfactors : int, optional
        Number of factors to extract.

    Returns
    -------
    dict with observed, simulated_95, n_factors, and optionally fa results.
    """
    s_json = json.dumps(S.tolist())
    v_json = json.dumps(V.tolist())

    result_json = parallel_analysis(s_json, v_json, r, p, diag=diag)
    result = json.loads(result_json)

    # Factor analysis mode
    if fa or fm is not None:
        result["fa"] = _run_factor_analysis(S, fm or "minres", nfactors or result["n_factors"])

    # Scree plot
    if save_pdf:
        pdf_path = save_pdf if isinstance(save_pdf, str) else "paLDSC_scree.pdf"
        _save_scree_plot(result["observed"], result["simulated_95"], pdf_path)
        print(f"Scree plot saved to: {pdf_path}")

    return result


def _run_factor_analysis(S: np.ndarray, fm: str, nfactors: int) -> dict:
    """Run factor analysis using available Python packages."""
    # Convert to correlation matrix
    sds = np.sqrt(np.maximum(np.diag(S), 1e-10))
    cor_mat = S / np.outer(sds, sds)
    np.fill_diagonal(cor_mat, 1.0)

    try:
        from factor_analyzer import FactorAnalyzer

        method = {"minres": "minres", "ml": "ml", "pa": "principal"}.get(fm, "minres")
        fa_obj = FactorAnalyzer(n_factors=nfactors, method=method, rotation=None)
        fa_obj.fit(cor_mat, is_corr_matrix=True)
        return {
            "loadings": fa_obj.loadings_.tolist(),
            "uniquenesses": fa_obj.get_uniquenesses().tolist(),
            "fm": fm,
            "nfactors": nfactors,
        }
    except ImportError:
        pass

    try:
        from sklearn.decomposition import FactorAnalysis

        fa_obj = FactorAnalysis(n_components=nfactors)
        fa_obj.fit(cor_mat)
        return {
            "loadings": fa_obj.components_.T.tolist(),
            "uniquenesses": fa_obj.noise_variance_.tolist(),
            "fm": fm,
            "nfactors": nfactors,
        }
    except ImportError:
        pass

    raise ImportError(
        "Factor analysis requires 'factor_analyzer' or 'scikit-learn'. "
        "Install with: pip install factor-analyzer"
    )


def _save_scree_plot(observed: list, simulated: list, path: str) -> None:
    """Save scree plot to PDF."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        k = len(observed)
        x = list(range(1, k + 1))
        plt.figure(figsize=(7, 5))
        plt.plot(x, observed, "b-o", label="Observed")
        plt.plot(x, simulated, "r--^", label="Simulated")
        plt.xlabel("Factor")
        plt.ylabel("Eigenvalue")
        plt.title("Parallel Analysis Scree Plot")
        plt.legend()
        plt.tight_layout()
        plt.savefig(path)
        plt.close()
    except ImportError:
        raise ImportError(
            "Scree plot requires matplotlib. Install with: pip install matplotlib"
        )
