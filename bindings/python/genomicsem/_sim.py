"""Extended simulation with rPheno, intercept, and sample overlap support."""

from __future__ import annotations

from typing import Optional

import numpy as np

from genomicsem.genomicsem import sim_ldsc


def sim_ldsc_extended(
    s_matrix,
    n_per_trait: list[float],
    ld_scores: list[float],
    m: float,
    r_pheno: Optional[np.ndarray] = None,
    intercepts: Optional[np.ndarray] = None,
    n_overlap: float = 0.0,
) -> np.ndarray:
    """Simulate GWAS summary statistics with optional environmental correlation.

    Parameters
    ----------
    s_matrix : np.ndarray or list of lists
        Genetic covariance matrix.
    n_per_trait : list of float
        Per-trait sample sizes.
    ld_scores : list of float
        LD scores per SNP.
    m : float
        Total number of SNPs.
    r_pheno : np.ndarray, optional
        Phenotypic correlation matrix.
    intercepts : np.ndarray, optional
        LDSC intercept matrix (default: identity).
    n_overlap : float
        Sample overlap proportion (0 to 1).

    Returns
    -------
    np.ndarray of shape (k traits, n_snps) with simulated Z-statistics.
    """
    S = np.ascontiguousarray(np.asarray(s_matrix, dtype=np.float64))
    int_arr = (
        np.ascontiguousarray(np.asarray(intercepts, dtype=np.float64))
        if intercepts is not None
        else None
    )
    r_arr = (
        np.ascontiguousarray(np.asarray(r_pheno, dtype=np.float64))
        if r_pheno is not None
        else None
    )
    return sim_ldsc(
        S,
        list(n_per_trait),
        list(ld_scores),
        float(m),
        int_arr,
        r_arr,
        float(n_overlap),
    )
