"""Extended simulation with rPheno, intercept, and sample overlap support."""

from __future__ import annotations

import json
from typing import Optional

import numpy as np

from genomicsem.genomicsem import sim_ldsc


def sim_ldsc_extended(
    s_matrix: list[list[float]],
    n_per_trait: list[float],
    ld_scores: list[float],
    m: float,
    r_pheno: Optional[np.ndarray] = None,
    intercepts: Optional[np.ndarray] = None,
    n_overlap: float = 0.0,
) -> list[list[float]]:
    """Simulate GWAS summary statistics with optional environmental correlation.

    Parameters
    ----------
    s_matrix : list of lists
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
    list of lists: k traits x n_snps Z-statistics.
    """
    # The Rust sim_ldsc currently takes the basic args.
    # For extended params, we call the Rust function with the extended config
    # via the JSON-based interface if available, otherwise fall back to Python.

    # Use Rust function directly for basic case
    if r_pheno is None and intercepts is None:
        return sim_ldsc(s_matrix, n_per_trait, ld_scores, m)

    # Extended simulation in Python using numpy
    k = len(s_matrix)
    n_snps = len(ld_scores)
    S = np.array(s_matrix)
    N = np.array(n_per_trait)
    sqrt_nn = np.sqrt(np.outer(N, N))

    int_mat = np.array(intercepts) if intercepts is not None else np.eye(k)
    env_cov = np.zeros((k, k))
    if r_pheno is not None and n_overlap > 0:
        env_cov = np.array(r_pheno) * n_overlap * sqrt_nn / n_snps

    z_all = np.zeros((k, n_snps))
    for s_idx in range(n_snps):
        ld = ld_scores[s_idx]
        sigma_z = int_mat + (S / m * ld * sqrt_nn) + env_cov

        # Ensure PD
        try:
            L = np.linalg.cholesky(sigma_z)
        except np.linalg.LinAlgError:
            eigvals, eigvecs = np.linalg.eigh(sigma_z)
            eigvals = np.maximum(eigvals, 0)
            L = eigvecs @ np.diag(np.sqrt(eigvals))

        z_all[:, s_idx] = L @ np.random.randn(k)

    return z_all.tolist()
