"""
CopBET_degree_distribution_entropy
=====================================
Copenhagen Brain Entropy Toolbox: Degree distribution entropy.

Evaluates temporal entropy as in Viol et al., 2017.
For each session, the Pearson correlation matrix is computed and entries
with non-significant p-values (p > 0.05) are set to zero. The matrix is
then binarized across a range of thresholds to yield a sweep of mean
degrees. For a set of desired mean degrees, the degree distribution is
computed and its Shannon entropy returned.

Input
-----
sessions : list of np.ndarray, each shape (T, N)

Returns
-------
entropy : list of np.ndarray, each shape (n_degrees,)
    Degree-distribution entropy across target degrees, per session.

Reference
---------
Viol et al., 2017.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from scipy import stats
import bct


def CopBET_degree_distribution_entropy(sessions,
                                       thresholds=None,
                                       chosen_degrees=None):
    """
    Compute degree-distribution Shannon entropy.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    thresholds : array-like or None
        Absolute value thresholds to sweep (default: 0.01 to 0.99, step 0.001).
    chosen_degrees : array-like or None
        Target mean degrees to evaluate (default: 1 to 100).

    Returns
    -------
    entropy : list of np.ndarray, each shape (n_chosen_degrees,)
    """
    if thresholds is None:
        thresholds = np.arange(0.01, 0.991, 0.001)
    if chosen_degrees is None:
        chosen_degrees = np.arange(1, 101)

    entropy = []

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)
        T, N = ts.shape

        # Pearson correlation + p-values
        R, pval = _corrcoef_pval(ts)

        # Zero out non-significant correlations
        R[pval > 0.05] = 0.0

        ent = _calc_viol17_entropy(R, thresholds, chosen_degrees)
        entropy.append(ent)

    return entropy


def _corrcoef_pval(ts):
    """Compute Pearson correlation matrix and matrix of p-values."""
    T, N = ts.shape
    R = np.corrcoef(ts.T)
    pval = np.ones((N, N))

    for i in range(N):
        for j in range(i + 1, N):
            r = R[i, j]
            # Two-tailed t-test: t = r * sqrt(T-2) / sqrt(1-r^2)
            r_clamped = np.clip(r, -1 + 1e-10, 1 - 1e-10)
            t_stat = r_clamped * np.sqrt(T - 2) / np.sqrt(1 - r_clamped ** 2)
            p = 2 * stats.t.sf(np.abs(t_stat), df=T - 2)
            pval[i, j] = p
            pval[j, i] = p

    np.fill_diagonal(pval, 0.0)
    return R, pval


def _calc_viol17_entropy(R, thresholds, chosen_degrees):
    """
    Compute degree-distribution Shannon entropy for a range of mean degrees.

    Parameters
    ----------
    R : np.ndarray, shape (N, N)
        (Thresholded) correlation matrix.
    thresholds : array-like
        Absolute value thresholds to sweep.
    chosen_degrees : array-like
        Target mean degrees.

    Returns
    -------
    entropy_out : np.ndarray, shape (len(chosen_degrees),)
    """
    N = R.shape[0]
    bins = np.arange(0.5, N + 1.5)  # degree bins 1..N

    # Binarize at each threshold and record mean degree
    bin_matrices = []
    mean_degrees = np.zeros(len(thresholds))
    for tt, thr in enumerate(thresholds):
        B = (np.abs(R) > thr).astype(float)
        np.fill_diagonal(B, 0)
        deg = bct.degrees_und(B)
        mean_degrees[tt] = np.mean(deg)
        bin_matrices.append(B)

    # For each desired mean degree, find nearest threshold and compute entropy
    entropy_out = np.zeros(len(chosen_degrees))
    for i, target_deg in enumerate(chosen_degrees):
        idx = np.argmin(np.abs(mean_degrees - target_deg))
        B = bin_matrices[idx]
        deg = bct.degrees_und(B)
        prob, _ = np.histogram(deg, bins=bins, density=False)
        prob = prob / prob.sum() if prob.sum() > 0 else prob
        mask = prob > 0
        entropy_out[i] = -np.sum(prob[mask] * np.log(prob[mask]))

    return entropy_out
