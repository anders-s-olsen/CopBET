"""
CopBET_diversity_coefficient
==============================
Copenhagen Brain Entropy Toolbox: Diversity coefficient.

Evaluates diversity coefficient as in Lebedev et al., 2015.
The correlation matrix is computed for each session and an average
correlation matrix is formed. Modularity is run 1000 times on the average
matrix; the assignment with the highest modularity score is used. For
each session, the positive diversity coefficient (Hpos from BCT's
diversity_coef_sign) is returned — a measure of diversity of out-network
connectivity.

Input
-----
sessions : list of np.ndarray, each shape (T, N)

Returns
-------
entropy : list of np.ndarray, each shape (N,)
    Per-ROI diversity coefficient for each session.

Reference
---------
Lebedev et al., 2015.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
import bct


def CopBET_diversity_coefficient(sessions, n_modularity=1000):
    """
    Compute diversity coefficient.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    n_modularity : int
        Number of modularity optimisation runs on average matrix.

    Returns
    -------
    entropy : list of np.ndarray, each shape (N,)
        Hpos (positive diversity coefficient) per ROI per session.
    """
    # Build per-session correlation matrices
    corr_matrices = []
    for ts in sessions:
        ts = np.asarray(ts, dtype=float)
        R = np.corrcoef(ts.T)
        corr_matrices.append(R)

    # Average correlation matrix
    R_avg = np.mean(np.stack(corr_matrices, axis=2), axis=2)

    # Find optimal module assignment via repeated Louvain modularity
    optimal_assignment = _optimal_modularity(R_avg, n_modularity)

    # Compute diversity coefficient per session
    entropy = []
    for R in corr_matrices:
        Hpos, _ = bct.diversity_coef_sign(R, optimal_assignment)
        entropy.append(Hpos)

    return entropy


def _optimal_modularity(W, n_iterations=1000):
    """
    Run Louvain modularity n_iterations times and return the best assignment.

    Parameters
    ----------
    W : np.ndarray, shape (N, N)
        Weighted, signed connectivity matrix.
    n_iterations : int

    Returns
    -------
    optimal_assignment : np.ndarray, shape (N,)
        Community labels (1-indexed, as expected by BCT functions).
    """
    best_q = -np.inf
    best_ci = None

    for _ in range(n_iterations):
        ci, q = bct.community_louvain(W, gamma=1.0, B='negative_sym')
        if q > best_q:
            best_q = q
            best_ci = ci.copy()

    return best_ci
