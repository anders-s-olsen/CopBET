"""
CopBET_geodesic_entropy
========================
Copenhagen Brain Entropy Toolbox: Geodesic entropy.

Calculates Geodesic entropy as in Viol et al., 2019.
For each session, the correlation matrix is binarized across a range of
thresholds. For a set of desired mean degrees, the shortest-path length
distribution is computed per ROI (ignoring infinite paths and
self-connections). The Shannon entropy of this distribution is the
geodesic entropy, reported as the mean across ROIs.

Input
-----
sessions : list of np.ndarray, each shape (T, N)

Returns
-------
entropy : list of np.ndarray, each shape (n_degrees,)
    Geodesic entropy across target degrees, per session.

Reference
---------
Viol et al., 2019.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
import bct
import networkx as nx


def CopBET_geodesic_entropy(sessions, thresholds=None, chosen_degrees=None):
    """
    Compute geodesic path-length Shannon entropy.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    thresholds : array-like or None
        Default: 0.01 to 0.99, step 0.001.
    chosen_degrees : array-like or None
        Default: 1 to 100.

    Returns
    -------
    entropy : list of np.ndarray, each shape (n_chosen_degrees,)
    """
    if thresholds is None:
        thresholds = np.arange(0.01, 1.0, 0.001)
    if chosen_degrees is None:
        chosen_degrees = np.arange(1, 101)

    entropy = []

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)
        R = np.corrcoef(ts.T)
        ent = _calc_viol19_entropy(R, thresholds, chosen_degrees)
        entropy.append(ent)

    return entropy


def _calc_viol19_entropy(R, thresholds, chosen_degrees):
    """
    Compute geodesic entropy across a range of mean degrees.
    """
    N = R.shape[0]
    bins = np.arange(0.5, N + 1.5)

    bin_matrices = []
    mean_degrees = np.zeros(len(thresholds))

    for tt, thr in enumerate(thresholds):
        B = (np.abs(R) > thr).astype(float)
        # np.fill_diagonal(B, 0)
        deg = bct.degrees_und(B)
        mean_degrees[tt] = np.mean(deg)
        bin_matrices.append(B)

    entropy_out = np.zeros(len(chosen_degrees))

    for i, target_deg in enumerate(chosen_degrees):
        idx = np.argmin(np.abs(mean_degrees - target_deg))
        B = bin_matrices[idx]

        # Compute shortest paths via networkx
        G = nx.from_numpy_array(B)
        paths = dict(nx.all_pairs_shortest_path_length(G))

        # Build distance matrix; inf where disconnected
        D = np.full((N, N), np.inf)
        for src, targets in paths.items():
            for tgt, length in targets.items():
                D[src, tgt] = length
        np.fill_diagonal(D, np.nan)  # exclude self

        entropy_n = np.full(N, np.nan)
        for n in range(N):
            row = D[n, :]
            row = row[~np.isnan(row)]
            row = row[~np.isinf(row)]
            if len(row) == 0:
                continue
            prob, _ = np.histogram(row, bins=bins, density=False)
            total = prob.sum()
            if total == 0:
                continue
            prob = prob / total
            mask = prob > 0
            entropy_n[n] = -np.sum(prob[mask] * np.log(prob[mask]))

        # Zero-entropy nodes (all paths inf) → NaN
        entropy_n[entropy_n == 0] = np.nan
        entropy_out[i] = np.nanmean(entropy_n)

    return entropy_out
