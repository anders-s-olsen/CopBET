"""
CopBET_motif_connectivity_entropy
===================================
Copenhagen Brain Entropy Toolbox: Motif-connectivity entropy.

Evaluates motif-connectivity entropy as in Tagliazucchi et al., 2014.
Non-overlapping sliding windows of varying lengths are slid across the
data. In each window, the partial correlation between the 4 ROIs
(conditioned on the remaining ROIs and a motion confound) is binarized
by p-value. The resulting 6-bit connectivity pattern is matched to one
of 64 possible graphs. Shannon entropy of the graph distribution is
returned per window length. The mean across window lengths is also returned.

NOTE: This function requires exactly 4 ROIs. If more are provided, the
first 4 are used unless roi_indices is specified.

Input
-----
sessions : list of np.ndarray, each shape (T, N) with N >= 4
    ROI mean time series per session (demeaned, standardized recommended).
motion : list of np.ndarray or None
    Motion confound time series, each shape (T, n_confounds).
    If None, partial correlations are replaced by Pearson correlations.
TR : float
    Repetition time in seconds (default 2.0).
roi_indices : list of int or None
    Indices of the 4 ROIs to use. Default: [0, 1, 2, 3].

Returns
-------
entropy : list of dict with keys:
    'per_window_length' : np.ndarray, shape (n_wl,) — entropy per window length
    'mean' : float — mean entropy across window lengths

Reference
---------
Tagliazucchi et al., 2014.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from scipy import stats
from itertools import product


# All 64 possible binary graphs on 6 edges
_POSSIBLE_GRAPHS = np.array(list(product([0, 1], repeat=6)), dtype=float).T  # (6, 64)
_ROI_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def CopBET_motif_connectivity_entropy(sessions, motion=None, TR=2.0,
                                      roi_indices=None,
                                      window_lengths_sec=None):
    """
    Compute motif-connectivity Shannon entropy.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N), N >= 4
    motion : list of np.ndarray or None
        Motion confound series per session, shape (T, n_confounds).
    TR : float
    roi_indices : list of int or None
        Indices for the 4 ROIs. Default [0, 1, 2, 3].
    window_lengths_sec : array-like or None
        Window lengths in seconds. Default: 15 to 150 seconds step 1.

    Returns
    -------
    entropy : list of dict
    """
    if roi_indices is None:
        roi_indices = [0, 1, 2, 3]
    if window_lengths_sec is None:
        window_lengths_sec = np.arange(15, 151)

    if motion is None:
        motion = [None] * len(sessions)

    entropy = []

    for ses, ts in enumerate(sessions):
        ts = np.asarray(ts, dtype=float)

        # Extract 4 ROIs, demean, standardize
        data4 = ts[:, roi_indices]  # (T, 4)
        data4 = (data4 - data4.mean(axis=0)) / (data4.std(axis=0) + 1e-10)
        data4 = data4.T  # (4, T)

        rp = motion[ses]
        if rp is not None:
            rp = np.asarray(rp, dtype=float)
            rp = _fwd_calc(rp)[:, np.newaxis]  # (T, 1)
        else:
            rp = np.zeros((ts.shape[0], 1))

        word_counts = _state_distribution(data4, rp, TR, window_lengths_sec)
        ent = _shannon_entropy(word_counts)

        entropy.append({
            'per_window_length': ent,
            'mean': float(np.nanmean(ent)),
            'window_lengths_sec': window_lengths_sec,
        })

    return entropy


def _fwd_calc(rp):
    """
    Compute framewise displacement from 6 motion parameters (3 trans, 3 rot).
    Rotation parameters are converted from radians to mm (r=50 mm).
    """
    radius = 50.0
    ts = rp.copy()
    ts[:, 3:6] = (2 * radius * np.pi / 360) * ts[:, 3:6] * (180 / np.pi)
    dts = np.diff(ts, axis=0)
    fwd = np.concatenate([[0], np.sum(np.abs(dts), axis=1)])
    fwd = (fwd - fwd.mean()) / (fwd.std() + 1e-10)
    return fwd


def _state_distribution(data, rp, TR, window_lengths_sec):
    """
    For each window length, count occurrences of each 64-bit graph pattern.

    Parameters
    ----------
    data : np.ndarray, shape (4, T)
    rp : np.ndarray, shape (T, n_confounds)
    TR : float
    window_lengths_sec : array-like

    Returns
    -------
    word_counts : np.ndarray, shape (n_wl, 64)
    """
    T = data.shape[1]
    eff_wl = np.round(np.asarray(window_lengths_sec) / TR).astype(int)
    eff_wl = np.maximum(eff_wl, 5)  # minimum window of 5 TRs

    word_counts = np.zeros((len(window_lengths_sec), 64), dtype=int)

    for wl_idx, wl in enumerate(eff_wl):
        if wl_idx > 0 and eff_wl[wl_idx] == eff_wl[wl_idx - 1]:
            word_counts[wl_idx] = word_counts[wl_idx - 1]
            continue

        # Non-overlapping windows
        window_starts = np.arange(0, T - wl + 1, wl)

        for ws in window_starts:
            data_win = data[:, ws:ws + wl]    # (4, wl)
            rp_win = np.abs(rp[ws:ws + wl])   # (wl, n_confounds)

            graph_bits = np.zeros(6, dtype=int)
            for pair_idx, (i, j) in enumerate(_ROI_PAIRS):
                # Confound: all other ROIs + motion
                other_rois = [k for k in range(4) if k != i and k != j]
                confounds = np.hstack([
                    data_win[other_rois, :].T,   # (wl, 2)
                    rp_win                        # (wl, n_confounds)
                ])  # (wl, n_confounds+2)

                xi = data_win[i, :]
                xj = data_win[j, :]

                try:
                    r, p = _partial_corr(xi, xj, confounds)
                    graph_bits[pair_idx] = int(p < 0.05 / 6)  # Bonferroni
                except Exception:
                    graph_bits[pair_idx] = 0

            # Find matching graph pattern
            diff = np.sum(np.abs(_POSSIBLE_GRAPHS - graph_bits[:, np.newaxis]), axis=0)
            word_idx = np.argmin(diff)
            word_counts[wl_idx, word_idx] += 1

    return word_counts


def _partial_corr(x, y, z):
    """
    Compute partial correlation of x and y controlling for z.

    Parameters
    ----------
    x, y : np.ndarray, shape (T,)
    z : np.ndarray, shape (T, k)

    Returns
    -------
    r : float, p : float
    """
    T = len(x)
    if z.shape[1] == 0:
        return stats.pearsonr(x, y)

    # Residualize x and y on z
    z_with_intercept = np.column_stack([np.ones(T), z])
    px = np.linalg.lstsq(z_with_intercept, x, rcond=None)[0]
    py = np.linalg.lstsq(z_with_intercept, y, rcond=None)[0]
    res_x = x - z_with_intercept @ px
    res_y = y - z_with_intercept @ py

    return stats.pearsonr(res_x, res_y)


def _shannon_entropy(word_counts):
    """
    Compute Shannon entropy of graph distribution for each window length.

    Parameters
    ----------
    word_counts : np.ndarray, shape (n_wl, 64)

    Returns
    -------
    ent : np.ndarray, shape (n_wl,)
    """
    ent = np.zeros(word_counts.shape[0])
    for wl in range(word_counts.shape[0]):
        total = word_counts[wl].sum()
        if total == 0:
            ent[wl] = np.nan
            continue
        prob = word_counts[wl] / total
        mask = prob > 0
        ent[wl] = np.sum(prob[mask] * np.log2(1.0 / prob[mask]))
    return ent
