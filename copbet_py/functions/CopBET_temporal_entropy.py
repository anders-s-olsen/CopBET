"""
CopBET_temporal_entropy
========================
Copenhagen Brain Entropy Toolbox: Temporal entropy.

Evaluates temporal entropy as in Luppi et al., 2021.
Tapered sliding-window connectivity matrices are constructed; the
Louvain community-detection algorithm is run on each window; module-degree
z-score and participation coefficient are computed. A 2D Cartographic
Profile (CP) histogram is built from these, then K=2 K-means clusters the
CP across windows. The entropy of the resulting binary state sequence
is returned.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
TR : float
    Repetition time in seconds (used to construct the tapered window).

Returns
-------
entropy : np.ndarray, shape (n_sessions,)

Reference
---------
Luppi et al., 2021.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from scipy.signal import windows as sig_windows
from sklearn.cluster import KMeans
import bct


def CopBET_temporal_entropy(sessions, TR=2.0, n_kmeans=500, random_state=42):
    """
    Compute temporal entropy.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    TR : float
        Repetition time in seconds.
    n_kmeans : int
        Number of K-means restarts for the CP clustering step.
    random_state : int

    Returns
    -------
    entropy : np.ndarray, shape (n_sessions,)
    """
    window = _construct_tapered_window(TR)
    entropy = np.full(len(sessions), np.nan)

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)
        ts = ts - ts.mean(axis=0)
        entropy[ses] = _luppi21_entropy(ts, window, n_kmeans, random_state)

    return entropy


def _construct_tapered_window(TR):
    """
    Build a tapered window: convolution of a rectangular and Gaussian window.

    Matches the MATLAB code: 22 TR window, FWHM 3 TRs.
    """
    # 22-TR time axis, odd length
    n_pts = 23  # ensures odd number (22+1 or trimmed)
    time_axis = np.arange(0, 23) * TR
    if len(time_axis) % 2 == 0:
        time_axis = time_axis[1:]

    n = len(time_axis)
    rect_win = np.ones(n)
    # Gaussian window: std ≈ 3 TRs (parameter controls sharpness)
    alpha = 3*TR
    std = (n - 1) / (2 * alpha)
    gauss_win = sig_windows.gaussian(n, std=std)

    window = np.convolve(rect_win, gauss_win, mode='same')
    window /= window.sum()
    return window


def _luppi21_entropy(data, window, n_kmeans, random_state):
    """
    Core Luppi 2021 temporal entropy computation for a single session.

    Parameters
    ----------
    data : np.ndarray, shape (T, N)
    window : np.ndarray, 1D tapered window
    n_kmeans : int
    random_state : int

    Returns
    -------
    entropy : float
    """
    T, N = data.shape
    win_size_half = (len(window) - 1) // 2
    end_pt = T - win_size_half
    midpoints = np.arange(win_size_half, end_pt)
    n_windows = len(midpoints)

    if n_windows < 2:
        return np.nan

    # Cartographic profile histogram bins
    x_bins = np.linspace(0, 1, 101)     # participation coef [0,1]
    y_bins = np.linspace(-5, 5, 101)    # module degree z-score [-5,5]
    xN = len(x_bins)
    yN = len(y_bins)

    CP = np.zeros((yN, xN, n_windows))

    for i, mid in enumerate(midpoints):
        print('Running window', i + 1, 'of', n_windows, end='\r')
        full_win = np.zeros(T)
        start = mid - win_size_half
        end = mid + win_size_half + 1
        full_win[start:end] = window[:end - start]

        # Weighted windowed data
        data_win = data * full_win[:, np.newaxis]

        # Remove all-zero rows (outside window)
        nonzero_rows = np.any(data_win != 0, axis=1)
        data_win_trim = data_win[nonzero_rows]

        if data_win_trim.shape[0] < 2:
            continue

        dFC = np.corrcoef(data_win_trim.T)
        dFC = np.nan_to_num(dFC)

        # Louvain community detection (100 iterations, keep best)
        best_q = -np.inf
        best_ci = None
        for _ in range(100):
            ci_init = np.arange(1, N + 1)
            ci, q = bct.community_louvain(dFC, gamma=1.0, ci=ci_init, B='negative_asym')
            if q > best_q:
                best_q = q
                best_ci = ci.copy()

        # Module degree z-score and participation coefficient
        WT = bct.module_degree_zscore(dFC, best_ci, flag=0)
        BT = bct.participation_coef_sign(dFC, best_ci)[0]

        # Bin into 2D CP histogram
        xi = np.round(np.interp(BT, x_bins, np.arange(xN))).astype(int)
        yi = np.round(np.interp(WT, y_bins, np.arange(yN))).astype(int)
        xi = np.clip(xi, 0, xN - 1)
        yi = np.clip(yi, 0, yN - 1)

        for n in range(N):
            CP[yi[n], xi[n], i] += 1

    # K-means on flattened CPs (each window is a feature vector)
    CP_flat = CP.reshape(xN * yN, n_windows).T  # (n_windows, xN*yN)

    best_inertia = np.inf
    best_labels = None
    for seed in range(n_kmeans):
        km = KMeans(n_clusters=2, n_init=1, random_state=random_state + seed)
        km.fit(CP_flat)
        if km.inertia_ < best_inertia:
            best_inertia = km.inertia_
            best_labels = km.labels_

    # Shannon entropy of the 2-state sequence
    prob, _ = np.histogram(best_labels, bins=[-0.5, 0.5, 1.5], density=False)
    prob = prob / prob.sum()
    mask = prob > 0
    H = -np.sum(prob[mask] * np.log2(prob[mask]))
    return float(H)
