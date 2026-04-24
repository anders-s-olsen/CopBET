"""
CopBET_metastate_series_complexity
====================================
Copenhagen Brain Entropy Toolbox: Metastate series complexity.

Evaluates Lempel-Ziv complexity (LZ76 exhaustive) for a binary series of
metastates. K-means (K=4) is run on all sessions concatenated, using
correlation distance. The four centroids are paired into two opposing
metastates (most anti-correlated pairs). The binary metastate activation
sequence is then subjected to LZ76 complexity.

NOTE: requires multiple sessions to be meaningful (clustering is run
across the pooled dataset).

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    ROI time series per session.

Returns
-------
entropy : np.ndarray, shape (n_sessions,)
    Normalized LZ76 complexity of the metastate series per session.

Reference
---------
Timmermann et al., 2019 / Lebedev et al., 2016 (metastate approach).

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from sklearn.cluster import KMeans
from .helper_functions.lempel_ziv import lz76_complexity
from .phase_coherence_kmeans import diametrical_clustering
from scipy.signal import butter, filtfilt
from scipy.signal import hilbert

def CopBET_metastate_series_complexity(sessions, n_replicates=200, random_state=42, function='kmeans', tr=None):
    """
    Compute metastate series LZ76 complexity.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
        Must have at least 2 sessions.
    n_replicates : int
        Number of K-means restarts.
    random_state : int
        Random seed.
    function : str
        Clustering function to use. Currently only 'kmeans' or 'complex_diametrical_clustering' is implemented.
    tr : list of float or None
        Repetition times for each session, required if function='complex_diametrical_clustering

    Returns
    -------
    entropy : np.ndarray, shape (n_sessions,)
    """
    if len(sessions) < 2:
        raise ValueError("Metastate complexity requires at least 2 sessions.")

    # Demean and concatenate across sessions
    datasizes = []
    data_list = []
    for ts in sessions:
        ts = np.asarray(ts, dtype=float)
        ts = ts - ts.mean(axis=0)
        data_list.append(ts)

    if function == 'kmeans':
        entropy = concatenate_data_run_kmeans(data_list, n_replicates, random_state)
    elif function == 'complex_diametrical_clustering':
        if tr is None:
            raise ValueError("TRs must be provided for complex diametrical clustering.")
        entropy = run_complex_diametrical_clustering_per_session(data_list, n_replicates, random_state, tr)

    return entropy

def run_complex_diametrical_clustering_per_session(data_list, n_replicates, random_state, tr):
    entropy = np.full(len(data_list), np.nan)
    for ses, ts in enumerate(data_list):
        if np.any(np.isnan(ts)):
            print(f"Warning: NaN values found in session {ses+1}. Skipping this session.")
            continue
        print(f"Processing session {ses + 1}/{len(data_list)}...", end='\r')

        # filter ts between 0.03 and 0.07 Hz
        fs = 1 / tr[ses]
        nyq = fs / 2
        low = 0.03 / nyq
        high = 0.07 / nyq
        b, a = butter(4, [low, high], btype='band')
        ts_filt = filtfilt(b, a, ts, axis=0)

        # do hilbert and extract phase
        ts_hilbert = hilbert(ts_filt, axis=0)
        ts_phase = np.angle(ts_hilbert)

        u = 1/np.sqrt(ts_phase.shape[1]) * np.exp(1j * ts_phase)  # (T, N)

        _,seg,_ = diametrical_clustering(u, K=2, num_repl=1)
        entropy[ses] = lz76_complexity(seg, normalize=True)
    return entropy

def concatenate_data_run_kmeans(data_list, n_replicates, random_state):
    datasizes = []
    for ts in data_list:
        datasizes.append(ts.shape[0])

    data_all = np.vstack(data_list)  # (T_total, N)

    # K-means with correlation distance (= 1 - Pearson correlation)
    # sklearn doesn't have correlation distance directly, but we can normalize
    # rows to unit norm, then cosine distance ≈ 1 - correlation
    data_norm = data_all / (np.linalg.norm(data_all, axis=1, keepdims=True) + 1e-10)

    best_inertia = np.inf
    best_partition = None
    best_centroids = None

    for seed in range(n_replicates):
        print(f"K-means replicate {seed+1}/{n_replicates}...", end='\r')
        km = KMeans(n_clusters=4, n_init=1, max_iter=1000, random_state=random_state + seed)
        km.fit(data_norm)
        if km.inertia_ < best_inertia:
            best_inertia = km.inertia_
            best_partition = km.labels_
            best_centroids = km.cluster_centers_

    # Group states into 2 metastates based on centroid anti-correlation
    cen_corr = np.corrcoef(best_centroids)
    min_idx = []
    for k in range(4):
        min_idx.append(np.argmin(cen_corr[k]))
    for k in range(4):
        if min_idx[min_idx[k]] != k:
            raise ValueError("Centroid correlation structure is not symmetric. Check K-means results.")

    metastates = np.zeros(4, dtype=int)
    assigned = set()
    meta_id = 1
    for k in range(4):
        if k not in assigned:
            partner = min_idx[k]
            metastates[k] = meta_id
            metastates[partner] = meta_id
            assigned.add(k)
            assigned.add(partner)
            meta_id += 1

    # Binary metastate series (False = metastate 0, True = metastate 1)
    binary_partition = np.array([metastates[label] == 1 for label in best_partition])

    # Compute LZ76 per session
    entropy = np.full(len(data_list), np.nan)
    idx = 0
    for ses, size in enumerate(datasizes):
        seg = binary_partition[idx:idx + size]
        entropy[ses] = lz76_complexity(seg, normalize=True)
        idx += size
    return entropy