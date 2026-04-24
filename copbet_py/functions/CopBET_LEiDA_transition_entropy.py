"""
CopBET_LEiDA_transition_entropy
================================
Copenhagen Brain Entropy Toolbox: LEiDA transition entropy.

Evaluates the leading eigenvectors of instantaneous phase-coherence matrices
(via Hilbert transform), clusters them with K-means (squared Euclidean
distance), and computes the Markov transition rate entropy from the
resulting state sequence, following Kringelbach et al., 2020.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    ROI time series per session.
K : int
    Number of clusters for K-means.

Returns
-------
entropy : np.ndarray, shape (n_sessions,)
    Markov entropy per session (normalized by log2(K)).

Reference
---------
Kringelbach et al., 2020. "Dynamic coupling of whole-brain neuronal and
neurotransmitter systems."

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from scipy.signal import hilbert
from sklearn.cluster import KMeans


def CopBET_LEiDA_transition_entropy(sessions, K=3, n_replicates=200, random_state=42):
    """
    Compute LEiDA Markov transition entropy.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    K : int
        Number of LEiDA states.
    n_replicates : int
        Number of K-means restarts.
    random_state : int

    Returns
    -------
    entropy : np.ndarray, shape (n_sessions,)
    """
    if K < 2:
        raise ValueError("K must be at least 2 for meaningful transition entropy.")
    # Compute leading eigenvectors for all sessions
    datasizes = []
    leading_eigs_list = []

    for ts in sessions:
        ts = np.asarray(ts, dtype=float)
        ts = ts - ts.mean(axis=0)
        eigs = _LEiDA(ts)          # (T, N)
        leading_eigs_list.append(eigs)
        datasizes.append(eigs.shape[0])

    Lead_eigs_all = np.vstack(leading_eigs_list)  # (T_total, N)

    # K-means on concatenated leading eigenvectors
    best_inertia = np.inf
    best_labels = None

    for seed in range(n_replicates):
        print(f"K-means replicate {seed + 1}/{n_replicates}...", end='\r')
        km = KMeans(n_clusters=K, n_init=1, max_iter=1000,
                    random_state=random_state + seed)
        km.fit(Lead_eigs_all)
        if km.inertia_ < best_inertia:
            best_inertia = km.inertia_
            best_labels = km.labels_

    # Compute Markov transition entropy per session
    entropy = np.full(len(sessions), np.nan)
    idx = 0
    for ses, size in enumerate(datasizes):
        seg = best_labels[idx:idx + size]
        P = _transition_matrix(seg, K)
        entropy[ses] = _markov_entropy(P, K)
        idx += size

    return entropy


def _LEiDA(ts):
    """
    Compute leading eigenvectors of instantaneous phase-coherence matrices.

    Parameters
    ----------
    ts : np.ndarray, shape (T, N)

    Returns
    -------
    leading_eigs : np.ndarray, shape (T, N)
    """
    T, N = ts.shape

    # Hilbert phase per region
    phase = np.angle(hilbert(ts, axis=0))  # (T, N)

    # via svd
    phase_stack = np.stack((np.cos(phase), np.sin(phase)), axis=-1)
    U, S, Vt = np.linalg.svd(phase_stack, full_matrices=False)
    leading_eigs = U[:,:,0]  # (N,)

    return leading_eigs


def _transition_matrix(labels, K):
    """
    Compute K×K Markov transition probability matrix from a state sequence.

    Parameters
    ----------
    labels : array-like, shape (T,)
        Integer state labels 0..K-1.
    K : int

    Returns
    -------
    P : np.ndarray, shape (K, K)
    """
    P = np.zeros((K, K))
    for t in range(len(labels) - 1):
        P[labels[t], labels[t + 1]] += 1

    # Normalize rows
    row_sums = P.sum(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        P = np.where(row_sums > 0, P / row_sums, 0.0)

    return P


def _markov_entropy(P, K):
    """
    Compute Markov entropy rate from a transition matrix.

    Uses the stationary distribution as the weighting over rows.

    Parameters
    ----------
    P : np.ndarray, shape (K, K)
    K : int

    Returns
    -------
    H : float
        Normalized Markov entropy rate.
    """
    # Stationary distribution via left eigenvector (eigenvalue = 1)
    eigvals, eigvecs = np.linalg.eig(P.T)
    # Find eigenvector corresponding to eigenvalue closest to 1
    idx = np.argmin(np.abs(eigvals - 1.0))
    stat = np.abs(eigvecs[:, idx])
    stat = stat / stat.sum()

    eps = np.finfo(float).eps
    Hi = np.zeros(K)
    for row in range(K):
        row_p = P[row, :] + eps
        Hi[row] = -stat[row] * np.sum(row_p * np.log2(row_p))

    H = np.sum(Hi) / np.log2(K)

    return float(H)
