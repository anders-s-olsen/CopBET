"""
CopBET_NGSC
============
Copenhagen Brain Entropy Toolbox: Normalized Geodesic Spectral Clustering
(NGSC) entropy.

Calculates NGSC as defined by Siegel et al., 2024. PCA is performed on
the time series of each ROI (or the whole brain). The normalized entropy
of the eigenvalue spectrum is returned:

    H = -sum(lambda_i * log(lambda_i)) / log(m)

where m is the number of non-zero eigenvalues.

This Python version accepts 2D ROI time series (T x N). The original
MATLAB code operates on voxelwise NIfTI data and computes both whole-brain
and per-parcel NGSC. Here the function computes:
  - Whole-brain NGSC on the full (T x N) matrix.
  - Per-parcel NGSC if a parcel_labels vector is provided.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
parcel_labels : array-like of int or None
    Label for each column (ROI/voxel). If None, only whole-brain NGSC
    is computed.

Returns
-------
results : list of dict with keys:
    'global' : float — whole-brain NGSC entropy
    'regional' : np.ndarray or None — per-parcel NGSC entropy

Reference
---------
Siegel et al., 2024.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np


def CopBET_NGSC(ts, parcel_labels=None):
    """
    Compute NGSC entropy.

    Parameters
    ----------
    ts : np.ndarray, each shape (T, N)
    parcel_labels : array-like of int or None
        Integer parcel label for each of the N columns. Labels of 0 are
        excluded. If None, no regional NGSC is computed.

    Returns
    -------
    results : dict with keys 'global' (float) and 'regional' (array or None).
    """
    if np.any(parcel_labels == 0):
        # warn if there are zero labels, since these will be ignored in regional NGSC
        print("Warning: parcel_labels contains zeros, which will be ignored in regional NGSC.")

    ts = np.asarray(ts, dtype=float)

    # # Standardize (zero mean, unit variance per ROI)
    # mu = ts.mean(axis=0)
    # sd = ts.std(axis=0)
    # sd[sd < 1e-10] = 1e-10
    # ts_std = (ts - mu) / sd
    ts_std = ts - ts.mean(axis=0) #don't divide by std

    # Whole-brain NGSC
    global_ngsc = _ngsc_entropy(ts_std)

    # Regional NGSC
    if parcel_labels is not None:
        labels = np.asarray(parcel_labels)
        unique_labels = np.unique(labels[labels != 0])
        regional_ngsc = np.zeros(len(unique_labels))
        for li, lbl in enumerate(unique_labels):
            region_data = ts_std[:, labels == lbl]
            regional_ngsc[li] = _ngsc_entropy(region_data)
    else:
        regional_ngsc = None

    results = {'global': global_ngsc, 'regional': regional_ngsc}

    return results


def _ngsc_entropy(data):
    """
    Compute normalized spectral entropy of the PCA eigenvalue spectrum.

    Parameters
    ----------
    data : np.ndarray, shape (T, K)

    Returns
    -------
    H : float
    """
    T, K = data.shape
    if K == 0 or T < 2:
        return np.nan

    # PCA via SVD of data matrix (equivalent to eigendecomposition of cov)
    # # We need eigenvalues of the covariance matrix
    # cov = np.cov(data.T)  # (K, K)
    # if cov.ndim == 0:
    #     # Single column: variance is the only eigenvalue
    #     eigvals = np.array([float(cov)])
    # else:
    #     eigvals = np.linalg.eigvalsh(cov)
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    eigvals = (S ** 2) / (T - 1)  # eigenvalues of covariance matrix

    # Keep only positive eigenvalues
    eigvals = eigvals[eigvals > 0]
    m = len(eigvals)

    if m == 0 or m == 1:
        return 0.0

    # Normalize
    lam = eigvals / eigvals.sum()
    lam[lam == 0] = np.finfo(float).eps

    H = -np.sum(lam * np.log(lam)) / np.log(m)
    return float(H)
