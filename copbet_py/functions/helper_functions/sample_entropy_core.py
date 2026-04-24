"""
Sample entropy and multiscale sample entropy.

Translated from:
  - sample_entropy.m (LOFT Complexity Toolbox) included in CopBET
  - Used in CopBET_sample_entropy.py (Lebedev et al., 2016)
"""
import numpy as np


def sample_entropy(u, m=2, r=None, scale=1):
    """
    Compute sample entropy of a 1D time series.

    Parameters
    ----------
    u : array-like, shape (T,)
        1D time series.
    m : int
        Template length (embedding dimension). Default 2.
    r : float or None
        Tolerance. If None, defaults to 0.3 * std(u).
    scale : int
        Coarse-graining scale for multiscale entropy.
        scale=1 gives standard sample entropy.

    Returns
    -------
    sampEn : float
        Sample entropy value (0 if undefined).
    """
    u = np.asarray(u, dtype=float).ravel()

    if r is None:
        r = 0.3 * np.std(u,ddof=1)

    # Coarse-graining for multiscale entropy
    if scale > 1:
        n = len(u) // scale
        uu = np.array([np.mean(u[i * scale:(i + 1) * scale]) for i in range(n)])
    else:
        n = len(u)
        uu = u.copy()

    if n <= m + 1:
        return 0.0

    # Build template matrix: each column is a window of length m+1
    # data_mat[i, j] = uu[j + i] for i in 0..m, j in 0..n-m-1
    data_mat = np.zeros((m + 1, n - m))
    for i in range(m + 1):
        data_mat[i, :] = uu[i:n - m + i]

    cor = np.zeros(2)
    for d in range(m, m + 2):
        count = np.zeros(n - m)
        temp_mat = data_mat[:d, :]
        for i in range(n - d):
            # Chebyshev (max) distance between template i and all others
            dist = np.max(np.abs(temp_mat[:, i + 1:n - m] - temp_mat[:, i:i + 1]), axis=0)
            count[i] = np.sum(dist < r) / (n - m - 1)
        cor[d - m] = np.sum(count) / (n - m)

    if cor[0] > 0 and cor[1] > 0:
        sampEn = -np.log(cor[1] / cor[0])
    else:
        sampEn = 0.0

    return sampEn


def multiscale_sample_entropy(u, m=2, r=None, scales=None):
    """
    Compute multiscale sample entropy across a range of scales.

    Parameters
    ----------
    u : array-like, shape (T,)
        1D time series.
    m : int
        Template length.
    r : float or None
        Tolerance. If None, defaults to 0.3 * std(u).
    scales : list of int or None
        Scales to compute. Defaults to [1, 2, 3, 4, 5].

    Returns
    -------
    mse : np.ndarray, shape (n_scales,)
        Sample entropy at each scale.
    """
    if scales is None:
        scales = [1, 2, 3, 4, 5]

    u = np.asarray(u, dtype=float).ravel()
    if r is None:
        r = 0.3 * np.std(u)

    mse = np.zeros(len(scales))
    for i, scale in enumerate(scales):
        mse[i] = sample_entropy(u, m=m, r=r, scale=scale)

    return mse


if __name__ == "__main__":
    import glob
    datapath = '/mrdata/np2/p3/denoised_parcellated_np2p3/parcellated_timeseries/acompcor_spikeregression_lowpass_gsr/'
    sessions_glob = glob.glob(datapath+'sub-53888/ses-PSI/func/*task-rest*_bold_parcellated_schaefer200.txt')
    TR = 2.0
    sessions = []
    for i in range(1):
        sessions.append(np.loadtxt(sessions_glob[i]))
        print(f"Loaded session {i+1} with shape {sessions[i].shape}")
    N = sessions[0].shape[1]  # number of ROIs

    x = sessions[0][:, 0]  # example ROI time series
    se = sample_entropy(x, scale=2)
    print(se)