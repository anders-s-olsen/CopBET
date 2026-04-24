"""
CopBET_sample_entropy
======================
Copenhagen Brain Entropy Toolbox: Sample entropy.

Evaluates multiscale sample entropy as in Lebedev et al., 2016.
For each session, sample entropy is computed across 5 coarse-graining
scales (m=2, r=0.3*std per region). This function accepts ROI time
series directly (rather than voxelwise NIfTI files as in the original).

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    ROI time series per session.

Returns
-------
entropy : list of np.ndarray, each shape (n_scales, N)
    Multi-scale sample entropy per scale and ROI per session.

Reference
---------
Lebedev et al., 2016.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from .helper_functions.sample_entropy_core import sample_entropy


def CopBET_sample_entropy(sessions, m=2, r_factor=0.3, scales=None):
    """
    Compute multiscale sample entropy per ROI.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    m : int
        Template length (default 2).
    r_factor : float
        Tolerance = r_factor * std(roi_ts). Default 0.3.
    scales : list of int or None
        Coarse-graining scales. Default [1, 2, 3, 4, 5].

    Returns
    -------
    entropy : list of np.ndarray, each shape (n_scales, N)
    """
    if scales is None:
        scales = [1, 2, 3, 4, 5]

    entropy = []

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)
        T, N = ts.shape
        ts = ts - ts.mean(axis=0)

        mse = np.zeros((len(scales), N))

        for roi in range(N):
            roi_ts = ts[:, roi]
            r = r_factor * np.std(roi_ts, ddof=1)
            for si, scale in enumerate(scales):
                mse[si, roi] = sample_entropy(roi_ts, m=m, r=r, scale=scale)

        entropy.append(mse)

    return entropy
