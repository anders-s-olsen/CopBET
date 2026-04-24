"""
CopBET_DCC_entropy
===================
Copenhagen Brain Entropy Toolbox: Dynamic Conditional Correlation entropy.

Evaluates DCC entropy as in Barrett et al., 2020.
A DCC(1,1)-GARCH(1,1) model is fit to each session's demeaned time series.
This produces a time-varying correlation matrix per volume. For each ROI
pair, the Shannon entropy of the histogram of correlation coefficients
across volumes is computed. A variance matrix is also returned.

NOTE: This function is computationally expensive. For large datasets or
many ROIs consider reducing n_regions or using precomputed DCC matrices.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    Demeaned ROI time series per session.

Returns
-------
results : list of dict with keys:
    'entropy'  : np.ndarray, shape (N, N) — Shannon entropy per ROI pair
    'variance' : np.ndarray, shape (N, N) — variance of DCC values per pair
    'Ct'       : np.ndarray, shape (N, N, T) — time-varying correlation matrices

Reference
---------
Barrett et al., 2020.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import sys
import os
import numpy as np

# Add external DCC-GARCH package to path
_EXT_PATH = os.path.join(os.path.dirname(__file__), '..', 'external', 'DCC_GARCH')
sys.path.insert(0, os.path.abspath(_EXT_PATH))

from DCC_GARCH.GARCH.GARCH import GARCH
from DCC_GARCH.GARCH.GARCH_loss import garch_loss_gen
from DCC_GARCH.DCC.DCC import DCC
from DCC_GARCH.DCC.DCC_loss import R_gen


def CopBET_DCC_entropy(files, n_bins=None):
    """
    Compute DCC entropy for each session.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    n_bins : int or None
        Number of histogram bins for entropy estimation.
        If None, uses Sturges' rule (ceil(log2(T)) + 1).

    Returns
    -------
    results : list of dict
    """
    results = []

    for ses, file in enumerate(files):
        DCC_file = file.replace('parcellated_timeseries','DCC_tmp').replace('.txt','_DCC.txt')
        Ct = np.loadtxt(DCC_file)
        N,_,T = Ct.shape
        entropy_out = np.zeros((N, N))

        if n_bins is None:
            n_bins = int(np.ceil(np.log2(T))) + 1

        for i in range(N):
            for j in range(i, N):
                hist, _ = np.histogram(Ct[i, j, :], bins=n_bins,
                                       range=(-1, 1), density=False)
                prob = hist / hist.sum()
                mask = prob > 0
                entropy_out[i, j] = -np.sum(prob[mask] * np.log(prob[mask]))

        # Symmetrize
        entropy_out = entropy_out + entropy_out.T

        results.append(entropy_out)

    return results
