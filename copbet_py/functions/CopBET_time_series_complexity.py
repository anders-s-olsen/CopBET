"""
CopBET_time_series_complexity
==============================
Copenhagen Brain Entropy Toolbox: Time series complexity.

Evaluates time-series complexity as in Varley et al., 2020.
The input data is Hilbert-transformed; the amplitude envelope is binarized
around its mean. The resulting binary matrix is flattened (spatially or
temporally) and Lempel-Ziv complexity is computed (LZ76 or LZ78), then
normalized by the complexity of a random permutation of the same series.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    ROI time series per session.
LZtype : str
    One of:
      'LZ78temporal'  - LZ78 flattened across time first (T x N)
      'LZ78spatial'   - LZ78 flattened across space first (N x T)
      'LZ76temporal'  - LZ76 flattened across time first
      'LZ76spatial'   - LZ76 flattened across space first

Returns
-------
entropy : np.ndarray, shape (n_sessions,)
    Normalized LZ complexity per session.

Reference
---------
Varley et al., 2020.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np
from scipy.signal import hilbert
from .helper_functions.lempel_ziv import lz76_complexity, lz78_complexity


def CopBET_time_series_complexity(sessions, LZtype='LZ78regional'):
    """
    Compute LZ complexity of Hilbert-transformed, binarized time series.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)
    LZtype : str
        One of 'LZ78temporal', 'LZ78spatial', 'LZ76temporal', 'LZ76spatial', 'LZ76regional','LZ78regional'.

    Returns
    -------
    entropy : np.ndarray, shape (n_sessions,)
    """
    valid_types = ('LZ78temporal', 'LZ78spatial', 'LZ76temporal', 'LZ76spatial', 'LZ76regional', 'LZ78regional')
    if LZtype not in valid_types:
        raise ValueError(f"LZtype must be one of {valid_types}, got '{LZtype}'")

    if 'regional' in LZtype:
        entropy = np.full((len(sessions),sessions[0].shape[1]), np.nan)
    else:
        entropy = np.full(len(sessions), np.nan)

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)  # (T, N)

        # Hilbert transform: take amplitude envelope
        abs_hts = np.abs(hilbert(ts, axis=0))  # (T, N)

        # Binarize: 1 where amplitude > column mean
        mean_bold = np.mean(abs_hts, axis=0)   # (N,)
        bin_abs_hts = (abs_hts > mean_bold).astype(bool)  # (T, N)

        # Random surrogate: shuffle each region independently
        random_ts = np.zeros_like(abs_hts)
        rng = np.random.default_rng()
        for roi in range(ts.shape[1]):
            perm = rng.permutation(ts.shape[0])
            # reverse the array
            # perm = np.arange(ts.shape[0])[::-1]
            random_ts[:, roi] = abs_hts[perm, roi]
        M_rand = (random_ts > mean_bold).astype(bool)

        # Flatten according to LZtype
        if 'regional' in LZtype:
            for roi in range(ts.shape[1]):
                # Compute LZ complexity
                if 'LZ78' in LZtype:
                    C = lz78_complexity(bin_abs_hts[:, roi])
                    C_rand = lz78_complexity(M_rand[:, roi])
                else:  # LZ76
                    C = lz76_complexity(bin_abs_hts[:, roi], normalize=True)
                    C_rand = lz76_complexity(M_rand[:, roi], normalize=True)

                entropy[ses,roi] = C / C_rand if C_rand > 0 else np.nan
        else:
            if 'temporal' in LZtype:
                long_ts = bin_abs_hts.ravel(order='F')    # column-major: time varies fastest
                long_rand = M_rand.ravel(order='F')
            elif 'spatial' in LZtype:  # spatial
                long_ts = bin_abs_hts.T.ravel(order='F')  # N x T, then flatten
                long_rand = M_rand.T.ravel(order='F')

            # Compute LZ complexity
            if 'LZ78' in LZtype:
                C = lz78_complexity(long_ts)
                C_rand = lz78_complexity(long_rand)
            else:  # LZ76
                C = lz76_complexity(long_ts, normalize=True)
                C_rand = lz76_complexity(long_rand, normalize=True)

            entropy[ses] = C / C_rand if C_rand > 0 else np.nan

    return entropy
