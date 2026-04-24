"""
CopBET_von_Neumann_entropy
==========================
Copenhagen Brain Entropy Toolbox: Von-Neumann entropy.

Evaluates Von-Neumann entropy as introduced in Felippe et al., 2021.
The Von-Neumann entropy is computed from the eigenvalues of the
density matrix rho = R / N, where R is the Pearson correlation matrix
and N is the number of regions.

Input
-----
sessions : list of np.ndarray, each shape (T, N)
    Time series for each session (T timepoints, N regions).

Returns
-------
entropy : np.ndarray, shape (n_sessions,)
    Normalized Von-Neumann entropy per session.

Reference
---------
Felippe et al., 2021.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np


def CopBET_von_Neumann_entropy(sessions):
    """
    Compute Von-Neumann entropy for a list of fMRI sessions.

    Parameters
    ----------
    sessions : list of np.ndarray, each shape (T, N)

    Returns
    -------
    entropy : np.ndarray, shape (n_sessions,)
    """
    entropy = np.full(len(sessions), np.nan)

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts, dtype=float)
        N = ts.shape[1]

        R = np.corrcoef(ts.T)  # (N, N) correlation matrix

        # Density matrix: rho = R / N
        rho = R / N

        # Eigenvalues of rho
        eigvals = np.linalg.eigvalsh(rho)

        # # Clamp tiny negatives from floating point
        # eigvals = np.where(eigvals < 0, np.finfo(float).eps, eigvals)

        # Von-Neumann entropy: S = -Tr(rho * log(rho)) = -sum(lambda * log(lambda))
        # Avoid log(0) by masking near-zero eigenvalues
        mask = eigvals > 0
        H = -np.sum(eigvals[mask] * np.log(eigvals[mask]))

        # Normalize by maximum entropy log(N)
        entropy[ses] = H / np.log(N)

    return entropy
