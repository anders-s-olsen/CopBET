"""
CopBET_intranetwork_synchrony
==============================
Copenhagen Brain Entropy Toolbox: Intranetwork synchrony.

Evaluates intranetwork synchrony as in Carhart-Harris et al., 2014.
For each time point and each ROI, the variance across voxels within the
ROI is computed. The Shannon entropy of the variance distribution over
time is returned per ROI.

This Python version accepts per-ROI voxelwise data rather than 4D NIfTI
files. Input is a list of sessions; each session is a list of 2D arrays
(one per ROI), where each array has shape (n_voxels_roi, T).

If ROI-level (single voxel per ROI) time series are provided, the
function automatically adds a small amount of noise to simulate voxel
spread, as this measure is only meaningful with multiple voxels per ROI.

Input
-----
sessions : list of list of np.ndarray
    sessions[i][roi] is shape (n_voxels, T) for ROI `roi` in session `i`.
    OR: list of np.ndarray each shape (T, N) — in this case each "voxel"
    is a single ROI signal and the function returns zeros (see note).

Returns
-------
entropy : list of np.ndarray, each shape (N_rois,)
    Per-ROI intranetwork synchrony entropy per session.

Reference
---------
Carhart-Harris et al., 2014.

Please cite McCulloch, Olsen et al., 2023 if you use CopBET.
"""

import numpy as np


def CopBET_intranetwork_synchrony(sessions, parcel_labels):
    """
    Compute intranetwork synchrony entropy.

    Parameters
    ----------
    sessions : list of (list of np.ndarray | np.ndarray)
        Each session is either:
        - A list of ROI arrays, where sessions[i][roi].shape = (n_voxels, T)
        - A 2D array of shape (T, N) — treated as N single-voxel ROIs

    Returns
    -------
    entropy : list of np.ndarray, each shape (N_rois,)
    """
    entropy = []

    labels = np.asarray(parcel_labels)
    unique_labels = np.unique(labels[labels != 0])
    N_rois = len(unique_labels)

    for ses, ts in enumerate(sessions):
        print(f"Processing session {ses + 1}/{len(sessions)}...", end='\r')
        ts = np.asarray(ts)

        T, N = ts.shape

        # Demean each voxel's time series
        ts = ts - ts.mean(axis=0)

        # Compute per-ROI, per-timepoint variance across voxels
        variances = np.zeros((T,N_rois))
        for i, roi in enumerate(unique_labels):
            roi_mask = labels == roi
            if np.sum(roi_mask) == 0:
                raise ValueError(f"ROI label {roi} has no voxels.")
            roi_ts = ts[:, roi_mask]  # (T, n_voxels_in_roi)
            variances[:,i] = np.nanvar(roi_ts, axis=1)  # (T,)

        # Shannon entropy of the variance distribution per ROI
        ent_roi = np.zeros(N_rois)
        for roi in range(N_rois):
            v = variances[:,roi]
            if np.all(v == 0):
                ent_roi[roi] = 0.0
                continue
            counts, _ = np.histogram(v, bins='auto')
            prob = counts / counts.sum()
            mask = prob > 0
            ent_roi[roi] = -np.sum(prob[mask] * np.log(prob[mask]))

        entropy.append(ent_roi)

    return entropy
