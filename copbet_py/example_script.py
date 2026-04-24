"""
CopBET_py Example Script
=========================
Demonstrates all 13 entropy measures in CopBET_py on simulated fMRI data.

Simulated data:
  - 5 sessions
  - 300 timepoints x 200 ROIs each
  - AR(1) structure with moderate spatial correlation

Produces:
  - A summary bar chart of scalar entropy values
  - Detailed subplots for measures that return arrays (e.g., per scale, per degree)
  - DCC correlation matrix heatmap

Run from the CopBET_py directory:
    python example_script.py
"""

import sys
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob 

# Add CopBET_py to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions import *

warnings.filterwarnings('ignore')
N_SESSIONS = 2

datapath = '/mrdata/np2/p3/denoised_parcellated_np2p3/parcellated_timeseries/acompcor_spikeregression_lowpass_gsr/'
sessions_glob = glob.glob(datapath+'sub-53888/ses-PSI/func/*task-rest*_bold_parcellated_schaefer200.txt')
TR = 2.0
sessions = []
for i in range(N_SESSIONS):
    sessions.append(np.loadtxt(sessions_glob[i]))
    print(f"Loaded session {i+1} with shape {sessions[i].shape}")
N = sessions[0].shape[1]  # number of ROIs
# ─────────────────────────────────────────────────────────────────
# 2. RUN ALL ENTROPY MEASURES
# ─────────────────────────────────────────────────────────────────

print("Running entropy measures:\n")

# ---- Von Neumann entropy (scalar per session) ---- OK
out = CopBET_von_Neumann_entropy(sessions)

# ---- Time series complexity (scalar per session) ----
out = CopBET_time_series_complexity(sessions, LZtype='LZ76temporal')

out = CopBET_time_series_complexity(sessions, LZtype='LZ78temporal') #OK

# ---- Metastate series complexity (scalar per session) ----
out = CopBET_metastate_series_complexity(sessions, n_replicates=20)

# ---- LEiDA transition entropy (scalar per session) ----
out = CopBET_LEiDA_transition_entropy(sessions, K=5, n_replicates=20)

# ---- NGSC (scalar per session) ----
# Use 8 pseudo-parcels for regional computation
parcel_labels = np.repeat(np.arange(1, 9), N // 8 + 1)[:N]
out = CopBET_NGSC(sessions, parcel_labels=parcel_labels)

# ---- Sample entropy (5 scales × 200 ROIs per session — use subset for speed) ----
# Use 10 ROIs only to keep example fast
small_sessions = [ts[:, :10] for ts in sessions]
out = CopBET_sample_entropy(sessions, scales=[1, 2, 3])

# ---- Intranetwork synchrony ----
# supposed to be run on voxelwise data
out = CopBET_intranetwork_synchrony(sessions, parcel_labels=parcel_labels)

# ---- Temporal entropy (use 10-ROI subset — full is slow) ----
#### temporal entropy is difficult to check for correctness
out = CopBET_temporal_entropy(small_sessions, TR=TR, n_kmeans=10)

# ---- Degree distribution entropy (use 20-ROI subset for speed) ----
# small20 = [ts[:, :20] for ts in sessions]
# thresholds_fast = np.arange(0.1, 1.0, 0.05)
# degrees_fast = np.arange(1, 16)
out = CopBET_degree_distribution_entropy(small_sessions, thresholds=None, chosen_degrees=None)

# ---- Geodesic entropy (use 20-ROI subset) ----
out = CopBET_geodesic_entropy(small_sessions, thresholds=None, chosen_degrees=None)

# ---- Diversity coefficient (use 30-ROI subset, fewer iterations) ----
# small30 = [ts[:, :30] for ts in sessions]
out = CopBET_diversity_coefficient(small_sessions)

# # ---- Motif connectivity entropy (4 ROIs, subset of windows) ----
motif_sessions = [ts[:, :4] for ts in sessions]
motion = np.random.rand(len(sessions))  # placeholder motion values
out = CopBET_motif_connectivity_entropy(motif_sessions, motion=motion, TR=TR, window_lengths_sec=np.arange(20, 61, 5))

