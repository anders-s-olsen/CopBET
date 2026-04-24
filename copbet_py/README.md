# CopBET_py

Python translation of the **Copenhagen Brain Entropy Toolbox (CopBET)**,
originally published in MATLAB at:

> https://github.com/anders-s-olsen/CopBET

Please cite:

> McCulloch, Olsen et al., 2023. "Navigating Chaos in Psychedelic Neuroimaging:
> A Rigorous Empirical Evaluation of the Entropic Brain Hypothesis."

---

## Overview

CopBET implements 13 entropy/complexity measures for resting-state fMRI data.
All functions take a **list of numpy arrays**, one per session, each of shape
**(T × N)** — T timepoints × N ROIs — and return entropy values as numpy arrays
or lists.

---

## Measures implemented

| Function | Method | Reference |
|---|---|---|
| `CopBET_von_Neumann_entropy` | Von-Neumann entropy of correlation matrix | Felippe et al., 2021 |
| `CopBET_time_series_complexity` | LZ76/LZ78 on Hilbert-transformed binary series | Varley et al., 2020 |
| `CopBET_metastate_series_complexity` | LZ76 on binary metastate sequence (K=4 k-means) | Timmermann et al., 2019 |
| `CopBET_LEiDA_transition_entropy` | Markov entropy of LEiDA phase-coherence states | Kringelbach et al., 2020 |
| `CopBET_degree_distribution_entropy` | Shannon entropy of degree distribution across thresholds | Viol et al., 2017 |
| `CopBET_geodesic_entropy` | Shannon entropy of geodesic path-length distributions | Viol et al., 2019 |
| `CopBET_diversity_coefficient` | BCT diversity coefficient (Hpos) | Lebedev et al., 2015 |
| `CopBET_temporal_entropy` | Cartographic-profile K=2 state entropy | Luppi et al., 2021 |
| `CopBET_sample_entropy` | Multi-scale sample entropy per ROI | Lebedev et al., 2016 |
| `CopBET_intranetwork_synchrony` | Entropy of voxel-variance distribution per ROI | Carhart-Harris et al., 2014 |
| `CopBET_motif_connectivity_entropy` | 4-ROI partial-correlation motif entropy | Tagliazucchi et al., 2014 |
| `CopBET_NGSC` | Normalized geodesic spectral clustering (PCA entropy) | Siegel et al., 2024 |
| `CopBET_DCC_entropy` | DCC-GARCH time-varying correlation entropy | Barrett et al., 2020 |

---

## Installation

```bash
pip install numpy scipy scikit-learn matplotlib networkx bctpy
```

Then clone or copy the `CopBET_py/` directory and run scripts from within it.

---

## Usage

```python
import numpy as np
from functions import CopBET_von_Neumann_entropy, CopBET_time_series_complexity

# List of sessions: each (T, N) array
sessions = [np.random.randn(300, 200) for _ in range(5)]

# Von-Neumann entropy
vne = CopBET_von_Neumann_entropy(sessions)
print(vne)  # shape (5,)

# LZ76 temporal complexity
lz = CopBET_time_series_complexity(sessions, LZtype='LZ76temporal')
print(lz)   # shape (5,)
```

See `example_script.py` for a full demonstration of all measures with plots.

---

## Run the example

```bash
cd CopBET_py
python example_script.py
```

This simulates 5 sessions of 300 × 200 fMRI data, runs all 13 measures, prints
a summary, and saves `CopBET_py_example_results.png`.

---

## Directory structure

```
CopBET_py/
├── README.md
├── requirements.txt
├── example_script.py          ← run this
├── functions/
│   ├── __init__.py
│   ├── CopBET_DCC_entropy.py
│   ├── CopBET_LEiDA_transition_entropy.py
│   ├── CopBET_NGSC.py
│   ├── CopBET_degree_distribution_entropy.py
│   ├── CopBET_diversity_coefficient.py
│   ├── CopBET_geodesic_entropy.py
│   ├── CopBET_intranetwork_synchrony.py
│   ├── CopBET_metastate_series_complexity.py
│   ├── CopBET_motif_connectivity_entropy.py
│   ├── CopBET_sample_entropy.py
│   ├── CopBET_temporal_entropy.py
│   ├── CopBET_time_series_complexity.py
│   ├── CopBET_von_Neumann_entropy.py
│   └── helper_functions/
│       ├── lempel_ziv.py          ← LZ76 & LZ78 complexity
│       └── sample_entropy_core.py ← sample entropy & MSE
└── external/
    └── DCC_GARCH/                 ← adapted from Topaceminem/DCC-GARCH
        ├── README.md              ← source & license info
        └── DCC_GARCH/
            ├── DCC/
            │   ├── DCC.py
            │   └── DCC_loss.py
            └── GARCH/
                ├── GARCH.py
                └── GARCH_loss.py
```

---

## Key design differences from MATLAB

| Aspect | MATLAB CopBET | CopBET_py |
|---|---|---|
| Input format | MATLAB `table` with cell-wrapped matrices | Python `list` of `np.ndarray` |
| Output | Table with added entropy column | `np.ndarray` or `list` |
| BCT functions | Included as `.m` files | Uses `bctpy` package |
| DCC | Custom MATLAB DCC toolbox | Adapted from `Topaceminem/DCC-GARCH` |
| Sample entropy | LOFT Complexity Toolbox | Translated inline |
| LZ complexity | External `calc_lz_complexity.m` | Translated in `helper_functions/lempel_ziv.py` |
| NIfTI voxelwise functions | Reads 4D NIfTI | Accepts 2D (T × N) arrays directly |
| Parallelism | MATLAB `parfor` | Single-threaded (add `joblib` to parallelise) |

---

## Notes on specific measures

- **Intranetwork synchrony**: requires voxelwise data per ROI. Pass a list of
  lists where `session[roi]` is `(n_voxels, T)`. Alternatively pass a `(T, N)`
  array with one voxel per ROI.

- **Motif connectivity entropy**: designed for exactly 4 ROIs. Specify
  `roi_indices` to pick which 4 from a larger matrix.

- **DCC entropy**: very slow for large N. The original paper ran for several days.
  Use `dcc_sessions = [ts[:, :K] for ts in sessions]` to reduce to K ROIs.

- **Degree distribution** and **Geodesic entropy**: return arrays of length
  `n_chosen_degrees` per session, not single scalars.

- **Sample entropy**: returns `(n_scales, N_rois)` per session.
