"""Interim stand-in for ROAST's real landmark detection (checkLandmarks.m
and the SPM/NiftyReg-to-MNI registration it relies on to map standard
landmark coordinates onto an individual head -- not yet ported, see
roast_py/README.md and the tracked follow-up task).

`heuristic_landmarks` derives nasion/inion/left/right/neck points directly
from the segmented scalp's bounding box instead: good enough to place
electrodes at roughly the right anatomical locations for smoke-testing the
pipeline (verified in tests/test_placement_integration.py and
tests/test_cap_fitting.py to produce anatomically sane placement -- Cz at
the vertex, Fpz frontal, Oz occipital, etc.), but NOT a substitute for real
landmark detection: it will be noticeably less accurate on any head that
isn't roughly upright and axis-aligned in the scan, and has no sub-voxel
precision. Replace this with real landmark detection before trusting
results for anything beyond pipeline smoke-testing.
"""

from __future__ import annotations

import numpy as np


def heuristic_landmarks(tissue_labels: np.ndarray) -> np.ndarray:
    """Derives the 6 landmarks electrode_placement() needs (nasion, inion,
    right, left, front_neck, back_neck) from a segmented head's scalp
    bounding box. Assumes RAS orientation (x=L-R, y=P-A, z=I-S), matching
    this package's convention throughout.
    """
    scalp_idx = np.argwhere(tissue_labels > 0)
    mid_x = int(round((scalp_idx[:, 0].min() + scalp_idx[:, 0].max()) / 2))

    slab = scalp_idx[np.abs(scalp_idx[:, 0] - mid_x) <= 3]
    nasion = slab[np.argmax(slab[:, 1])].astype(float)
    inion = slab[np.argmin(slab[:, 1])].astype(float)

    mid_y = int(round((scalp_idx[:, 1].min() + scalp_idx[:, 1].max()) / 2))
    mid_z = int(round((scalp_idx[:, 2].min() + scalp_idx[:, 2].max()) / 2))
    band = scalp_idx[(np.abs(scalp_idx[:, 1] - mid_y) <= 5) & (np.abs(scalp_idx[:, 2] - mid_z) <= 5)]
    right = band[np.argmax(band[:, 0])].astype(float)
    left = band[np.argmin(band[:, 0])].astype(float)

    # Neck landmarks: lowest-z scalp points near the midline, front/back
    # split by y. Only used if the montage includes neck electrodes.
    low_z = scalp_idx[scalp_idx[:, 2] <= scalp_idx[:, 2].min() + 3]
    low_z_mid = low_z[np.abs(low_z[:, 0] - mid_x) <= 5]
    if low_z_mid.shape[0] > 0:
        front_neck = low_z_mid[np.argmax(low_z_mid[:, 1])].astype(float)
        back_neck = low_z_mid[np.argmin(low_z_mid[:, 1])].astype(float)
    else:
        front_neck = nasion.copy()
        back_neck = inion.copy()

    return np.array([nasion, inion, right, left, front_neck, back_neck])
