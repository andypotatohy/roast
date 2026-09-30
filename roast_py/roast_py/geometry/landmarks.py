"""Head landmarks from the TPM, as roast.m computes them.

ROAST defines its landmarks once, as voxel coordinates in the eTPM.nii
atlas (`landmarksInTPM` in roast.m), and carries them onto each head
through the subject-to-MNI registration (runNiftyReg.m, or SPM's
_seg8.mat):

    tpm2mri   = inv(image(1).mat) * inv(Affine) * tpm(1).mat
    landmarks = round(tpm2mri * [landmarksInTPM 1]')

This replaces the scalp-bounding-box heuristic roast_py used before
registration was ported.

Rows, in roast.m's order: nasion, inion, right, left (ear points; eTPM is
LAS, hence "right" at the low x index there), front_neck, back_neck, the
scalp center, and nine 10-10 electrodes on the central sagittal line.
electrode_placement() uses the first six; the rest are what roast.m's
manual landmark correction (checkLandmarks, not ported) refits the
registration with.
"""

from __future__ import annotations

import numpy as np

LANDMARK_NAMES = (
    "nasion",
    "inion",
    "right",
    "left",
    "front_neck",
    "back_neck",
    "scalp_center",
    *(f"sagittal_10-10_{i}" for i in range(1, 10)),
)

# roast.m's landmarksInTPM, verbatim: 1-based (MATLAB) voxel coordinates in eTPM.nii.
LANDMARKS_IN_TPM = np.array(
    [
        [61, 139, 98],  # nasion
        [61, 9, 100],  # inion
        [11, 62, 93],  # right
        [111, 63, 93],  # left; note here because eTPM is LAS orientation
        [61, 113, 7],  # front_neck
        [61, 7, 20],  # back_neck
        [61.1698, 74.8445, 128.6539],  # scalp center
        [61.1266, 6.7765, 128.3046],  # nine 10-10 electrodes on the central sagittal line
        [61.0821, 13.1721, 153.1825],
        [61.0449, 26.7791, 176.6905],
        [61.0294, 49.4432, 192.0928],
        [61.0444, 76.4435, 193.3487],
        [61.0751, 101.3033, 185.8634],
        [61.1154, 122.7582, 172.3329],
        [61.1658, 137.4929, 151.4222],
        [61.2171, 141.6240, 126.5092],
    ]
)


def _matlab_round(x: np.ndarray) -> np.ndarray:
    """MATLAB's round(): halves go away from zero (numpy's go to even)."""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def tpm_landmarks_to_subject(tpm2mri: np.ndarray) -> np.ndarray:
    """Maps LANDMARKS_IN_TPM onto the subject: (16, 3) int, 0-based voxels.

    `tpm2mri` is the 0-based eTPM-voxel -> subject-voxel matrix
    (roast_py.registration.Registration.tpm2mri). Rounding happens on the
    1-based coordinates, as in roast.m, so the result is exactly MATLAB's
    landmarks minus one.
    """
    tpm_zero_based = LANDMARKS_IN_TPM - 1.0
    homogeneous = np.column_stack([tpm_zero_based, np.ones(len(tpm_zero_based))])
    subject_zero_based = (np.asarray(tpm2mri, dtype=float) @ homogeneous.T).T[:, :3]
    return (_matlab_round(subject_zero_based + 1.0) - 1.0).astype(int)
