"""Reads the registration SPM's segmentation saves in `<subj>_seg8.mat`.

This is roast.m's other source of `Affine` (and `image`, `tpm`): with SPM
segmentation it loads them from `_seg8.mat` rather than running NiftyReg.
`Affine` maps subject world -> MNI/TPM world (mm), exactly what
Registration.affine holds. SPM's `.mat` voxel-to-world matrices take
1-based voxel indices; Registration's are nibabel affines on 0-based ones,
so they are shifted by one voxel: affine = mat @ T(+1).
"""

from __future__ import annotations

import os

import numpy as np

from .niftyreg import Registration

# 0-based voxel -> 1-based voxel
_ONE_BASED = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1]], dtype=float)


def spm_mat_to_affine(mat) -> np.ndarray:
    """SPM's 1-based voxel-to-world `.mat` -> a nibabel (0-based) affine."""
    return np.asarray(mat, dtype=float) @ _ONE_BASED


def registration_from_seg8(path: str | os.PathLike) -> Registration:
    """Registration from an SPM `_seg8.mat` (variables Affine, image, tpm)."""
    from scipy.io import loadmat

    data = loadmat(str(path), struct_as_record=False, squeeze_me=True)
    image = data["image"]
    image = image[0] if np.ndim(image) else image  # struct array when T1 + T2
    tpm = data["tpm"]
    tpm = tpm[0] if np.ndim(tpm) else tpm  # one struct per tissue class
    return Registration(
        affine=np.asarray(data["Affine"], dtype=float),
        image_affine=spm_mat_to_affine(image.mat),
        tpm_affine=spm_mat_to_affine(tpm.mat),
    )
