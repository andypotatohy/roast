"""Subject-to-MNI registration: NiftyReg (ports runNiftyReg.m), or the one
SPM's segmentation saves in `_seg8.mat`."""

from .niftyreg import Registration, find_reg_aladin, registration_path, run_niftyreg
from .seg8 import registration_from_seg8, spm_mat_to_affine

__all__ = [
    "Registration",
    "find_reg_aladin",
    "registration_from_seg8",
    "registration_path",
    "run_niftyreg",
    "spm_mat_to_affine",
]
