"""NIfTI I/O and header/geometry utilities.

Ports the pure-numerical preprocessing steps that ROAST performs in MATLAB
via the bundled ``lib/NIFTI_20110921`` toolbox and SPM header calls:

- ``convertToRAS.m``    -> :func:`convert_to_ras`
- ``zeroPadding.m``     -> :func:`zero_pad`
- ``alignHeader2mni.m`` -> :func:`align_header_to_mni`

``resampToOneMM.m`` and ``realignT2.m`` are NOT ported here: both call into
SPM (``spm_reslice`` / ``spm_jobman`` coreg-estimate), which do real
registration/resampling optimization rather than pure header math. They are
deferred to the segmentation phase of the Python port, where the SPM
dependency as a whole gets resolved (see /roast_py/README.md).
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import nibabel as nib
import numpy as np


@dataclass
class RasResult:
    path: str
    was_reoriented: bool


def _split(path: str) -> tuple[str, str, str]:
    dirname, filename = os.path.split(path)
    if not dirname:
        dirname = os.getcwd()
    base, ext = os.path.splitext(filename)
    if ext == ".gz":
        base2, ext2 = os.path.splitext(base)
        base, ext = base2, ext2 + ext
    return dirname, base, ext


def convert_to_ras(mri_path: str) -> RasResult:
    """Re-orient ``mri_path`` to RAS+ if it isn't already, mirroring convertToRAS.m.

    Uses nibabel's ``as_closest_canonical``, which performs the same
    axis-permutation + flip + affine-origin update as the hand-rolled
    sform/qform math in convertToRAS.m, but is the well-tested nibabel
    primitive for this operation rather than a reimplementation of it.
    Returns the (possibly unchanged) path to use downstream, matching
    ROAST's behavior of caching the reoriented file with a ``_ras`` suffix
    and reusing it if already present.
    """
    img = nib.load(mri_path)
    canonical = nib.as_closest_canonical(img)

    # as_closest_canonical returns a new image even when no reorientation
    # was needed; detect that case by comparing affines/orientation.
    orig_ornt = nib.io_orientation(img.affine)
    identity_ornt = np.array([[0, 1], [1, 1], [2, 1]])
    if np.array_equal(orig_ornt, identity_ornt):
        return RasResult(path=mri_path, was_reoriented=False)

    dirname, base, ext = _split(mri_path)
    ras_path = os.path.join(dirname, f"{base}_ras{ext}")
    if not os.path.exists(ras_path):
        nib.save(canonical, ras_path)
    return RasResult(path=ras_path, was_reoriented=True)


def zero_pad(mri_path: str, pad_num: int) -> str:
    """Pad ``mri_path`` with ``pad_num`` empty voxel slices on all six faces.

    Mirrors zeroPadding.m: the image content is unchanged, only surrounded
    by zeros, and the affine's translation is updated so world coordinates
    of the original content are preserved (a voxel that used to be at index
    (0,0,0) is now at (pad_num,pad_num,pad_num)).
    """
    if "_padded" in mri_path:
        raise ValueError(
            f"{mri_path} looks already zero-padded; pass the MRI name "
            "without the _padded suffix."
        )

    dirname, base, ext = _split(mri_path)
    pad_path = os.path.join(dirname, f"{base}_padded{pad_num}{ext}")
    if os.path.exists(pad_path):
        return pad_path

    img = nib.load(mri_path)
    data = np.asarray(img.dataobj)
    p = pad_num
    padded = np.zeros(tuple(s + 2 * p for s in data.shape), dtype=data.dtype)
    padded[p : p + data.shape[0], p : p + data.shape[1], p : p + data.shape[2]] = data

    affine = img.affine.copy()
    # New voxel (p,p,p) must map to the same world point as old voxel (0,0,0):
    # t_new = t_old - R @ [p, p, p]
    affine[:3, 3] = affine[:3, 3] - affine[:3, :3] @ np.array([p, p, p], dtype=float)

    padded_img = nib.Nifti1Image(padded, affine, img.header)
    nib.save(padded_img, pad_path)
    return pad_path


def align_header_to_mni(mri_path: str, mri2mni: np.ndarray, out_name: str) -> str:
    """Rewrite ``mri_path``'s header so voxel-to-world maps directly to MNI space.

    Mirrors alignHeader2mni.m's ``update_affine`` helper: replaces the
    image's affine with the supplied 4x4 ``mri2mni`` matrix (voxel -> MNI
    world), forces it to be recorded via sform, and writes the result with
    an ``_MNI`` suffix.
    """
    if mri2mni.shape != (4, 4):
        raise ValueError(f"mri2mni must be a 4x4 matrix, got shape {mri2mni.shape}")

    img = nib.load(mri_path)
    out_img = nib.Nifti1Image(np.asarray(img.dataobj), mri2mni, img.header)
    out_img.header.set_qform(mri2mni, code=0)
    out_img.header.set_sform(mri2mni, code=1)

    dirname, _, ext = _split(mri_path)
    out_path = os.path.join(dirname, f"{out_name}_MNI{ext}")
    nib.save(out_img, out_path)
    return out_path
