"""Ports binaryMaskGenerate.m, sizeOfObject.m, and the core cleanup pipeline
of segTouchup.m to numpy/scipy.

This is the post-processing ROAST applies to SPM's raw per-tissue
probability maps (c1..c6) to get a single clean 6-label volume. It's needed
regardless of which segmenter produced those probability maps, so it's kept
separate from roast_py.segmentation.spm_standalone.

Known gap vs. segTouchup.m: the three final "patching" passes (gray-matter,
CSF, and bone patching) rely on `eyes_vol`, `holes_vol`, and
`WMexclude_vol` -- three extra probability volumes that ROAST's own patched
copy of SPM's `spm_preproc_write8.m` (lib/spm12/spm_preproc_write8.m)
computes and saves alongside the normal c1..c6 outputs, from extra tissue
classes in the extended eTPM.nii atlas beyond the standard 6. Reproducing
that requires understanding those SPM-internal modifications in detail,
which is out of scope for this pass -- see the tracked follow-up task.
Everything else in segTouchup.m's pipeline (smoothing, binarization, CSF
continuity fix, disconnected-voxel removal, empty-voxel relabeling, final
tissue-label assembly) is ported here.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

# Tissue label order in the final mask (matches ROAST's 'conductivities'
# option order, and roast_py.segmentation.multiaxial's output labels).
WHITE, GRAY, CSF, BONE, SKIN, AIR = 1, 2, 3, 4, 5, 6


def _gaussian_kernel_5x5(sigma: float) -> np.ndarray:
    """Matches MATLAB's fspecial('gaussian', 5, sigma)."""
    ax = np.arange(-2, 3)
    xx, yy = np.meshgrid(ax, ax, indexing="ij")
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    return kernel / kernel.sum()


def _smooth_slices(volume: np.ndarray, sigma: float) -> np.ndarray:
    """Matches segTouchup.m's per-slice `imfilter(slice, fspecial('gaussian',5,sigma))`:
    correlation (== convolution here, kernel is symmetric) with zero-padded
    'same'-size output, applied independently to each axial slice."""
    kernel = _gaussian_kernel_5x5(sigma)
    out = np.empty_like(volume, dtype=float)
    for i in range(volume.shape[2]):
        out[:, :, i] = ndimage.convolve(volume[:, :, i], kernel, mode="constant", cval=0.0)
    return out


def binary_mask_generate(*arrays: np.ndarray):
    """Ports binaryMaskGenerate.m: stacks the inputs behind an implicit
    all-zero "background" volume and assigns each voxel to whichever input
    is largest there. Returns (empty_mask, mask_for_each_input...), where
    empty_mask marks voxels where the background "won" (i.e. none of the
    inputs claimed that voxel).
    """
    if not arrays:
        raise ValueError("binary_mask_generate needs at least one array")
    stacked = np.stack([np.zeros_like(arrays[0]), *arrays], axis=-1)
    maxind = np.argmax(stacked, axis=-1)
    empty = maxind == 0
    masks = tuple(maxind == (i + 1) for i in range(len(arrays)))
    return (empty, *masks)


def _connectivity_structure(ndim: int, conn: int) -> np.ndarray:
    """Maps MATLAB bwconncomp connectivity values to scipy's connectivity
    parameter for ndimage.generate_binary_structure: MATLAB's 4/6 (face
    neighbors only) -> scipy 1, MATLAB's 8/18 (+edges) -> scipy 2, MATLAB's
    26 (+corners) -> scipy 3."""
    mapping = {4: 1, 8: 2, 6: 1, 18: 2, 26: 3}
    if conn not in mapping:
        raise ValueError(f"Unsupported connectivity {conn}")
    return ndimage.generate_binary_structure(ndim, mapping[conn])


def size_of_object(mask: np.ndarray, conn: int | None = None):
    """Ports sizeOfObject.m: sizes of connected components in a binary mask,
    sorted descending, plus the (1-based, to match MATLAB) label index each
    size corresponds to."""
    if conn is None:
        conn = 8 if mask.ndim == 2 else 26
    structure = _connectivity_structure(mask.ndim, conn)
    labeled, num = ndimage.label(mask, structure=structure)
    if num == 0:
        return np.array([]), np.array([], dtype=int)
    sizes = ndimage.sum(mask, labeled, index=np.arange(1, num + 1))
    order = np.argsort(sizes)[::-1]
    return sizes[order], order + 1


def bwareaopen(mask: np.ndarray, thres: float, conn: int) -> np.ndarray:
    """Ports MATLAB's bwareaopen: removes connected components with fewer
    than `thres` voxels."""
    structure = _connectivity_structure(mask.ndim, conn)
    labeled, num = ndimage.label(mask, structure=structure)
    if num == 0:
        return mask.copy()
    sizes = ndimage.sum(mask, labeled, index=np.arange(1, num + 1))
    keep = np.zeros(num + 1, dtype=bool)
    keep[1:] = sizes >= thres
    return keep[labeled]


def touchup(
    gray: np.ndarray,
    white: np.ndarray,
    csf: np.ndarray,
    bone: np.ndarray,
    skin: np.ndarray,
    air: np.ndarray,
    is_smooth: bool = True,
    conn: int = 18,
) -> np.ndarray:
    """Ports segTouchup.m's cleanup pipeline (minus the eyes/holes/WM-exclude
    patching passes -- see module docstring) from six per-tissue probability
    volumes to a single uint8 label volume using WHITE/GRAY/CSF/BONE/SKIN/AIR.
    """
    gray_t, white_t, csf_t, bone_t, skin_t, air_t = (
        a.astype(float) for a in (gray, white, csf, bone, skin, air)
    )

    if is_smooth:
        gray_t = _smooth_slices(gray_t, 0.2)
        white_t = _smooth_slices(white_t, 0.1)
        csf_t = _smooth_slices(csf_t, 0.1)
        bone_t = _smooth_slices(bone_t, 0.4)
        skin_t = _smooth_slices(skin_t, 1.0)
        air_t = _smooth_slices(air_t, 1.0)

    empt_t, gray_t, white_t, csf_t, bone_t, skin_t, air_t = binary_mask_generate(
        gray_t, white_t, csf_t, bone_t, skin_t, air_t
    )

    # Fix CSF continuity: CSF-adjacent empty voxels, and bone-adjacent GM
    # voxels, both get folded into CSF before the final assignment.
    se = np.ones((3, 3, 3), dtype=bool)
    dcsf = ndimage.binary_dilation(csf_t, structure=se)
    dbone = ndimage.binary_dilation(bone_t, structure=se)
    contin = (empt_t & dcsf) | (dbone & gray_t)
    csf_t = csf_t | contin

    _, csf_t, bone_t, gray_t = binary_mask_generate(csf_t, bone_t, gray_t)

    # Remove disconnected voxels, keeping (for most tissues) only the
    # largest connected component -- see segTouchup.m for the per-tissue
    # rationale on why CSF keeps its top-3 components.
    def _keep_largest(mask, n_keep=1):
        sizes, _ = size_of_object(mask, conn)
        if len(sizes) >= n_keep + 1:
            thres = sizes[n_keep] + 1
            return bwareaopen(mask, thres, conn)
        return mask

    gray_t = _keep_largest(gray_t)
    white_t = _keep_largest(white_t)
    csf_t = _keep_largest(csf_t, n_keep=3)
    bone_t = _keep_largest(bone_t)
    skin_t = _keep_largest(skin_t)
    air_t = bwareaopen(air_t, 20, 26 if air_t.ndim == 3 else 8)

    # Iteratively relabel empty voxels to whichever tissue's smoothed,
    # binarized field reaches them first (nearest tissue by Gaussian-blurred
    # distance), same approach as segTouchup.m's while loop.
    empt_t = binary_mask_generate(gray_t, white_t, csf_t, bone_t, skin_t, air_t)[0]
    tissues = [gray_t, white_t, csf_t, bone_t, skin_t, air_t]
    sigmas = [1, 1, 1, 1, 1, 1]
    guard = 0
    while np.any(empt_t):
        guard += 1
        if guard > 50:
            raise RuntimeError("touchup(): empty-voxel relabeling did not converge")
        filled = [
            _smooth_slices((t.astype(np.uint8) * 255).astype(float), s)
            for t, s in zip(tissues, sigmas)
        ]
        _, *fil_masks = binary_mask_generate(*filled)
        new_tissues = []
        newly_assigned = np.zeros_like(empt_t)
        for t, fil in zip(tissues, fil_masks):
            assigned = empt_t & fil
            newly_assigned |= assigned
            new_tissues.append(assigned | t)
        tissues = new_tissues
        empt_t = empt_t & ~newly_assigned

    gray_t, white_t, csf_t, bone_t, skin_t, air_t = tissues

    all_mask = np.zeros_like(gray_t, dtype=np.uint8)
    all_mask[white_t] = WHITE
    all_mask[gray_t] = GRAY
    all_mask[csf_t] = CSF
    all_mask[bone_t] = BONE
    all_mask[skin_t] = SKIN
    all_mask[air_t] = AIR
    return all_mask
