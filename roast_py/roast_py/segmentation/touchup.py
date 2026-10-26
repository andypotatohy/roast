"""Ports segTouchup.m (with binaryMaskGenerate.m and sizeOfObject.m) to
numpy/scipy: the cleanup ROAST applies to SPM's six per-tissue probability
maps (c1..c6) to get one clean 6-label volume.

All of segTouchup.m is here: smoothing, binarization, CSF continuity,
removal of disconnected pieces, relabeling of empty voxels, removal of
outside air, and the three patching passes (gray matter, CSF, bone). The
patching passes use three atlas masks ROAST's modified SPM12 warped into
the subject (`_rmask.mat`), which roast_py.spm.routines.rmask_volumes()
reproduces from official SPM's deformation field.

Details kept from MATLAB because they change the result:

* SPM's maps are uint8 (0-255) and MATLAB smooths them as integers: each
  `imfilter` result is rounded, which decides ties in the argmax that
  follows. Pass the raw uint8 maps to match MATLAB exactly.
* Ties between tissues go to the first in each call's argument order,
  which is not the same in every step (air comes first when relabeling
  empty voxels).
* The even-sized boxes (10x10x10, 12x12x12) centre where MATLAB's strel
  does, and erosion treats the outside of the volume as foreground, as
  MATLAB's imerode does.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from ..morphology import imdilate_box, imerode_box

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
    """segTouchup.m's per-slice `imfilter(slice, fspecial('gaussian',5,sigma))`:
    correlation (== convolution here, the kernel is symmetric) with
    zero-padded 'same'-size output, on each axial slice.

    Like imfilter, an integer input gives an integer output: rounded (half
    away from zero) and saturated to its type's range.
    """
    kernel = _gaussian_kernel_5x5(sigma)
    out = np.empty(volume.shape, dtype=float)
    data = volume.astype(float)
    for i in range(volume.shape[2]):
        out[:, :, i] = ndimage.convolve(data[:, :, i], kernel, mode="constant", cval=0.0)
    if np.issubdtype(volume.dtype, np.integer):
        info = np.iinfo(volume.dtype)
        out = np.clip(np.floor(out + 0.5), info.min, info.max).astype(volume.dtype)
    return out


# segTouchup.m's vect_g: the 18 neighbours whose offset has L1 norm 1 or 2
# (faces and edges of the 3x3x3 cube), in ndgrid order.
_NEIGHBOURS_18 = [
    (i, j, k)
    for k in (-1, 0, 1) for j in (-1, 0, 1) for i in (-1, 0, 1)
    if abs(i) + abs(j) + abs(k) in (1, 2)
]


def _interior(shape) -> np.ndarray:
    inside = np.zeros(shape, dtype=bool)
    inside[1:-1, 1:-1, 1:-1] = True
    return inside


def _shifted(volume: np.ndarray, offset) -> np.ndarray:
    """volume[p + offset] for interior p (wraps at the border, which the
    callers never read)."""
    return np.roll(volume, shift=tuple(-o for o in offset), axis=(0, 1, 2))


def _neighbour_has(labels: np.ndarray, wanted) -> np.ndarray:
    """Voxels with at least one of their 18 neighbours labelled in `wanted`."""
    hit = np.zeros(labels.shape, dtype=bool)
    target = np.isin(labels, wanted)
    for offset in _NEIGHBOURS_18:
        hit |= _shifted(target, offset)
    return hit


def _neighbours_of(source: np.ndarray) -> np.ndarray:
    """Voxels that are one of the 18 neighbours of a `source` voxel."""
    hit = np.zeros(source.shape, dtype=bool)
    for offset in _NEIGHBOURS_18:
        hit |= np.roll(source, shift=offset, axis=(0, 1, 2))
    return hit


def patch_labels(all_mask: np.ndarray, holes, eyes, wm_exclude) -> np.ndarray:
    """segTouchup.m's three patching passes, then skull holes to skin.

    * gray matter: white matter touching CSF/bone/skin/air becomes gray,
      except near the marked white-matter surface (`WMexclude_vol`);
    * CSF: bone/skin/air touching white or gray matter becomes CSF;
    * bone: skin/air touching CSF (outside the eye regions) becomes bone;
    * bone inside the skull holes (optic canals etc.) becomes skin.

    "Touching" is the 18-neighbourhood, and only voxels off the volume's
    border are sources, as in MATLAB. Each pass reads the labels the
    previous one left; within a pass the result doesn't depend on order.
    """
    labels = all_mask.copy()
    interior = _interior(labels.shape)

    white = labels == WHITE
    keep = white & imdilate_box(np.asarray(wm_exclude) > 1e-4, 10)
    source = white & ~keep & interior
    labels[source & _neighbour_has(labels, [CSF, BONE, SKIN, AIR])] = GRAY

    source = np.isin(labels, [WHITE, GRAY]) & interior
    labels[_neighbours_of(source) & np.isin(labels, [BONE, SKIN, AIR])] = CSF

    source = (labels == CSF) & ~(np.asarray(eyes) > 1e-4) & interior
    labels[_neighbours_of(source) & np.isin(labels, [SKIN, AIR])] = BONE

    labels[(labels == BONE) & (np.asarray(holes) > 1e-4)] = SKIN
    return labels


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
    holes=None,
    eyes=None,
    wm_exclude=None,
) -> np.ndarray:
    """Ports segTouchup.m: six per-tissue probability maps (SPM's c1..c6
    order: gray, white, CSF, bone, skin, air) to one uint8 label volume
    using WHITE/GRAY/CSF/BONE/SKIN/AIR.

    `holes`, `eyes` and `wm_exclude` are the `_rmask.mat` volumes for the
    patching passes (see patch_labels()); leave them out to skip those
    passes.
    """
    gray_t, white_t, csf_t, bone_t, skin_t, air_t = (
        a if np.issubdtype(np.asarray(a).dtype, np.integer) else np.asarray(a, dtype=float)
        for a in (gray, white, csf, bone, skin, air)
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
    guard = 0
    while np.any(empt_t):
        guard += 1
        if guard > 50:
            raise RuntimeError("touchup(): empty-voxel relabeling did not converge")
        # uint8(mask)*255, smoothed as integers (sigma 1 for every tissue)
        filled = [_smooth_slices(t.astype(np.uint8) * np.uint8(255), 1.0) for t in tissues]
        # MATLAB's call order here puts air first: binaryMaskGenerate(air_fil,
        # gray_fil, white_fil, csf_fil, bone_fil, skin_fil); ties go to air.
        _, air_fil, *rest = binary_mask_generate(filled[5], *filled[:5])
        fil_masks = [*rest, air_fil]
        new_tissues = []
        newly_assigned = np.zeros_like(empt_t)
        for t, fil in zip(tissues, fil_masks):
            assigned = empt_t & fil
            newly_assigned |= assigned
            new_tissues.append(assigned | t)
        tissues = new_tissues
        empt_t = empt_t & ~newly_assigned

    gray_t, white_t, csf_t, bone_t, skin_t, air_t = tissues

    # Remove outside air: keep only air enclosed by the head.
    not_air = ~air_t
    closed = imerode_box(imdilate_box(not_air, 10), 10)
    enclosed = imerode_box(ndimage.binary_fill_holes(closed), 12)
    air_t = air_t & enclosed

    all_mask = np.zeros_like(gray_t, dtype=np.uint8)
    all_mask[white_t] = WHITE
    all_mask[gray_t] = GRAY
    all_mask[csf_t] = CSF
    all_mask[bone_t] = BONE
    all_mask[skin_t] = SKIN
    all_mask[air_t] = AIR

    if holes is not None or eyes is not None or wm_exclude is not None:
        zeros = np.zeros(all_mask.shape)
        all_mask = patch_labels(
            all_mask,
            zeros if holes is None else holes,
            zeros if eyes is None else eyes,
            zeros if wm_exclude is None else wm_exclude,
        )
    return all_mask
