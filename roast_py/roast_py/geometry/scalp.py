"""Ports cleanScalp.m: morphological clean-up of the scalp mask so it has
a smooth, closed outer surface suitable for computing local normal vectors
when placing electrodes."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

_CUBE3 = np.ones((3, 3, 3), dtype=bool)


def _close_cube(mask: np.ndarray, se: int) -> np.ndarray:
    """binary_closing(mask, structure=ones((se,se,se))), computed as `se`
    iterations of dilation/erosion with a 3x3x3 cube instead of a single
    call with an se^3 structuring element -- mathematically identical for
    an all-ones cube (dilating k times with a 3x3x3 cube == dilating once
    with a (2k+1)^3 cube), but tractable at full head-volume resolution
    where a literal se=30-ish cube would be far too slow.
    """
    dilated = ndimage.binary_dilation(mask, structure=_CUBE3, iterations=se, border_value=0)
    return ndimage.binary_erosion(dilated, structure=_CUBE3, iterations=se, border_value=1)


def _open_cube(mask: np.ndarray, se: int) -> np.ndarray:
    """binary_opening(mask, structure=ones((se,se,se))), same iterated-cube
    trick as _close_cube."""
    eroded = ndimage.binary_erosion(mask, structure=_CUBE3, iterations=se, border_value=1)
    return ndimage.binary_dilation(eroded, structure=_CUBE3, iterations=se, border_value=0)


def clean_scalp(scalp_in: np.ndarray, scalp_surf: np.ndarray):
    """Returns (scalp_out, scalp_filled): scalp_filled is the scalp mask
    closed and filled solid (no interior holes, e.g. from hair/sinuses at
    the mask level); scalp_out is scalp_filled further opened to smooth
    out the outer surface for stable normal-vector estimation.
    """
    scalp_in = np.array(scalp_in, dtype=bool, copy=True)

    lo = scalp_surf.min(axis=0).astype(int)
    hi = scalp_surf.max(axis=0).astype(int)
    # Force-close the bottom-most slice: some MRIs cut off mid-face, which
    # would otherwise leave the scalp shell open there.
    scalp_in[lo[0] : hi[0] + 1, lo[1] : hi[1] + 1, lo[2]] = True

    xtemp, ytemp, ztemp = np.nonzero(scalp_in)
    centroid = (
        int(round(xtemp.mean())),
        int(round(ytemp.mean())),
        int(round(ztemp.mean())),
    )

    se = 0
    is_filled = False
    scalp_mid = scalp_in
    while not is_filled:
        se += 10
        closed = _close_cube(scalp_in, se)
        scalp_mid = ndimage.binary_fill_holes(closed)
        is_filled = bool(scalp_mid[centroid])

    nx, ny, nz = scalp_in.shape
    se = 0
    is_open = np.array([True])
    scalp_out = scalp_mid
    while np.any(is_open):
        se += 10
        if se > 30:
            scalp_mid = scalp_mid.copy()
            scalp_mid[:, :, 0] = False
            scalp_mid[:, :, nz - 1] = False
            scalp_mid[:, 0, :] = False
            scalp_mid[:, ny - 1, :] = False
            scalp_mid[0, :, :] = False
            scalp_mid[nx - 1, :, :] = False

        scalp_out = _open_cube(scalp_mid, se)
        is_open = np.concatenate(
            [
                scalp_out[:, :, nz - 1].ravel(),
                scalp_out[:, 0, :].ravel(),
                scalp_out[0, :, :].ravel(),
                scalp_out[nx - 1, :, :].ravel(),
            ]
        )

    return scalp_out, scalp_mid
