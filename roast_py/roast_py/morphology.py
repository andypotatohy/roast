"""MATLAB's binary morphology with a box structuring element,
`ones(n,n,n)`, as ROAST's segTouchup.m and cleanScalp.m use it.

Each box is applied as three 1-D filters (it is separable), so any size
costs about the same. Even sizes are centred where MATLAB centres them
(floor((n+1)/2)), and the border is handled MATLAB's way: outside the
image counts as background for dilation and as foreground for erosion.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def imdilate_box(mask: np.ndarray, n: int) -> np.ndarray:
    """MATLAB's imdilate(mask, ones(n,n,n)), even n included."""
    out = mask.astype(np.uint8)
    for axis in range(mask.ndim):
        out = ndimage.maximum_filter1d(out, n, axis=axis, mode="constant", cval=0, origin=0)
    return out.astype(bool)


def imerode_box(mask: np.ndarray, n: int) -> np.ndarray:
    """MATLAB's imerode(mask, ones(n,n,n)): for even n the box's centre is
    MATLAB's floor((n+1)/2), and the outside counts as foreground."""
    out = mask.astype(np.uint8)
    origin = -1 if n % 2 == 0 else 0
    for axis in range(mask.ndim):
        out = ndimage.minimum_filter1d(out, n, axis=axis, mode="constant", cval=1, origin=origin)
    return out.astype(bool)


def imclose_box(mask: np.ndarray, n: int) -> np.ndarray:
    """MATLAB's imclose(mask, ones(n,n,n)): dilation, then erosion."""
    return imerode_box(imdilate_box(mask, n), n)


def imopen_box(mask: np.ndarray, n: int) -> np.ndarray:
    """MATLAB's imopen(mask, ones(n,n,n)): erosion, then dilation."""
    return imdilate_box(imerode_box(mask, n), n)
