"""Ports lib/ncs2daprox/ncs2dapprox.m (plus its MergeSortRemoveDuplicates /
MaxSqDistAndRowIndexbw2Mat helpers, folded in directly since Python's `set`
+ `argmax` already do exactly what those hand-rolled routines did).

Used by fitCap2individual.m to get a smooth, compact spline representation
of the scalp's central-sagittal-slice boundary, for measuring the
nasion-inion arc length without being thrown off by pixel-level jaggedness
in the raw edge-detected contour.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline


def ncs2d_approx(x: np.ndarray, y: np.ndarray, max_allowed_sq_dist: float = 1.0):
    """Adaptive natural-cubic-spline approximation: starts with just the
    first and last points as knots, then repeatedly adds whichever data
    point the current spline fits worst, until every point is within
    `max_allowed_sq_dist` of the spline (evaluated at each data point's own
    parameter value 1..n, matching how fitCap2individual.m calls this with
    `spline(t, [bx,by]')` and `ppval` at `1:length(x)`).

    Note MATLAB's `spline()` uses not-a-knot end conditions by default
    (not natural/zero-curvature ones, despite this function's "ncs" name),
    so this uses `scipy.interpolate.CubicSpline(..., bc_type='not-a-knot')`
    to match.

    Returns (bx, by, breaks): the knot x/y coordinates and their 1-based
    indices into (x, y) -- kept 1-based (not roast_py's usual 0-based
    voxel convention) since these are spline parameter values passed
    straight back into CubicSpline, not voxel coordinates.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 2:
        raise ValueError("ncs2d_approx needs at least two points")

    breaks = {1, n}

    def fit():
        t = np.array(sorted(breaks))
        bx = x[t - 1]
        by = y[t - 1]
        if len(t) < 2:
            raise ValueError("need at least two distinct breakpoints")
        cs = CubicSpline(t, np.stack([bx, by], axis=1), bc_type="not-a-knot")
        return cs(np.arange(1, n + 1)), bx, by, t

    xi_yi, bx, by, t = fit()
    sq_dist = (x - xi_yi[:, 0]) ** 2 + (y - xi_yi[:, 1]) ** 2
    worst = int(np.argmax(sq_dist))

    guard = 0
    while sq_dist[worst] > max_allowed_sq_dist:
        guard += 1
        if guard > n:
            break  # every point has become a breakpoint; can't improve further
        breaks.add(worst + 1)
        xi_yi, bx, by, t = fit()
        sq_dist = (x - xi_yi[:, 0]) ** 2 + (y - xi_yi[:, 1]) ** 2
        worst = int(np.argmax(sq_dist))

    return bx, by, t
