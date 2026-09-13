"""Ports postGetDP.m: parses getDP's NodeTable .pos output and interpolates
node-wise scalar/vector fields onto the regular MRI voxel grid.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator


def read_pos_node_table(path: str, n_components: int) -> tuple[np.ndarray, np.ndarray]:
    """Reads a getDP NodeTable .pos file: a header line (ignored, matching
    postGetDP.m's `fgetl(fid)`), then one row per node of
    `node_id value_1 .. value_n`, until a non-numeric (footer) line.

    Returns (node_ids, values): node_ids are 1-based (mesh node numbering,
    not roast_py's usual 0-based convention -- see roast_py.fem's package
    docstring), values has shape (n_rows, n_components).
    """
    ids = []
    vals = []
    with open(path) as f:
        f.readline()  # header line, e.g. `View "v" {`
        for line in f:
            parts = line.split()
            if len(parts) != 1 + n_components:
                break  # footer line (e.g. `};`), matching textscan's implicit stop
            try:
                row = [float(p) for p in parts]
            except ValueError:
                break
            ids.append(int(row[0]))
            vals.append(row[1:])
    return np.array(ids, dtype=np.int64), np.array(vals, dtype=float)


def interpolate_to_grid(points: np.ndarray, values: np.ndarray, grid_shape: tuple[int, int, int]) -> np.ndarray:
    """Ports postGetDP.m's `TriScatteredInterp(...)` + evaluation on
    `ndgrid(1:dim(1),1:dim(2),1:dim(3))`: linear interpolation (barycentric,
    on the scattered points' own Delaunay triangulation) of `values`
    (shape (n_points,) or (n_points, k)) from `points` (mesh node physical
    coordinates) onto every voxel-center point of a `grid_shape` grid.
    Voxels outside the data's convex hull come back as NaN, matching
    MATLAB's TriScatteredInterp default (no extrapolation).
    """
    xi, yi, zi = np.meshgrid(
        np.arange(1, grid_shape[0] + 1),
        np.arange(1, grid_shape[1] + 1),
        np.arange(1, grid_shape[2] + 1),
        indexing="ij",
    )
    query = np.stack([xi.ravel(), yi.ravel(), zi.ravel()], axis=1)

    interp = LinearNDInterpolator(points, values)
    out = interp(query)
    return out.reshape(grid_shape) if values.ndim == 1 else out.reshape(*grid_shape, values.shape[1])
