"""Ports cylinder2P.m, drawCylinder.m, drawCuboid.m, and drawLine.m --
point-cloud generators used to voxelize electrode/gel shapes."""

from __future__ import annotations

import numpy as np


def cylinder2p(radii: np.ndarray, n_theta: int, r1: np.ndarray, r2: np.ndarray):
    """Ports cylinder2P.m: an N-sided cylindrical (or frustum-like, if
    `radii` varies) surface between points r1 and r2.

    Deviates from the original in one way: MATLAB picks an arbitrary
    vector perpendicular to the cylinder axis via a *random* vector
    (`rand(1,3)`) and Gram-Schmidt. That arbitrary in-plane rotation
    doesn't affect the final result anywhere it's used here (the
    circumference is always sampled densely and the points get voxelized,
    i.e. rounded to integer coordinates and de-duplicated), so this uses a
    deterministic perpendicular basis (the standard basis vector least
    parallel to the axis) instead, for reproducibility.

    Returns X, Y, Z, each shape (m, n_theta) where m = max(len(radii), 2).
    """
    radii = np.atleast_1d(np.asarray(radii, dtype=float))
    if len(radii) == 1:
        radii = np.array([radii[0], radii[0]])
    m = len(radii)

    theta = np.linspace(0, 2 * np.pi, n_theta)
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)
    v = r2 - r1
    v = v / np.linalg.norm(v)

    e = np.eye(3)[np.argmin(np.abs(v))]
    x2 = e - v * np.dot(e, v)
    x2 = x2 / np.linalg.norm(x2)
    x3 = np.cross(v, x2)
    x3 = x3 / np.linalg.norm(x3)

    t = np.linspace(0, 1, m)
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    X = np.empty((m, n_theta))
    Y = np.empty((m, n_theta))
    Z = np.empty((m, n_theta))
    for j in range(m):
        center = r1 + (r2 - r1) * t[j]
        ring = center[:, None] + radii[j] * cos_t[None, :] * x2[:, None] + radii[j] * sin_t[None, :] * x3[:, None]
        X[j], Y[j], Z[j] = ring
    return X, Y, Z


def draw_line(p0, direction, length: float, density: float) -> np.ndarray:
    """Ports drawLine.m: a point cloud of a straight line segment."""
    p0 = np.asarray(p0, dtype=float)
    direction = np.asarray(direction, dtype=float)
    n_samp = int(round(length * density))
    j = np.arange(n_samp + 1)
    coords = p0[None, :] + (1.0 / density) * direction[None, :] * j[:, None]
    return np.unique(coords, axis=0)


def draw_cylinder(inner_radius: float, outer_radius: float, top, bottom, density: float) -> np.ndarray:
    """Ports drawCylinder.m: a point cloud of a (possibly hollow, for a
    ring electrode) cylinder between `top` and `bottom`."""
    top = np.asarray(top, dtype=float)
    bottom = np.asarray(bottom, dtype=float)
    height = np.linalg.norm(top - bottom)

    n_rings = int(round(height * density))
    radii = np.arange(inner_radius + 1 / density, outer_radius + 1e-9, 1 / density)

    coords = []
    for r in radii:
        n_theta = int(round(2 * np.pi * r * density))
        if n_rings == 0 or n_theta == 0:
            continue
        X, Y, Z = cylinder2p(np.full(n_rings, r), n_theta, top, bottom)
        coords.append(np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1))
    return np.concatenate(coords, axis=0) if coords else np.empty((0, 3))


def draw_cuboid(center, dim, long_axis, short_axis, normal_axis, density: float) -> np.ndarray:
    """Ports drawCuboid.m: a point cloud of a cuboid (used for pad
    electrodes), built as a grid of line segments along `long_axis`."""
    center = np.asarray(center, dtype=float)
    long_axis = np.asarray(long_axis, dtype=float)
    short_axis = np.asarray(short_axis, dtype=float)
    normal_axis = np.asarray(normal_axis, dtype=float)

    corner = center - long_axis * dim[0] / 2 - short_axis * dim[1] / 2 - normal_axis * dim[2] / 2
    n_short = int(round(dim[1] * density + 1))
    n_normal = int(round(dim[2] * density + 1))

    coords = []
    for s in range(n_short):
        for n in range(n_normal):
            start = corner + (1 / density) * short_axis * s + (1 / density) * normal_axis * n
            coords.append(draw_line(start, long_axis, dim[0], density))
    return np.concatenate(coords, axis=0) if coords else np.empty((0, 3))
