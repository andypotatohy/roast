"""Ports mask2EdgePointCloud.m, project2ClosestSurfacePoints.m,
map2Points.m, and convertToRASpointCloud.m."""

from __future__ import annotations

import nibabel as nib
import numpy as np
from scipy import ndimage


def mask2edge_point_cloud(mask: np.ndarray, method: str, structure: np.ndarray):
    """Ports mask2EdgePointCloud.m: the voxel shell just inside (erode) or
    just outside (dilate) `mask`, as a point cloud of 0-based (i, j, k)
    voxel coordinates.

    Border handling matches MATLAB's imerode/imdilate defaults: erosion
    treats voxels outside the volume as foreground (so the true edge isn't
    artificially eaten away at the image boundary), dilation treats them as
    background.
    """
    mask = np.asarray(mask, dtype=bool)
    if method == "erode":
        processed = ndimage.binary_erosion(mask, structure=structure, border_value=1)
        edge = mask & ~processed
    elif method == "dilate":
        processed = ndimage.binary_dilation(mask, structure=structure, border_value=0)
        edge = processed & ~mask
    else:
        raise ValueError("method must be 'erode' or 'dilate'")
    edge_point_cloud = np.argwhere(edge).astype(float)
    return edge_point_cloud, processed


def project2closest_surface_points(points: np.ndarray, surf: np.ndarray, surf_center: np.ndarray):
    """Ports project2ClosestSurfacePoints.m: ranks every surface point by
    how closely its direction from `surf_center` aligns (cosine similarity)
    with each query point's direction from `surf_center`.

    Simplified from the original's per-point-pair repmat/dot broadcasting
    to a single matrix product -- same cosine-similarity computation,
    vectorized directly instead of manually.

    Returns (cosine_sorted, index_on_surf), each of shape
    (n_surface_points, n_points), sorted descending (most-aligned first)
    along axis 0.
    """
    points = np.atleast_2d(np.asarray(points, dtype=np.float32))
    surf = np.asarray(surf, dtype=np.float32)
    surf_center = np.asarray(surf_center, dtype=np.float32)

    vec_p = points - surf_center
    vec_p = vec_p / np.linalg.norm(vec_p, axis=1, keepdims=True)

    vec_s = surf - surf_center
    vec_s = vec_s / np.linalg.norm(vec_s, axis=1, keepdims=True)

    cosine = vec_s @ vec_p.T  # (n_surface, n_points)
    order = np.argsort(-cosine, axis=0)
    cosine_sorted = np.take_along_axis(cosine, order, axis=0)
    return cosine_sorted, order


def map2points(input_points: np.ndarray, goal_points: np.ndarray, criterion: str, num_of_pts: int | None = None):
    """Ports map2Points.m: for each input point, ranks goal points by
    Euclidean distance and returns the closest/farthest (or the closest/
    farthest `num_of_pts`).
    """
    input_points = np.atleast_2d(np.asarray(input_points, dtype=np.float32))
    goal_points = np.asarray(goal_points, dtype=np.float32)

    diff = goal_points[:, np.newaxis, :] - input_points[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=2))  # (n_goal, n_input)

    order = np.argsort(dist, axis=0)
    dist_sorted = np.take_along_axis(dist, order, axis=0)

    if criterion == "closest":
        return dist_sorted[0], order[0]
    if criterion == "farthest":
        return dist_sorted[-1], order[-1]
    if criterion == "closer":
        if num_of_pts is None:
            raise ValueError("num_of_pts is required for criterion='closer'")
        if num_of_pts > goal_points.shape[0]:
            raise ValueError("num_of_pts exceeds size of goal point cloud")
        return dist_sorted[:num_of_pts], order[:num_of_pts]
    if criterion == "farther":
        if num_of_pts is None:
            raise ValueError("num_of_pts is required for criterion='farther'")
        if num_of_pts > goal_points.shape[0]:
            raise ValueError("num_of_pts exceeds size of goal point cloud")
        return dist_sorted[-num_of_pts:], order[-num_of_pts:]
    raise ValueError("criterion must be one of closest/farthest/closer/farther")


def convert_to_ras_point_cloud(affine: np.ndarray, data: np.ndarray, shape: tuple[int, int, int]):
    """Ports convertToRASpointCloud.m: flips/permutes a point cloud of
    0-based voxel-index coordinates the same way roast_py.io.nifti.convert_to_ras
    reorients the volume itself, using nibabel's own orientation machinery
    (nib.io_orientation) rather than re-deriving the sform/qform flip logic
    by hand, matching io/nifti.py's approach.

    Returns (data_ras, perm_order).
    """
    ornt = nib.io_orientation(affine)
    perm = ornt[:, 0].astype(int)
    flip = ornt[:, 1]

    data = np.asarray(data, dtype=float)
    out = np.empty_like(data)
    for out_axis in range(3):
        src_axis = perm[out_axis]
        vals = data[:, src_axis]
        if flip[out_axis] < 0:
            vals = shape[src_axis] - 1 - vals
        out[:, out_axis] = vals
    return out, perm
