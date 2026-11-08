"""Ports mask2EdgePointCloud.m, project2ClosestSurfacePoints.m,
map2Points.m, and convertToRASpointCloud.m."""

from __future__ import annotations

import nibabel as nib
import numpy as np
from scipy import ndimage

from ..morphology import imdilate_box, imerode_box


def mask2edge_point_cloud(mask: np.ndarray, method: str, structure: np.ndarray):
    """Ports mask2EdgePointCloud.m: the voxel shell just inside (erode) or
    just outside (dilate) `mask`, as a point cloud of 0-based (i, j, k)
    voxel coordinates.

    Border handling matches MATLAB's imerode/imdilate defaults: erosion
    treats voxels outside the volume as foreground (so the true edge isn't
    artificially eaten away at the image boundary), dilation treats them as
    background.

    Two details follow MATLAB exactly, as electrode placement depends on
    them. A box `ones(n,n,n)` goes through roast_py.morphology, which
    centres even sizes where MATLAB does: scipy centres them one voxel
    over, which moved the 4 mm pad of runExamples cmd 25 (POz) a voxel
    outward. And the points come in find()'s column-major order, which is
    how MATLAB breaks ties when it ranks them.
    """
    mask = np.asarray(mask, dtype=bool)
    structure = np.asarray(structure, dtype=bool)
    is_box = structure.ndim == mask.ndim and structure.all() and len(set(structure.shape)) == 1
    if method == "erode":
        if is_box:
            processed = imerode_box(mask, structure.shape[0])
        else:
            processed = ndimage.binary_erosion(mask, structure=structure, border_value=1)
        edge = mask & ~processed
    elif method == "dilate":
        if is_box:
            processed = imdilate_box(mask, structure.shape[0])
        else:
            processed = ndimage.binary_dilation(mask, structure=structure, border_value=0)
        edge = processed & ~mask
    else:
        raise ValueError("method must be 'erode' or 'dilate'")
    index = np.flatnonzero(edge.ravel(order="F"))  # find(): column-major
    edge_point_cloud = np.stack(np.unravel_index(index, edge.shape, order="F"), axis=1).astype(float)
    return edge_point_cloud, processed


def project2closest_surface_points(points: np.ndarray, surf: np.ndarray, surf_center: np.ndarray):
    """Ports project2ClosestSurfacePoints.m: ranks every surface point by
    how closely its direction from `surf_center` aligns (cosine similarity)
    with each query point's direction from `surf_center`.

    As in MATLAB, the directions are normalised in double precision and then
    rounded to single, and the ranking is a stable sort, so equal cosines
    keep the surface points' order (see mask2edge_point_cloud).

    Returns (cosine_sorted, index_on_surf), each of shape
    (n_surface_points, n_points), sorted descending (most-aligned first)
    along axis 0.
    """
    points = np.atleast_2d(np.asarray(points, dtype=float))
    surf = np.asarray(surf, dtype=float)
    surf_center = np.asarray(surf_center, dtype=float)

    vec_p = points - surf_center
    vec_p = (vec_p / np.linalg.norm(vec_p, axis=1, keepdims=True)).astype(np.float32)  # vecP = single(vecP)

    vec_s = surf - surf_center
    vec_s = (vec_s / np.linalg.norm(vec_s, axis=1, keepdims=True)).astype(np.float32)

    # dot(vecP, vecS, 2) in single precision, summed per pair as MATLAB does (not a BLAS product)
    cosine = (vec_s[:, None, :] * vec_p[None, :, :]).sum(axis=2, dtype=np.float32)  # (n_surface, n_points)
    order = np.argsort(-cosine, axis=0, kind="stable")  # sort(distMetric,'descend')
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

    order = np.argsort(dist, axis=0, kind="stable")  # MATLAB's sort() is stable: ties keep point order
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
