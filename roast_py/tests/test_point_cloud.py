import nibabel as nib
import numpy as np

from roast_py.geometry.point_cloud import (
    convert_to_ras_point_cloud,
    map2points,
    mask2edge_point_cloud,
    project2closest_surface_points,
)

from .test_nifti import MNI152


def test_mask2edge_point_cloud_erode_returns_shell_of_solid_cube():
    mask = np.zeros((10, 10, 10), dtype=bool)
    mask[2:8, 2:8, 2:8] = True  # 6x6x6 solid cube, 216 voxels

    edge, eroded = mask2edge_point_cloud(mask, "erode", np.ones((3, 3, 3), dtype=bool))

    # Erosion with a 3x3x3 structuring element removes exactly the outer
    # 1-voxel shell, leaving a 4x4x4 interior (64 voxels) -> edge = 216-64.
    assert eroded.sum() == 4**3
    assert edge.shape[0] == 6**3 - 4**3
    # Every edge point must be inside the original cube.
    idx = edge.astype(int)
    assert np.all(mask[idx[:, 0], idx[:, 1], idx[:, 2]])


def test_mask2edge_point_cloud_dilate_returns_shell_just_outside():
    mask = np.zeros((10, 10, 10), dtype=bool)
    mask[4:6, 4:6, 4:6] = True  # 2x2x2 cube, 8 voxels

    edge, dilated = mask2edge_point_cloud(mask, "dilate", np.ones((3, 3, 3), dtype=bool))
    idx = edge.astype(int)
    assert np.all(~mask[idx[:, 0], idx[:, 1], idx[:, 2]])
    assert np.all(dilated[idx[:, 0], idx[:, 1], idx[:, 2]])


def test_project2closest_surface_points_picks_aligned_direction():
    # A ring of surface points around the origin in the xy-plane.
    angles = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    surf = np.stack([np.cos(angles), np.sin(angles), np.zeros_like(angles)], axis=1)
    center = np.array([0.0, 0.0, 0.0])

    # A query point straight out along +x should best align with the
    # surface point at angle 0 (also +x).
    query = np.array([[5.0, 0.0, 0.0]])
    cosine_sorted, order = project2closest_surface_points(query, surf, center)
    best_surface_point = surf[order[0, 0]]
    np.testing.assert_allclose(best_surface_point, [1.0, 0.0, 0.0], atol=1e-6)
    assert cosine_sorted[0, 0] > cosine_sorted[-1, 0]


def test_map2points_closest_and_farthest():
    goal = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]])
    query = np.array([[1.0, 0.0]])

    dist, idx = map2points(query, goal, "closest")
    assert idx[0] == 0

    dist, idx = map2points(query, goal, "farthest")
    assert idx[0] == 2

    dist, idx = map2points(query, goal, "closer", num_of_pts=2)
    np.testing.assert_array_equal(idx[:, 0], [0, 1])


def test_convert_to_ras_point_cloud_matches_convert_to_ras_on_real_mni152():
    # Real regression check: a voxel-index point in the original (LAS)
    # MNI152 volume, converted with convert_to_ras_point_cloud(), should
    # land on the same *world* coordinate as before -- i.e. applying the
    # RAS affine to the transformed point must equal applying the original
    # affine to the original point.
    img = nib.load(MNI152)
    shape = img.shape

    rng = np.random.default_rng(0)
    voxel_pts = rng.integers(0, np.array(shape), size=(20, 3)).astype(float)

    ras_pts, perm = convert_to_ras_point_cloud(img.affine, voxel_pts, shape)

    # Derive the RAS affine the same way roast_py.io.nifti.convert_to_ras does.
    ras_img = nib.as_closest_canonical(img)

    world_before = (img.affine @ np.c_[voxel_pts, np.ones(len(voxel_pts))].T).T
    world_after = (ras_img.affine @ np.c_[ras_pts, np.ones(len(ras_pts))].T).T
    np.testing.assert_allclose(world_before, world_after, atol=1e-4)
