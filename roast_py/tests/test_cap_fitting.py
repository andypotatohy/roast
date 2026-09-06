"""Integration test for fit_cap_to_individual against a synthetic
idealized head (a solid ellipsoid), since there's no real MATLAB output to
compare against in this environment. Landmarks are hand-placed on the
ellipsoid's surface. This checks the algorithm converges and produces
geometrically sane output (points on the scalp surface, roughly
front/back/left/right in the right places), not bit-parity with MATLAB --
that's Phase 5's job once a real MATLAB run is available to compare
against.
"""

import numpy as np
import pytest

from roast_py.geometry.cap_fitting import fit_cap_to_individual
from roast_py.geometry.cap_info import load_cap_info
from roast_py.geometry.point_cloud import mask2edge_point_cloud


@pytest.fixture(scope="module")
def synthetic_head():
    # Axes: 0=left-right (x), 1=posterior-anterior (y), 2=inferior-superior (z).
    shape = (60, 70, 70)
    cx, cy, cz = 30, 35, 35
    rx, ry, rz = 25, 30, 30

    xx, yy, zz = np.meshgrid(
        np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing="ij"
    )
    scalp = (
        ((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2 + ((zz - cz) / rz) ** 2
    ) <= 1.0

    scalp_surface, _ = mask2edge_point_cloud(scalp, "erode", np.ones((3, 3, 3), dtype=bool))

    # Landmarks placed on the ellipsoid surface at (approximately) anterior,
    # posterior, right, and left, each at a distinct height so the
    # sagittal-contour z-matching in fit_cap_to_individual's nasion/inion
    # search can uniquely locate them.
    nasion = np.array([cx, cy + ry, cz - 5], dtype=float)
    inion = np.array([cx, cy - ry, cz + 5], dtype=float)
    right = np.array([cx + rx, cy, cz], dtype=float)
    left = np.array([cx - rx, cy, cz], dtype=float)
    landmarks = np.stack([nasion, inion, right, left, nasion, inion])  # neck landmarks unused here

    return scalp, scalp_surface, landmarks, (cx, cy, cz), (rx, ry, rz)


def test_fit_cap_to_individual_places_points_near_ellipsoid_surface(synthetic_head):
    scalp, scalp_surface, landmarks, cxyz, rxyz = synthetic_head
    names, template = load_cap_info("1010")
    ind_need = np.array([names.index("Fp1"), names.index("Fp2"), names.index("Cz"), names.index("Oz")])

    electrode_coord, center = fit_cap_to_individual(
        scalp,
        scalp_surface,
        landmarks,
        voxel_size=np.array([1.0, 1.0, 1.0]),
        cap_names=names,
        cap_template=template,
        ind_need=ind_need,
        is_biosemi=False,
        is_egi=False,
    )

    assert electrode_coord.shape == (4, 3)

    cx, cy, cz = cxyz
    rx, ry, rz = rxyz
    # Every placed electrode must land (approximately) ON the ellipsoid
    # surface -- i.e. the normalized quadratic form should be close to 1 --
    # since fit_cap_to_individual only ever returns points taken directly
    # from scalp_surface.
    q = (
        ((electrode_coord[:, 0] - cx) / rx) ** 2
        + ((electrode_coord[:, 1] - cy) / ry) ** 2
        + ((electrode_coord[:, 2] - cz) / rz) ** 2
    )
    np.testing.assert_allclose(q, 1.0, atol=0.15)

    # The optimized cap center should be within the head, roughly near the
    # ellipsoid's own center.
    assert np.linalg.norm(center - np.array(cxyz)) < max(rxyz)


def test_fit_cap_to_individual_fp1_left_of_fp2(synthetic_head):
    # Fp1 is over the left hemisphere, Fp2 over the right, in the 10-05
    # system -- a basic sanity check that left/right didn't get flipped.
    scalp, scalp_surface, landmarks, cxyz, rxyz = synthetic_head
    names, template = load_cap_info("1010")
    ind_need = np.array([names.index("Fp1"), names.index("Fp2")])

    electrode_coord, _ = fit_cap_to_individual(
        scalp,
        scalp_surface,
        landmarks,
        voxel_size=np.array([1.0, 1.0, 1.0]),
        cap_names=names,
        cap_template=template,
        ind_need=ind_need,
    )
    fp1, fp2 = electrode_coord
    assert fp1[0] < fp2[0]
