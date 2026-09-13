"""Tests roast_py.fem.boundary against a small, hand-built conforming
tetrahedral mesh: two unit cubes stacked along z (a "gel" cube z in [0,1]
and an "electrode" cube z in [1,2]), sharing their interface nodes exactly
-- the same conforming-mesh structure a real gel/electrode mesh region
pair has.
"""

import numpy as np
import pytest

from roast_py.fem.boundary import extract_electrode_outer_surface, free_boundary_faces, triangle_areas

# Fan decomposition of a cube (corners in the order 000,100,110,010,001,101,111,011)
# into 6 tets around the v0-v6 diagonal.
_CUBE_TETS_LOCAL = [
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
]


def _cube_corners(x0, y0, z0):
    return np.array(
        [
            [x0, y0, z0],
            [x0 + 1, y0, z0],
            [x0 + 1, y0 + 1, z0],
            [x0, y0 + 1, z0],
            [x0, y0, z0 + 1],
            [x0 + 1, y0, z0 + 1],
            [x0 + 1, y0 + 1, z0 + 1],
            [x0, y0 + 1, z0 + 1],
        ]
    )


@pytest.fixture
def stacked_cubes():
    # Gel cube: global nodes 1-8 (1-based). Electrode cube: shares its
    # bottom face (nodes 5,6,7,8, the gel cube's top face) and adds new
    # nodes 9-12 for its own top face.
    gel_corners = _cube_corners(0, 0, 0)
    elec_top = _cube_corners(0, 0, 1)[4:]  # just the new top-face corners

    node = np.concatenate([gel_corners, elec_top], axis=0)  # 12 points, 0-based rows

    gel_local = np.arange(8)  # nodes 1-8
    elec_local = np.array([4, 5, 6, 7, 8, 9, 10, 11])  # nodes 5-12

    gel_tets = np.array([[gel_local[i] for i in tet] for tet in _CUBE_TETS_LOCAL]) + 1
    elec_tets = np.array([[elec_local[i] for i in tet] for tet in _CUBE_TETS_LOCAL]) + 1

    return node, gel_tets, elec_tets


def test_free_boundary_faces_of_single_cube_is_12_triangles(stacked_cubes):
    node, gel_tets, _ = stacked_cubes
    faces = free_boundary_faces(gel_tets)
    assert faces.shape == (12, 3)  # 6 faces x 2 triangles, all free (no internal tet touches another tet here... )
    total_area = triangle_areas(faces, node).sum()
    assert np.isclose(total_area, 6.0)  # unit cube surface area


def test_extract_electrode_outer_surface_excludes_gel_interface(stacked_cubes):
    node, gel_tets, elec_tets = stacked_cubes
    outer_faces, area = extract_electrode_outer_surface(gel_tets, elec_tets, node)

    # 5 of the electrode cube's 6 faces (all but the gel-facing bottom
    # face) survive -> 10 triangles, area 5.
    assert outer_faces.shape == (10, 3)
    assert np.isclose(area, 5.0)

    # None of the surviving faces should be entirely on the z=1 interface plane.
    pts = node[outer_faces - 1]  # (n, 3, 3)
    on_interface = np.all(np.isclose(pts[:, :, 2], 1.0), axis=1)
    assert not np.any(on_interface)


def test_extract_electrode_outer_surface_raises_on_empty_gel(stacked_cubes):
    node, gel_tets, elec_tets = stacked_cubes
    with pytest.raises(ValueError, match="Gel under this electrode"):
        extract_electrode_outer_surface(np.empty((0, 4), dtype=int), elec_tets, node)


def test_extract_electrode_outer_surface_raises_on_empty_electrode(stacked_cubes):
    node, gel_tets, elec_tets = stacked_cubes
    with pytest.raises(ValueError, match="Electrode was not meshed"):
        extract_electrode_outer_surface(gel_tets, np.empty((0, 4), dtype=int), node)
