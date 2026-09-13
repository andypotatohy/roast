"""Tests roast_py.meshing.cgal_mesher against the REAL bundled cgalmesh
binary (lib/iso2mesh/bin/cgalmesh.mexa64) -- this is the validation of the
finding that it's a genuine standalone executable, not a MATLAB MEX file
that would need MATLAB to load, despite the filename.
"""

import numpy as np
import pytest

from roast_py.meshing.cgal_mesher import find_cgalmesh_binary, mesh_by_iso2mesh, run_cgal_mesher


def test_find_cgalmesh_binary_locates_bundled_binary():
    path = find_cgalmesh_binary()
    assert path.exists()
    assert path.name in ("cgalmesh.mexa64", "cgalmesh.mexmaci64", "cgalmesh.exe", "cgalmesh")


@pytest.fixture(scope="module")
def solid_sphere_labels():
    shape = (40, 40, 40)
    c = 20
    xx, yy, zz = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing="ij")
    r2 = (xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2
    vol = np.zeros(shape, dtype=np.uint8)
    vol[r2 <= 15**2] = 1  # outer region
    vol[r2 <= 8**2] = 2  # inner "core" region
    return vol


@pytest.mark.slow
def test_run_cgal_mesher_produces_valid_tetrahedral_mesh(tmp_path, solid_sphere_labels):
    node, elem, face = run_cgal_mesher(solid_sphere_labels, tmp_path, maxvol=20)

    assert node.shape[0] > 0
    assert elem.shape[0] > 0
    assert node.shape[1] == 4  # x,y,z,ref
    assert elem.shape[1] == 5  # 4 node refs + region

    # Every element's node references must be valid 1-based indices into node.
    refs = elem[:, :4].astype(int)
    assert refs.min() >= 1
    assert refs.max() <= node.shape[0]

    # Both labeled regions (1 and 2) should show up as element regions.
    regions = set(np.unique(elem[:, 4]).astype(int).tolist())
    assert regions <= {1, 2}
    assert len(regions) == 2

    # Every tetrahedron should have positive volume (no degenerate elements).
    pts = node[refs - 1, :3]  # (n_elem, 4, 3)
    v0, v1, v2, v3 = pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3]
    vol6 = np.abs(np.einsum("ij,ij->i", v1 - v0, np.cross(v2 - v0, v3 - v0)))
    assert np.all(vol6 > 0)


@pytest.mark.slow
def test_mesh_by_iso2mesh_end_to_end_with_electrode_and_gel_regions(tmp_path):
    # cgalmesh assigns consecutive region ids to *sorted distinct nonzero
    # input labels* (see cgal_mesher.py's mesh_by_iso2mesh docstring) --
    # so, to actually exercise the GEL/ELEC numbering meshByIso2mesh.m
    # relies on, this uses a realistic, gap-free label set: all 6 tissue
    # labels present (1-6), then gel at 7, electrode at 8, contiguous.
    shape = (40, 40, 40)
    c = 20
    xx, yy, zz = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing="ij")
    r2 = (xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2

    # 6 concentric "onion" shells for the 6 tissue labels, innermost first.
    tissue = np.zeros(shape, dtype=np.uint8)
    radii = [16, 14, 12, 10, 8, 6]  # outer -> inner, label 1..6
    for label, r in enumerate(radii, start=1):
        tissue[r2 <= r**2] = label
    tissue = 7 - tissue  # innermost shell (smallest r) should be label 6 (air), not 1
    tissue[tissue == 7] = 0  # voxels outside all shells stay background

    gel_mask = np.zeros(shape, dtype=np.uint8)
    gel_mask[r2 <= 4**2] = 1  # one gel region, id 1

    elec_mask = np.zeros(shape, dtype=np.uint8)
    elec_mask[r2 <= 2**2] = 1  # one electrode region, id 1

    tissue[(gel_mask > 0) | (elec_mask > 0)] = 0  # no overlap, mirrors electrode_placement's own cleanup

    out_path = tmp_path / "test.msh"
    node, elem, face = mesh_by_iso2mesh(
        tissue, elec_mask, gel_mask, voxel_size=np.array([1.0, 1.0, 1.0]), work_dir=tmp_path, out_path=out_path, maxvol=20
    )

    assert out_path.exists()
    regions = set(np.unique(elem[:, 4]).astype(int).tolist())
    assert regions == {1, 2, 3, 4, 5, 6, 7, 8}  # all 6 tissues + GEL1 (7) + ELEC1 (8)

    text = out_path.read_text()
    assert "$MeshFormat" in text
    assert "$Nodes" in text
    assert "$Elements" in text
