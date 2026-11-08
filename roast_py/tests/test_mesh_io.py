import numpy as np

from roast_py.meshing.mesh_io import read_medit, save_inr, save_msh, sort_mesh


def test_save_inr_header_is_exactly_256_bytes_and_round_trips_shape(tmp_path):
    vol = (np.random.default_rng(0).integers(0, 7, size=(5, 6, 7))).astype(np.uint8)
    path = tmp_path / "test.inr"
    save_inr(vol, str(path))

    data = path.read_bytes()
    assert len(data) == 256 + vol.size
    header = data[:256].decode("ascii")
    assert header.startswith("#INRIMAGE-4#{\n")
    assert "XDIM=5\nYDIM=6\nZDIM=7\n" in header
    assert header.endswith("##}\n")

    # Fortran-order (x fastest) payload, matching MATLAB's native
    # column-major fwrite.
    payload = np.frombuffer(data[256:], dtype=np.uint8).reshape(vol.shape, order="F")
    np.testing.assert_array_equal(payload, vol)


def test_read_medit_parses_vertices_triangles_tetrahedra(tmp_path):
    # A minimal, hand-written MEDIT file mimicking cgalmesh's output shape,
    # including the MeshVersionFormatted/Dimension header lines that
    # read_medit must skip via its generic single-int consumption.
    content = """MeshVersionFormatted 1
Dimension 3
Vertices
2
0.0 0.0 0.0 1
1.0 2.0 3.0 2
Triangles
1
1 2 2 5
Tetrahedra
1
1 2 2 2 9
End
"""
    path = tmp_path / "test.mesh"
    path.write_text(content)

    node, elem, face = read_medit(str(path))
    np.testing.assert_allclose(node, [[0, 0, 0, 1], [1, 2, 3, 2]])
    np.testing.assert_allclose(face, [[1, 2, 2, 5]])
    np.testing.assert_allclose(elem, [[1, 2, 2, 2, 9]])


def test_sort_mesh_preserves_node_set_and_remaps_elem_consistently():
    node = np.array(
        [
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 0.0, 5.0],
        ]
    )
    # One tetrahedron using all 4 nodes (1-based references).
    elem = np.array([[1, 2, 3, 4, 1]])

    no, el, _ = sort_mesh(None, node, elem, [0, 1, 2, 3])

    # Same set of node coordinates, just possibly reordered.
    orig_set = {tuple(row) for row in node}
    new_set = {tuple(row) for row in no}
    assert orig_set == new_set

    # The single tetrahedron's remapped node references, read back against
    # the *sorted* node array, must reproduce the same 4 physical points.
    referenced_pts = {tuple(no[idx - 1]) for idx in el[0, :4].astype(int)}
    assert referenced_pts == orig_set


def test_save_msh_writes_expected_gmsh_sections(tmp_path):
    node = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    elem = np.array([[1, 2, 3, 4, 1]])
    path = tmp_path / "test.msh"
    save_msh(node, elem, str(path), ["WHITE"])

    text = path.read_text()
    assert "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n" in text
    assert "$Nodes\n4\n" in text
    assert "1 0.0000000000 0.0000000000 0.0000000000\n" in text
    assert "$Elements\n1\n" in text
    assert "1 4 2 1 1 1 2 3 4\n" in text
