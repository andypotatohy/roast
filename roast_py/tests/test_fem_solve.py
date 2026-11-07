"""End-to-end test of the FEM solve pipeline against the REAL getdp
binary: mesh -> boundary extraction -> .pro generation -> getdp solve ->
.pos parsing -> voxel-grid interpolation, on a small synthetic two-
electrode "head" (concentric tissue shells + a small air pocket, so all 6
tissue labels are present with no gaps -- see cgal_mesher.py's docstring
on why that matters -- plus disc-shaped gel/electrode pads poking outward
on opposite sides, closer to the real electrode geometry
placeAndModelElectrodes.m produces than a fully gel-enclosed electrode).

Not just a plumbing check: asserts the solved potential is actually higher
near the current-injecting anode than near the current-sinking cathode,
which only holds if the boundary conditions, region wiring, and area
computation all have the right sign/magnitude/association.
"""

import numpy as np
import pytest

from roast_py.fem.getdp_runner import find_getdp_binary, run_getdp
from roast_py.fem.pos_parser import interpolate_to_grid, read_pos_node_table
from roast_py.fem.prepare import attached_to_electrodes, prepare_for_getdp
from roast_py.fem.pro_writer import Conductivities, write_pro_file
from roast_py.meshing.cgal_mesher import mesh_by_iso2mesh
from roast_py.meshing.mesh_io import save_msh


def test_find_getdp_binary_locates_bundled_binary():
    path = find_getdp_binary()
    assert path.exists()


@pytest.fixture(scope="module")
def two_electrode_sphere_mesh(tmp_path_factory):
    shape = (56, 56, 56)
    c = 27
    xx, yy, zz = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing="ij")
    r2 = (xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2

    tissue = np.zeros(shape, dtype=np.uint8)
    tissue[r2 <= 20**2] = 5  # skin (outermost, touches gel)
    tissue[r2 <= 18**2] = 4  # bone
    tissue[r2 <= 16**2] = 3  # csf
    tissue[r2 <= 14**2] = 2  # gray
    tissue[r2 <= 12**2] = 1  # white
    tissue[r2 <= 3**2] = 6  # tiny air pocket at the very center

    def cylinder_mask(x_lo, x_hi, y_radius):
        m = np.zeros(shape, dtype=bool)
        for xv in range(x_lo, x_hi + 1):
            m[xv] = (yy[xv] - c) ** 2 + (zz[xv] - c) ** 2 <= y_radius**2
        return m

    gel_mask = np.zeros(shape, dtype=np.uint8)
    elec_mask = np.zeros(shape, dtype=np.uint8)
    # anode at +x: gel pad just outside skin, electrode metal further out.
    gel_mask[cylinder_mask(c + 21, c + 22, 5)] = 1
    elec_mask[cylinder_mask(c + 23, c + 25, 5)] = 1
    # cathode at -x, mirrored.
    gel_mask[cylinder_mask(c - 22, c - 21, 5)] = 2
    elec_mask[cylinder_mask(c - 25, c - 23, 5)] = 2

    tissue[(gel_mask > 0) | (elec_mask > 0)] = 0

    work_dir = tmp_path_factory.mktemp("fem_solve")
    node, elem, face = mesh_by_iso2mesh(
        tissue, elec_mask, gel_mask, voxel_size=np.array([1.0, 1.0, 1.0]),
        work_dir=work_dir, out_path=work_dir / "subj_test.msh", maxvol=20,
    )
    return work_dir, node, elem, shape, c


def test_attached_to_electrodes_finds_floating_pieces():
    # Two tetrahedra sharing a face (skin + gel), one sharing a single node
    # with them (still attached), and one sharing nothing (floating skin).
    elem = np.array([
        [1, 2, 3, 4, 5],
        [2, 3, 4, 5, 7],
        [5, 6, 8, 9, 5],
        [10, 11, 12, 13, 5],
    ])
    assert attached_to_electrodes(elem).tolist() == [True, True, True, False]


@pytest.mark.slow
def test_floating_pieces_left_out_of_solver_mesh(two_electrode_sphere_mesh, tmp_path):
    # runExamples cmd 12: specks of skin floating free of the head made
    # MUMPS stop with "numerically singular matrix". They are left out.
    work_dir, node, elem, shape, c = two_electrode_sphere_mesh
    far = node[:, :3].max(axis=0) + 10
    speck = far + np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2], [2, 2, 2]], dtype=float)
    n0 = node.shape[0]
    node2 = np.vstack([node[:, :3], speck])
    elem2 = np.vstack([elem[:, :5], [[n0 + 1, n0 + 2, n0 + 3, n0 + 4, 5], [n0 + 2, n0 + 3, n0 + 4, n0 + 5, 5]]])
    save_msh(node2, elem2, str(tmp_path / "s.msh"))

    faces, area = prepare_for_getdp(str(tmp_path / "s.msh"), str(tmp_path / "s_ready.msh"), node2, elem2, ["a", "b"])
    text = (tmp_path / "s_ready.msh").read_text()
    block = text.split("$Elements\n", 1)[1].split("$EndElements", 1)[0].splitlines()
    assert int(block[0]) == len(elem) + sum(f.shape[0] for f in faces) == len(block) - 1
    assert [int(line.split()[0]) for line in block[1:]] == list(range(1, len(block)))
    used = {int(v) for line in block[1:] for v in line.split()[5:]}
    assert max(used) == n0  # no element uses the speck's nodes

    sigma = Conductivities(gel=[0.3, 0.3], electrode=[5.9e7, 5.9e7])
    write_pro_file(str(tmp_path / "s.pro"), current=[1.0, -1.0], sigma=sigma, area_elec_needed=area, ind_use=[1, 2])
    run_getdp(str(tmp_path / "s.pro"), str(tmp_path / "s_ready.msh"))
    v_ids, v_vals = read_pos_node_table(str(tmp_path / "s_v.pos"), 1)
    assert v_ids.size > 0 and not np.any(np.isnan(v_vals))


@pytest.mark.slow
def test_prepare_for_getdp_and_getdp_solve_end_to_end(two_electrode_sphere_mesh):
    work_dir, node, elem, shape, c = two_electrode_sphere_mesh
    elec_names = ["anode", "cathode"]

    element_elec_needed, area_elec_needed = prepare_for_getdp(
        str(work_dir / "subj_test.msh"), str(work_dir / "subj_test_ready.msh"), node, elem, elec_names
    )
    assert all(a > 0 for a in area_elec_needed)
    assert len(element_elec_needed) == 2

    sigma = Conductivities(gel=[0.3, 0.3], electrode=[5.9e7, 5.9e7])
    write_pro_file(
        str(work_dir / "subj_test.pro"),
        current=[1.0, -1.0],
        sigma=sigma,
        area_elec_needed=area_elec_needed,
        ind_use=[1, 2],
    )

    run_getdp(str(work_dir / "subj_test.pro"), str(work_dir / "subj_test_ready.msh"))

    assert (work_dir / "subj_test_v.pos").exists()
    assert (work_dir / "subj_test_e.pos").exists()
    # getdp's own intermediate files should have been cleaned up.
    assert not (work_dir / "subj_test.pre").exists()
    assert not (work_dir / "subj_test.res").exists()

    v_ids, v_vals = read_pos_node_table(str(work_dir / "subj_test_v.pos"), 1)
    e_ids, e_vals = read_pos_node_table(str(work_dir / "subj_test_e.pos"), 3)
    assert v_ids.shape[0] == node.shape[0]
    assert e_ids.shape[0] == node.shape[0]
    assert not np.any(np.isnan(v_vals))
    assert not np.any(np.isnan(e_vals))

    v_vals = v_vals[:, 0] - v_vals[:, 0].min()
    vol_v = interpolate_to_grid(node[v_ids - 1, :3], v_vals, shape)

    # Physical sanity: higher potential near the current-injecting anode
    # (+x) than near the current-sinking cathode (-x) -- only holds if
    # the boundary condition sign, region association, and area
    # computation all line up correctly.
    v_near_anode = vol_v[c + 18, c, c]
    v_near_cathode = vol_v[c - 18, c, c]
    assert not np.isnan(v_near_anode)
    assert not np.isnan(v_near_cathode)
    assert v_near_anode > v_near_cathode

    ef_mag = np.sqrt(np.sum(interpolate_to_grid(node[e_ids - 1, :3], e_vals, shape) ** 2, axis=-1))
    assert np.nanmax(ef_mag) > 0


def test_getdp_killed_for_memory_says_so(tmp_path):
    """A SIGKILLed getDP (the OOM killer's signature) gets an explanation,
    not the generic 'cannot work properly on your system'."""
    import sys

    if sys.platform.startswith("win"):
        pytest.skip("POSIX signals")
    fake = tmp_path / "getdp"
    fake.write_text("#!/bin/sh\necho 'Info    : Solve[Sys_Ele]'\nkill -9 $$\n")
    fake.chmod(0o755)
    (tmp_path / "x.pro").write_text("")
    (tmp_path / "x_ready.msh").write_text("")
    with pytest.raises(RuntimeError, match=r"(?s)killed while solving.*out of memory.*Solve\[Sys_Ele\]"):
        run_getdp(tmp_path / "x.pro", tmp_path / "x_ready.msh", bin_path=fake)


def test_getdp_out_of_memory_in_mumps_says_so(tmp_path):
    """MUMPS failing to allocate the factors (INFOG(1)=-13, what subject1
    hits on a 14 GB machine) is reported as running out of memory. The solve
    passes -ksp_error_if_not_converged: without it getdp writes NaN and exits 0."""
    import sys

    if sys.platform.startswith("win"):
        pytest.skip("shell script stand-in for getdp")
    fake = tmp_path / "getdp"
    fake.write_text(
        "#!/bin/sh\necho \"$@\" > args.txt\n"
        "echo '[0]PETSC ERROR: Error reported by MUMPS in numerical factorization phase: "
        "INFOG(1)=-13, INFO(2)=1166537700' >&2\nexit 1\n"
    )
    fake.chmod(0o755)
    (tmp_path / "x.pro").write_text("")
    (tmp_path / "x_ready.msh").write_text("")
    with pytest.raises(RuntimeError, match=r"(?s)ran out of memory.*INFOG\(1\)=-13"):
        run_getdp(tmp_path / "x.pro", tmp_path / "x_ready.msh", bin_path=fake)
    assert "-ksp_error_if_not_converged" in (tmp_path / "args.txt").read_text()
