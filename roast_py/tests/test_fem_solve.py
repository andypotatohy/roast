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
from roast_py.fem.prepare import prepare_for_getdp
from roast_py.fem.pro_writer import Conductivities, write_pro_file
from roast_py.meshing.cgal_mesher import mesh_by_iso2mesh


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
