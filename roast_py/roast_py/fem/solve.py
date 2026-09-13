"""Top-level orchestrator porting solveByGetDP.m + postGetDP.m's roast()
(non-lead-field) path: prepares boundary elements, writes the .pro file,
runs getdp, and interpolates the result onto the MRI voxel grid.

roast_target()'s lead-field path (looping solveByGetDP.m once per
electrode with LFtag set) is deferred to Phase 6 (targeting) -- the
building blocks here (write_pro_file's ind_use/lf_tag, run_getdp,
read_pos_node_table) are already general enough to reuse for it without
rework.
"""

from __future__ import annotations

import os

import numpy as np

from .getdp_runner import run_getdp
from .pos_parser import interpolate_to_grid, read_pos_node_table
from .prepare import prepare_for_getdp
from .pro_writer import Conductivities, write_pro_file


def solve_and_postprocess(
    work_dir: str,
    subj_tag: str,
    node: np.ndarray,
    elem: np.ndarray,
    elec_names: list[str],
    current: list[float],
    sigma: Conductivities,
    voxel_size: np.ndarray,
    grid_shape: tuple[int, int, int],
    bin_path: str | os.PathLike | None = None,
):
    """Ports the roast() (non-lead-field) path of solveByGetDP.m +
    postGetDP.m. `node`/`elem` are roast_py.meshing.cgal_mesher.mesh_by_iso2mesh's
    output (node already in physical/mm space); `work_dir` must contain
    `{subj_tag}.msh` from that same call. `current` is in mA, one value
    per electrode in `elec_names`' order (must sum to ~0).

    Returns (vol_v, vol_e, ef_mag): voltage (mV) and E-field (V/m, 3
    components) volumes on the `grid_shape` voxel grid, plus E-field
    magnitude -- matching postGetDP.m's `vol_all`/`ef_all`/`ef_mag`.
    """
    msh_path = os.path.join(work_dir, f"{subj_tag}.msh")
    ready_msh_path = os.path.join(work_dir, f"{subj_tag}_ready.msh")
    pro_path = os.path.join(work_dir, f"{subj_tag}.pro")

    _element_elec_needed, area_elec_needed = prepare_for_getdp(msh_path, ready_msh_path, node, elem, elec_names)

    ind_use = list(range(1, len(elec_names) + 1))
    write_pro_file(pro_path, current, sigma, area_elec_needed, ind_use, lf_tag="")

    run_getdp(pro_path, ready_msh_path, bin_path=bin_path)

    v_path = os.path.join(work_dir, f"{subj_tag}_v.pos")
    e_path = os.path.join(work_dir, f"{subj_tag}_e.pos")
    v_ids, v_vals = read_pos_node_table(v_path, 1)
    e_ids, e_vals = read_pos_node_table(e_path, 3)

    if v_ids.size == 0 or e_ids.size == 0:
        raise RuntimeError("getDP did not converge. Please check getDP before proceeding.")

    v_vals = v_vals[:, 0] - v_vals[:, 0].min()  # re-reference voltage, matching postGetDP.m

    # Mesh node coordinates are in physical (mm) space (mesh_by_iso2mesh's
    # output); convert back to voxel-index space for interpolation onto
    # the MRI grid, matching postGetDP.m's `node(:,i)/imgHdr.mat(i,i)`.
    node_voxel = node[:, :3] / np.asarray(voxel_size)

    vol_v = interpolate_to_grid(node_voxel[v_ids - 1], v_vals, grid_shape)
    if np.all(np.isnan(vol_v)):
        raise RuntimeError("getDP did not converge. Please check getDP before proceeding.")

    vol_e = interpolate_to_grid(node_voxel[e_ids - 1], e_vals, grid_shape)
    if np.all(np.isnan(vol_e)):
        raise RuntimeError("getDP did not converge. Please check getDP before proceeding.")

    ef_mag = np.sqrt(np.sum(vol_e**2, axis=-1))
    return vol_v, vol_e, ef_mag
