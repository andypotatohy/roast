"""Top-level orchestrator: roast_py's equivalent of ROAST's `roast()` entry
point, chaining Phases 1-4 (segmentation -> electrode placement -> meshing
-> FEM solve) into a single call.

Not yet a full replacement for MATLAB's roast(subj, recipe, varargin): it
covers the default disc-electrode, 10-05-cap, single-montage case that
exercises the whole pipeline, not every option MATLAB's roast() accepts
(pad/ring electrodes work via ElectrodeParams but aren't wired into this
function's own keyword arguments yet; neck/custom electrodes, T2-assisted
segmentation, and zero-padding aren't wired in either). Landmarks are
estimated with the interim heuristic in geometry/landmarks.py, not real
landmark detection (not yet ported -- see README).
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import nibabel as nib
import numpy as np

from .fem.pro_writer import Conductivities
from .fem.solve import solve_and_postprocess
from .geometry.cap_info import load_cap_info
from .geometry.landmarks import heuristic_landmarks
from .geometry.placement import ElectrodeParams, electrode_placement
from .meshing.cgal_mesher import mesh_by_iso2mesh
from .segmentation.multiaxial import segment

# Matches ROAST's own default recipe (anode Fp1 1 mA, cathode P4 -1 mA).
DEFAULT_RECIPE = {"Fp1": 1.0, "P4": -1.0}


@dataclass
class RoastResult:
    subj: str
    work_dir: str
    tissue_labels: np.ndarray
    elec_mask: np.ndarray
    gel_mask: np.ndarray
    vol_v: np.ndarray
    vol_e: np.ndarray
    ef_mag: np.ndarray
    voxel_size: np.ndarray
    affine: np.ndarray


def roast(
    subj: str,
    recipe: dict[str, float] | None = None,
    *,
    cap_type: str = "1010",
    elec_type: str = "disc",
    elec_size=(6.0, 2.0),
    conductivities: Conductivities | None = None,
    work_dir: str | None = None,
    maxvol: float = 10.0,
    model_dir=None,
    cgalmesh_bin=None,
    getdp_bin=None,
) -> RoastResult:
    """roast_py's equivalent of `roast(subj, recipe, ...)`.

    `subj` is a path to a T1 NIfTI (must already be RAS-oriented -- run it
    through roast_py.io.nifti.convert_to_ras first if not; see README).
    `recipe` is an electrode-name -> injected-current(mA) dict; currents
    must sum to ~0. Defaults to ROAST's own default recipe.

    Saves voltage/E-field/E-field-magnitude NIfTI outputs next to `subj`
    (or in `work_dir` if given), matching postGetDP.m's outputs, and
    returns them directly (along with the intermediate tissue/electrode/
    gel masks) for inspection.
    """
    if recipe is None:
        recipe = DEFAULT_RECIPE
    total_current = sum(recipe.values())
    if abs(total_current) > 1e-6:
        raise ValueError(f"recipe currents must sum to ~0, got {total_current}")
    if conductivities is None:
        conductivities = Conductivities()

    subj = os.path.abspath(subj)
    base = os.path.splitext(os.path.basename(subj))[0]
    work_dir = os.path.abspath(work_dir) if work_dir else os.path.dirname(subj)
    os.makedirs(work_dir, exist_ok=True)

    t1_img = nib.load(subj)
    if nib.aff2axcodes(t1_img.affine) != ("R", "A", "S"):
        raise ValueError(
            f"{subj} is not RAS-oriented (got {nib.aff2axcodes(t1_img.affine)}). "
            "Run it through roast_py.io.nifti.convert_to_ras first."
        )
    voxel_size = np.abs(np.diag(t1_img.affine)[:3])

    print(f"[1/5] Segmenting {subj} ...")
    mask_path = segment(subj, model_dir=model_dir)
    tissue_img = nib.load(mask_path)
    tissue_labels = np.asarray(tissue_img.dataobj, dtype=np.uint8)

    print("[2/5] Estimating landmarks (interim heuristic -- see geometry/landmarks.py) ...")
    landmarks = heuristic_landmarks(tissue_labels)

    print("[3/5] Placing electrodes ...")
    elec_names = list(recipe.keys())
    current = list(recipe.values())
    cap_names, cap_template = load_cap_info(cap_type)
    elec_paras = [
        ElectrodeParams(elec_type=elec_type, elec_size=np.array(elec_size, dtype=float), cap_type=cap_type)
        for _ in elec_names
    ]
    elec_mask, gel_mask = electrode_placement(
        tissue_labels, landmarks, elec_names, elec_paras, cap_names, cap_template, voxel_size
    )

    print("[4/5] Meshing ...")
    msh_path = os.path.join(work_dir, f"{base}.msh")
    node, elem, _face = mesh_by_iso2mesh(
        tissue_labels, elec_mask, gel_mask, voxel_size, work_dir=work_dir, out_path=msh_path,
        maxvol=maxvol, bin_path=cgalmesh_bin,
    )

    print("[5/5] Solving FEM ...")
    n_elec = len(elec_names)
    if not conductivities.gel:
        conductivities.gel = [0.3] * n_elec
    if not conductivities.electrode:
        conductivities.electrode = [5.9e7] * n_elec
    vol_v, vol_e, ef_mag = solve_and_postprocess(
        work_dir, base, node, elem, elec_names, current, conductivities, voxel_size,
        tissue_labels.shape, bin_path=getdp_bin,
    )

    print("Saving NIfTI outputs ...")
    nib.save(nib.Nifti1Image(vol_v.astype(np.float32), t1_img.affine), os.path.join(work_dir, f"{base}_v.nii"))
    nib.save(nib.Nifti1Image(ef_mag.astype(np.float32), t1_img.affine), os.path.join(work_dir, f"{base}_emag.nii"))
    nib.save(nib.Nifti1Image(vol_e.astype(np.float32), t1_img.affine), os.path.join(work_dir, f"{base}_e.nii"))

    print(f"Done. Outputs saved in {work_dir}")
    return RoastResult(
        subj=subj,
        work_dir=work_dir,
        tissue_labels=tissue_labels,
        elec_mask=elec_mask,
        gel_mask=gel_mask,
        vol_v=vol_v,
        vol_e=vol_e,
        ef_mag=ef_mag,
        voxel_size=voxel_size,
        affine=t1_img.affine,
    )
