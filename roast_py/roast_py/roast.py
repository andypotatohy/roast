"""Top-level orchestrator: roast_py's equivalent of ROAST's `roast()` entry
point, chaining Phases 1-4 (segmentation -> electrode placement -> meshing
-> FEM solve) into a single call.

Not yet a full replacement for MATLAB's roast(subj, recipe, varargin): it
covers the default disc-electrode, 10-05-cap, single-montage case that
exercises the whole pipeline, not every option MATLAB's roast() accepts
(pad/ring electrodes work via ElectrodeParams but aren't wired into this
function's own keyword arguments yet; neck/custom electrodes, T2-assisted
segmentation, and zero-padding aren't wired in either). Landmarks come
from roast.m's TPM landmarks, carried onto the head by a NiftyReg
registration to MNI space (runNiftyReg.m); manual landmark correction
(checkLandmarks.m) isn't ported.

Like MATLAB's roast(), it finishes by visualizing the results (see
roast_py.viz: MRI, segmentation and electrode-placement views plus the
voltage/E-field renderings of visualizeRes.m), and it saves everything
roast_py.viz.review_res() needs to redraw them later.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING

from .dependencies import ensure_dependencies

# Every third-party and pipeline import lives inside roast() rather than
# here, so that importing this module -- and therefore getting hold of the
# roast function at all -- needs nothing but the standard library. That is
# what lets roast() install its own missing dependencies when it is
# called: an import at module scope would fail first and leave the user
# with no way to reach the installer. (It also keeps TensorFlow out of the
# process until it is actually needed, which matters beyond speed:
# segmentation/_keras_compat.py has to set TF_USE_LEGACY_KERAS before
# anything imports tensorflow.)
if TYPE_CHECKING:  # annotations only; `from __future__ import annotations` keeps these lazy
    import numpy as np

    from .fem.pro_writer import Conductivities

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
    recipe: dict[str, float] | None = None
    landmarks: np.ndarray | None = None  # 0-based voxel coords, geometry/landmarks.py order
    mri2mni: np.ndarray | None = None  # 0-based voxel -> MNI mm (roast.m's Affine*image(1).mat)
    mesh_node: np.ndarray | None = None  # mesh node coordinates (mm, see meshing/cgal_mesher.py)
    mesh_elem: np.ndarray | None = None  # tetrahedra: 1-based node ids + region label


def output_paths(work_dir: str, base: str) -> dict[str, str]:
    """Where roast() writes each output for subject `base` in `work_dir`.

    One place for the naming, shared with roast_py.viz.review_res(), which
    reads these back. The MATLAB equivalents carry a simulation tag
    (`<subj>_<simTag>_...`); roast_py doesn't tag simulations yet, so a
    new run on the same subject and work_dir replaces the previous one.
    """
    stem = os.path.join(work_dir, base)
    return {
        "options": f"{stem}_roastOptions.json",
        "elec_mask": f"{stem}_mask_elec.nii",
        "gel_mask": f"{stem}_mask_gel.nii",
        "mesh": f"{stem}_mesh.npz",
        "v_pos": f"{stem}_v.pos",
        "e_pos": f"{stem}_e.pos",
        "v": f"{stem}_v.nii",
        "e": f"{stem}_e.nii",
        "emag": f"{stem}_emag.nii",
    }


def roast(
    subj: str,
    recipe: dict[str, float] | None = None,
    *,
    cap_type: str = "1010",
    elec_type: str = "disc",
    elec_size=(6.0, 2.0),
    conductivities: Conductivities | None = None,
    work_dir: str | None = None,
    mesh_options: dict[str, float] | None = None,
    model_dir=None,
    cgalmesh_bin=None,
    getdp_bin=None,
    niftyreg_bin=None,
    install_missing: bool = True,
    visualize: bool = True,
) -> RoastResult:
    """roast_py's equivalent of `roast(subj, recipe, ...)`.

    `subj` is a path to a T1 NIfTI (must already be RAS-oriented -- run it
    through roast_py.io.nifti.convert_to_ras first if not; see README).
    `recipe` is an electrode-name -> injected-current(mA) dict; currents
    must sum to ~0. Defaults to ROAST's own default recipe.

    Head landmarks (nasion, inion, ears, neck) are found as in MATLAB with
    Multiaxial: the head is registered to the MNI152 template with
    NiftyReg's reg_aladin (a couple of minutes; saved as
    `<subj>_niftyReg.json`), and ROAST's landmarks defined in eTPM.nii are
    mapped through that registration onto the head.

    `mesh_options` is MATLAB's 'meshOptions': any of 'radbound',
    'angbound', 'distbound', 'reratio' and 'maxvol' (iso2mesh's meaning),
    defaulting to roast.m's {radbound: 5, angbound: 30, distbound: 0.3,
    reratio: 3, maxvol: 10}.

    Saves voltage/E-field/E-field-magnitude NIfTI outputs next to `subj`
    (or in `work_dir` if given), matching postGetDP.m's outputs, and
    returns them directly (along with the intermediate tissue/electrode/
    gel masks) for inspection.

    With `visualize` (the default) the results are displayed at the end,
    as MATLAB's roast() does: slice views of the MRI and segmentation, a
    3D view of the electrode placement, and the voltage and E-field on the
    gray matter in 3D and in slices (see roast_py.viz). Windows open once
    the simulation has finished; with no display available the figures
    are saved as PNG files in the work directory instead. Redraw them any
    time later with roast_py.viz.review_res(subj).

    Before any work starts, the running interpreter is checked against
    roast_py's tested environment (exact package versions that ran this
    pipeline end to end, see roast_py.dependencies) and anything missing or
    at a different version is installed with pip. Needs Python 3.11-3.13.
    Set `install_missing=False`, or the ROAST_PY_NO_AUTO_INSTALL
    environment variable, to install nothing: missing packages then raise
    an actionable error and other version differences only warn.
    """
    # Match the tested environment before starting several minutes of
    # work, rather than failing partway through. Pass install_missing=False
    # (or set ROAST_PY_NO_AUTO_INSTALL) to install nothing.
    ensure_dependencies(install_missing=install_missing)

    # Imported here, not at module scope -- see the note at the top of this
    # module. By this point ensure_dependencies() has guaranteed they exist.
    import nibabel as nib
    import numpy as np

    from .fem.pro_writer import Conductivities
    from .fem.solve import solve_and_postprocess
    from .geometry.cap_info import load_cap_info
    from .geometry.landmarks import tpm_landmarks_to_subject
    from .geometry.placement import ElectrodeParams, electrode_placement
    from .meshing.cgal_mesher import mesh_by_iso2mesh, resolve_mesh_options
    from .registration import run_niftyreg

    if recipe is None:
        recipe = DEFAULT_RECIPE
    total_current = sum(recipe.values())
    if abs(total_current) > 1e-6:
        raise ValueError(f"recipe currents must sum to ~0, got {total_current}")
    if conductivities is None:
        conductivities = Conductivities()
    mesh_options = resolve_mesh_options(mesh_options)

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
    from .segmentation.multiaxial import segment  # imports TensorFlow; see note at top

    mask_path = segment(subj, model_dir=model_dir)
    tissue_img = nib.load(mask_path)
    tissue_labels = np.asarray(tissue_img.dataobj, dtype=np.uint8)

    print("[2/5] Registering to MNI space (NiftyReg) and mapping the TPM landmarks ...")
    registration = run_niftyreg(subj, out_dir=work_dir, bin_path=niftyreg_bin)
    landmarks = tpm_landmarks_to_subject(registration.tpm2mri)
    mri2mni = registration.mri2mni

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
        **mesh_options, bin_path=cgalmesh_bin,
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
    paths = output_paths(work_dir, base)
    nib.save(nib.Nifti1Image(vol_v.astype(np.float32), t1_img.affine), paths["v"])
    nib.save(nib.Nifti1Image(ef_mag.astype(np.float32), t1_img.affine), paths["emag"])
    nib.save(nib.Nifti1Image(vol_e.astype(np.float32), t1_img.affine), paths["e"])

    # What review_res() needs to redraw everything later without re-running
    # the pipeline -- the counterparts of MATLAB's _mask_elec.nii,
    # _mask_gel.nii, <subj>_<tag>.mat (mesh) and _roastOptions.mat.
    nib.save(nib.Nifti1Image(elec_mask.astype(np.uint8), t1_img.affine), paths["elec_mask"])
    nib.save(nib.Nifti1Image(gel_mask.astype(np.uint8), t1_img.affine), paths["gel_mask"])
    np.savez_compressed(paths["mesh"], node=node, elem=elem)
    options = {
        "subj": subj,
        "work_dir": work_dir,
        "recipe": {name: float(i) for name, i in recipe.items()},
        "cap_type": cap_type,
        "elec_type": elec_type,
        "elec_size": [float(x) for x in np.ravel(elec_size)],
        "conductivities": asdict(conductivities),
        "mesh_options": mesh_options,
        "masks": os.path.abspath(mask_path),
        "landmarks": np.asarray(landmarks, dtype=float).tolist(),
        # 0-based voxel -> MNI mm; MATLAB's opt.mri2mni is the 1-based equivalent
        "mri2mni": mri2mni.tolist(),
        "Affine": registration.affine.tolist(),
    }
    with open(paths["options"], "w") as f:
        json.dump(options, f, indent=2)

    print(f"Done. Outputs saved in {work_dir}")
    result = RoastResult(
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
        recipe=dict(recipe),
        landmarks=np.asarray(landmarks),
        mri2mni=mri2mni,
        mesh_node=node,
        mesh_elem=elem,
    )

    if visualize:
        # The results are already on disk, so a display problem must not
        # cost the user the simulation: report it and return normally.
        try:
            from .viz.results import show_roast_results

            show_roast_results(result)
        except Exception as e:  # noqa: BLE001
            import warnings

            warnings.warn(
                f"Visualization failed ({type(e).__name__}: {e}). The simulation "
                f"results are saved in {work_dir}; retry the display with "
                f"roast_py.viz.review_res({subj!r}).",
                stacklevel=2,
            )
    return result
