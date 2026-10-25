"""Top-level orchestrator: roast_py's equivalent of ROAST's `roast()` entry
point, chaining segmentation -> registration/landmarks -> electrode
placement -> meshing -> FEM solve into a single call.

Takes MATLAB roast()'s options as keyword arguments (see roast()'s
docstring and roast_py/options.py for the Python spellings). Implemented:
capType, elecType, elecSize, elecOri, T2, multiaxial, meshOptions,
conductivities, simulationTag, resampling and zeroPadding, plus neck and
custom electrodes, and automatic conversion to RAS. Not ported yet, and
refused with NotImplementedError: the New York head, lead-field generation
and the manual landmark GUI.

Segmentation is Multiaxial or SPM; SPM, resampling and T2 alignment run
through SPM's official Python interface (roast_py.spm), which installs the
MATLAB Runtime it needs on first use, after asking to accept its license.
Landmarks come from roast.m's TPM landmarks, carried onto the head by the
registration to MNI space: NiftyReg's (runNiftyReg.m) with Multiaxial,
SPM's own (`_seg8.mat`) with SPM, as in MATLAB.

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
    sim_tag: str | None = None  # simulationTag; None = untagged output names
    model_mri: str | None = None  # the MRI actually modeled (after RAS conversion / zero-padding)


def simulation_stem(base: str, sim_tag: str | None = None) -> str:
    """`<subj>_<simTag>` (MATLAB's naming), or just `<subj>` untagged."""
    return f"{base}_{sim_tag}" if sim_tag else base


def output_paths(work_dir: str, base: str, sim_tag: str | None = None) -> dict[str, str]:
    """Where roast() writes each output for subject `base` in `work_dir`.

    One place for the naming, shared with roast_py.viz.review_res(), which
    reads these back: MATLAB's `<subj>_<simTag>_...`. roast() always has a
    tag (the user's, or a date-time stamp; see roast_py.simulation_log);
    `sim_tag=None` gives the untagged `<subj>_...` names that roast_py
    versions before tags wrote.
    """
    stem = os.path.join(work_dir, simulation_stem(base, sim_tag))
    return {
        "options": f"{stem}_roastOptions.json",
        "elec_mask": f"{stem}_mask_elec.nii",
        "gel_mask": f"{stem}_mask_gel.nii",
        "mesh": f"{stem}_mesh.npz",
        "msh": f"{stem}.msh",
        "v_pos": f"{stem}_v.pos",
        "e_pos": f"{stem}_e.pos",
        "v": f"{stem}_v.nii",
        "e": f"{stem}_e.nii",
        "emag": f"{stem}_emag.nii",
    }


def _multiaxial_masks_path(model_mri: str) -> str:
    """`<mri>_multiaxial_masks.nii` next to the modeled MRI, as
    segmentation.multiaxial.segment() names it (computed here so a cached
    segmentation doesn't need TensorFlow imported)."""
    base = os.path.splitext(os.path.basename(model_mri))[0]
    return os.path.join(os.path.dirname(model_mri), f"{base}_multiaxial_masks.nii")


def _write_atomically(path: str, write) -> None:
    """`write(tmp)` to a temporary name next to `path`, then rename it into
    place: an interrupted run never leaves a half-written file that the
    next run's "already done, skipping" check would trust."""
    tmp = os.path.join(os.path.dirname(path), f".partial_{os.path.basename(path)}")
    try:
        write(tmp)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _mri_base(path: str) -> str:
    name = os.path.basename(path)
    for ext in (".gz", ".nii", ".hdr", ".img"):
        name = name.removesuffix(ext)
    return name


def _default_subject() -> str:
    """example/MNI152_T1_1mm.nii, MATLAB roast()'s default subject."""
    from pathlib import Path

    for parent in Path(__file__).resolve().parents:
        candidate = parent / "example" / "MNI152_T1_1mm.nii"
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("Could not locate example/MNI152_T1_1mm.nii, roast()'s default subject.")


def roast(
    subj: str | None = None,
    recipe=None,
    *,
    cap_type: str = "1010",
    elec_type="disc",
    elec_size=None,
    elec_ori=None,
    t2: str | None = None,
    multiaxial=None,
    manual_gui=False,
    mesh_options: dict[str, float] | None = None,
    conductivities=None,
    simulation_tag: str | None = None,
    resampling=False,
    zero_padding: int = 0,
    work_dir: str | None = None,
    model_dir=None,
    cgalmesh_bin=None,
    getdp_bin=None,
    niftyreg_bin=None,
    install_missing: bool = True,
    accept_matlab_runtime_license: bool | None = None,
    visualize: bool = True,
) -> RoastResult:
    """roast_py's equivalent of MATLAB's `roast(subj, recipe, ...)`.

    Arguments follow MATLAB's, as Python keywords: `roast(subj, recipe,
    'elecType', 'pad', 'zeroPadding', 20)` is `roast(subj, recipe,
    elec_type="pad", zero_padding=20)`. See roast_py/options.py for how
    MATLAB's cell arrays, matrices and structs are spelled.

    `subj` is a T1 NIfTI; default example/MNI152_T1_1mm.nii, as in MATLAB.
    A non-RAS image is converted to RAS first (saved as `<subj>_ras.nii`).
    `recipe` is a dict of electrode name -> current (mA), summing to 0;
    default {'Fp1': 1, 'P4': -1}. Neck electrodes are 'Nk1'..'Nk4';
    custom ones ('custom...') are read from `<subj>_customLocations`.

    Segmentation: `multiaxial=False` (or 'off') segments with SPM, as do
    calls with a `t2` image (SPM uses both; the T2 is first aligned to the
    T1 if needed). Otherwise roast_py segments with Multiaxial -- its
    default, unlike MATLAB's, whose default is SPM. `resampling=True`
    reslices the MRI to 1 mm voxels with SPM first. With Multiaxial,
    landmarks come from a NiftyReg registration; with SPM, from SPM's own.

    SPM needs the MATLAB Runtime R2025b (free, several GB). If it isn't
    installed, roast() installs it before doing any other work, once its
    license is accepted: at an interactive prompt, or with
    `accept_matlab_runtime_license=True` or
    ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE=yes (see roast_py.spm.runtime).

    Not ported, and raising NotImplementedError: `manual_gui`, the New
    York head and recipe='leadfield'.

    `mesh_options` is MATLAB's 'meshOptions': any of 'radbound',
    'angbound', 'distbound', 'reratio' and 'maxvol' (iso2mesh's meaning),
    defaulting to roast.m's {radbound: 5, angbound: 30, distbound: 0.3,
    reratio: 3, maxvol: 10}.

    Head landmarks (nasion, inion, ears, neck) are found as in MATLAB with
    Multiaxial: the head is registered to the MNI152 template with
    NiftyReg's reg_aladin (a couple of minutes; saved as
    `<subj>_niftyReg.json`), and ROAST's landmarks defined in eTPM.nii are
    mapped through that registration onto the head.

    Saves voltage/E-field/E-field-magnitude NIfTI outputs next to `subj`
    (or in `work_dir` if given), named `<subj>_<tag>_v.nii` etc., and
    returns them directly (along with the intermediate tissue/electrode/gel
    masks) for inspection. Results are in the space of the modeled MRI
    (`RoastResult.model_mri`), i.e. after RAS conversion and zero-padding,
    as in MATLAB.

    Simulations and their tags work as in MATLAB. `<tag>` is
    `simulation_tag`, or by default the date and time
    (`20261025T143012`, MATLAB's `datestr(now,30)`). The options are saved
    as `<subj>_<tag>_roastOptions.json` and logged in `<subj>_roastLog`.
    If a simulation of the same MRI was already saved with the same
    options, roast() carries on under its tag instead, warning if you gave
    a different one. Reusing a tag for different options is an error. Each
    of the six steps (segmentation, registration, electrode placement,
    meshing, solving, saving the results) is skipped if its output already
    exists. A repeated call therefore just loads the results, and an
    interrupted one resumes where it stopped. Delete a step's output to
    redo it; the steps after a step that ran are redone too, which
    MATLAB doesn't do. `RoastResult.sim_tag` is the tag used.

    With `visualize` (the default) the results are displayed at the end,
    as MATLAB's roast() does: slice views of the MRI and segmentation, a
    3D view of the electrode placement, and the voltage and E-field on the
    gray matter in 3D and in slices (see roast_py.viz). Windows open once
    the simulation has finished; with no display available the figures
    are saved as PNG files in the work directory instead. Redraw them any
    time later with roast_py.review_res(subj, sim_tag=...).

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
    import warnings

    import nibabel as nib
    import numpy as np

    from .fem.solve import postprocess, solution_exists, solve
    from .geometry.cap_info import load_cap_info
    from .geometry.landmarks import tpm_landmarks_to_subject
    from .geometry.placement import electrode_placement
    from .io.nifti import clear_padding, convert_to_ras, zero_pad
    from .meshing.cgal_mesher import mesh_by_iso2mesh, resolve_mesh_options
    from .options import (
        on_off,
        read_custom_locations,
        resolve_conductivities,
        resolve_electrodes,
        resolve_recipe,
    )
    from .registration import Registration, registration_path, run_niftyreg
    from .simulation_log import check_tag, register_simulation
    from .simulation_log import log_path as sim_log_path

    # -- check inputs, in roast.m's order --------------------------------------
    if not subj:
        subj = _default_subject()
    if str(subj).lower() in ("nyhead", "example/nyhead.nii"):
        raise NotImplementedError(
            "The New York head is not ported to roast_py yet: it is modeled from its bundled SPM "
            "segmentation and _seg8.mat registration rather than an MRI (and needs >50 GB of memory)."
        )
    subj = os.path.abspath(subj)
    if not os.path.exists(subj):
        raise FileNotFoundError(f"The subject MRI you provided {subj} does not exist.")

    if recipe is None or (not isinstance(recipe, str) and len(recipe) == 0):
        recipe = DEFAULT_RECIPE
    elec_names, current = resolve_recipe(recipe)
    elec_paras = resolve_electrodes(len(elec_names), cap_type, elec_type, elec_size, elec_ori)
    cap_type = cap_type.lower()

    if t2 is not None:
        if not os.path.exists(t2):
            raise FileNotFoundError(f"The T2 MRI you provided {t2} does not exist.")
    multiaxial_on = None if multiaxial is None else on_off(multiaxial, "multiaxial")
    if t2 is not None and multiaxial_on:
        raise ValueError(
            "Multiaxial cannot be run with both T1 and T2 images. If you meant to use Multiaxial, please "
            "only provide T1 image with option 'multiaxial' turned on."
        )
    # Explicitly off, or a T2 given: SPM. Otherwise roast_py's default, Multiaxial.
    use_spm = multiaxial_on is False or t2 is not None
    if on_off(manual_gui, "manualGui"):
        raise NotImplementedError("The manual landmark GUI (checkLandmarks.m) is not ported to roast_py yet.")

    if mesh_options is not None:
        warnings.warn(
            "You're changing the advanced options of ROAST. Unless you know what you're doing, please keep "
            "mesh options default.",
            stacklevel=2,
        )
    mesh_options = resolve_mesh_options(mesh_options)
    conductivities = resolve_conductivities(conductivities, len(elec_names))
    check_tag(simulation_tag)

    do_resample = on_off(resampling, "resampling")
    if zero_padding is None:
        zero_padding = 0
    if zero_padding < 0 or zero_padding != int(zero_padding):
        raise ValueError(
            "Unrecognized option value. Please enter positive integer value for option 'zeroPadding'. "
            "A recommended value is 10."
        )
    zero_padding = int(zero_padding)

    # -- preprocess the MRI: RAS, then zero-padding (roast.m) ----------------------
    base = _mri_base(subj)
    work_dir = os.path.abspath(work_dir) if work_dir else os.path.dirname(subj)
    os.makedirs(work_dir, exist_ok=True)

    original_img = nib.load(subj)
    pixdim = np.asarray(original_img.header.get_zooms()[:3], dtype=float)
    if np.any(pixdim < 0.8) and not do_resample:
        warnings.warn(
            "The MRI has higher resolution (<0.8mm) in at least one direction. This will make the modeling "
            "process more computationally expensive and thus slower. If you wish to run faster using just "
            "1-mm model, you can ask ROAST to re-sample the MRI into 1 mm first, by turning on the "
            "'resampling' option.",
            stacklevel=2,
        )
    if len(np.unique(pixdim)) > 1 and not do_resample:
        warnings.warn(
            "The MRI has anisotropic resolution. It is highly recommended that you turn on the 'resampling' "
            "option, as the electrode size will not be exact if the model is built from an MRI with "
            "anisotropic resolution.",
            stacklevel=2,
        )

    # SPM (segmentation, T2 alignment, resampling) needs its runtime: get it
    # in place -- and the license question asked -- before any long work.
    spm_runtime = None
    if use_spm or do_resample:
        from .spm.runtime import ensure_spm

        spm_runtime = ensure_spm(accept_license=accept_matlab_runtime_license, install_missing=install_missing)

    # RAS, then resampling, then zero-padding, then T2 alignment (roast.m).
    ras = convert_to_ras(subj)
    model_mri = ras.path
    if do_resample:
        from .spm.routines import resample_to_1mm

        model_mri = resample_to_1mm(model_mri, runtime=spm_runtime)
    if zero_padding > 0:
        model_mri = zero_pad(model_mri, zero_padding)
    model_mri = os.path.abspath(model_mri)
    if t2 is not None:
        from .spm.routines import realign_t2

        t2 = realign_t2(os.path.abspath(t2), model_mri, runtime=spm_runtime)
    t1_img = nib.load(model_mri)
    voxel_size = np.abs(np.diag(t1_img.affine)[:3])

    recipe_text = ", ".join(f"{n} ({c:g} mA)" for n, c in zip(elec_names, current))
    print(f"ROAST {subj}\nUSING RECIPE: {recipe_text}")

    # Each step is skipped when its output already exists, as in roast.m:
    # steps 1-2 are shared by every simulation of the modeled MRI, steps
    # 3-6 belong to one simulation tag.
    if use_spm:
        from .registration import registration_from_seg8
        from .spm.routines import seg8_path, spm_masks_path, spm_segment

        with_t2 = t2 is not None
        redo = not (spm_masks_path(model_mri, with_t2).exists() and seg8_path(model_mri, with_t2).exists())
        print(f"[1/6] SPM segmentation of {model_mri} ...")
        mask_path = spm_segment(model_mri, t2, runtime=spm_runtime)  # skips what's done itself
        print("[2/6] Mapping the TPM landmarks through SPM's registration (_seg8.mat) ...")
        registration = registration_from_seg8(seg8_path(model_mri, with_t2=t2 is not None))
    else:
        mask_path = _multiaxial_masks_path(model_mri)
        redo = not os.path.exists(mask_path)
        if not redo:
            print(f"[1/6] Multiaxial segmentation done, skipping ({mask_path})")
        else:
            print(f"[1/6] Segmenting {model_mri} with Multiaxial ...")
            from .segmentation.multiaxial import segment  # imports TensorFlow; see note at top

            _write_atomically(mask_path, lambda tmp: segment(model_mri, out_path=tmp, model_dir=model_dir))
        reg_path = registration_path(model_mri, work_dir)
        if reg_path.exists():  # depends only on the MRI, not on the segmentation
            print(f"[2/6] NiftyReg registration done, skipping ({reg_path})")
            registration = Registration.load(reg_path)
        else:
            print("[2/6] Registering to MNI space (NiftyReg) and mapping the TPM landmarks ...")
            registration = run_niftyreg(model_mri, out_dir=work_dir, bin_path=niftyreg_bin)
            redo = True
    if zero_padding > 0:
        tissue_labels = clear_padding(mask_path, zero_padding)  # not in MATLAB; see clear_padding
    else:
        tissue_labels = np.asarray(nib.load(mask_path).dataobj, dtype=np.uint8)

    landmarks = tpm_landmarks_to_subject(registration.tpm2mri)
    mri2mni = registration.mri2mni

    # -- the simulation's tag: a saved simulation with the same options, or a new one
    def plain(value):
        return value.tolist() if isinstance(value, np.ndarray) else value

    options = {
        "subj": subj,
        "model_mri": model_mri,
        "work_dir": work_dir,
        "simulation_tag": simulation_tag,
        "recipe": {name: c for name, c in zip(elec_names, current)},
        "config_txt": recipe_text,
        "cap_type": cap_type,
        "electrodes": [
            {"elec_type": p.elec_type, "elec_size": plain(p.elec_size), "elec_ori": plain(p.elec_ori)}
            for p in elec_paras
        ],
        "conductivities": asdict(conductivities),
        "mesh_options": mesh_options,
        "zero_padding": zero_padding,
        "resampling": do_resample,
        "segmentation": "spm" if use_spm else "multiaxial",
        "t2": t2,
        "is_non_ras": ras.was_reoriented,
        "masks": os.path.abspath(mask_path),
        "landmarks": np.asarray(landmarks, dtype=float).tolist(),
        # 0-based voxel -> MNI mm; MATLAB's opt.mri2mni is the 1-based equivalent
        "mri2mni": mri2mni.tolist(),
        "Affine": registration.affine.tolist(),
    }
    options, reused = register_simulation(options, work_dir, base)
    simulation_tag = options["simulation_tag"]
    if reused:
        print(f"Same options as simulation '{simulation_tag}' run before: carrying on under its tag")
    else:
        print(f"Simulation options saved in {sim_log_path(work_dir, base)}, under tag: {simulation_tag}")
    stem = simulation_stem(base, simulation_tag)
    paths = output_paths(work_dir, base, simulation_tag)

    # `redo`: an earlier step just ran, so whatever later steps saved
    # before was built from something else and is redone too.
    if not redo and os.path.exists(paths["elec_mask"]) and os.path.exists(paths["gel_mask"]):
        print("[3/6] Electrodes already placed, skipping")
        elec_mask = np.asarray(nib.load(paths["elec_mask"]).dataobj, dtype=np.uint8)
        gel_mask = np.asarray(nib.load(paths["gel_mask"]).dataobj, dtype=np.uint8)
    else:
        print("[3/6] Placing electrodes ...")
        cap_names, cap_template = load_cap_info(cap_type)
        custom_names = custom_coords = None
        if any("custom" in name.lower() for name in elec_names):
            custom_names, custom_coords = read_custom_locations(subj)
            # MATLAB's coordinates are 1-based voxels of the input MRI; carry
            # them into the modeled MRI's 0-based voxels (RAS + padding).
            homogeneous = np.column_stack([custom_coords - 1.0, np.ones(len(custom_coords))])
            to_model = np.linalg.inv(t1_img.affine) @ original_img.affine
            custom_coords = (to_model @ homogeneous.T).T[:, :3]
        elec_mask, gel_mask = electrode_placement(
            tissue_labels, landmarks, elec_names, elec_paras, cap_names, cap_template, voxel_size,
            is_biosemi=cap_type == "biosemi", is_egi=cap_type == "egi",
            custom_names=custom_names, custom_coords=custom_coords,
        )
        for key, mask in (("gel_mask", gel_mask), ("elec_mask", elec_mask)):
            image = nib.Nifti1Image(mask.astype(np.uint8), t1_img.affine)
            _write_atomically(paths[key], lambda tmp, image=image: nib.save(image, tmp))
        redo = True

    if not redo and os.path.exists(paths["mesh"]) and os.path.exists(paths["msh"]):
        print("[4/6] Mesh already generated, skipping")
        with np.load(paths["mesh"]) as mesh:
            node, elem = mesh["node"], mesh["elem"]
    else:
        print("[4/6] Meshing ...")
        node, elem, _face = mesh_by_iso2mesh(
            tissue_labels, elec_mask, gel_mask, voxel_size, work_dir=work_dir, out_path=paths["msh"],
            **mesh_options, bin_path=cgalmesh_bin,
        )
        # Written after the .msh, so its presence means both are complete.
        _write_atomically(paths["mesh"], lambda tmp: np.savez_compressed(tmp, node=node, elem=elem))
        redo = True

    if not redo and solution_exists(work_dir, stem):
        print("[5/6] Model already solved, skipping")
    else:
        print("[5/6] Solving FEM ...")
        solve(work_dir, stem, node, elem, elec_names, current, conductivities, bin_path=getdp_bin)
        redo = True

    if not redo and all(os.path.exists(paths[key]) for key in ("v", "e", "emag")):
        print("[6/6] All steps done, loading the results")
        vol_v, vol_e, ef_mag = (
            np.asarray(nib.load(paths[key]).dataobj, dtype=np.float32) for key in ("v", "e", "emag")
        )
    else:
        print("[6/6] Saving the results ...")
        # float32, as saved: a later run that loads them returns the same arrays
        vol_v, vol_e, ef_mag = (
            vol.astype(np.float32) for vol in postprocess(work_dir, stem, node, voxel_size, tissue_labels.shape)
        )
        for key, vol in (("v", vol_v), ("e", vol_e), ("emag", ef_mag)):
            image = nib.Nifti1Image(vol, t1_img.affine)
            _write_atomically(paths[key], lambda tmp, image=image: nib.save(image, tmp))

    print(f"Done. Outputs saved in {work_dir}, under tag: {simulation_tag}")
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
        recipe=dict(zip(elec_names, current)),
        landmarks=np.asarray(landmarks),
        mri2mni=mri2mni,
        mesh_node=node,
        mesh_elem=elem,
        sim_tag=simulation_tag,
        model_mri=model_mri,
    )

    if visualize:
        # The results are already on disk, so a display problem must not
        # cost the user the simulation: report it and return normally.
        try:
            from .viz.review import review_res

            # Draws in a child process: TensorFlow is loaded in this one,
            # and VTK can crash next to it (see review_res).
            review_res(subj, work_dir=work_dir, sim_tag=simulation_tag)
        except Exception as e:  # noqa: BLE001
            warnings.warn(
                f"Visualization failed ({type(e).__name__}: {e}). The simulation "
                f"results are saved in {work_dir}; retry the display with "
                f"roast_py.review_res({subj!r}, sim_tag={simulation_tag!r}).",
                stacklevel=2,
            )
    return result
