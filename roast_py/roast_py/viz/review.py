"""Ports reviewRes.m: redraw a finished simulation from its saved outputs.

    from roast_py.viz import review_res
    review_res("example/subject1.nii")                  # latest simulation, brain
    review_res("example/subject1.nii", tissue="skin")   # results on another tissue
    review_res("example/subject1.nii", sim_tag="awesomeSimulation", tissue="all")

Differences from MATLAB:

* `sim_tag` is a keyword (MATLAB's 2nd argument). Left out, it is the
  subject's most recent simulation, where MATLAB asks for a tag.
* Only roast() results. roast_target() hasn't been ported, so the
  targeting branch (`tarTag`, the montage topoplot) isn't either.

Third-party imports live inside review_res() -- like roast() -- so this
module imports with nothing installed and review_res() can install the
dependencies it needs.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from ..dependencies import ensure_dependencies


def review_res(
    subj: str,
    tissue: str = "brain",
    fast_render: bool = True,
    work_dir: str | None = None,
    *,
    sim_tag: str | None = None,
    show: bool = True,
    save_dir: str | None = None,
    install_missing: bool = True,
    in_child: bool | None = None,
):
    """reviewRes(subj, simTag, tissue, fastRender) for a simulation roast() ran.

    `subj` is the same MRI path given to roast() (default
    example/MNI152_T1_1mm.nii, as in MATLAB), and `sim_tag` the
    simulation's tag (RoastResult.sim_tag; default: the subject's most
    recent simulation in `work_dir`). `work_dir` is where roast() wrote
    its outputs, if not next to `subj`. `tissue` is one of
    'white', 'gray', 'csf', 'bone', 'skin', 'air', 'brain' (default) or
    'all'. `fast_render=False` smooths the 3D surface for display (slower).

    Shows the figures (blocking until their windows are closed); with
    `save_dir` they're also saved there as PNGs, and with `show=False`
    only saved. Returns the FigureSet.

    If TensorFlow is already loaded in this process -- as it is right
    after roast() -- the figures are drawn in a child Python process
    instead, and this returns None once it finishes. With TensorFlow
    loaded, VTK can crash the whole process (a segmentation fault) while
    setting up OpenGL for the 3D views; seen with Mesa's software
    renderer on Linux. The child never imports TensorFlow. `in_child`
    forces either way (default: decide from whether TensorFlow is loaded).
    """
    from ..roast import _default_subject, _mri_base, output_paths, simulation_stem
    from ..simulation_log import latest_simulation, saved_simulations

    if not subj:
        subj = _default_subject()
    if str(subj).lower() == "nyhead":
        raise NotImplementedError("The New York head is not ported to roast_py yet.")
    subj = os.path.abspath(subj)
    base = _mri_base(subj)
    work_dir = os.path.abspath(work_dir) if work_dir else os.path.dirname(subj)
    if sim_tag is None:
        latest = latest_simulation(work_dir, base, subj)
        if latest is None:
            raise FileNotFoundError(
                f"No simulation of {subj} found in {work_dir} (no {base}_<tag>_roastOptions.json). "
                "Please run roast() first."
            )
        sim_tag = latest["simulation_tag"]
        print(f"No simulation tag given: showing the latest simulation of {base}, '{sim_tag}'")
    stem = simulation_stem(base, sim_tag)
    paths = output_paths(work_dir, base, sim_tag)

    if not os.path.exists(paths["options"]):
        tags = [opt["simulation_tag"] for opt in saved_simulations(work_dir, base, subj)]
        raise FileNotFoundError(
            f"Option file {paths['options']} not found. Simulation {sim_tag} may never have been run. "
            "Please run roast() first."
            + (f" Simulations of {base} saved in {work_dir}: {', '.join(tags)}." if tags else "")
        )

    if in_child is None:
        in_child = "tensorflow" in sys.modules
    if in_child:
        _review_res_in_child(
            subj, tissue=tissue, fast_render=fast_render, work_dir=work_dir,
            sim_tag=sim_tag, show=show, save_dir=save_dir,
        )
        return None

    ensure_dependencies(install_missing=install_missing)

    import nibabel as nib
    import numpy as np

    from .results import all_views, figures_dir

    with open(paths["options"]) as f:
        options = json.load(f)
    print(f"Showing results for Simulation {sim_tag} ...")

    def need(key: str, what: str, step: str) -> str:
        path = options["masks"] if key == "masks" else paths[key]
        if not os.path.exists(path):
            raise FileNotFoundError(f"{what} {path} not found. Check if you ran through {step}.")
        return path

    # The MRI that was modeled: after RAS conversion and zero-padding.
    model_mri = options.get("model_mri", subj)
    if not os.path.exists(model_mri):
        raise FileNotFoundError(
            f"The subject MRI {model_mri} does not exist. Check if you ran through RAS conversion or "
            "zero-padding if you tried to do that."
        )
    t1 = nib.load(model_mri)
    tissue_labels = np.asarray(nib.load(need("masks", "Segmentation masks", "MRI segmentation")).dataobj)
    elec = np.asarray(nib.load(need("elec_mask", "Electrode mask", "electrode placement")).dataobj)
    gel = np.asarray(nib.load(need("gel_mask", "Gel mask", "electrode placement")).dataobj)
    mesh = np.load(need("mesh", "Mesh file", "meshing"))
    need("v_pos", "Solution file", "solving")
    need("e_pos", "Solution file", "solving")
    vol_v = np.asarray(nib.load(need("v", "Result file", "post processing after solving")).dataobj)
    vol_e = np.asarray(nib.load(need("e", "Result file", "post processing after solving")).dataobj)
    ef_mag = np.asarray(nib.load(need("emag", "Result file", "post processing after solving")).dataobj)

    mri2mni = options.get("mri2mni")
    figs = all_views(
        model_mri,
        tissue_labels,
        elec,
        gel,
        np.asarray(options["landmarks"]),
        mesh["node"],
        mesh["elem"],
        list(options["recipe"].values()),
        t1.affine,
        np.abs(np.diag(t1.affine)[:3]),
        vol_v,
        ef_mag,
        vol_e,
        work_dir,
        tissue=tissue,
        fast_render=fast_render,
        mri2mni=None if mri2mni is None else np.asarray(mri2mni),
        t2=_existing_t2(options.get("t2")),
        stem=stem,
        sim_tag=sim_tag,
    )
    if save_dir:
        saved = figs.save(save_dir)
        print(f"Saved {len(saved)} figure(s) in {save_dir}")
    if show:
        figs.show(fallback_dir=None if save_dir else figures_dir(work_dir, stem))
    return figs


def _existing_t2(path):
    """The (aligned) T2 the simulation used, for viewMRI; None if none or gone."""
    if path and not os.path.exists(path):
        raise FileNotFoundError(
            f"T2 file {path} does not exist. You used that to run ROAST but maybe later deleted it."
        )
    return path or None


def _review_res_in_child(subj, *, tissue, fast_render, work_dir, sim_tag, show, save_dir) -> None:
    """Runs review_res() in a fresh interpreter (see review_res's docstring)."""
    cmd = [sys.executable, "-m", "roast_py.viz.review", str(subj) if subj else "", "--tissue", tissue]
    if not fast_render:
        cmd.append("--slow-render")
    if work_dir:
        cmd += ["--work-dir", str(work_dir)]
    if sim_tag:
        cmd += ["--sim-tag", sim_tag]
    if save_dir:
        cmd += ["--save-dir", str(save_dir)]
    if not show:
        cmd.append("--no-show")
    # Make roast_py importable in the child even when it's used from a git
    # clone through a sys.path entry rather than installed.
    env = dict(os.environ)
    package_parent = str(Path(__file__).resolve().parents[2])
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [package_parent, env.get("PYTHONPATH")]))
    result = subprocess.run(cmd, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Drawing the figures failed (child process exit status {result.returncode}).")


def _main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m roast_py.viz.review", description="Redraw a finished roast() simulation (reviewRes)."
    )
    parser.add_argument("subj", nargs="?", default="", help="the MRI given to roast() (default: MNI152)")
    parser.add_argument("--sim-tag", default=None)
    parser.add_argument("--tissue", default="brain")
    parser.add_argument("--slow-render", action="store_true", help="smooth the 3D surface (fastRender off)")
    parser.add_argument("--work-dir", default=None)
    parser.add_argument("--save-dir", default=None, help="also save the figures here as PNGs")
    parser.add_argument("--no-show", action="store_true", help="only save, don't open windows")
    args = parser.parse_args(argv)
    review_res(
        args.subj or None, args.tissue, not args.slow_render, args.work_dir,
        sim_tag=args.sim_tag, show=not args.no_show, save_dir=args.save_dir, install_missing=False,
        in_child=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
