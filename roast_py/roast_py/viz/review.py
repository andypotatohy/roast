"""Ports reviewRes.m: redraw a finished simulation from its saved outputs.

    from roast_py.viz import review_res
    review_res("example/subject1.nii")                  # brain, as roast() shows it
    review_res("example/subject1.nii", tissue="skin")   # results on another tissue

Differences from MATLAB:

* No simulation tag: roast_py doesn't tag simulations yet, so a subject's
  outputs are found by its name alone (see roast.output_paths).
* Only roast() results. roast_target() hasn't been ported, so the
  targeting branch (`tarTag`, the montage topoplot) isn't either.

Third-party imports live inside review_res() -- like roast() -- so this
module imports with nothing installed and review_res() can install the
dependencies it needs.
"""

from __future__ import annotations

import json
import os

from ..dependencies import ensure_dependencies


def review_res(
    subj: str,
    tissue: str = "brain",
    fast_render: bool = True,
    work_dir: str | None = None,
    *,
    show: bool = True,
    save_dir: str | None = None,
    install_missing: bool = True,
):
    """reviewRes(subj, [], tissue, fastRender) for a simulation roast() ran.

    `subj` is the same MRI path given to roast(); `work_dir` is where
    roast() wrote its outputs, if not next to `subj`. `tissue` is one of
    'white', 'gray', 'csf', 'bone', 'skin', 'air', 'brain' (default) or
    'all'. `fast_render=False` smooths the 3D surface for display (slower).

    Shows the figures (blocking until their windows are closed); with
    `save_dir` they're also saved there as PNGs, and with `show=False`
    only saved. Returns the FigureSet.
    """
    ensure_dependencies(install_missing=install_missing)

    import nibabel as nib
    import numpy as np

    from ..roast import output_paths
    from .results import all_views, figures_dir

    subj = os.path.abspath(subj)
    base = os.path.splitext(os.path.basename(subj))[0]
    work_dir = os.path.abspath(work_dir) if work_dir else os.path.dirname(subj)
    paths = output_paths(work_dir, base)

    if not os.path.exists(paths["options"]):
        raise FileNotFoundError(
            f"Option file {paths['options']} not found. The simulation for {subj} may never "
            "have been run (or ran with an older roast_py). Please run roast() first."
        )
    with open(paths["options"]) as f:
        options = json.load(f)
    print(f"Showing results for {base} ...")

    def need(key: str, what: str, step: str) -> str:
        path = options["masks"] if key == "masks" else paths[key]
        if not os.path.exists(path):
            raise FileNotFoundError(f"{what} {path} not found. Check if you ran through {step}.")
        return path

    if not os.path.exists(subj):
        raise FileNotFoundError(f"The subject MRI you provided {subj} does not exist.")
    t1 = nib.load(subj)
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
        subj,
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
    )
    if save_dir:
        saved = figs.save(save_dir)
        print(f"Saved {len(saved)} figure(s) in {save_dir}")
    if show:
        figs.show(fallback_dir=None if save_dir else figures_dir(work_dir, subj))
    return figs
