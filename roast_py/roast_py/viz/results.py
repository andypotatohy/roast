"""Ports visualizeRes.m (and the result-drawing half of reviewRes.m).

For a finished simulation this draws, like MATLAB:

* voltage and E-field magnitude on a tissue surface in 3D (gray matter by
  default), with each electrode colored by its injected current and two
  color bars -- the field on the right, the injected current on the left;
* voltage and E-field slice views restricted to that tissue, the E-field
  one with arrows for the field direction, both cropped to the brain.

The 3D values come from the node-wise .pos files getDP wrote (as in
MATLAB), and the mesh is taken to world coordinates, so left and right in
the 3D view are the subject's own. The slice views are in voxel space.
"""

from __future__ import annotations

import os

import numpy as np

from ._display import FigureSet, Scene3D
from .sliceshow import SliceViewer
from .views import brain_crop, view_electrodes, view_mri, view_seg

NUM_TISSUE = 6  # hard coded across ROAST

# reviewRes.m's `tissue` option -> (surface label for 3D, labels for slices).
TISSUES = {
    "white": (1, (1,)),
    "gray": (2, (2,)),
    "csf": (3, (3,)),
    "bone": (4, (4,)),
    "skin": (5, (5,)),
    "air": (6, (6,)),
    "brain": (2, (1, 2)),
    "all": (5, (1, 2, 3, 4, 5, 6)),
}


def _tissue_choice(tissue: str):
    try:
        return TISSUES[tissue.lower()]
    except KeyError:
        raise ValueError(
            "Supported tissues to be displayed are: 'white', 'gray', 'CSF', 'bone', "
            "'skin', 'air', 'brain' and 'all'."
        ) from None


def mesh_to_world(node, voxel_size, affine) -> np.ndarray:
    """Mesh node coordinates (mm, as meshing/cgal_mesher.py returns them) ->
    world coordinates. The mm coordinates are 1-based voxel positions
    scaled by the voxel size (see fem/solve.py), so undo the scaling, shift
    to 0-based and apply the image affine -- visualizeRes.m's
    `node/mat(i,i)` followed by `mat*[node 1]'`."""
    voxel = np.asarray(node, dtype=float)[:, :3] / np.asarray(voxel_size, dtype=float) - 1.0
    affine = np.asarray(affine, dtype=float)
    return voxel @ affine[:3, :3].T + affine[:3, 3]


def _tet_grid(points, elem, labels):
    """A PyVista grid of the tetrahedra whose region label is in `labels`."""
    import pyvista as pv

    tets = elem[np.isin(elem[:, 4], labels), :4].astype(np.int64) - 1  # 1-based node ids
    cells = np.hstack([np.full((tets.shape[0], 1), 4), tets]).ravel()
    celltypes = np.full(tets.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, celltypes, points), np.unique(tets)


def field_scene(
    points,
    elem,
    node_values,
    currents,
    surface_label: int,
    title: str,
    bar_title: str,
    upper: str = "max",
    fast_render: bool = True,
) -> Scene3D:
    """The plotmesh() rendering in visualizeRes/reviewRes: `node_values`
    (one per mesh node, NaN where unknown) on the surface of the
    `surface_label` tissue, plus the electrodes colored by their current.

    `upper` is "max" (voltage) or "p95" (E-field: the 95th percentile, so
    a few hot spots at the electrode edges don't wash out the brain).
    """
    currents = np.asarray(currents, dtype=float)
    num_gel = len(currents)

    tissue_grid, tissue_nodes = _tet_grid(points, elem, [surface_label])
    if tissue_grid.n_cells == 0:
        raise ValueError(f"The mesh has no elements of tissue {surface_label} to display.")
    tissue_grid.point_data["value"] = node_values
    surface = tissue_grid.extract_surface(algorithm="dataset_surface")
    if not fast_render:
        # reviewRes.m's sms(): smooth the displayed surface only.
        surface = surface.smooth(n_iter=10, relaxation_factor=0.5)

    shown = node_values[tissue_nodes]
    shown = shown[np.isfinite(shown)]
    lo = float(shown.min())
    hi = float(shown.max() if upper == "max" else np.percentile(shown, 95))

    electrodes = []
    for i, current in enumerate(currents, start=1):
        grid, _ = _tet_grid(points, elem, [NUM_TISSUE + num_gel + i])
        if grid.n_cells:
            elec_surface = grid.extract_surface(algorithm="dataset_surface")
            elec_surface.point_data["current"] = np.full(elec_surface.n_points, current)
            electrodes.append(elec_surface)
    current_range = (float(currents.min()), float(currents.max()))

    def draw(plotter):
        # Two color bars, as in MATLAB: the field, and the injected current
        # the electrodes are colored by. PyVista merges bars that share a
        # title across the whole window, so each panel's current bar gets
        # an (invisible) unique title of its own.
        current_title = "Injected current (mA)" + " " * len(plotter.scalar_bars)
        # Only the end values are labeled: middle labels would sit under
        # the centered title.
        bar_style = dict(vertical=False, position_y=0.03, width=0.42, height=0.1, color="black",
                         title_font_size=14, label_font_size=12, n_labels=2, fmt="%.3g")
        # Split normals at sharp edges: plain smooth shading draws dark
        # cracks along the mesh's non-manifold edges in the sulci.
        shading = dict(smooth_shading=True, split_sharp_edges=True)
        plotter.add_mesh(
            surface,
            scalars="value",
            cmap="jet",
            clim=(lo, hi),
            **shading,
            scalar_bar_args=dict(title=bar_title, position_x=0.54, **bar_style),
        )
        for elec in electrodes:
            plotter.add_mesh(
                elec,
                scalars="current",
                cmap="jet",
                clim=current_range,
                **shading,
                scalar_bar_args=dict(title=current_title, position_x=0.04, **bar_style),
            )

    return Scene3D(title, draw)


def read_node_field(pos_path: str, n_nodes: int, vector: bool) -> np.ndarray:
    """A getDP NodeTable .pos file -> one value per mesh node (NaN if absent):
    the voltage re-referenced to its minimum, or the E-field magnitude."""
    from ..fem.pos_parser import read_pos_node_table

    ids, vals = read_pos_node_table(pos_path, 3 if vector else 1)
    out = np.full(n_nodes, np.nan)
    if vector:
        out[ids - 1] = np.sqrt(np.sum(vals**2, axis=1))
    else:
        out[ids - 1] = vals[:, 0] - vals[:, 0].min()
    return out


def result_views(
    tissue_labels,
    node,
    elem,
    currents,
    affine,
    voxel_size,
    vol_v,
    ef_mag,
    vol_e,
    v_pos: str,
    e_pos: str,
    tag: str = "",
    tissue: str = "brain",
    fast_render: bool = True,
    mri2mni=None,
) -> FigureSet:
    """The result figures shared by visualize_res() and review_res()."""
    surface_label, slice_labels = _tissue_choice(tissue)
    tag_suffix = f": {tag}" if tag else ""
    figs = FigureSet()

    print("generating 3D renderings...")
    points = mesh_to_world(node, voxel_size, affine)
    voltage = read_node_field(v_pos, len(points), vector=False)
    efield = read_node_field(e_pos, len(points), vector=True)
    figs.scenes.append(
        field_scene(points, elem, voltage, currents, surface_label,
                    f"Voltage in Simulation{tag_suffix}", "Voltage (mV)", "max", fast_render)
    )
    figs.scenes.append(
        field_scene(points, elem, efield, currents, surface_label,
                    f"Electric field in Simulation{tag_suffix}", "Electric field (V/m)", "p95", fast_render)
    )

    print("generating slice views...")
    tissue_labels = np.asarray(tissue_labels)
    in_tissue = np.isin(tissue_labels, slice_labels)
    nan_mask = np.where(in_tissue, 1.0, np.nan)
    if tissue.lower() in ("white", "gray", "brain"):
        bbox = brain_crop(tissue_labels)
        pos = None if bbox is None else np.round(bbox.mean(axis=0)).astype(int)
    else:
        bbox, pos = None, None

    viewer = SliceViewer(
        np.asarray(vol_v) * nan_mask, pos=pos, cmap="jet", label="Voltage (mV)",
        fig_name=f"Voltage in Simulation{tag_suffix}. Click anywhere to navigate.",
        mri2mni=mri2mni, bbox=bbox,
    )
    figs.figures.append(viewer.fig)

    ef_masked = np.asarray(ef_mag) * nan_mask
    values = ef_masked[np.isfinite(ef_masked)]
    viewer = SliceViewer(
        ef_masked, pos=pos, cmap="jet", clim=(float(values.min()), float(np.percentile(values, 95))),
        label="Electric field (V/m)",
        fig_name=f"Electric field in Simulation{tag_suffix}. Click anywhere to navigate.",
        vec_img=np.asarray(vol_e) * nan_mask[..., None], mri2mni=mri2mni, bbox=bbox,
    )
    figs.figures.append(viewer.fig)
    return figs


def visualize_res(
    subj: str,
    mask,
    mri2mni,
    node,
    elem,
    currents,
    affine,
    voxel_size,
    vol_v,
    ef_mag,
    vol_e,
    work_dir: str | None = None,
    tag: str = "",
) -> FigureSet:
    """visualizeRes(subj, mask, mri2mni, node, elem, face, inCurrent, imgHdr,
    uniTag, vol_all, ef_mag, ef_all) for roast() results, in roast_py's
    terms: the image header becomes `affine` + `voxel_size`, the unused
    `face` is dropped, and the .pos files are read from `work_dir`
    (default: next to `subj`) as `<subj>_v.pos` / `<subj>_e.pos`.

    Returns the figures; display them with ``.show()``.
    """
    base = os.path.splitext(os.path.basename(subj))[0]
    work_dir = work_dir or os.path.dirname(os.path.abspath(subj))
    return result_views(
        mask, node, elem, currents, affine, voxel_size, vol_v, ef_mag, vol_e,
        os.path.join(work_dir, f"{base}_v.pos"), os.path.join(work_dir, f"{base}_e.pos"),
        tag=tag or base, mri2mni=mri2mni,
    )


def all_views(
    subj: str,
    tissue_labels,
    elec_mask,
    gel_mask,
    landmarks,
    node,
    elem,
    currents,
    affine,
    voxel_size,
    vol_v,
    ef_mag,
    vol_e,
    work_dir: str,
    tissue: str = "brain",
    fast_render: bool = True,
    mri2mni=None,
    t2=None,
) -> FigureSet:
    """Every figure MATLAB's roast() draws for a simulation, in its order:
    viewMRI, viewSeg, viewElectrodes, then visualizeRes."""
    base = os.path.splitext(os.path.basename(subj))[0]
    figs = FigureSet(title=f"ROAST: {base}")
    print("showing MRI...")
    figs.extend(view_mri(subj, t2, mri2mni))
    print("showing segmentations...")
    figs.extend(view_seg(tissue_labels, mri2mni))
    print("showing electrode placement...")
    figs.extend(view_electrodes(tissue_labels, elec_mask, gel_mask, landmarks, affine, tag=base))
    figs.extend(
        result_views(
            tissue_labels, node, elem, currents, affine, voxel_size, vol_v, ef_mag, vol_e,
            os.path.join(work_dir, f"{base}_v.pos"), os.path.join(work_dir, f"{base}_e.pos"),
            tag=base, tissue=tissue, fast_render=fast_render, mri2mni=mri2mni,
        )
    )
    return figs


def figures_dir(work_dir: str, subj: str) -> str:
    base = os.path.splitext(os.path.basename(subj))[0]
    return os.path.join(work_dir, f"{base}_figures")


def show_roast_results(result, block: bool = True) -> FigureSet:
    """What roast() calls at the end: draws everything and shows it (or,
    without a display, saves PNGs to `<work_dir>/<subj>_figures/`)."""
    figs = all_views(
        result.subj, result.tissue_labels, result.elec_mask, result.gel_mask, result.landmarks,
        result.mesh_node, result.mesh_elem, list(result.recipe.values()), result.affine,
        result.voxel_size, result.vol_v, result.ef_mag, result.vol_e, result.work_dir,
        mri2mni=result.mri2mni,
    )
    figs.show(fallback_dir=figures_dir(result.work_dir, result.subj), block=block)
    return figs
