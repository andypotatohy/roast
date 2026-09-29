"""Tests for roast_py.viz (the reviewRes/visualizeRes/sliceshow port).

Everything runs on matplotlib's Agg backend with synthetic data. The 3D
scenes are exercised through a fake plotter that records what would be
drawn, since test machines usually have no OpenGL; actual rendering and
the interactive windows were verified under Xvfb (see README).
"""

import json

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.backend_bases import MouseEvent  # noqa: E402

from roast_py.viz import _display  # noqa: E402
from roast_py.viz.results import (  # noqa: E402
    NUM_TISSUE,
    field_scene,
    mesh_to_world,
    read_node_field,
    result_views,
)
from roast_py.viz.sliceshow import SliceViewer, sliceshow  # noqa: E402
from roast_py.viz.views import brain_crop, view_electrodes, view_seg, voxel_to_world  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


class FakePlotter:
    """Records add_mesh/add_point_labels calls instead of rendering."""

    def __init__(self):
        self.meshes = []
        self.labels = []
        self.scalar_bars = {}

    def add_mesh(self, mesh, **kwargs):
        self.meshes.append((mesh, kwargs))
        title = kwargs.get("scalar_bar_args", {}).get("title")
        if title is not None:
            self.scalar_bars[title] = kwargs.get("clim")

    def add_point_labels(self, points, labels, **kwargs):
        self.labels.append((np.asarray(points), list(labels)))


# --------------------------------------------------------------------------
# brainCrop
# --------------------------------------------------------------------------


def test_brain_crop_is_the_white_matter_bounding_box():
    mask = np.zeros((30, 40, 50), dtype=np.uint8)
    mask[5:20, 10:30, 10:40] = 1  # 15 x 20 = 300 voxels per axial slice
    bbox = brain_crop(mask)
    assert bbox.tolist() == [[5, 10, 10], [19, 29, 39]]


def test_brain_crop_ignores_thin_axial_slices_like_matlab():
    mask = np.zeros((30, 40, 50), dtype=np.uint8)
    mask[5:20, 10:30, 10:40] = 1
    mask[10, 15, 2:10] = 1  # a brainstem-like strand: 1 voxel per axial slice
    bbox = brain_crop(mask)
    assert bbox[0, 2] == 10  # S/I extent needs > 100 voxels per slice
    assert bbox[0, 0] == 5  # R/L and A/P have no threshold


def test_brain_crop_without_white_matter_is_none():
    assert brain_crop(np.zeros((5, 5, 5))) is None


# --------------------------------------------------------------------------
# sliceshow
# --------------------------------------------------------------------------


def _volume(shape=(20, 24, 28)):
    i, j, k = np.indices(shape)
    return (i + 100 * j + 10000 * k).astype(float)  # every voxel distinct


def _click(viewer, panel, x, y):
    ax = viewer.axes[panel]
    xd, yd = ax.transData.transform((x, y))
    event = MouseEvent("button_press_event", viewer.fig.canvas, xd, yd, button=1)
    viewer.fig.canvas.callbacks.process("button_press_event", event)


def test_sliceshow_starts_at_the_center_and_reports_the_value():
    img = _volume()
    viewer = sliceshow(img)
    assert viewer.voxel.tolist() == [9, 11, 13]  # (shape - 1) // 2
    assert viewer.value == img[9, 11, 13]
    assert viewer.clim == (img.min(), img.max())
    # Coronal / sagittal / axial panels show the right planes.
    assert np.array_equal(viewer._images[0].get_array(), img[:, 11, :].T)
    assert np.array_equal(viewer._images[1].get_array(), img[9, :, :].T)
    assert np.array_equal(viewer._images[2].get_array(), img[:, :, 13].T)


@pytest.mark.parametrize(
    "panel, xy, expected",
    [
        (0, (3, 7), [3, 11, 7]),  # coronal: x, z
        (1, (4, 5), [9, 4, 5]),  # sagittal: y, z
        (2, (2, 6), [2, 6, 13]),  # axial: x, y
    ],
)
def test_clicking_a_panel_moves_to_that_voxel(panel, xy, expected):
    img = _volume()
    viewer = SliceViewer(img)
    _click(viewer, panel, *xy)
    assert viewer.voxel.tolist() == expected
    assert viewer.value == img[tuple(expected)]
    assert viewer._value_text.get_text() == f"{img[tuple(expected)]:.2f}"
    assert [viewer._boxes[("Voxel", j)].text for j in range(3)] == [str(v) for v in expected]


def test_clicking_outside_the_volume_does_nothing():
    viewer = SliceViewer(_volume())
    _click(viewer, 2, 25, 26)  # beyond x=19 (axes span the largest dimension)
    assert viewer.voxel.tolist() == [9, 11, 13]


def test_typing_voxel_coordinates_navigates():
    viewer = SliceViewer(_volume())
    viewer._boxes[("Voxel", 0)].set_val("4")  # set_val submits, like pressing enter
    viewer._on_submit("Voxel", 0, "4")
    assert viewer.voxel.tolist() == [4, 11, 13]
    viewer._on_submit("Voxel", 2, "not a number")
    assert viewer.voxel.tolist() == [4, 11, 13]
    viewer._on_submit("Voxel", 1, "999")  # outside: stays put
    assert viewer.voxel.tolist() == [4, 11, 13]


def test_mni_coordinates_follow_the_mapping_and_can_be_typed():
    mri2mni = np.array([[1, 0, 0, -10], [0, 1, 0, -20], [0, 0, 1, -30], [0, 0, 0, 1]], float)
    viewer = SliceViewer(_volume(), mri2mni=mri2mni)
    assert viewer.mni.tolist() == [9 - 10, 11 - 20, 13 - 30]
    viewer._on_submit("MNI", 0, "-5")
    assert viewer.voxel.tolist() == [5, 11, 13]
    assert viewer._boxes[("MNI", 0)].text == "-5"


def test_bounding_box_crops_but_coordinates_stay_in_the_full_volume():
    img = _volume()
    bbox = np.array([[2, 3, 4], [12, 15, 20]])
    viewer = SliceViewer(img, pos=(5, 6, 7), bbox=bbox)
    assert viewer.img.shape == (11, 13, 17)
    assert viewer.voxel.tolist() == [5, 6, 7]
    assert viewer.value == img[5, 6, 7]
    _click(viewer, 2, 0, 0)  # cropped origin
    assert viewer.voxel.tolist() == [2, 3, 7]
    with pytest.raises(ValueError, match="outside of the bounding box"):
        SliceViewer(img, pos=(1, 6, 7), bbox=bbox)


def test_nan_voxels_are_white_and_shown_as_nan():
    img = _volume()
    img[:, :, :10] = np.nan
    viewer = SliceViewer(img, pos=(3, 3, 3))
    assert viewer._value_text.get_text() == "nan"
    assert viewer.cmap.get_bad().tolist() == [1.0, 1.0, 1.0, 1.0]


def test_vector_field_arrows_are_drawn_every_5_voxels():
    img = _volume()
    vec = np.zeros((*img.shape, 3))
    vec[..., 0] = 1.0
    viewer = SliceViewer(img, vec_img=vec)
    quiver = viewer._quivers[2]  # axial: 20 x 24 -> 4 x 5 arrows
    assert quiver.N == 4 * 5
    assert np.allclose(quiver.U, 1.0) and np.allclose(quiver.V, 0.0)


def test_sliceshow_validates_its_inputs():
    img = _volume()
    with pytest.raises(ValueError, match="meaningful values"):
        SliceViewer(np.full((4, 4, 4), np.nan))
    with pytest.raises(ValueError, match="Vector field"):
        SliceViewer(img, vec_img=np.zeros((*img.shape, 2)))
    with pytest.raises(ValueError, match="voxel-to-MNI"):
        SliceViewer(img, mri2mni=np.eye(3))
    with pytest.raises(ValueError, match="3D"):
        SliceViewer(np.zeros((4, 4)))


def test_segmentation_view_uses_one_color_per_tissue():
    mask = np.zeros((10, 10, 10), dtype=np.uint8)
    for label in range(7):
        mask[label, :, :] = label
    viewer = view_seg(mask).figures[0]._roast_viewer
    colors = [viewer.cmap(viewer._images[0].norm(label)) for label in range(7)]
    assert len({tuple(np.round(c, 3)) for c in colors}) == 7
    assert np.allclose(colors[0][:3], (0, 0, 0)) and np.allclose(colors[1][:3], (1, 1, 1))


# --------------------------------------------------------------------------
# coordinates and .pos reading
# --------------------------------------------------------------------------


def test_mesh_nodes_map_to_the_same_world_point_as_their_voxel():
    affine = np.array([[-1.0, 0, 0, 90], [0, 1.2, 0, -126], [0, 0, 1.2, -72], [0, 0, 0, 1]])
    voxel_size = np.abs(np.diag(affine)[:3])
    voxels = np.array([[0, 0, 0], [10, 20, 30]], float)
    # Mesh coordinates: 1-based voxel positions scaled by voxel size (fem/solve.py).
    node = (voxels + 1) * voxel_size
    assert np.allclose(mesh_to_world(node, voxel_size, affine), voxel_to_world(voxels, affine))


def test_read_node_field_rereferences_voltage_and_takes_field_magnitude(tmp_path):
    v = tmp_path / "v.pos"
    v.write_text('View "v" {\n1 5.0\n3 2.0\n};\n')
    e = tmp_path / "e.pos"
    e.write_text('View "e" {\n2 3.0 4.0 0.0\n};\n')
    volts = read_node_field(str(v), 4, vector=False)
    assert np.allclose(volts[[0, 2]], [3.0, 0.0]) and np.isnan(volts[[1, 3]]).all()
    field = read_node_field(str(e), 3, vector=True)
    assert field[1] == 5.0 and np.isnan(field[[0, 2]]).all()


# --------------------------------------------------------------------------
# 3D scenes (drawn into a fake plotter)
# --------------------------------------------------------------------------


def _tiny_mesh(n_elec=2):
    """Tetrahedra on a small grid: gray matter below z=2, electrodes above."""
    from scipy.spatial import Delaunay

    grid = np.stack(np.meshgrid(*[np.arange(4.0)] * 3, indexing="ij"), -1).reshape(-1, 3)
    tets = Delaunay(grid).simplices
    centroid = grid[tets].mean(axis=1)
    labels = np.full(len(tets), 1)  # white matter by default
    labels[centroid[:, 2] < 2] = 2  # gray
    labels[(centroid[:, 2] < 2) & (centroid[:, 1] < 1)] = 5  # skin
    top = centroid[:, 2] >= 2
    labels[top & (centroid[:, 0] < 1.5)] = NUM_TISSUE + n_elec + 1
    labels[top & (centroid[:, 0] > 1.5)] = NUM_TISSUE + n_elec + 2
    elem = np.column_stack([tets + 1, labels])  # 1-based node ids
    return grid + 1.0, elem  # 1-based voxel positions, 1 mm voxels


def test_field_scene_colors_the_tissue_and_electrodes_like_visualizeres():
    node, elem = _tiny_mesh()
    values = node[:, 2] * 10.0  # voltage rising with z
    scene = field_scene(node, elem, values, [1.0, -1.0], 2, "Voltage", "Voltage (mV)", "max")
    plotter = FakePlotter()
    scene.draw(plotter)

    surface, surface_kw = plotter.meshes[0]
    gray_nodes = np.unique(elem[elem[:, 4] == 2, :4]) - 1
    assert surface_kw["clim"] == (values[gray_nodes].min(), values[gray_nodes].max())
    assert surface.n_cells > 0

    electrodes = plotter.meshes[1:]
    assert len(electrodes) == 2
    assert [set(m.point_data["current"]) for m, _ in electrodes] == [{1.0}, {-1.0}]
    assert all(kw["clim"] == (-1.0, 1.0) for _, kw in electrodes)
    assert set(plotter.scalar_bars) == {"Voltage (mV)", "Injected current (mA)"}


def test_each_panels_current_bar_gets_its_own_title():
    """PyVista merges same-titled color bars across a window's panels."""
    node, elem = _tiny_mesh()
    plotter = FakePlotter()
    for title in ("Voltage (mV)", "Electric field (V/m)"):
        field_scene(node, elem, node[:, 2], [1.0, -1.0], 2, "x", title).draw(plotter)
    current_bars = [t for t in plotter.scalar_bars if t.startswith("Injected current")]
    assert len(current_bars) == 2


def test_field_scene_caps_the_efield_at_the_95th_percentile():
    node, elem = _tiny_mesh()
    rng = np.random.default_rng(0)
    values = rng.random(len(node))
    scene = field_scene(node, elem, values, [1.0, -1.0], 2, "E", "Electric field (V/m)", "p95")
    plotter = FakePlotter()
    scene.draw(plotter)
    gray_nodes = np.unique(elem[elem[:, 4] == 2, :4]) - 1
    assert plotter.meshes[0][1]["clim"][1] == pytest.approx(np.percentile(values[gray_nodes], 95))


def test_view_electrodes_draws_each_structure_and_the_four_landmarks():
    shape = (30, 30, 30)
    i, j, k = np.indices(shape)
    r = np.sqrt((i - 15) ** 2 + (j - 15) ** 2 + (k - 15) ** 2)
    mask = np.zeros(shape, dtype=np.uint8)
    mask[r < 12] = 5
    mask[r < 7] = 2
    elec = ((r >= 12) & (r < 15) & (k > 24)).astype(np.uint8)
    gel = ((r >= 12) & (r < 15) & (k > 20) & (k <= 24)).astype(np.uint8)
    landmarks = np.array([[15, 27, 15], [15, 3, 15], [27, 15, 15], [3, 15, 15], [15, 20, 3], [15, 10, 3]])
    affine = np.diag([2.0, 2.0, 2.0, 1.0])

    scene = view_electrodes(mask, elec, gel, landmarks, affine, tag="demo").scenes[0]
    assert scene.title == "Electrode placement in Simulation: demo"
    plotter = FakePlotter()
    scene.draw(plotter)
    assert [kw["color"] for _, kw in plotter.meshes] == [
        (229 / 255, 181 / 255, 161 / 255), (1.0, 0.6, 0.8), "blue", "green"
    ]
    assert [kw["opacity"] for _, kw in plotter.meshes] == [0.2, 1.0, 0.8, 0.8]
    points, names = plotter.labels[0]
    assert names == ["Nasion", "Inion", "Right Ear", "Left Ear"]
    assert np.allclose(points, landmarks[:4] * 2.0)  # world coordinates
    # Isosurfaces land in world space too: the scalp sphere is ~24 mm across.
    scalp = plotter.meshes[0][0]
    assert 40 < scalp.bounds[1] - scalp.bounds[0] < 56


# --------------------------------------------------------------------------
# result views and showing/saving
# --------------------------------------------------------------------------


def _result_inputs(tmp_path):
    node, elem = _tiny_mesh()
    shape = (6, 6, 6)
    labels = np.full(shape, 2, dtype=np.uint8)
    labels[1:4, 1:4, 1:4] = 1  # some white matter
    vol_v = np.random.default_rng(1).random(shape)
    vol_e = np.random.default_rng(2).random((*shape, 3))
    ef_mag = np.linalg.norm(vol_e, axis=-1)
    (tmp_path / "v.pos").write_text("x\n" + "".join(f"{n + 1} {n * 1.0}\n" for n in range(len(node))) + "};\n")
    (tmp_path / "e.pos").write_text("x\n" + "".join(f"{n + 1} 1 0 0\n" for n in range(len(node))) + "};\n")
    return labels, node, elem, vol_v, ef_mag, vol_e


def test_result_views_builds_two_3d_scenes_and_two_slice_viewers(tmp_path):
    labels, node, elem, vol_v, ef_mag, vol_e = _result_inputs(tmp_path)
    figs = result_views(
        labels, node, elem, [1.0, -1.0], np.eye(4), np.ones(3), vol_v, ef_mag, vol_e,
        str(tmp_path / "v.pos"), str(tmp_path / "e.pos"), tag="t", tissue="all",
    )
    assert [s.title for s in figs.scenes] == ["Voltage in Simulation: t", "Electric field in Simulation: t"]
    voltage, efield = (f._roast_viewer for f in figs.figures)
    assert np.allclose(voltage.img, vol_v)  # tissue='all': nothing masked out
    assert efield.clim == pytest.approx((ef_mag.min(), np.percentile(ef_mag, 95)))
    assert efield.vec is not None


def test_result_views_masks_slices_to_the_chosen_tissue(tmp_path):
    labels, node, elem, vol_v, ef_mag, vol_e = _result_inputs(tmp_path)
    labels[0] = 5  # some skin
    figs = result_views(
        labels, node, elem, [1.0, -1.0], np.eye(4), np.ones(3), vol_v, ef_mag, vol_e,
        str(tmp_path / "v.pos"), str(tmp_path / "e.pos"), tissue="skin",
    )
    voltage = figs.figures[0]._roast_viewer
    assert np.isfinite(voltage.img[0]).all() and np.isnan(voltage.img[1:]).all()
    with pytest.raises(ValueError, match="Supported tissues"):
        result_views(
            labels, node, elem, [1.0, -1.0], np.eye(4), np.ones(3), vol_v, ef_mag, vol_e,
            str(tmp_path / "v.pos"), str(tmp_path / "e.pos"), tissue="liver",
        )


def test_without_a_display_figures_are_saved_instead(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(_display, "display_available", lambda: False)
    monkeypatch.setattr(_display, "can_render_3d", lambda: False)
    figs = _display.FigureSet()
    figs.figures.append(SliceViewer(_volume(), fig_name="MRI: Click anywhere to navigate.").fig)
    figs.scenes.append(_display.Scene3D("x", lambda p: None))

    saved = figs.show(fallback_dir=tmp_path / "figs")

    assert [p.name for p in saved] == ["01_mri.png"]
    assert saved[0].stat().st_size > 0
    out = capsys.readouterr().out
    assert "No display available" in out and "3D views skipped" in out


def test_review_res_explains_a_missing_simulation(tmp_path):
    from roast_py.viz.review import review_res

    with pytest.raises(FileNotFoundError, match="Please run roast"):
        review_res(str(tmp_path / "subject1.nii"), show=False, install_missing=False)


def test_review_res_reads_back_what_roast_saved(tmp_path, monkeypatch):
    """A miniature roast() output directory, then review_res() on it."""
    import nibabel as nib

    from roast_py.roast import output_paths
    from roast_py.viz import review as review_module

    labels, node, elem, vol_v, ef_mag, vol_e = _result_inputs(tmp_path)
    affine = np.eye(4)
    subj = tmp_path / "subj.nii"
    nib.save(nib.Nifti1Image(np.ones(labels.shape, np.float32), affine), subj)
    masks = tmp_path / "subj_multiaxial_masks.nii"
    nib.save(nib.Nifti1Image(labels, affine), masks)
    paths = output_paths(str(tmp_path), "subj")
    for key, data in (("elec_mask", np.zeros(labels.shape, np.uint8)), ("gel_mask", np.zeros(labels.shape, np.uint8)),
                      ("v", vol_v), ("e", vol_e), ("emag", ef_mag)):
        nib.save(nib.Nifti1Image(np.asarray(data, np.float32 if key in "v e emag" else np.uint8), affine), paths[key])
    np.savez(paths["mesh"], node=node, elem=elem)
    (tmp_path / "v.pos").rename(paths["v_pos"])
    (tmp_path / "e.pos").rename(paths["e_pos"])
    with open(paths["options"], "w") as f:
        json.dump({"recipe": {"Fp1": 1.0, "P4": -1.0}, "masks": str(masks), "landmarks": [[1, 1, 1]] * 6,
                   "mri2mni": None}, f)

    figs = review_module.review_res(str(subj), show=False, install_missing=False)
    assert [s.title for s in figs.scenes] == [
        "Electrode placement in Simulation: subj", "Voltage in Simulation: subj", "Electric field in Simulation: subj",
    ]
    assert len(figs.figures) == 4  # MRI, segmentation, voltage, E-field

    saved = figs.save(tmp_path / "out", include_3d=False)
    assert [p.name for p in saved] == [
        "01_mri.png", "02_segmentation.png", "03_voltage_in_simulation_subj.png",
        "04_electric_field_in_simulation_subj.png",
    ]
