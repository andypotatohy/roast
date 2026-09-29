"""Ports viewMRI.m, viewSeg.m, viewElectrodes.m and brainCrop.m.

The first two are slice viewers (sliceshow); viewElectrodes is a 3D
rendering of the scalp, gray matter, electrodes, gel and landmarks, drawn
here with PyVista in world (scanner RAS, mm) coordinates -- the same space
as the voltage/E-field renderings, so all 3D views share one camera.
"""

from __future__ import annotations

import numpy as np

from ._display import FigureSet, Scene3D
from .sliceshow import SliceViewer

TISSUE_NAMES = ["background", "white", "gray", "CSF", "bone", "skin", "air"]

# viewSeg.m's "more anatomical looking colormap" (Andrew Birnbaum).
SEG_COLORS = [
    (0, 0, 0),  # background: black
    (1, 1, 1),  # white matter: white
    (0.7, 0.7, 0.7),  # gray matter: gray
    (105 / 255, 175 / 255, 255 / 255),  # CSF: blue
    (241 / 255, 214 / 255, 145 / 255),  # bone
    (177 / 255, 122 / 255, 101 / 255),  # skin
    (0.6863, 0.8824, 0.6863),  # air cavities
]

LANDMARK_NAMES = ["Nasion", "Inion", "Right Ear", "Left Ear"]


def brain_crop(mask) -> np.ndarray | None:
    """Ports brainCrop.m: the white matter's bounding box, rows (min, max),
    columns R/L, A/P, S/I -- as 0-based inclusive voxel indices.

    Like the original, the S/I extent only counts axial slices with more
    than 100 white-matter voxels, which trims stray voxels in the brainstem
    and cerebellum. Returns None if there is no white matter.
    """
    wm = np.asarray(mask) == 1
    thresholds = (0, 0, 100)
    bbox = np.zeros((2, 3), dtype=int)
    for axis in range(3):
        other = tuple(a for a in range(3) if a != axis)
        idx = np.nonzero(wm.sum(axis=other) > thresholds[axis])[0]
        if idx.size == 0:
            return None
        bbox[:, axis] = idx.min(), idx.max()
    return bbox


def _load(img):
    """A path (via nibabel) or an array -> ndarray."""
    if isinstance(img, (str, bytes)) or hasattr(img, "__fspath__"):
        import nibabel as nib

        return np.asarray(nib.load(img).dataobj)
    return np.asarray(img)


def view_mri(t1, t2=None, mri2mni=None) -> FigureSet:
    """Ports viewMRI.m: grayscale slice viewers of the T1 (and T2, if given)."""
    figs = FigureSet()
    viewer = SliceViewer(_load(t1), cmap="gray", fig_name="MRI: Click anywhere to navigate.", mri2mni=mri2mni)
    figs.figures.append(viewer.fig)
    if t2 is not None:
        viewer = SliceViewer(
            _load(t2), cmap="gray", fig_name="MRI: T2. Click anywhere to navigate.", mri2mni=mri2mni
        )
        figs.figures.append(viewer.fig)
    return figs


def view_seg(mask, mri2mni=None) -> FigureSet:
    """Ports viewSeg.m: the 6-tissue segmentation with its anatomical colormap."""
    from matplotlib.colors import ListedColormap

    viewer = SliceViewer(
        _load(mask),
        cmap=ListedColormap(SEG_COLORS, name="roast_seg"),
        clim=(-0.5, len(SEG_COLORS) - 0.5),  # one color per label
        label="Tissue index",
        fig_name="Segmentation. Click anywhere to navigate.",
        mri2mni=mri2mni,
        ticks=range(len(SEG_COLORS)),
        ticklabels=TISSUE_NAMES,
    )
    return FigureSet(figures=[viewer.fig])


def voxel_to_world(points, affine) -> np.ndarray:
    """0-based voxel coordinates (n, 3) -> world coordinates via a nibabel affine."""
    points = np.asarray(points, dtype=float)
    return points @ np.asarray(affine)[:3, :3].T + np.asarray(affine)[:3, 3]


def isosurface(volume, affine, level=0.5, sigma=1.0):
    """imgaussfilt3 + isosurface: a smoothed binary volume's level surface,
    as a PyVista mesh in world coordinates (None if the volume is empty)."""
    import pyvista as pv
    from scipy.ndimage import gaussian_filter
    from skimage.measure import marching_cubes

    volume = np.asarray(volume, dtype=np.float32)
    if not volume.any():
        return None
    # imgaussfilt3's default kernel is 2*ceil(2*sigma)+1 wide: truncate at 2 sigma.
    smooth = gaussian_filter(volume, sigma=sigma, truncate=2.0)
    if smooth.max() <= level:
        return None
    # Pad so surfaces touching the volume edge still close.
    verts, faces, _normals, _values = marching_cubes(np.pad(smooth, 1), level=level)
    verts = voxel_to_world(verts - 1, affine)
    return pv.PolyData(verts, np.hstack([np.full((faces.shape[0], 1), 3), faces]).ravel())


def view_electrodes(mask, elec, gel, landmarks, affine, tag="") -> FigureSet:
    """Ports viewElectrodes.m: scalp (translucent), gray matter (pink),
    electrodes (blue), gel (green) and the four head landmarks (red).

    `landmarks` are 0-based voxel coordinates in geometry/landmarks.py's
    order (nasion, inion, right, left, ...); only the first four are drawn,
    as in MATLAB. Returns a FigureSet holding one 3D scene.
    """
    mask = _load(mask)
    surfaces = [
        (isosurface(mask == 5, affine), dict(color=(229 / 255, 181 / 255, 161 / 255), opacity=0.2)),
        (isosurface(mask == 2, affine), dict(color=(1.0, 0.6, 0.8), opacity=1.0)),
        (isosurface(_load(elec) > 0, affine), dict(color="blue", opacity=0.8)),
        (isosurface(_load(gel) > 0, affine), dict(color="green", opacity=0.8)),
    ]
    points = None
    if landmarks is not None and len(landmarks) >= 4:
        points = voxel_to_world(np.asarray(landmarks)[:4], affine)

    def draw(plotter):
        for surface, style in surfaces:
            if surface is not None:
                plotter.add_mesh(surface, smooth_shading=True, **style)
        if points is not None:
            plotter.add_point_labels(
                points,
                LANDMARK_NAMES,
                point_color="red",
                point_size=14,
                render_points_as_spheres=True,
                text_color="red",
                font_size=14,
                bold=True,
                shape=None,
                always_visible=True,
            )

    title = f"Electrode placement in Simulation: {tag}" if tag else "Electrode placement"
    return FigureSet(scenes=[Scene3D(title, draw)])
