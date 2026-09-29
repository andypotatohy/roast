"""Ports sliceshow.m: an interactive three-plane slice viewer.

Layout matches the MATLAB version: coronal (top left), sagittal (top
right) and axial (bottom left) slices through the current voxel, with a
crosshair and marker at it, plus a panel (bottom right) showing the value
there and editable voxel/MNI coordinates. Click in any slice to move
there, or type coordinates into the boxes. An optional vector field (the
E-field) is overlaid as arrows every 5 voxels, and a bounding box (see
views.brain_crop) restricts the display to the brain.

Differences from MATLAB, all deliberate:

* Voxel coordinates are **0-based**, like everything else in roast_py (and
  nibabel): MATLAB's voxel (129, 129, 129) is (128, 128, 128) here.
* `mri2mni` likewise maps 0-based voxel coordinates to MNI -- i.e. it
  composes like a nibabel affine, not like SPM's 1-based `.mat`.
* NaN voxels (outside the displayed tissue) are drawn white via the
  colormap's "bad" color, rather than by prepending white to the
  colormap, so the lowest data value keeps its own color.
"""

from __future__ import annotations

import numpy as np

from ._display import new_figure

_QUIVER_STEP = 5  # MATLAB: 1:5:end
_QUIVER_SCALE = 2  # MATLAB: quiver3(..., 2, ...)


class SliceViewer:
    """One sliceshow window. Keep a reference (FigureSet does) or its
    widgets stop responding once garbage-collected."""

    def __init__(
        self,
        img,
        pos=None,
        cmap="jet",
        clim=None,
        label=None,
        fig_name="",
        vec_img=None,
        mri2mni=None,
        bbox=None,
        ticks=None,
        ticklabels=None,
    ):
        import matplotlib
        import matplotlib.pyplot as plt  # noqa: F401 -- figure management

        img = np.asarray(img)
        if img.ndim != 3:
            raise ValueError("At least give us a volume to display (expected a 3D array).")
        img = img.astype(float) if not np.issubdtype(img.dtype, np.floating) else img
        full_shape = np.array(img.shape)
        finite = img[np.isfinite(img)]
        if finite.size == 0:
            raise ValueError("The image volume you provided does not have any meaningful values.")

        if clim is None:
            lo, hi = float(finite.min()), float(finite.max())
            clim = None if lo == hi else (lo, hi)
        self.clim = clim

        if vec_img is not None:
            vec_img = np.asarray(vec_img)
            if vec_img.shape != (*img.shape, 3):
                raise ValueError("Vector field does not have correct size.")

        if mri2mni is not None:
            mri2mni = np.asarray(mri2mni, dtype=float)
            if mri2mni.shape != (4, 4) or np.any(np.round(mri2mni[3]) != [0, 0, 0, 1]):
                raise ValueError("Unrecognized format of the voxel-to-MNI mapping.")
        self.mri2mni = mri2mni

        if bbox is None:
            bbox = np.array([[0, 0, 0], full_shape - 1])
        self.bbox = np.asarray(bbox, dtype=int)
        lo_, hi_ = self.bbox
        self.img = img[lo_[0] : hi_[0] + 1, lo_[1] : hi_[1] + 1, lo_[2] : hi_[2] + 1]
        self.vec = (
            None
            if vec_img is None
            else vec_img[lo_[0] : hi_[0] + 1, lo_[1] : hi_[1] + 1, lo_[2] : hi_[2] + 1]
        )

        if pos is None:
            pos = (full_shape - 1) // 2
        pos = np.round(np.asarray(pos, dtype=float)).astype(int) - lo_
        if np.any(pos < 0) or np.any(pos >= self.img.shape):
            raise ValueError("Voxel selected falls outside of the bounding box.")
        self.pos = pos  # in cropped coordinates

        if isinstance(cmap, str):
            cmap = matplotlib.colormaps[cmap]
        self.cmap = cmap.with_extremes(bad="white")
        self.label = label or ""
        self._updating = False

        self._build(fig_name, ticks, ticklabels)
        self._redraw()

    # -- coordinates ----------------------------------------------------------

    @property
    def voxel(self) -> np.ndarray:
        """Current position in full-volume, 0-based voxel coordinates."""
        return self.pos + self.bbox[0]

    @property
    def mni(self) -> np.ndarray | None:
        if self.mri2mni is None:
            return None
        return np.round(self.mri2mni @ np.append(self.voxel, 1.0))[:3].astype(int)

    @property
    def value(self) -> float:
        return float(self.img[tuple(self.pos)])

    def goto(self, voxel) -> bool:
        """Moves to a full-volume voxel; returns False (and stays) if outside."""
        new = np.round(np.asarray(voxel, dtype=float)).astype(int) - self.bbox[0]
        if np.any(new < 0) or np.any(new >= self.img.shape):
            return False
        self.pos = new
        self._redraw()
        return True

    def goto_mni(self, mni) -> bool:
        if self.mri2mni is None:
            raise ValueError("No voxel-to-MNI mapping for this viewer.")
        voxel = np.linalg.solve(self.mri2mni, np.append(np.asarray(mni, dtype=float), 1.0))[:3]
        return self.goto(voxel)

    # -- drawing --------------------------------------------------------------

    def _build(self, fig_name, ticks, ticklabels) -> None:
        from matplotlib.widgets import TextBox

        # Width-to-height ratio of sliceshow.m's window.
        self.fig = fig = new_figure(fig_name or "sliceshow", figsize=(8.5, 8.5 / 1.0187))
        fig._roast_viewer = self  # keep widgets alive with the figure
        if fig_name:
            fig.suptitle(fig_name.split(". ")[0], fontsize=11)
        gs = fig.add_gridspec(2, 2, left=0.02, right=0.98, bottom=0.02, top=0.95, wspace=0.02, hspace=0.06)
        self.axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0])]
        info = fig.add_subplot(gs[1, 1])
        info.axis("off")

        # (horizontal, vertical) volume axes for each panel: coronal, sagittal, axial.
        self._dims = [(0, 2), (1, 2), (0, 1)]
        side = max(self.img.shape)
        self._images, self._markers, self._vlines, self._hlines = [], [], [], []
        self._quivers = [None, None, None]
        for ax in self.axes:
            im = ax.imshow(
                np.zeros((2, 2)), origin="lower", cmap=self.cmap, interpolation="nearest", aspect="equal"
            )
            if self.clim is not None:
                im.set_clim(*self.clim)
            self._images.append(im)
            self._vlines.append(ax.axvline(0, color="0.15", lw=0.8))
            self._hlines.append(ax.axhline(0, color="0.15", lw=0.8))
            (marker,) = ax.plot([], [], "o", mfc="none", mec="m", mew=3, ms=12)
            self._markers.append(marker)
            ax.set_xlim(-0.5, side - 0.5)
            ax.set_ylim(-0.5, side - 0.5)
            ax.set_axis_off()

        cbar = fig.colorbar(
            self._images[1], ax=self.axes[1], location="bottom", fraction=0.05, pad=0.02, ticks=ticks
        )
        if ticklabels is not None:
            cbar.ax.set_xticklabels(ticklabels, fontsize=8)
        if self.label:
            cbar.set_label(self.label, fontsize=11)

        # Info panel: value readout and editable coordinates.
        box = info.get_position()
        x0, y0, w, h = box.x0, box.y0, box.width, box.height
        self._label_text = fig.text(x0 + w / 2, y0 + 0.40 * h, self.label, ha="center", fontsize=15, weight="bold")
        self._value_text = fig.text(x0 + w / 2, y0 + 0.28 * h, "", ha="center", fontsize=15, weight="bold")

        self._boxes = {}
        rows = [("Voxel", 0.78)]
        if self.mri2mni is not None:
            rows.append(("MNI", 0.62))
        for name, row_y in rows:
            fig.text(x0 + 0.02 * w, y0 + (row_y + 0.03) * h, name, fontsize=13, weight="bold")
            for j, axis_name in enumerate("XYZ"):
                bx = fig.add_axes([x0 + (0.30 + 0.23 * j) * w, y0 + row_y * h, 0.15 * w, 0.10 * h])
                tb = TextBox(bx, f"{axis_name} ", initial="", textalignment="right")
                tb.label.set_fontsize(12)
                tb.text_disp.set_fontsize(12)
                tb.on_submit(lambda text, n=name, k=j: self._on_submit(n, k, text))
                self._boxes[(name, j)] = tb

        fig.canvas.mpl_connect("button_press_event", self._on_click)

    def _slice(self, panel: int) -> np.ndarray:
        x, y, z = self.pos
        return [self.img[:, y, :], self.img[x, :, :], self.img[:, :, z]][panel]

    def _redraw(self) -> None:
        for panel, ax in enumerate(self.axes):
            h, v = self._dims[panel]
            self._images[panel].set_data(self._slice(panel).T)
            # set_data doesn't resize the image; one pixel per voxel.
            self._images[panel].set_extent(
                (-0.5, self.img.shape[h] - 0.5, -0.5, self.img.shape[v] - 0.5)
            )
            if self.clim is None:
                self._images[panel].autoscale()
            self._vlines[panel].set_xdata([self.pos[h]] * 2)
            self._hlines[panel].set_ydata([self.pos[v]] * 2)
            self._markers[panel].set_data([self.pos[h]], [self.pos[v]])
            if self.vec is not None:
                self._draw_arrows(panel, ax)

        value = self.value
        self._value_text.set_text("nan" if np.isnan(value) else f"{value:.2f}")
        self._updating = True
        try:
            for j in range(3):
                self._boxes[("Voxel", j)].set_val(str(self.voxel[j]))
                if self.mri2mni is not None:
                    self._boxes[("MNI", j)].set_val(str(self.mni[j]))
        finally:
            self._updating = False
        self.fig.canvas.draw_idle()

    def _draw_arrows(self, panel: int, ax) -> None:
        if self._quivers[panel] is not None:
            self._quivers[panel].remove()
        h, v = self._dims[panel]
        sl = [slice(None)] * 3
        fixed = ({0, 1, 2} - {h, v}).pop()
        for d in (h, v):
            sl[d] = slice(0, None, _QUIVER_STEP)
        sl[fixed] = self.pos[fixed]
        vec = self.vec[tuple(sl)]  # (n_h, n_v, 3)
        uu, vv = vec[..., h], vec[..., v]
        gh, gv = np.meshgrid(
            np.arange(0, self.img.shape[h], _QUIVER_STEP),
            np.arange(0, self.img.shape[v], _QUIVER_STEP),
            indexing="ij",
        )
        ok = np.isfinite(uu) & np.isfinite(vv)
        if not ok.any():
            self._quivers[panel] = None
            return
        # MATLAB's autoscale fits the longest arrow to ~0.9 grid spacing,
        # then quiver3's scale argument multiplies that.
        longest = float(np.hypot(uu[ok], vv[ok]).max()) or 1.0
        scale = longest / (0.9 * _QUIVER_STEP * _QUIVER_SCALE)
        self._quivers[panel] = ax.quiver(
            gh[ok], gv[ok], uu[ok], vv[ok],
            angles="xy", scale_units="xy", scale=scale, color="k", width=0.002,
        )

    # -- interaction ------------------------------------------------------------

    def _on_click(self, event) -> None:
        if event.inaxes not in self.axes or event.button != 1 or event.xdata is None:
            return
        toolbar = getattr(self.fig.canvas, "toolbar", None)
        if toolbar is not None and getattr(toolbar, "mode", ""):
            return  # zooming/panning, not navigating
        panel = self.axes.index(event.inaxes)
        h, v = self._dims[panel]
        new = self.pos.copy()
        new[h], new[v] = int(round(event.xdata)), int(round(event.ydata))
        if np.all(new >= 0) and np.all(new < self.img.shape):
            self.pos = new
            self._redraw()

    def _on_submit(self, row: str, axis: int, text: str) -> None:
        if self._updating:
            return
        try:
            number = float(text)
        except ValueError:
            self._redraw()  # restore the box
            return
        if row == "Voxel":
            target = self.voxel.astype(float)
            target[axis] = number
            ok = self.goto(target)
        else:
            target = self.mni.astype(float)
            target[axis] = number
            ok = self.goto_mni(target)
        if not ok:
            self._redraw()


def sliceshow(
    img,
    pos=None,
    color="jet",
    clim=None,
    label=None,
    fig_name="",
    vec_img=None,
    mri2mni=None,
    bbox=None,
    **kwargs,
) -> SliceViewer:
    """sliceshow(img, pos, color, clim, label, figName, vecImg, mri2mni, bbox),
    argument for argument (see the module docstring for the 0-based
    coordinate convention). Returns the SliceViewer; show it with
    roast_py.viz.show(...) or matplotlib's plt.show()."""
    return SliceViewer(img, pos, color, clim, label, fig_name, vec_img, mri2mni, bbox, **kwargs)
