"""Showing and saving roast_py's figures.

MATLAB draws every ROAST figure in its own non-blocking window. Python has
two separate GUI toolkits in play -- matplotlib for the slice viewers
(sliceshow) and VTK, through PyVista, for the 3D renderings -- each with its
own event loop, so this module is what makes them behave like one session:

* all 3D views go into a single PyVista window, one panel each, with
  linked cameras (they share world coordinates);
* while that window is open, a VTK timer pumps matplotlib's events, so the
  slice viewers stay clickable at the same time;
* once the 3D window is closed, matplotlib keeps the slice viewers open
  until they are closed too.

Without a display (a remote server, CI) nothing is shown and every figure
is saved as a PNG instead. The 3D views additionally need OpenGL, which a
headless machine often lacks; that is probed in a throwaway subprocess,
because VTK failing to get a context can abort the whole process rather
than raise.
"""

from __future__ import annotations

import functools
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

# matplotlib backends that can only write files, never open a window.
_FILE_ONLY_BACKENDS = {"agg", "pdf", "ps", "svg", "pgf", "cairo", "template"}

# MATLAB's view(3): azimuth -37.5 deg, elevation 30 deg, z up.
MATLAB_VIEW3 = (-0.5299, -0.6906, 0.5)


def display_available() -> bool:
    """Whether a window could be opened at all."""
    if sys.platform.startswith(("linux", "freebsd")):
        return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    return True  # macOS and Windows always have a window server for a desktop user


def matplotlib_can_show() -> bool:
    """Whether matplotlib picked a backend that displays figures.

    With a display but no GUI toolkit (no tkinter, no Qt), matplotlib
    silently falls back to Agg, which can only save files.
    """
    import matplotlib

    return matplotlib.get_backend().lower() not in _FILE_ONLY_BACKENDS


_GL_PROBE = (
    "import pyvista as pv\n"
    "p = pv.Plotter(off_screen=True, window_size=(64, 64))\n"
    "p.add_mesh(pv.Sphere())\n"
    "p.screenshot(return_img=True)\n"
    "p.close()\n"
)


@functools.cache
def can_render_3d() -> bool:
    """Whether VTK can get an OpenGL context here (probed once, out of process)."""
    try:
        result = subprocess.run(
            [sys.executable, "-c", _GL_PROBE], capture_output=True, timeout=120
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


@dataclass
class Scene3D:
    """One 3D view: a title plus a function that draws into the current
    PyVista renderer (``draw(plotter)``)."""

    title: str
    draw: Callable


@dataclass
class FigureSet:
    """The figures from one roast()/review_res() call, shown or saved together.

    ``figures`` are matplotlib figures (mostly slice viewers, each kept
    alive through its figure), ``scenes`` the 3D views.
    """

    title: str = "ROAST"
    figures: list = field(default_factory=list)
    scenes: list[Scene3D] = field(default_factory=list)

    def extend(self, other: FigureSet) -> None:
        self.figures += other.figures
        self.scenes += other.scenes

    # -- 3D ---------------------------------------------------------------

    def build_plotter(self, off_screen: bool):
        import pyvista as pv

        n = len(self.scenes)
        plotter = pv.Plotter(
            shape=(1, n),
            window_size=(620 * n, 680),
            off_screen=off_screen,
            title=f"{self.title} -- 3D views. Drag to rotate.",
        )
        plotter.set_background("white")
        for i, scene in enumerate(self.scenes):
            plotter.subplot(0, i)
            scene.draw(plotter)
            plotter.add_text(scene.title, position="upper_edge", font_size=10, color="black")
            plotter.view_vector(MATLAB_VIEW3, viewup=(0, 0, 1))
            plotter.reset_camera()
        if n > 1:
            plotter.link_views()
        return plotter

    # -- saving -------------------------------------------------------------

    def save(self, directory: str | os.PathLike, include_3d: bool = True) -> list[Path]:
        """Saves every figure (and, if OpenGL is available, the 3D views) as PNG."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        saved = []
        for i, fig in enumerate(self.figures, start=1):
            path = directory / f"{i:02d}_{_slug(_figure_name(fig))}.png"
            fig.savefig(path, dpi=100, facecolor="white")
            saved.append(path)
        if include_3d and self.scenes and can_render_3d():
            path = directory / f"{len(self.figures) + 1:02d}_3d_views.png"
            plotter = self.build_plotter(off_screen=True)
            plotter.screenshot(str(path))
            plotter.close()
            saved.append(path)
        return saved

    # -- showing ------------------------------------------------------------

    def show(self, fallback_dir: str | os.PathLike | None = None, block: bool = True) -> list[Path]:
        """Shows everything that can be shown; saves the rest to `fallback_dir`.

        Returns the paths of any files saved. With `block` (the default,
        unless matplotlib is in interactive mode) this returns once every
        window has been closed -- a script that exits would take its
        windows with it.
        """
        import matplotlib.pyplot as plt

        has_display = display_available()
        show_2d = has_display and matplotlib_can_show() and bool(self.figures)
        show_3d = has_display and bool(self.scenes) and can_render_3d()
        block = block and not plt.isinteractive()

        saved: list[Path] = []
        unshown_2d = self.figures if not show_2d else []
        unshown_3d = bool(self.scenes) and not show_3d
        if (unshown_2d or unshown_3d) and fallback_dir is not None:
            partial = FigureSet(self.title, list(unshown_2d), self.scenes if unshown_3d else [])
            saved = partial.save(fallback_dir)
            _explain_fallback(saved, fallback_dir, has_display, unshown_2d, unshown_3d)
        elif unshown_3d and not can_render_3d():
            print("3D views skipped: VTK could not get an OpenGL context on this machine.")

        if show_2d:
            plt.show(block=False)
            plt.pause(0.001)

        if show_3d:
            plotter = self.build_plotter(off_screen=False)
            if show_2d:
                # Keep the slice viewers responsive while VTK's event loop
                # owns the main thread.
                plotter.add_timer_event(max_steps=10**9, duration=50, callback=_pump_matplotlib)
            plotter.show(auto_close=True)

        if show_2d and block and plt.get_fignums():
            plt.show(block=True)

        if not show_2d:
            for fig in self.figures:
                plt.close(fig)
        return saved


def _pump_matplotlib(_step=None) -> None:
    import matplotlib.pyplot as plt

    for num in plt.get_fignums():
        try:
            plt.figure(num).canvas.flush_events()
        except Exception:  # noqa: BLE001 -- a closing window must not kill the 3D view
            pass


def _figure_name(fig) -> str:
    manager = getattr(fig.canvas, "manager", None)
    name = getattr(fig, "_roast_name", None) or (manager.get_window_title() if manager else "")
    return name or "figure"


def _slug(text: str) -> str:
    text = re.sub(r"[.:]?\s*Click anywhere to navigate\.?\s*$", "", text)
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()[:60] or "figure"


def _explain_fallback(saved, directory, has_display, unshown_2d, unshown_3d) -> None:
    if not has_display:
        reason = "No display available"
    elif unshown_2d:
        reason = (
            "matplotlib has no GUI backend here (install tkinter, or `pip install "
            "PyQt6`, for interactive slice viewers)"
        )
    else:
        reason = "3D windows unavailable"
    what = f"{len(saved)} figure(s)" if saved else "no figures"
    print(f"{reason}; saved {what} as PNG in {directory}")
    if unshown_3d and not can_render_3d():
        print("  (3D views skipped: VTK could not get an OpenGL context on this machine)")


def new_figure(name: str, figsize: tuple[float, float]):
    """A pyplot-managed figure whose window title (and saved filename) is `name`."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize)
    fig._roast_name = name
    manager = getattr(fig.canvas, "manager", None)
    if manager is not None:
        manager.set_window_title(name)
    return fig
