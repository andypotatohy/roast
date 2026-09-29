"""Ports ROAST's visualization: reviewRes.m, visualizeRes.m, viewMRI.m,
viewSeg.m, viewElectrodes.m, sliceshow.m and brainCrop.m.

Slice viewers use matplotlib; 3D renderings use PyVista (VTK). See
_display.py for how the two are shown together.

    from roast_py.viz import review_res
    review_res("example/subject1.nii")

Names resolve lazily so `import roast_py.viz` (and review_res itself)
works before matplotlib/PyVista are installed -- review_res() installs
them, as roast() does.
"""

import importlib

_LAZY = {
    "review_res": ".review",
    "visualize_res": ".results",
    "show_roast_results": ".results",
    "sliceshow": ".sliceshow",
    "SliceViewer": ".sliceshow",
    "view_mri": ".views",
    "view_seg": ".views",
    "view_electrodes": ".views",
    "brain_crop": ".views",
    "FigureSet": "._display",
}

__all__ = sorted(_LAZY)


def __getattr__(name: str):
    module_name = _LAZY.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return __all__
