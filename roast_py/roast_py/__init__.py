"""Python port of ROAST's core simulation pipeline (roast/roast_target/reviewRes),
built to remove the MATLAB dependency. See /roast_py/README.md for status and scope.
"""

import importlib

__all__ = ["roast", "RoastResult", "DEFAULT_RECIPE", "check_dependencies", "install_dependencies"]

_LAZY = {
    "roast": ".roast",
    "RoastResult": ".roast",
    "DEFAULT_RECIPE": ".roast",
    "check_dependencies": ".dependencies",
    "install_dependencies": ".dependencies",
}


def __getattr__(name: str):
    """Resolves the public API lazily (PEP 562).

    Importing roast_py must not drag in TensorFlow, pandas and friends --
    partly because it's slow, but mainly because roast_py.dependencies has
    to stay reachable in an environment that is missing them. An eager
    `from .roast import roast` here would raise a raw ImportError from deep
    inside the import chain before the user could call the very helper
    meant to fix it.
    """
    module_name = _LAZY.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    if module_name == ".roast":
        # One clear message listing every missing dependency, instead of a
        # raw ImportError for whichever happens to be imported first.
        from .dependencies import check_dependencies

        check_dependencies()

    module = importlib.import_module(module_name, __name__)
    value = getattr(module, name)
    # Cache in the package namespace. This also has to *overwrite* the
    # submodule attribute that importing `.roast` just bound here --
    # otherwise `from roast_py import roast` would return the function the
    # first time and the module every time after.
    globals()[name] = value
    return value


def __dir__():
    return sorted(__all__)
