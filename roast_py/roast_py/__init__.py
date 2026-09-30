"""Python port of ROAST's core simulation pipeline (roast/roast_target/reviewRes),
built to remove the MATLAB dependency. See /roast_py/README.md for status and scope.
"""

import importlib
import sys
import types

__all__ = [
    "roast",
    "RoastResult",
    "DEFAULT_RECIPE",
    "check_dependencies",
    "install_dependencies",
    "ensure_dependencies",
    "review_res",
]

_LAZY = {
    "roast": ".roast",
    "RoastResult": ".roast",
    "DEFAULT_RECIPE": ".roast",
    "check_dependencies": ".dependencies",
    "install_dependencies": ".dependencies",
    "ensure_dependencies": ".dependencies",
    "review_res": ".viz.review",
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

    # No dependency check here on purpose: `from roast_py import roast`
    # must succeed in an environment with nothing installed, because
    # roast() is what installs the missing dependencies. roast.py keeps
    # its own imports inside the function body to make that possible.
    module = importlib.import_module(module_name, __name__)
    # Cache every name this module provides, so later lookups skip this hook.
    for attr, source in _LAZY.items():
        if source == module_name:
            globals()[attr] = getattr(module, attr)
    return globals()[name]


def __dir__():
    return sorted(__all__)


class _Package(types.ModuleType):
    """Keeps `roast_py.roast` pointing at the roast() function.

    The function lives in roast.py, and whenever anything imports that
    submodule, Python binds the *module* as the package attribute `roast`,
    replacing the function -- so `from roast_py import roast` would start
    returning a module depending on what had been imported before. This
    catches that binding and stores the module's function instead. The
    module itself stays importable as usual (`from roast_py.roast import
    output_paths`, via sys.modules).
    """

    def __setattr__(self, name, value):
        if isinstance(value, types.ModuleType) and name in _SAME_NAME_AS_MODULE:
            value = getattr(value, name)
        super().__setattr__(name, value)


# Public functions that share their name with the submodule defining them.
_SAME_NAME_AS_MODULE = {"roast"}
sys.modules[__name__].__class__ = _Package
