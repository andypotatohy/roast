"""Dependency preflight for roast_py.

Deliberately imports nothing but the standard library: this module has to
stay importable in exactly the situation it exists to fix (a fresh
environment where roast_py's third-party dependencies aren't installed
yet), so it can't depend on any of them.

Run a check, or install what's missing, from the command line::

    python -m roast_py.dependencies            # report what's missing
    python -m roast_py.dependencies --install  # pip-install what's missing

or from Python::

    from roast_py.dependencies import check_dependencies, install_dependencies
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class Dependency:
    import_name: str  # what `import x` uses
    pip_name: str  # what `pip install x` uses
    needed_for: str


# Everything roast() needs at runtime. Kept in sync with pyproject.toml's
# [project] dependencies -- see test_dependencies.py, which fails if the
# two drift apart.
REQUIRED: tuple[Dependency, ...] = (
    Dependency("numpy", "numpy", "arrays, used everywhere"),
    Dependency("scipy", "scipy", "morphology, interpolation, spline fitting"),
    Dependency("nibabel", "nibabel", "reading/writing NIfTI volumes"),
    Dependency("pandas", "pandas", "reading capInfo.xlsx electrode templates"),
    Dependency("openpyxl", "openpyxl", "pandas' .xlsx engine, for capInfo.xlsx"),
    Dependency("skimage", "scikit-image", "resampling in multiaxial segmentation"),
    Dependency("tensorflow", "tensorflow", "running the bundled multiaxial segmentation models"),
    # The bundled lib/multiaxial/*.h5 models were saved under Keras 2 and
    # cannot be loaded by Keras 3 (TF >= 2.16's default) -- see
    # segmentation/_keras_compat.py. tf-keras is what makes them loadable.
    Dependency("tf_keras", "tf-keras", "legacy Keras 2 runtime for the bundled .h5 models"),
)


def missing_dependencies(required: tuple[Dependency, ...] | None = None) -> list[Dependency]:
    """Which of `required` can't be imported in this interpreter.

    Uses importlib.util.find_spec rather than actually importing, so this
    stays fast and side-effect free -- importing tensorflow here would
    both cost seconds and resolve Keras before
    segmentation/_keras_compat.py gets its chance to set
    TF_USE_LEGACY_KERAS.

    `required` defaults to REQUIRED, resolved at call time rather than as
    a default argument value, so the module-level list stays overridable.
    """
    required = REQUIRED if required is None else required
    missing = []
    for dep in required:
        try:
            found = importlib.util.find_spec(dep.import_name) is not None
        except (ImportError, ValueError):
            found = False
        if not found:
            missing.append(dep)
    return missing


def install_command(deps: list[Dependency]) -> list[str]:
    """The pip command that would install `deps` into *this* interpreter."""
    return [sys.executable, "-m", "pip", "install", *(d.pip_name for d in deps)]


def format_missing(deps: list[Dependency]) -> str:
    lines = ["roast_py is missing these dependencies:", ""]
    width = max(len(d.pip_name) for d in deps)
    for dep in deps:
        lines.append(f"  {dep.pip_name:<{width}}  ({dep.needed_for})")
    lines += [
        "",
        "Install them with:",
        "",
        f"    {' '.join(install_command(deps))}",
        "",
        "or let roast_py do it for you:",
        "",
        "    python -m roast_py.dependencies --install",
    ]
    return "\n".join(lines)


def check_dependencies(raise_on_missing: bool = True) -> list[Dependency]:
    """Reports every missing dependency at once.

    Without this, a fresh environment surfaces them one at a time: you
    install the package the first ImportError named, re-run the several-
    minute pipeline, and hit the next one.
    """
    missing = missing_dependencies()
    if missing and raise_on_missing:
        raise ImportError(format_missing(missing))
    return missing


def install_dependencies(deps: list[Dependency] | None = None, quiet: bool = False) -> None:
    """pip-installs the missing dependencies into the running interpreter.

    Deliberately NOT called automatically on import or from roast(): this
    downloads on the order of a gigabyte (tensorflow alone), needs network
    access, and pip-installing into a conda environment can shadow or
    break conda-managed packages. It runs only when you ask it to.

    In a conda environment, prefer conda/mamba for numpy/scipy/pandas and
    keep pip for the rest.
    """
    deps = missing_dependencies() if deps is None else deps
    if not deps:
        if not quiet:
            print("All roast_py dependencies are already installed.")
        return

    cmd = install_command(deps)
    if not quiet:
        print("Installing: " + ", ".join(d.pip_name for d in deps))
        print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)

    still_missing = missing_dependencies()
    if still_missing:
        raise RuntimeError(
            "Some dependencies are still missing after install:\n\n"
            + format_missing(still_missing)
        )
    if not quiet:
        print("Done. All roast_py dependencies are installed.")


def _main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m roast_py.dependencies",
        description="Check (or install) roast_py's runtime dependencies.",
    )
    parser.add_argument("--install", action="store_true", help="pip-install whatever is missing")
    args = parser.parse_args(argv)

    missing = missing_dependencies()
    if not missing:
        print("All roast_py dependencies are installed.")
        return 0

    if args.install:
        install_dependencies(missing)
        return 0

    print(format_missing(missing))
    return 1


if __name__ == "__main__":
    raise SystemExit(_main())
