"""Dependency preflight and installer for roast_py.

Deliberately imports nothing but the standard library: this module has to
stay importable in exactly the situation it exists to fix (a fresh
environment where roast_py's third-party dependencies aren't installed
yet), so it can't depend on any of them.

Check or install from the command line::

    python -m roast_py.dependencies            # report what's missing
    python -m roast_py.dependencies --install  # install what's missing

or from Python::

    from roast_py.dependencies import check_dependencies, install_dependencies

roast() calls install_dependencies() itself when something is missing (see
its `install_missing` argument), so an end user normally never has to.

Conda vs pip
------------
In a conda environment the installer uses `conda install` for the
scientific stack and pip for the rest, which is the recommended ordering
when the two are mixed (conda first, pip last, so conda's solver sees the
environment before pip writes into it). Outside conda it uses pip for
everything. Either way it targets the *running interpreter's* environment
explicitly -- `--prefix sys.prefix` for conda, `sys.executable -m pip` for
pip -- rather than whatever environment happens to be active, which are
not always the same thing.

TensorFlow and tf-keras are installed with pip even under conda: pip is
TensorFlow's official distribution channel, and tf-keras (the Keras 2
compatibility package the bundled .h5 models need) is a recent package
whose conda-forge availability is not something this code should assume.
Anything conda fails to install falls back to pip automatically, so a
wrong guess about channel contents self-corrects instead of dead-ending.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

# Set this to disable roast()'s automatic installation (for CI, locked
# environments, reproducible builds).
NO_AUTO_INSTALL_ENV_VAR = "ROAST_PY_NO_AUTO_INSTALL"

CONDA_CHANNEL = "conda-forge"


@dataclass(frozen=True)
class Dependency:
    import_name: str  # what `import x` uses
    pip_name: str  # what `pip install x` uses
    needed_for: str
    conda_name: str | None = None  # None => install with pip even under conda


# Everything roast() needs at runtime. Kept in sync with pyproject.toml's
# [project] dependencies -- see test_dependencies.py, which fails if the
# two drift apart.
REQUIRED: tuple[Dependency, ...] = (
    Dependency("numpy", "numpy", "arrays, used everywhere", "numpy"),
    Dependency("scipy", "scipy", "morphology, interpolation, spline fitting", "scipy"),
    Dependency("nibabel", "nibabel", "reading/writing NIfTI volumes", "nibabel"),
    Dependency("pandas", "pandas", "reading capInfo.xlsx electrode templates", "pandas"),
    Dependency("openpyxl", "openpyxl", "pandas' .xlsx engine, for capInfo.xlsx", "openpyxl"),
    Dependency("skimage", "scikit-image", "resampling in multiaxial segmentation", "scikit-image"),
    # pip is TensorFlow's official distribution channel; conda builds lag
    # and solve slowly.
    Dependency("tensorflow", "tensorflow", "running the bundled multiaxial segmentation models"),
    # The bundled lib/multiaxial/*.h5 models were saved under Keras 2 and
    # cannot be loaded by Keras 3 (TF >= 2.16's default) -- see
    # segmentation/_keras_compat.py.
    Dependency("tf_keras", "tf-keras", "legacy Keras 2 runtime for the bundled .h5 models"),
)


# --------------------------------------------------------------------------
# detection
# --------------------------------------------------------------------------


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


def in_conda_environment() -> bool:
    """Whether the *running interpreter* lives in a conda environment.

    Checks for the conda-meta directory next to sys.prefix, which conda
    creates in every environment it manages. Deliberately not based on
    CONDA_PREFIX/CONDA_DEFAULT_ENV: those describe the shell's activated
    environment, which isn't necessarily the one this interpreter belongs
    to (e.g. an activated env shelling out to a different python).
    """
    return (Path(sys.prefix) / "conda-meta").is_dir()


def find_conda() -> str | None:
    """Path to a usable conda-family executable, or None.

    Prefers $CONDA_EXE (set by conda's own shell integration, and points
    at the installation this environment came from), then mamba/micromamba
    /conda on PATH.
    """
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe and Path(conda_exe).exists():
        return conda_exe
    for name in ("mamba", "conda", "micromamba"):
        found = shutil.which(name)
        if found:
            return found
    return None


def use_conda_for(deps: list[Dependency]) -> tuple[list[Dependency], list[Dependency]]:
    """Splits `deps` into (install with conda, install with pip).

    Everything goes to pip unless the running interpreter is in a conda
    environment, a conda executable is available, and the package declares
    a conda_name.
    """
    if not (in_conda_environment() and find_conda()):
        return [], list(deps)
    conda_deps = [d for d in deps if d.conda_name]
    pip_deps = [d for d in deps if not d.conda_name]
    return conda_deps, pip_deps


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------


def pip_install_command(deps: list[Dependency]) -> list[str]:
    """The pip command that installs `deps` into *this* interpreter."""
    return [sys.executable, "-m", "pip", "install", *(d.pip_name for d in deps)]


def conda_install_command(deps: list[Dependency], conda_exe: str | None = None) -> list[str]:
    """The conda command that installs `deps` into *this* interpreter's env.

    Uses `--prefix sys.prefix` rather than relying on the activated
    environment, so it writes into the environment the running interpreter
    actually belongs to -- the conda equivalent of `sys.executable -m pip`.
    """
    conda_exe = conda_exe or find_conda() or "conda"
    return [
        conda_exe,
        "install",
        "--prefix",
        sys.prefix,
        "-c",
        CONDA_CHANNEL,
        "-y",
        *(d.conda_name or d.pip_name for d in deps),
    ]


def install_commands(deps: list[Dependency]) -> list[list[str]]:
    """Every command needed to install `deps`, conda before pip."""
    conda_deps, pip_deps = use_conda_for(deps)
    commands = []
    if conda_deps:
        commands.append(conda_install_command(conda_deps))
    if pip_deps:
        commands.append(pip_install_command(pip_deps))
    return commands


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def format_missing(deps: list[Dependency]) -> str:
    lines = ["roast_py is missing these dependencies:", ""]
    width = max(len(d.pip_name) for d in deps)
    for dep in deps:
        lines.append(f"  {dep.pip_name:<{width}}  ({dep.needed_for})")
    lines += ["", "Install them with:", ""]
    for cmd in install_commands(deps):
        lines.append(f"    {' '.join(cmd)}")
    lines += [
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


# --------------------------------------------------------------------------
# installation
# --------------------------------------------------------------------------


def auto_install_disabled() -> bool:
    return os.environ.get(NO_AUTO_INSTALL_ENV_VAR, "").strip() not in ("", "0", "false", "False")


def install_dependencies(
    deps: list[Dependency] | None = None, quiet: bool = False
) -> None:
    """Installs the missing dependencies into the running interpreter's env.

    Uses conda for the scientific stack when the interpreter is in a conda
    environment (see module docstring), pip otherwise. Anything conda
    fails on is retried with pip, so a package that isn't on conda-forge
    doesn't dead-end the install.

    This downloads on the order of a gigabyte (TensorFlow alone), so it
    reports what it's about to do before doing it.
    """
    deps = missing_dependencies() if deps is None else deps
    if not deps:
        if not quiet:
            print("All roast_py dependencies are already installed.")
        return

    conda_deps, pip_deps = use_conda_for(deps)

    if not quiet:
        print("Installing missing roast_py dependencies: " + ", ".join(d.pip_name for d in deps))
        print("  target environment: " + sys.prefix)
        if conda_deps:
            print("  via conda: " + ", ".join(d.conda_name or d.pip_name for d in conda_deps))
        if pip_deps:
            print("  via pip:   " + ", ".join(d.pip_name for d in pip_deps))
        print("  (TensorFlow is a large download; this can take several minutes)")

    if conda_deps:
        cmd = conda_install_command(conda_deps)
        if not quiet:
            print("  " + " ".join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            # Most likely one of these isn't on the channel. Rather than
            # dead-end, hand the whole conda set to pip.
            if not quiet:
                print("  conda install failed; falling back to pip for those packages")
            pip_deps = conda_deps + pip_deps

    if pip_deps:
        cmd = pip_install_command(pip_deps)
        if not quiet:
            print("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)

    still_missing = missing_dependencies(tuple(deps))
    if still_missing:
        raise RuntimeError(
            "Some dependencies are still missing after install:\n\n" + format_missing(still_missing)
        )
    if not quiet:
        print("Done. All roast_py dependencies are installed.")


def ensure_dependencies(install_missing: bool = True, quiet: bool = False) -> None:
    """Preflight used by roast(): check, and optionally install, in one call.

    With `install_missing` false (or ROAST_PY_NO_AUTO_INSTALL set) this
    raises the usual actionable ImportError instead of installing.
    """
    missing = missing_dependencies()
    if not missing:
        return

    if not install_missing or auto_install_disabled():
        raise ImportError(format_missing(missing))

    install_dependencies(missing, quiet=quiet)


def _main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m roast_py.dependencies",
        description="Check (or install) roast_py's runtime dependencies.",
    )
    parser.add_argument("--install", action="store_true", help="install whatever is missing")
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
