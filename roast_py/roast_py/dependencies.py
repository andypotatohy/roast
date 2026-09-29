"""Dependency preflight and installer for roast_py.

roast_py pins its runtime to one exact, verified environment rather than
to version ranges. That environment ran the full pipeline end to end --
`roast("../example/subject1.nii")`, segmentation through FEM solve --
and is listed in TESTED_ENVIRONMENT below (and, identically, in
requirements-lock.txt at the project root). roast() checks the running
interpreter against it before doing any work and installs the exact
versions with pip when anything differs.

Why exact pins: the bundled segmentation models need TensorFlow plus the
tf-keras compatibility package on *matching* versions, and "any recent
version" of each has proven fragile in practice (a conda-forge TensorFlow
newer than any tf-keras release, pip backtracking to a tf-keras that calls
TensorFlow internals removed in 2.16, ...). Pinning the whole stack to
what was actually run removes that whole class of problems.

This module deliberately imports nothing but the standard library: it has
to stay importable in exactly the situation it exists to fix (a fresh
environment with none of roast_py's dependencies yet).

Check or install from the command line::

    python -m roast_py.dependencies            # report what differs
    python -m roast_py.dependencies --install  # install the tested versions

The recommended setup is a fresh environment on the tested Python::

    conda create -n roast_py python=3.11 -y
    conda activate roast_py
    python -m roast_py.dependencies --install   # or just call roast()

Everything is installed with pip, into the running interpreter
(`sys.executable -m pip`), even inside a conda environment -- that is how
the tested environment was built. Conda only provides Python itself.
"""

from __future__ import annotations

import importlib.metadata
import importlib.util
import os
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# the tested environment
# --------------------------------------------------------------------------

# The interpreter the environment below was verified on.
TESTED_PYTHON = "3.11.15"

# Python versions the exact pins below can be installed on (inclusive).
# 3.10 is out because keras 3.15 / numpy 2.4 / scipy 1.17 / pandas 3.0 all
# require >= 3.11; 3.14 is out because TensorFlow 2.21 ships no wheels for
# it. Wheels exist for Linux (x86_64 and aarch64, glibc >= 2.27), macOS on
# Apple silicon, and Windows x86_64 -- not for Intel Macs, which TensorFlow
# stopped supporting after 2.16.
SUPPORTED_PYTHON = ((3, 11), (3, 13))

# Exactly what was installed (all by pip, none by conda) when roast()
# last ran end to end: roast_py's direct dependencies plus everything they
# pull in, as reported by importlib.metadata. setuptools, wheel, packaging
# and six are left out on purpose -- they're interpreter plumbing that
# came with the base Python and nothing here depends on their version. So
# is pyobjc-framework-Cocoa, which PyVista pulls in on macOS only; pip
# picks its version there.
#
# Kept identical to requirements-lock.txt; tests/test_dependencies.py
# fails if the two drift apart.
TESTED_ENVIRONMENT: dict[str, str] = {
    "absl-py": "2.5.0",
    "astunparse": "1.6.3",
    "attrs": "26.1.0",
    "certifi": "2026.2.25",
    "charset-normalizer": "3.4.6",
    "contourpy": "1.3.3",
    "cycler": "0.12.1",
    "cyclopts": "5.1.0",
    "docstring-parser": "0.18.0",
    "et-xmlfile": "2.0.0",
    "flatbuffers": "25.12.19",
    "fonttools": "4.66.1",
    "gast": "0.7.0",
    "google-pasta": "0.2.0",
    "grpcio": "1.83.1",
    "h5py": "3.14.0",
    "idna": "3.11",
    "imageio": "2.37.4",
    "importlib-resources": "7.1.0",
    "keras": "3.15.1",
    "kiwisolver": "1.5.1",
    "lazy-loader": "0.5",
    "libclang": "18.1.1",
    "markdown-it-py": "4.2.0",
    "matplotlib": "3.11.2",
    "mdurl": "0.1.2",
    "ml-dtypes": "0.6.0",
    "namex": "0.1.0",
    "networkx": "3.6.1",
    "nibabel": "5.4.2",
    "numpy": "2.4.6",
    "openpyxl": "3.1.5",
    "opt-einsum": "3.4.0",
    "optree": "0.20.0",
    "pandas": "3.0.5",
    "pillow": "12.3.0",
    "platformdirs": "4.12.2",
    "pooch": "1.9.0",
    "protobuf": "7.36.1",
    "pygments": "2.21.0",
    "pyparsing": "3.1.1",
    "python-dateutil": "2.9.0.post0",
    "pyvista": "0.49.0",
    "pyvista-validation": "0.2.2",
    "requests": "2.33.1",
    "rich": "15.0.0",
    "rich-rst": "2.2.0",
    "scikit-image": "0.26.0",
    "scipy": "1.17.1",
    "scooby": "0.12.0",
    "tensorflow": "2.21.0",
    "termcolor": "3.3.0",
    "tf-keras": "2.21.0",
    "tifffile": "2026.3.3",
    "typing-extensions": "4.16.0",
    "urllib3": "2.6.3",
    "vtk": "9.7.1",
    "wrapt": "2.4.0",
}

# Set this to stop roast() from installing anything (for CI, or an
# environment you manage yourself). Missing packages and a broken
# TensorFlow/tf-keras pairing still raise; other version differences only
# warn.
NO_AUTO_INSTALL_ENV_VAR = "ROAST_PY_NO_AUTO_INSTALL"


@dataclass(frozen=True)
class Dependency:
    import_name: str  # what `import x` uses
    pip_name: str  # what `pip install x` uses
    needed_for: str


# The packages roast_py itself imports. Kept in sync with pyproject.toml's
# [project] dependencies -- see test_dependencies.py.
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
    # segmentation/_keras_compat.py.
    Dependency("tf_keras", "tf-keras", "legacy Keras 2 runtime for the bundled .h5 models"),
    Dependency("matplotlib", "matplotlib", "slice viewers (sliceshow, viewMRI, viewSeg)"),
    Dependency("pyvista", "pyvista", "3D renderings (viewElectrodes, visualizeRes)"),
)


# --------------------------------------------------------------------------
# detection
# --------------------------------------------------------------------------


def python_version_problem(version_info: tuple[int, ...] | None = None) -> str | None:
    """Explains why this Python can't host the tested environment, or None."""
    version_info = tuple(version_info or sys.version_info)
    lo, hi = SUPPORTED_PYTHON
    if lo <= version_info[:2] <= hi:
        return None
    running = ".".join(str(p) for p in version_info[:3])
    return (
        f"roast_py needs Python {lo[0]}.{lo[1]}-{hi[0]}.{hi[1]} (tested on "
        f"{TESTED_PYTHON}), but this is Python {running} ({sys.executable}).\n\n"
        "TensorFlow and the rest of roast_py's pinned dependencies have no "
        "builds for this version. Create an environment on the tested Python "
        "and run roast() from there:\n\n" + fresh_environment_instructions()
    )


def missing_dependencies(required: tuple[Dependency, ...] | None = None) -> list[Dependency]:
    """Which of `required` (default REQUIRED) can't be imported here.

    Uses importlib.util.find_spec rather than importing, so it's fast and
    side-effect free -- importing tensorflow here would resolve Keras
    before segmentation/_keras_compat.py can set TF_USE_LEGACY_KERAS.
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


def installed_version(dist_name: str) -> str | None:
    """Version of an installed distribution (from its metadata), or None."""
    try:
        return importlib.metadata.version(dist_name)
    except importlib.metadata.PackageNotFoundError:
        return None


def installed_by(dist_name: str) -> str | None:
    """Which tool installed a distribution ('pip', 'conda', ...), if recorded."""
    try:
        text = importlib.metadata.distribution(dist_name).read_text("INSTALLER")
    except importlib.metadata.PackageNotFoundError:
        return None
    return text.strip() if text else None


def version_drift() -> list[tuple[str, str | None, str]]:
    """(package, installed version or None, tested version) for every
    package in TESTED_ENVIRONMENT that isn't at its tested version."""
    drift = []
    for name, wanted in TESTED_ENVIRONMENT.items():
        have = installed_version(name)
        if have != wanted:
            drift.append((name, have, wanted))
    return drift


def _major_minor(version: str | None) -> tuple[int, int] | None:
    if not version:
        return None
    parts = version.split(".")
    try:
        return int(parts[0]), int(parts[1])
    except (IndexError, ValueError):
        return None


def tensorflow_keras_mismatch() -> str | None:
    """Describes a TensorFlow/tf-keras version mismatch, or None if fine.

    tf-keras X.Y only works with tensorflow X.Y.*; anything else fails on
    import, e.g. tf-keras < 2.16 beside TensorFlow >= 2.16 raises
    `AttributeError: ... no attribute 'register_load_context_function'`.
    """
    tf_version = installed_version("tensorflow")
    keras_version = installed_version("tf-keras")
    tf_mm, keras_mm = _major_minor(tf_version), _major_minor(keras_version)
    if tf_mm is None or keras_mm is None or tf_mm == keras_mm:
        return None
    return (
        f"tensorflow {tf_version} and tf-keras {keras_version} are incompatible "
        f"(tf-keras {keras_mm[0]}.{keras_mm[1]} requires tensorflow "
        f"{keras_mm[0]}.{keras_mm[1]}.*).\n\n"
        "This is what produces errors like:\n"
        "    AttributeError: module 'tensorflow._api.v2.compat.v2.__internal__'\n"
        "    has no attribute 'register_load_context_function'\n\n"
        "Install the tested versions with:\n\n"
        "    python -m roast_py.dependencies --install"
    )


# --------------------------------------------------------------------------
# commands and messages
# --------------------------------------------------------------------------


def lock_specs() -> list[str]:
    """TESTED_ENVIRONMENT as pip requirement specs ('name==version')."""
    return [f"{name}=={version}" for name, version in TESTED_ENVIRONMENT.items()]


def install_command() -> list[str]:
    """The pip command that installs the tested environment into *this*
    interpreter -- all pins in one call, so pip resolves them together."""
    return [sys.executable, "-m", "pip", "install", *lock_specs()]


def lock_file() -> Path:
    """requirements-lock.txt at the project root (exists in a source checkout)."""
    return Path(__file__).resolve().parents[1] / "requirements-lock.txt"


def fresh_environment_instructions() -> str:
    tested_minor = ".".join(TESTED_PYTHON.split(".")[:2])
    return (
        f"    conda create -n roast_py python={tested_minor} -y\n"
        "    conda activate roast_py\n"
        "    python -m roast_py.dependencies --install   # or just call roast()"
    )


def format_missing(deps: list[Dependency]) -> str:
    lines = ["roast_py is missing these dependencies:", ""]
    width = max(len(d.pip_name) for d in deps)
    for dep in deps:
        lines.append(f"  {dep.pip_name:<{width}}  ({dep.needed_for})")
    lines += [
        "",
        "Install the tested versions with:",
        "",
        "    python -m roast_py.dependencies --install",
        "",
        "or let roast() do it for you (it installs them automatically unless",
        f"install_missing=False or {NO_AUTO_INSTALL_ENV_VAR} is set).",
    ]
    return "\n".join(lines)


def format_drift(drift: list[tuple[str, str | None, str]]) -> str:
    width = max(len(name) for name, _, _ in drift)
    lines = ["These packages differ from roast_py's tested environment:", ""]
    for name, have, wanted in drift:
        lines.append(f"  {name:<{width}}  installed {have or '(none)':<14} tested {wanted}")
    return "\n".join(lines)


def check_dependencies(raise_on_missing: bool = True) -> list[Dependency]:
    """Reports every missing dependency at once, plus a Python version or
    TensorFlow/tf-keras problem that would stop roast() from running."""
    missing = missing_dependencies()
    if raise_on_missing:
        problem = python_version_problem()
        if problem:
            raise RuntimeError(problem)
        if missing:
            raise ImportError(format_missing(missing))
        mismatch = tensorflow_keras_mismatch()
        if mismatch:
            raise ImportError(mismatch)
    return missing


# --------------------------------------------------------------------------
# installation
# --------------------------------------------------------------------------

_SMOKE_TEST = (
    "import os; os.environ['TF_USE_LEGACY_KERAS'] = '1'; "
    "import tensorflow, tf_keras, nibabel, scipy, skimage, pandas, openpyxl, "
    "matplotlib, pyvista"
)


def smoke_test() -> str | None:
    """Imports the heavy dependencies in a fresh interpreter, the way
    roast() will. Returns the error output, or None if it worked.

    A separate process because this one mustn't import tensorflow before
    _keras_compat does, and because pip may just have replaced packages
    this process already has loaded.
    """
    result = subprocess.run([sys.executable, "-c", _SMOKE_TEST], capture_output=True, text=True)
    if result.returncode == 0:
        return None
    return (result.stderr or result.stdout).strip()


def auto_install_disabled() -> bool:
    return os.environ.get(NO_AUTO_INSTALL_ENV_VAR, "").strip() not in ("", "0", "false", "False")


def install_dependencies(quiet: bool = False) -> None:
    """Installs the tested environment into the running interpreter with pip.

    One `pip install name==version ...` for every pinned package; pip
    skips the ones already at their tested version. TensorFlow alone is
    several hundred MB, so this says what it's about to do first.
    """
    problem = python_version_problem()
    if problem:
        raise RuntimeError(problem)

    drift = version_drift()
    if not drift:
        if not quiet:
            print("roast_py's tested environment is already installed.")
        return

    conda_owned = [name for name, have, _ in drift if have and installed_by(name) == "conda"]
    cmd = install_command()
    if not quiet:
        print(format_drift(drift))
        print()
        print(f"Installing roast_py's tested environment with pip into {sys.prefix}")
        if conda_owned:
            print(
                "  note: replacing conda-installed "
                + ", ".join(conda_owned)
                + " with the tested pip builds"
            )
        print("  (TensorFlow is a large download; this can take several minutes)")
        print(f"  {' '.join(cmd[:4])} <{len(cmd) - 4} pinned packages>")

    returncode = subprocess.run(cmd).returncode
    remaining = version_drift()
    broken = smoke_test() if not remaining else None

    if returncode != 0 or remaining or broken:
        details = []
        if returncode != 0:
            details.append(f"pip exited with status {returncode}.")
        if remaining:
            details.append(format_drift(remaining))
        if broken:
            details.append("Importing TensorFlow/tf-keras still fails:\n\n" + broken)
        raise RuntimeError(
            f"Could not install roast_py's tested environment into {sys.prefix}.\n\n"
            + "\n\n".join(details)
            + "\n\nThis usually means the environment already holds packages pip "
            "can't cleanly replace (a conda-installed TensorFlow, for instance). "
            "The reliable fix is a fresh environment:\n\n"
            + fresh_environment_instructions()
        )
    if not quiet:
        print("Done. roast_py's tested environment is installed.")


def ensure_dependencies(install_missing: bool = True, quiet: bool = False) -> None:
    """Preflight used by roast(): make this interpreter match the tested
    environment, installing it if needed.

    With `install_missing` false (or ROAST_PY_NO_AUTO_INSTALL set) nothing
    is installed: missing packages and a TensorFlow/tf-keras mismatch
    raise, and other version differences only warn.
    """
    problem = python_version_problem()
    if problem:
        raise RuntimeError(problem)

    drift = version_drift()
    if not drift:
        return

    if install_missing and not auto_install_disabled():
        install_dependencies(quiet=quiet)
        return

    missing = missing_dependencies()
    if missing:
        raise ImportError(format_missing(missing))
    mismatch = tensorflow_keras_mismatch()
    if mismatch:
        raise ImportError(mismatch)
    warnings.warn(
        format_drift(drift)
        + "\n\nContinuing anyway because automatic installation is disabled. "
        "Run `python -m roast_py.dependencies --install` to match it.",
        stacklevel=2,
    )


def _main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m roast_py.dependencies",
        description="Check (or install) roast_py's tested dependency environment.",
    )
    parser.add_argument(
        "--install", action="store_true", help="install the tested versions of everything"
    )
    args = parser.parse_args(argv)

    problem = python_version_problem()
    if problem:
        print(problem)
        return 1

    drift = version_drift()
    if not drift:
        print(
            f"roast_py's tested environment is installed (Python "
            f"{sys.version.split()[0]}, tested on {TESTED_PYTHON})."
        )
        return 0

    if args.install:
        install_dependencies()
        return 0

    print(format_drift(drift))
    print("\nInstall the tested versions with:\n\n    python -m roast_py.dependencies --install")
    return 1


if __name__ == "__main__":
    raise SystemExit(_main())
