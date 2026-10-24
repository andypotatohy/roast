"""Gets SPM running through its official Python interface, spm-python
(https://github.com/spm/spm-python), installing what's missing.

spm-python runs a compiled build of SPM on the free MATLAB Runtime, so two
things have to be present:

1. The Python packages: spm-python and the compiled SPM they call
   (spm-runtime-r2025b, built for MATLAB Runtime R2025b), plus
   mpython-core and matlab-runtime. They are pinned below and installed
   with pip the first time SPM is needed. They are not part of roast_py's
   always-installed environment: they are only needed for the SPM
   routines, and the compiled SPM alone is an 80 MB wheel.
2. The MATLAB Runtime R2025b itself, several GB from MathWorks. It is
   installed with the `matlab-runtime` package's installer (by one of
   spm-python's authors), which downloads MathWorks' official installer
   and runs it silently. An existing installation is used if one is found
   (MathWorks' default locations, $MATLAB_RUNTIME_PATH, or a MATLAB R2025b
   installation).

Installing the runtime means accepting the MATLAB Runtime license. roast_py
never does that silently: it asks once, at an interactive prompt, or takes
the answer from `accept_license=True` / the environment variable
ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE=yes (for scripts and batch jobs).

By default the runtime goes into a per-user folder (see default_prefix()),
so no administrator rights are needed on Linux or macOS; set
$MATLAB_RUNTIME_PATH to put it elsewhere.

SPM always runs in a separate Python process (roast_py.spm.worker), never
in the caller's: the MATLAB Runtime is large, slow to start and ships its
own copies of common libraries, and roast_py's own process usually already
has TensorFlow and VTK loaded. On macOS that process has to be started
through MathWorks' `mwpython`, which matlab-runtime's `mwpython2` launcher
finds and runs.

Everything at module level is standard library, so this imports in an
environment that has nothing installed yet.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from ..dependencies import NO_AUTO_INSTALL_ENV_VAR, auto_install_disabled, installed_version

# The MATLAB release the pinned compiled SPM was built for. The runtime
# must be exactly this release.
MATLAB_RELEASE = "R2025b"

# Exact versions, as with the rest of roast_py's environment.
SPM_PACKAGES: dict[str, str] = {
    "spm-python": "25.1.2.post2",
    "mpython-core": "25.4rc1",
    "spm-runtime-r2025b": "25.1.2.post4",
    "matlab-runtime": "0.0.6",
}

LICENSE_ENV_VAR = "ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE"
PREFIX_ENV_VAR = "MATLAB_RUNTIME_PATH"  # matlab-runtime's own variable
LICENSE_URL = "https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html"
RUNTIME_DEPS_URL = "https://github.com/mathworks-ref-arch/container-images/tree/main/matlab-runtime-deps"

# Approximate: MathWorks' installer download plus the unpacked runtime.
APPROX_INSTALL_GB = 10


class MatlabRuntimeLicenseNotAccepted(RuntimeError):
    """The MATLAB Runtime is needed but its license hasn't been accepted."""


# --------------------------------------------------------------------------
# where the runtime lives
# --------------------------------------------------------------------------


def default_prefix() -> Path:
    """Per-user install location for the MATLAB Runtime.

    $MATLAB_RUNTIME_PATH if set; otherwise
    ~/.local/share/roast_py/MATLAB_Runtime on Linux ($XDG_DATA_HOME if set),
    ~/Library/Application Support/roast_py/MATLAB_Runtime on macOS,
    %LOCALAPPDATA%\\roast_py\\MATLAB_Runtime on Windows.
    """
    if os.environ.get(PREFIX_ENV_VAR):
        return Path(os.environ[PREFIX_ENV_VAR])
    system = platform.system()
    if system == "Windows":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif system == "Darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
    return base / "roast_py" / "MATLAB_Runtime"


def find_runtime(prefix: str | os.PathLike | None = None) -> Path | None:
    """An installed MATLAB Runtime (or MATLAB) of MATLAB_RELEASE, or None.

    Looks under `prefix` (default: default_prefix()) and everywhere
    matlab-runtime looks: MathWorks' default install locations,
    $MATLAB_PATH, and a `matlab` on the PATH.
    """
    from matlab_runtime.utils import find_runtime as _find

    found = _find(MATLAB_RELEASE, str(prefix or default_prefix()))
    return Path(found) if found else None


# --------------------------------------------------------------------------
# python packages
# --------------------------------------------------------------------------


def package_drift() -> list[tuple[str, str | None, str]]:
    """(package, installed version or None, pinned version) for each SPM
    package that isn't at its pinned version."""
    return [
        (name, installed_version(name), wanted)
        for name, wanted in SPM_PACKAGES.items()
        if installed_version(name) != wanted
    ]


def install_command() -> list[str]:
    return [sys.executable, "-m", "pip", "install", *(f"{n}=={v}" for n, v in SPM_PACKAGES.items())]


def ensure_packages(install_missing: bool = True) -> None:
    drift = package_drift()
    if not drift:
        return
    listing = ", ".join(f"{n} {have or '(none)'} -> {want}" for n, have, want in drift)
    if not install_missing or auto_install_disabled():
        raise ImportError(
            f"SPM's Python packages are missing or at other versions ({listing}). Install them with:\n\n"
            f"    {' '.join(install_command())}\n\n"
            f"(automatic installation is off: install_missing=False or {NO_AUTO_INSTALL_ENV_VAR} is set)"
        )
    print(f"Installing SPM's Python interface ({listing}) ...")
    result = subprocess.run(install_command())
    if result.returncode != 0 or package_drift():
        raise RuntimeError(
            f"Could not install SPM's Python packages (pip exit status {result.returncode}). "
            f"Try by hand:\n\n    {' '.join(install_command())}"
        )


# --------------------------------------------------------------------------
# the MATLAB Runtime
# --------------------------------------------------------------------------


def license_accepted_by_env() -> bool:
    return os.environ.get(LICENSE_ENV_VAR, "").strip().lower() in ("1", "yes", "y", "true", "on")


def _license_text(prefix: Path) -> str:
    return (
        f"ROAST's SPM routines need the free MATLAB Runtime {MATLAB_RELEASE} from MathWorks "
        f"(about {APPROX_INSTALL_GB} GB of disk space while installing).\n"
        f"roast_py can download and install it into:\n    {prefix}\n"
        "Installing it means accepting the MATLAB Runtime license agreement:\n"
        f"    {LICENSE_URL}\n"
        "(the installer also leaves a copy, matlabruntime_license_agreement.pdf, next to the runtime)."
    )


def _ask_for_license(prefix: Path, accept_license: bool | None) -> None:
    if accept_license is True or license_accepted_by_env():
        return
    if accept_license is False:
        raise MatlabRuntimeLicenseNotAccepted("The MATLAB Runtime license was declined (accept_license=False).")
    if not sys.stdin or not sys.stdin.isatty():
        raise MatlabRuntimeLicenseNotAccepted(
            _license_text(prefix)
            + "\n\nThis isn't an interactive session, so roast_py can't ask. To accept the license and "
            f"let roast_py install the runtime, set {LICENSE_ENV_VAR}=yes (or pass accept_license=True), "
            "or install it yourself with\n\n"
            f"    {sys.executable} -m roast_py.spm install --accept-license"
        )
    print(_license_text(prefix))
    answer = input("Accept the MATLAB Runtime license and install it now? Type 'yes' to accept: ")
    if answer.strip().lower() not in ("y", "yes"):
        raise MatlabRuntimeLicenseNotAccepted("The MATLAB Runtime license was not accepted; nothing installed.")


def install_runtime(prefix: str | os.PathLike | None = None, accept_license: bool | None = None) -> Path:
    """Downloads and installs MATLAB Runtime MATLAB_RELEASE into `prefix`
    (default: default_prefix()), after the license is accepted (see
    module docstring). Returns the runtime's folder.
    """
    from matlab_runtime import install

    prefix = Path(prefix or default_prefix()).expanduser().resolve()
    _ask_for_license(prefix, accept_license)

    prefix.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(prefix).free / 1e9
    if free_gb < APPROX_INSTALL_GB:
        print(
            f"warning: only {free_gb:.1f} GB free at {prefix}; the MATLAB Runtime needs about "
            f"{APPROX_INSTALL_GB} GB while installing."
        )

    # The installer downloads into a temporary folder; keep it on the same
    # disk as the runtime rather than in a (often small) system /tmp.
    download_dir = prefix / ".download"
    download_dir.mkdir(exist_ok=True)
    saved_tempdir = tempfile.tempdir
    tempfile.tempdir = str(download_dir)
    try:
        print(f"Installing MATLAB Runtime {MATLAB_RELEASE} into {prefix} (a large download; this takes a while) ...")
        # auto_answer=True answers the installer's own questions, including
        # its license question -- which _ask_for_license() has just put to the user.
        install(MATLAB_RELEASE, prefix=str(prefix), auto_answer=True)
    except Exception as error:
        raise RuntimeError(
            f"Could not download or install MATLAB Runtime {MATLAB_RELEASE} ({type(error).__name__}: {error}).\n\n"
            "The installer comes from MathWorks (https://ssd.mathworks.com), so that site must be reachable. "
            "Alternatively, install the runtime yourself from "
            "https://www.mathworks.com/products/compiler/matlab-runtime.html and set "
            f"{PREFIX_ENV_VAR} to the folder that contains its {MATLAB_RELEASE} folder."
        ) from error
    finally:
        tempfile.tempdir = saved_tempdir
        shutil.rmtree(download_dir, ignore_errors=True)

    found = find_runtime(prefix)
    if found is None:
        raise RuntimeError(f"MATLAB Runtime {MATLAB_RELEASE} installation finished but isn't in {prefix}.")
    return found


def ensure_runtime(accept_license: bool | None = None, install_missing: bool = True) -> Path:
    """The MATLAB Runtime folder, installing it first if there is none."""
    found = find_runtime()
    if found is not None:
        return found
    if not install_missing or auto_install_disabled():
        raise RuntimeError(
            f"MATLAB Runtime {MATLAB_RELEASE} is not installed, and automatic installation is off. Install it "
            f"with:\n\n    {sys.executable} -m roast_py.spm install --accept-license"
        )
    return install_runtime(accept_license=accept_license)


def ensure_spm(accept_license: bool | None = None, install_missing: bool = True) -> Path:
    """Makes sure SPM can run: the pinned Python packages, then the MATLAB
    Runtime. Returns the runtime's folder. Call this before starting any
    long work, so missing pieces (and the license question) come up first.
    """
    ensure_packages(install_missing=install_missing)
    return ensure_runtime(accept_license=accept_license, install_missing=install_missing)


# --------------------------------------------------------------------------
# running SPM
# --------------------------------------------------------------------------


def worker_command(job_path: str | os.PathLike) -> list[str]:
    """How to start roast_py.spm.worker on this platform."""
    if platform.system() == "Darwin":
        # MathWorks only supports the MATLAB Runtime from Python on macOS
        # through its mwpython launcher.
        launcher = "import sys; from matlab_runtime.cli import mwpython2; sys.exit(mwpython2(sys.argv[1:]))"
        return [sys.executable, "-c", launcher, "-variant", MATLAB_RELEASE, "-m", "roast_py.spm.worker", str(job_path)]
    return [sys.executable, "-m", "roast_py.spm.worker", str(job_path)]


def worker_env(runtime: Path) -> dict[str, str]:
    env = dict(os.environ)
    # The folder *containing* <release>/, which is what matlab-runtime
    # expects in MATLAB_RUNTIME_PATH.
    env[PREFIX_ENV_VAR] = str(Path(runtime).parent)
    package_parent = str(Path(__file__).resolve().parents[2])
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [package_parent, env.get("PYTHONPATH")]))
    return env


def run_job(job: dict, runtime: Path | None = None, accept_license: bool | None = None) -> None:
    """Runs one roast_py.spm.worker job (see worker.py) in its own process."""
    runtime = runtime or ensure_spm(accept_license=accept_license)
    with tempfile.TemporaryDirectory(prefix="roast_spm_") as tmp:
        job_path = Path(tmp) / "job.json"
        job_path.write_text(json.dumps(job))
        result = subprocess.run(worker_command(job_path), env=worker_env(runtime))
    if result.returncode != 0:
        hint = ""
        if platform.system() == "Linux":
            hint = (
                "\nIf the MATLAB Runtime failed to start, it may be missing system libraries; MathWorks lists "
                f"them for each release at {RUNTIME_DEPS_URL}."
            )
        raise RuntimeError(f"SPM failed (worker exit status {result.returncode}).{hint}")
