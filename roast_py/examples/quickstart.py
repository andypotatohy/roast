"""Quick simulation on the bundled subject1.nii with ROAST's default
recipe (anode Fp1 1 mA, cathode P4 -1 mA) -- the Python equivalent of
MATLAB's `roast('example/subject1.nii')`.

Run from the roast_py/ directory:

    python examples/quickstart.py                    # install the tested env if needed, then run
    python examples/quickstart.py --check-deps       # just report what differs from it
    python examples/quickstart.py --no-install-deps  # don't install anything

Works straight from a git clone: no `pip install` step needed. roast()
installs roast_py's tested environment (exact versions, with pip) into the
running interpreter when anything is missing or at a different version.
Needs Python 3.11-3.13; `conda create -n roast_py python=3.11` gives you
the exact Python it was tested on.

Takes several minutes on CPU: ~2-3 min for segmentation, ~1-2 min for
meshing + the FEM solve, plus the dependency install on first run. See the
top-level README for what each phase does and how it's been verified.
"""

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]  # the roast_py/ project directory
REPO_ROOT = Path(__file__).resolve().parents[2]
SUBJECT1 = REPO_ROOT / "example" / "subject1.nii"

# Work straight from a git clone, with no `pip install -e .` step: running
# a script in examples/ puts examples/ on sys.path, not the project root,
# so roast_py wouldn't otherwise be importable.
if importlib.util.find_spec("roast_py") is None:
    sys.path.insert(0, str(PACKAGE_ROOT))

# Importing roast_py never requires the heavy dependencies -- that is what
# lets this script (and roast() itself) install them when they're missing.
from roast_py.dependencies import (  # noqa: E402
    format_drift,
    format_missing,
    missing_dependencies,
    python_version_problem,
    version_drift,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--no-install-deps",
        action="store_true",
        help="install nothing: fail if dependencies are missing, warn if versions differ",
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="report how this environment differs from the tested one, then exit",
    )
    parser.add_argument(
        "--work-dir",
        default="/tmp/roast_py_quickstart",
        help="where to copy the subject and write outputs (default: %(default)s)",
    )
    args = parser.parse_args()

    if args.check_deps:
        problem = python_version_problem()
        if problem:
            print(problem)
            return 1
        drift = version_drift()
        if drift:
            print(format_drift(drift))
            return 1
        print("roast_py's tested environment is installed.")
        return 0

    missing = missing_dependencies()
    if missing and args.no_install_deps:
        print(format_missing(missing))
        return 1

    from roast_py import roast  # importable even with nothing installed

    # roast() writes its outputs (.msh, .pro, .pos, and the final _v/_e/_emag
    # .nii files) next to the input, so copy subject1.nii out of the repo's
    # example/ folder first rather than writing into it.
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    subj = work_dir / "subject1.nii"
    shutil.copy(SUBJECT1, subj)

    # recipe defaults to {'Fp1': 1.0, 'P4': -1.0}
    result = roast(str(subj), install_missing=not args.no_install_deps)

    print(f"Voltage volume:  {result.vol_v.shape}")
    print(f"E-field volume:  {result.vol_e.shape}")
    print(f"Outputs saved under: {result.work_dir}")
    print(f"  {subj.stem}_v.nii     -- voltage (mV)")
    print(f"  {subj.stem}_e.nii     -- E-field, 3 components (V/m)")
    print(f"  {subj.stem}_emag.nii  -- E-field magnitude (V/m)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
