"""Quick simulation on ROAST's default head, the MNI152 template
(example/MNI152_T1_1mm.nii), with ROAST's default recipe (anode Fp1 1 mA,
cathode P4 -1 mA) -- the Python equivalent of MATLAB's `roast()`.

Run from the roast_py/ directory:

    python examples/quickstart.py                    # install the tested env if needed, then run
    python examples/quickstart.py --check-deps       # just report what differs from it
    python examples/quickstart.py --no-install-deps  # don't install anything
    python examples/quickstart.py --no-visualize     # skip the result figures

Works straight from a git clone: no `pip install` step needed. roast()
installs roast_py's tested environment (exact versions, with pip) into the
running interpreter when anything is missing or at a different version.
Needs Python 3.11-3.13; `conda create -n roast_py python=3.11` gives you
the exact Python it was tested on.

Takes about 10-15 minutes on 4 CPU cores, plus the dependency install on
the first run, and needs about 8 GB of memory (the FEM solve). See
README.md for what each step does and how it's been verified.
"""

import argparse
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]  # the roast_py/ project directory
REPO_ROOT = Path(__file__).resolve().parents[2]
MNI152 = REPO_ROOT / "example" / "MNI152_T1_1mm.nii"

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
        "--no-visualize",
        action="store_true",
        help="skip the result figures roast() shows at the end",
    )
    parser.add_argument(
        "--work-dir",
        default=str(Path(tempfile.gettempdir()) / "roast_py_quickstart"),
        help="where to copy the MRI and write outputs (default: %(default)s)",
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
    # .nii files) next to the input, so copy the MRI out of the repo's
    # example/ folder first rather than writing into it.
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    subj = work_dir / MNI152.name
    if not subj.exists():
        shutil.copy(MNI152, subj)

    # recipe defaults to {'Fp1': 1.0, 'P4': -1.0}
    result = roast(
        str(subj), install_missing=not args.no_install_deps, visualize=not args.no_visualize
    )

    print(f"Voltage volume:  {result.vol_v.shape}")
    print(f"E-field volume:  {result.vol_e.shape}")
    stem = f"{subj.stem}_{result.sim_tag}"
    print(f"Outputs saved under: {result.work_dir}")
    print(f"  {stem}_v.nii     -- voltage (mV)")
    print(f"  {stem}_e.nii     -- E-field, 3 components (V/m)")
    print(f"  {stem}_emag.nii  -- E-field magnitude (V/m)")
    print(f"Redraw the figures any time with: roast_py.review_res({str(subj)!r})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
