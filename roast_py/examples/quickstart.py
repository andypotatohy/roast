"""Quick simulation on the bundled subject1.nii with ROAST's default
recipe (anode Fp1 1 mA, cathode P4 -1 mA) -- the Python equivalent of
MATLAB's `roast('example/subject1.nii')`.

Run from the roast_py/ directory:

    python examples/quickstart.py                 # run the simulation
    python examples/quickstart.py --install-deps  # install missing deps first, then run
    python examples/quickstart.py --check-deps    # just report what's missing

Takes several minutes on CPU: ~2-3 min for segmentation, ~1-2 min for
meshing + the FEM solve. See the top-level README for what each phase
does and how it's been verified.
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

# Importing roast_py is cheap and never requires the heavy dependencies --
# only resolving roast_py.roast does. That's what lets this script check
# for (and install) them before touching the pipeline.
from roast_py.dependencies import (  # noqa: E402
    format_missing,
    install_dependencies,
    missing_dependencies,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--install-deps",
        action="store_true",
        help="pip-install any missing dependencies before running (downloads ~1GB: TensorFlow is large)",
    )
    parser.add_argument(
        "--check-deps", action="store_true", help="report missing dependencies and exit"
    )
    parser.add_argument(
        "--work-dir",
        default="/tmp/roast_py_quickstart",
        help="where to copy the subject and write outputs (default: %(default)s)",
    )
    args = parser.parse_args()

    missing = missing_dependencies()

    if args.check_deps:
        if missing:
            print(format_missing(missing))
            return 1
        print("All roast_py dependencies are installed.")
        return 0

    if missing:
        if args.install_deps:
            install_dependencies(missing)
        else:
            print(format_missing(missing))
            print("\nOr re-run this script with --install-deps to do that automatically.")
            return 1

    # Imported only after the dependency check, so a missing dependency
    # produces the message above rather than a raw ImportError.
    from roast_py import roast

    # roast() writes its outputs (.msh, .pro, .pos, and the final _v/_e/_emag
    # .nii files) next to the input, so copy subject1.nii out of the repo's
    # example/ folder first rather than writing into it.
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    subj = work_dir / "subject1.nii"
    shutil.copy(SUBJECT1, subj)

    result = roast(str(subj))  # recipe defaults to {'Fp1': 1.0, 'P4': -1.0}

    print(f"Voltage volume:  {result.vol_v.shape}")
    print(f"E-field volume:  {result.vol_e.shape}")
    print(f"Outputs saved under: {result.work_dir}")
    print(f"  {subj.stem}_v.nii     -- voltage (mV)")
    print(f"  {subj.stem}_e.nii     -- E-field, 3 components (V/m)")
    print(f"  {subj.stem}_emag.nii  -- E-field magnitude (V/m)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
