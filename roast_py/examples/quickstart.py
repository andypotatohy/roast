"""Quick simulation on the bundled subject1.nii with ROAST's default
recipe (anode Fp1 1 mA, cathode P4 -1 mA) -- the Python equivalent of
MATLAB's `roast('example/subject1.nii')`.

Run from the roast_py/ directory (after `pip install -e ".[multiaxial]"`):

    python examples/quickstart.py

Takes several minutes on CPU: ~2-3 min for segmentation, ~1-2 min for
meshing + the FEM solve. See the top-level README for what each phase
does and how it's been verified.
"""

import shutil
from pathlib import Path

from roast_py import roast

REPO_ROOT = Path(__file__).resolve().parents[2]
SUBJECT1 = REPO_ROOT / "example" / "subject1.nii"

# roast() writes its outputs (.msh, .pro, .pos, and the final _v/_e/_emag
# .nii files) next to the input by default -- copy subject1.nii out of the
# repo's example/ folder first rather than writing into it.
work_dir = Path("/tmp/roast_py_quickstart")
work_dir.mkdir(exist_ok=True)
subj = work_dir / "subject1.nii"
shutil.copy(SUBJECT1, subj)

result = roast(str(subj))  # recipe defaults to {'Fp1': 1.0, 'P4': -1.0}

print(f"Voltage volume:  {result.vol_v.shape}")
print(f"E-field volume:  {result.vol_e.shape}")
print(f"Outputs saved under: {result.work_dir}")
print(f"  {subj.stem}_v.nii     -- voltage (mV)")
print(f"  {subj.stem}_e.nii     -- E-field, 3 components (V/m)")
print(f"  {subj.stem}_emag.nii  -- E-field magnitude (V/m)")
