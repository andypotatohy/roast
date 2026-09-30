"""End-to-end test of roast_py's top-level roast() orchestrator -- the
Python equivalent of calling MATLAB's `roast()` with the default recipe
(anode Fp1 1 mA, cathode P4 -1 mA). Chains every phase (segmentation,
NiftyReg registration + TPM landmarks, electrode placement, meshing, FEM
solve) through the real bundled multiaxial models and getdp/cgalmesh/
reg_aladin binaries, on example/MNI152_T1_1mm.nii -- MATLAB roast()'s own
default subject.

MNI152 rather than subject1 because at roast.m's mesh settings getDP's
direct solver needs over 8 GB of RAM for subject1, more than some test
machines have; the MNI152 head fits comfortably. Takes ~10 minutes on CPU.
"""

import shutil

import numpy as np
import pytest

from roast_py import DEFAULT_RECIPE, roast

from .test_nifti import MNI152


@pytest.mark.slow
def test_roast_default_recipe_on_mni152(tmp_path):
    from roast_py.io.nifti import convert_to_ras

    shutil.copy(MNI152, tmp_path / "MNI152_T1_1mm.nii")
    # MNI152_T1_1mm.nii is stored LAS; roast() needs RAS input.
    subj = convert_to_ras(str(tmp_path / "MNI152_T1_1mm.nii")).path
    base = "MNI152_T1_1mm_ras"

    result = roast(subj, visualize=False)

    assert result.vol_v.shape == result.tissue_labels.shape
    assert result.vol_e.shape == (*result.tissue_labels.shape, 3)
    assert not np.all(np.isnan(result.vol_v))
    assert not np.all(np.isnan(result.ef_mag))

    # Physically sane magnitudes for a 1 mA montage (not the 10^10-scale
    # nonsense a badly-set-up model produces -- see fem/solve.py's tests).
    v_valid = result.vol_v[~np.isnan(result.vol_v)]
    assert 0 <= v_valid.min()
    assert v_valid.max() < 10000

    ef_valid = result.ef_mag[~np.isnan(result.ef_mag)]
    assert ef_valid.max() > 0
    assert ef_valid.max() < 10000

    for suffix in ("v.nii", "emag.nii", "e.nii"):
        assert (tmp_path / f"{base}_{suffix}").exists(), suffix

    # Everything review_res() needs to redraw the results later.
    for name in ("roastOptions.json", "mask_elec.nii", "mask_gel.nii", "mesh.npz", "v.pos", "e.pos",
                 "niftyReg.json"):
        assert (tmp_path / f"{base}_{name}").exists(), name

    # Both DEFAULT_RECIPE electrodes placed, and where the 10-10 system puts
    # them: this head *is* MNI space, so the registration is near identity.
    assert set(np.unique(result.elec_mask).tolist()) == {0, 1, 2}
    assert len(DEFAULT_RECIPE) == 2
    centers = [(result.mri2mni @ np.append(np.argwhere(result.elec_mask == i).mean(0), 1))[:3] for i in (1, 2)]
    fp1, p4 = centers
    assert fp1[0] < -10 and fp1[1] > 60  # left frontal
    assert p4[0] > 30 and p4[1] < -50 and p4[2] > 20  # right parietal
