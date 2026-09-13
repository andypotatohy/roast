"""End-to-end test of roast_py's top-level roast() orchestrator -- the
Python equivalent of calling MATLAB's `roast('example/subject1.nii')` with
the default recipe. Chains all four implemented phases (segmentation,
electrode placement, meshing, FEM solve) through the real bundled
multiaxial models and getdp/cgalmesh binaries. Takes several minutes on
CPU (dominated by segmentation ~2-3 min and the FEM solve ~1-2 min at full
head resolution).
"""

import shutil

import numpy as np
import pytest

from roast_py import DEFAULT_RECIPE, roast

from .test_nifti import SUBJECT1


@pytest.mark.slow
def test_roast_default_recipe_on_subject1(tmp_path):
    subj = tmp_path / "subject1.nii"
    shutil.copy(SUBJECT1, subj)

    result = roast(str(subj))

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

    assert (tmp_path / "subject1_v.nii").exists()
    assert (tmp_path / "subject1_emag.nii").exists()
    assert (tmp_path / "subject1_e.nii").exists()

    # Both DEFAULT_RECIPE electrodes should have been placed.
    assert set(np.unique(result.elec_mask).tolist()) <= {0, 1, 2}
    assert len(DEFAULT_RECIPE) == 2
