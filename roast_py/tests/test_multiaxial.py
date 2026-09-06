"""Tests for roast_py.segmentation.multiaxial.

The full end-to-end test runs real inference with the bundled pretrained
models (lib/multiaxial/*.h5) against the repo's example/subject1.nii --
takes ~2-3 minutes on CPU, so it's marked slow and skipped by default
(`pytest -m slow` to run it explicitly).
"""

import shutil

import numpy as np
import pytest

from roast_py.segmentation.multiaxial import _remap_to_roast_labels, find_model_dir

from .test_nifti import REPO_ROOT, SUBJECT1


def test_find_model_dir_locates_bundled_weights():
    model_dir = find_model_dir()
    for name in ("sagittal_model.h5", "axial_model.h5", "coronal_model.h5", "consensus_layer.h5"):
        assert (model_dir / name).exists()


def test_remap_to_roast_labels_matches_segment_py_mapping():
    # multiaxial labels: 0 background, 1 air, 2 gray, 3 white, 4 csf, 5 bone, 6 skin
    multiaxial_labels = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.uint8)
    roast_labels = _remap_to_roast_labels(multiaxial_labels)
    # ROAST labels: 0 background, 1 white, 2 gray, 3 csf, 4 bone, 5 skin, 6 air
    np.testing.assert_array_equal(roast_labels, [0, 6, 2, 1, 3, 4, 5])


@pytest.mark.slow
def test_segment_end_to_end_on_subject1(tmp_path):
    import nibabel as nib

    t1_copy = tmp_path / "subject1.nii"
    shutil.copy(SUBJECT1, t1_copy)

    from roast_py.segmentation.multiaxial import segment

    out_path = segment(str(t1_copy))
    assert out_path.endswith("subject1_multiaxial_masks.nii")

    out_img = nib.load(out_path)
    t1_img = nib.load(t1_copy)
    assert out_img.shape == t1_img.shape

    data = out_img.get_fdata()
    labels = set(np.unique(data).astype(int).tolist())
    assert labels <= {0, 1, 2, 3, 4, 5, 6}

    # Sanity ranges for a real head scan: brain tissue (white+gray) and skin
    # should each be a clearly non-trivial fraction of the volume, and the
    # rarer tissues (csf/bone/air) should be present but small.
    total = data.size
    frac = {label: np.sum(data == label) / total for label in range(7)}
    assert 0.02 < frac[1] < 0.20  # white
    assert 0.02 < frac[2] < 0.20  # gray
    assert 0.0 < frac[3] < 0.10  # csf
    assert 0.0 < frac[4] < 0.20  # bone
    assert 0.02 < frac[5] < 0.40  # skin
    assert 0.0 < frac[6] < 0.05  # air
