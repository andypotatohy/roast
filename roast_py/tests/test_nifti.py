"""Regression tests for roast_py.io.nifti against ROAST's MATLAB behavior.

example/MNI152_T1_1mm.nii is documented in the top-level README ("you'll
notice this if you just run ROAST on the default subject... as it's in LAS
orientation") as ROAST's own worked example of a non-RAS input, which makes
it a natural fixture here.
"""

import shutil
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from roast_py.io.nifti import align_header_to_mni, convert_to_ras, zero_pad

REPO_ROOT = Path(__file__).resolve().parents[2]
MNI152 = REPO_ROOT / "example" / "MNI152_T1_1mm.nii"
SUBJECT1 = REPO_ROOT / "example" / "subject1.nii"


@pytest.fixture
def mni152_copy(tmp_path):
    dst = tmp_path / "MNI152_T1_1mm.nii"
    shutil.copy(MNI152, dst)
    return dst


@pytest.fixture
def subject1_copy(tmp_path):
    dst = tmp_path / "subject1.nii"
    shutil.copy(SUBJECT1, dst)
    return dst


def test_convert_to_ras_flips_las_mni152(mni152_copy):
    img = nib.load(mni152_copy)
    assert nib.aff2axcodes(img.affine) == ("L", "A", "S")

    result = convert_to_ras(str(mni152_copy))
    assert result.was_reoriented is True
    assert result.path.endswith("_ras.nii")

    out = nib.load(result.path)
    assert nib.aff2axcodes(out.affine) == ("R", "A", "S")
    assert out.shape == img.shape

    # World coordinates of anatomical content must be preserved, i.e. this
    # is a re-labeling of the same volume, not a resample: derive the
    # expected flipped affine directly from ROAST's own translation-update
    # formula (srow(4) = -origin*srow(1:3)) and compare.
    orig = img.affine
    expected = orig.copy()
    expected[0, 0] *= -1
    expected[0, 3] = orig[0, 0] * (img.shape[0] - 1) + orig[0, 3]
    assert np.allclose(out.affine, expected)


def test_convert_to_ras_is_noop_for_already_ras_subject(subject1_copy):
    img = nib.load(subject1_copy)
    axcodes = nib.aff2axcodes(img.affine)

    result = convert_to_ras(str(subject1_copy))
    if axcodes == ("R", "A", "S"):
        assert result.was_reoriented is False
        assert result.path == str(subject1_copy)
    else:
        # If subject1.nii isn't already RAS, just check it round-trips to RAS.
        out = nib.load(result.path)
        assert nib.aff2axcodes(out.affine) == ("R", "A", "S")


def test_zero_pad_preserves_world_coordinates(mni152_copy):
    ras = convert_to_ras(str(mni152_copy))
    img = nib.load(ras.path)

    padded_path = zero_pad(ras.path, 10)
    padded = nib.load(padded_path)

    assert padded.shape == tuple(s + 20 for s in img.shape)

    data = np.asarray(img.dataobj)
    padded_data = np.asarray(padded.dataobj)
    np.testing.assert_array_equal(
        padded_data[10:-10, 10:-10, 10:-10], data
    )

    # A voxel and its shifted counterpart in the padded volume must map to
    # the same world coordinate.
    vox = np.array([50, 60, 70, 1])
    world_before = img.affine @ vox
    world_after = padded.affine @ (vox + [10, 10, 10, 0])
    assert np.allclose(world_before, world_after)


def test_zero_pad_rejects_already_padded_name(tmp_path):
    fake = tmp_path / "foo_padded10.nii"
    fake.touch()
    with pytest.raises(ValueError):
        zero_pad(str(fake), 10)


def test_align_header_to_mni_forces_sform(mni152_copy, tmp_path):
    target = np.eye(4)
    target[:3, 3] = [1.0, 2.0, 3.0]

    out_path = align_header_to_mni(str(mni152_copy), target, "myscan")
    out = nib.load(out_path)

    assert out_path.endswith("myscan_MNI.nii")
    assert np.allclose(out.affine, target)
    sform, code = out.header.get_sform(coded=True)
    assert code == 1
    assert np.allclose(sform, target)
