"""Tests for registration-based landmarks (roast.m's landmarksInTPM mapped
through runNiftyReg.m's registration).

The fast tests check the matrix algebra and conventions against MATLAB's
1-based formulas, with a stand-in for the reg_aladin binary. The slow test
runs the real bundled reg_aladin on example/subject1.nii and checks the
landmarks land on the scalp.
"""

import os
import stat
import sys

import nibabel as nib
import numpy as np
import pytest

from roast_py.geometry.landmarks import LANDMARK_NAMES, LANDMARKS_IN_TPM, tpm_landmarks_to_subject
from roast_py.registration import Registration, registration_path, run_niftyreg
from roast_py.registration.niftyreg import find_mni_template, find_reg_aladin, find_tpm

from .test_nifti import SUBJECT1

# SPM's .mat takes 1-based voxel indices: mat = affine @ SHIFT.
SHIFT = np.array([[1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1], [0, 0, 0, 1]], float)


def _registration():
    tpm = nib.load(str(find_tpm())).affine
    subject = nib.load(str(SUBJECT1)).affine
    # A plausible subject->MNI affine: small rotation, scaling, shift.
    angle = np.deg2rad(7)
    affine = np.array(
        [
            [0.97, 0, 0, 2.0],
            [0, 1.05 * np.cos(angle), -np.sin(angle), -20.0],
            [0, np.sin(angle), 1.05 * np.cos(angle), -60.0],
            [0, 0, 0, 1],
        ]
    )
    return Registration(affine, subject, tpm)


def test_landmark_table_is_roast_ms():
    assert LANDMARKS_IN_TPM.shape == (16, 3) == (len(LANDMARK_NAMES), 3)
    assert LANDMARKS_IN_TPM[0].tolist() == [61, 139, 98]  # nasion
    assert LANDMARKS_IN_TPM[3].tolist() == [111, 63, 93]  # left (eTPM is LAS)
    assert LANDMARKS_IN_TPM[-1].tolist() == [61.2171, 141.6240, 126.5092]


def test_mapping_equals_matlabs_1_based_formula_minus_one():
    reg = _registration()
    # roast.m, verbatim, with SPM's 1-based matrices:
    image_mat, tpm_mat = reg.image_affine @ SHIFT, reg.tpm_affine @ SHIFT
    tpm2mri = np.linalg.inv(image_mat) @ np.linalg.inv(reg.affine) @ tpm_mat
    matlab = (tpm2mri @ np.column_stack([LANDMARKS_IN_TPM, np.ones(16)]).T).T[:, :3]
    matlab = np.sign(matlab) * np.floor(np.abs(matlab) + 0.5)  # MATLAB round()

    assert np.array_equal(tpm_landmarks_to_subject(reg.tpm2mri), matlab.astype(int) - 1)


def test_mri2mni_is_matlabs_mapping_on_0_based_voxels():
    reg = _registration()
    voxel0 = np.array([10.0, 20.0, 30.0, 1.0])
    matlab_mri2mni = reg.affine @ (reg.image_affine @ SHIFT)  # Affine*image(1).mat
    assert np.allclose(reg.mri2mni @ voxel0, matlab_mri2mni @ (voxel0 + [1, 1, 1, 0]))


def test_rounding_goes_half_away_from_zero_like_matlab():
    # A pure translation that puts TPM landmark coordinates exactly on .5
    # (numpy's round() would send 2.5 to 2; MATLAB's to 3).
    tpm2mri = np.eye(4)
    tpm2mri[:3, 3] = 0.5
    out = tpm_landmarks_to_subject(tpm2mri)
    assert out[0].tolist() == [61, 139, 98]  # 0-based 60.5 -> 1-based 61.5 -> 62 -> 61


def test_registration_round_trips_through_json(tmp_path):
    reg = _registration()
    path = tmp_path / "x_niftyReg.json"
    reg.save(path)
    loaded = Registration.load(path)
    assert np.allclose(loaded.affine, reg.affine)
    assert np.allclose(loaded.tpm2mri, reg.tpm2mri)


def test_bundled_files_are_found():
    assert find_mni_template().name == "MNI152_T1_1mm.nii"
    assert find_tpm().name == "eTPM.nii"
    if sys.platform.startswith("linux"):
        assert os.access(find_reg_aladin(), os.X_OK)


def _fake_reg_aladin(tmp_path, matrix, exit_code=0):
    """A stand-in reg_aladin that writes `matrix` to the -aff path."""
    rows = "\\n".join(" ".join(str(v) for v in row) for row in matrix)
    script = tmp_path / "reg_aladin"
    script.write_text(
        "#!/bin/sh\n"
        'while [ "$#" -gt 0 ]; do\n'
        '  case "$1" in -aff) aff="$2"; shift;; -res) res="$2"; shift;; esac; shift\n'
        "done\n"
        f'[ {exit_code} -ne 0 ] && {{ echo "boom" >&2; exit {exit_code}; }}\n'
        f'printf "{rows}\\n" > "$aff"; : > "$res"\n'
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


@pytest.mark.skipif(sys.platform.startswith("win"), reason="shell-script stand-in")
def test_run_niftyreg_inverts_reg_aladins_matrix_and_saves_it(tmp_path):
    ref_to_flo = np.array([[1.03, 0.02, 0, -1.8], [-0.017, 0.92, -0.12, 25.5], [0.006, 0.13, 0.88, 62.8], [0, 0, 0, 1]])
    fake = _fake_reg_aladin(tmp_path, ref_to_flo)

    reg = run_niftyreg(SUBJECT1, out_dir=tmp_path, bin_path=fake)

    assert np.allclose(reg.affine, np.linalg.inv(ref_to_flo))  # "consistent with SPM format"
    assert np.allclose(reg.image_affine, nib.load(str(SUBJECT1)).affine)
    assert np.allclose(reg.tpm_affine, nib.load(str(find_tpm())).affine)
    saved = registration_path(SUBJECT1, tmp_path)
    assert saved.name == "subject1_niftyReg.json"
    assert np.allclose(Registration.load(saved).affine, reg.affine)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["reg_aladin", "subject1_niftyReg.json"]  # temp files gone


@pytest.mark.skipif(sys.platform.startswith("win"), reason="shell-script stand-in")
def test_run_niftyreg_reports_a_failure(tmp_path):
    fake = _fake_reg_aladin(tmp_path, np.eye(4), exit_code=3)
    with pytest.raises(RuntimeError, match=r"(?s)niftyReg failed.*boom"):
        run_niftyreg(SUBJECT1, out_dir=tmp_path, bin_path=fake)


@pytest.mark.slow
def test_real_registration_puts_the_landmarks_on_subject1s_scalp(tmp_path):
    """runs the bundled reg_aladin (~2 min) and checks the head landmarks
    against a scalp mask: the head's outline at an intensity threshold."""
    from scipy.ndimage import binary_fill_holes, distance_transform_edt

    reg = run_niftyreg(SUBJECT1, out_dir=tmp_path)
    landmarks = tpm_landmarks_to_subject(reg.tpm2mri)

    t1 = np.asarray(nib.load(str(SUBJECT1)).dataobj, dtype=float)
    head = binary_fill_holes(t1 > 0.1 * np.percentile(t1, 99))
    outside, inside = distance_transform_edt(~head), distance_transform_edt(head)
    for name, point in zip(LANDMARK_NAMES[:4], landmarks[:4]):
        p = tuple(point)
        distance = outside[p] if not head[p] else inside[p]
        assert distance <= 5, f"{name} at {point.tolist()} is {distance:.1f} voxels from the scalp"

    nasion, inion, right, left = landmarks[:4]
    assert nasion[1] > inion[1] + 150  # nasion well anterior of inion
    assert right[0] > left[0] + 120  # RAS: right ear at high x
    # Cz-ish: the vertex electrode on the sagittal line sits above the scalp center.
    assert landmarks[11][2] > landmarks[6][2] + 50
    # MNI coordinates of the nasion are near the template's (0, ~85, ~-40..0).
    mni = reg.mri2mni @ np.append(nasion, 1)
    assert abs(mni[0]) < 10 and 70 < mni[1] < 100
