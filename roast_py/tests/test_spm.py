"""Tests for roast_py.spm (SPM through spm-python) and its Python-side steps.

The MATLAB Runtime isn't available to these tests, so SPM itself never
runs here: the tests cover what roast_py does around it -- the batch jobs
(against start_seg.m / realignT2.m / resampToOneMM.m), their conversion to
spm-python objects, reading a real SPM `_seg8.mat`, the resampling grid,
the atlas-mask sampling through SPM's deformation, segTouchup's patching
passes, and the installer's license and location handling.
"""

import os
import subprocess
import sys
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from roast_py.spm import runtime
from roast_py.spm.batches import TISSUE_NGAUS, coregister_batch, reslice_batch, segment_batch, to_matlabbatch
from roast_py.spm.routines import (
    _same_space,
    one_mm_grid,
    sample_atlas_mask,
    spm_output_name,
)

from .test_nifti import MNI152, SUBJECT1

REPO = Path(__file__).resolve().parents[2]
SEG8 = REPO / "example" / "nyhead_T1orT2_seg8.mat"
NYHEAD_MASKS_ZIP = REPO / "example" / "nyhead_T1orT2_SPM_masks.nii.zip"


# --------------------------------------------------------------------------
# batches
# --------------------------------------------------------------------------


def test_segment_batch_is_start_seg_m():
    job = segment_batch("/data/t1.nii", tpm="/roast/eTPM.nii")["spm"]["spatial"]["preproc"]
    assert job["channel"] == [{"vols": ["/data/t1.nii,1"], "biasreg": 0.001, "biasfwhm": 60, "write": [0, 0]}]
    assert [t["ngaus"] for t in job["tissue"]] == list(TISSUE_NGAUS) == [1, 1, 2, 3, 4, 2]
    assert [t["tpm"] for t in job["tissue"]] == [[f"/roast/eTPM.nii,{i}"] for i in range(1, 7)]
    assert all(t["native"] == [1, 0] and t["warped"] == [0, 0] for t in job["tissue"])
    assert job["warp"] == {
        "reg": [0, 0.001, 0.5, 0.05, 0.2], "affreg": "mni", "samp": 3,
        "write": [1, 0],  # + inverse deformation iy_, for the rmask volumes (start_seg.m: [0 0])
        "mrf": 0, "cleanup": 0, "fwhm": 0,
    }


def test_segment_batch_with_t2_and_analyze_input():
    channels = segment_batch("/data/t1.hdr", "/data/t2.nii", tpm="/t.nii")["spm"]["spatial"]["preproc"]["channel"]
    assert [c["vols"] for c in channels] == [["/data/t1.img,1"], ["/data/t2.nii,1"]]  # .hdr -> .img
    assert channels[1]["biasreg"] == 0.001 and channels[1]["biasfwhm"] == 60


def test_coregister_batch_is_realign_t2_m():
    job = coregister_batch("/d/t1.nii", "/d/t1_T2.nii")["spm"]["spatial"]["coreg"]["estwrite"]
    assert job["ref"] == ["/d/t1.nii,1"] and job["source"] == ["/d/t1_T2.nii,1"] and job["other"] == [""]
    assert job["eoptions"]["cost_fun"] == "nmi" and job["eoptions"]["sep"] == [4, 2]
    assert job["eoptions"]["fwhm"] == [7, 7] and len(job["eoptions"]["tol"]) == 12
    assert job["roptions"] == {"interp": 7, "wrap": [0, 0, 0], "mask": 0, "prefix": "r"}


def test_reslice_batch_matches_resamp_to_one_mm():
    job = reslice_batch("/d/grid.nii", "/d/t1.nii")["spm"]["spatial"]["coreg"]["write"]
    assert job["roptions"]["interp"] == 7 and job["roptions"]["mask"] == 0  # spm_reslice defaults otherwise


def test_batches_convert_to_spm_python_objects():
    pytest.importorskip("spm")
    batch = to_matlabbatch([segment_batch("/d/t1.nii", "/d/t2.nii", tpm="/t.nii"), coregister_batch("/a", "/b")])
    entries = batch._as_runtime()
    assert len(entries) == 2
    preproc = entries[0]["spm"]["spatial"]["preproc"]
    # channel and tissue must be struct arrays, not cells of structs: SPM
    # rejects the latter ("job is not a struct")
    for field in ("channel", "tissue"):
        assert preproc[field]["type__"] == "structarray", field
    assert [c["vols"] for c in preproc["channel"]["data__"]] == [["/d/t1.nii,1"], ["/d/t2.nii,1"]]
    assert [int(t["ngaus"]) for t in preproc["tissue"]["data__"]] == [1, 1, 2, 3, 4, 2]
    assert preproc["warp"]["reg"].tolist() == [0, 0.001, 0.5, 0.05, 0.2]
    assert preproc["warp"]["affreg"] == "mni"
    assert entries[1]["spm"]["spatial"]["coreg"]["estwrite"]["eoptions"]["cost_fun"] == "nmi"
    # Every number is a double, as in MATLAB: SPM fails on int64 inputs.
    assert preproc["warp"]["reg"].dtype == np.float64
    assert np.asarray(preproc["tissue"]["data__"][0]["ngaus"]).dtype == np.float64
    assert np.asarray(entries[1]["spm"]["spatial"]["coreg"]["estwrite"]["roptions"]["interp"]).dtype == np.float64


# --------------------------------------------------------------------------
# _seg8.mat, from a real SPM run (the New York head's)
# --------------------------------------------------------------------------


def test_real_seg8_registration_matches_the_images(tmp_path):
    from roast_py.registration import registration_from_seg8

    reg = registration_from_seg8(SEG8)
    with zipfile.ZipFile(NYHEAD_MASKS_ZIP) as z:
        z.extractall(tmp_path)
    masks = nib.load(str(tmp_path / "nyhead_T1orT2_SPM_masks.nii"))
    # SPM's 1-based .mat matrices, shifted to nibabel's 0-based affines:
    assert np.allclose(reg.image_affine, masks.affine)
    assert np.allclose(reg.tpm_affine, nib.load(str(REPO / "eTPM.nii")).affine, atol=1e-6)


@pytest.mark.slow
def test_landmarks_from_real_seg8_land_on_the_scalp(tmp_path):
    from scipy.ndimage import distance_transform_edt

    from roast_py.geometry.landmarks import tpm_landmarks_to_subject
    from roast_py.registration import registration_from_seg8

    reg = registration_from_seg8(SEG8)
    with zipfile.ZipFile(NYHEAD_MASKS_ZIP) as z:
        z.extractall(tmp_path)
    labels = np.asarray(nib.load(str(tmp_path / "nyhead_T1orT2_SPM_masks.nii")).dataobj)
    head = labels > 0
    outside, inside = distance_transform_edt(~head), distance_transform_edt(head)
    for point in tpm_landmarks_to_subject(reg.tpm2mri)[:4]:  # nasion, inion, ears
        p = tuple(point)
        distance_mm = 0.5 * (outside[p] if not head[p] else inside[p])  # 0.5 mm voxels
        assert distance_mm < 3


# --------------------------------------------------------------------------
# resampling grid, T2 space check
# --------------------------------------------------------------------------


def test_one_mm_grid_covers_the_image_like_spm_get_bbox():
    img = nib.load(str(SUBJECT1))  # 1 x 0.988 x 0.988 mm
    affine, shape = one_mm_grid(img)
    assert np.allclose(np.diag(affine)[:3], 1)
    corners = np.array([[0, 0, 0], np.array(img.shape) - 1])
    world = corners @ img.affine[:3, :3].T + img.affine[:3, 3]
    assert np.allclose(affine[:3, 3], world.min(axis=0))
    extent = world.max(axis=0) - world.min(axis=0)
    assert shape == tuple(int(v) for v in np.ceil(extent + 1 - 0.1))
    assert shape[0] == img.shape[0] and shape[1] < img.shape[1]  # 256 x 0.988 mm -> 253 voxels


def test_same_space_compares_like_realign_t2_m():
    t1 = nib.load(str(SUBJECT1))
    assert _same_space(t1, nib.load(str(SUBJECT1)))
    assert not _same_space(t1, nib.load(str(MNI152)))


def test_spm_output_names_are_roast_ms():
    assert spm_output_name("/d/subject1_ras.nii", with_t2=False) == "subject1_ras_T1orT2"
    assert spm_output_name("/d/subject1.nii", with_t2=True) == "subject1_T1andT2"


# --------------------------------------------------------------------------
# rmask volumes (spm_load_mask.m + spm_sample_mask.m)
# --------------------------------------------------------------------------


def test_atlas_mask_sampling_is_matlabs_quadratic_bspline_of_the_log_mask():
    """At voxel centres, a degree-2 B-spline with the log mask as its
    coefficients is a separable [1 6 1]/8 blur of the log mask (SPM's
    mirror boundary at the edges)."""
    from scipy.ndimage import correlate1d

    rng = np.random.default_rng(0)
    mask = (rng.random((12, 13, 14)) > 0.5).astype(float)
    grid = np.stack(np.meshgrid(*[np.arange(n, dtype=float) for n in mask.shape], indexing="ij"), axis=-1)
    expected = np.log(mask + 1e-5)
    for axis in range(3):
        expected = correlate1d(expected, [1 / 8, 6 / 8, 1 / 8], axis=axis, mode="mirror")
    assert np.allclose(sample_atlas_mask(mask, grid), np.exp(expected))

    single = np.zeros((7, 7, 7))
    single[3, 3, 3] = 1
    sampled = sample_atlas_mask(single, grid[:7, :7, :7])
    # the voxel itself: exp((1 - (6/8)^3) * log(1e-5)) ~ 1.3e-3, above segTouchup's 1e-4 ...
    assert sampled[3, 3, 3] == pytest.approx(np.exp((1 - (6 / 8) ** 3) * np.log(1e-5)), rel=1e-3)
    # ... a face neighbour, with 1/8 * (6/8)^2 of its weight on it, ~2e-5: below
    assert 1e-5 < sampled[2, 3, 3] < 1e-4


def test_atlas_mask_sampling_is_zero_outside_the_atlas():
    mask = np.ones((5, 5, 5))
    coords = np.array([[[[-1.0, 2, 2], [2, 2, 2], [2, 2, 4.5]]]])
    assert sample_atlas_mask(mask, coords).ravel().tolist()[0] == 0
    assert sample_atlas_mask(mask, coords).ravel()[1] == pytest.approx(1 + 1e-5)
    assert sample_atlas_mask(mask, coords).ravel()[2] == 0


def test_rmask_volumes_through_an_identity_deformation(tmp_path):
    """A deformation that maps every eTPM voxel to itself gives back the
    atlas masks themselves."""
    from roast_py.spm.routines import RMASK_FILES, rmask_volumes

    tpm = nib.load(str(REPO / "eTPM.nii"))
    shape = tpm.shape[:3]
    voxels = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1).astype(float)
    mm = voxels @ tpm.affine[:3, :3].T + tpm.affine[:3, 3]
    y = tmp_path / "iy_identity.nii"
    nib.save(nib.Nifti1Image(mm[:, :, :, None, :].astype(np.float32), tpm.affine), str(y))
    volumes = rmask_volumes(y)
    for name, relative in RMASK_FILES.items():
        atlas = np.asarray(nib.load(str(REPO / relative)).dataobj, dtype=float).reshape(shape)
        inner = (slice(2, -2),) * 3
        agree = ((volumes[name] > 1e-4) == (atlas > 1e-4))[inner].mean()
        assert agree > 0.97, name


# --------------------------------------------------------------------------
# segTouchup's patching passes
# --------------------------------------------------------------------------


def test_patching_passes():
    from roast_py.segmentation.touchup import BONE, CSF, GRAY, SKIN, WHITE, patch_labels

    labels = np.full((9, 9, 9), SKIN, dtype=np.uint8)
    labels[2:7, 2:7, 2:7] = BONE
    labels[3:6, 3:6, 3:6] = CSF
    labels[4, 4, 4] = WHITE
    labels[4, 4, 3] = BONE  # bone touching white matter
    zeros = np.zeros(labels.shape)
    out = patch_labels(labels, zeros, zeros, zeros)
    assert out[4, 4, 4] == GRAY  # white touching bone -> gray
    assert out[4, 4, 3] == CSF  # ...and that bone, now touching gray matter, -> CSF
    assert out[2, 4, 4] == BONE  # skin touching CSF -> bone: none left untouched
    # The marked white-matter surface keeps white matter white.
    keep = zeros.copy()
    keep[4, 4, 4] = 1
    assert patch_labels(labels, zeros, zeros, keep)[4, 4, 4] == WHITE
    # Skull holes turn bone into skin.
    holes = zeros.copy()
    holes[2, 2, 2] = 1
    assert patch_labels(labels, holes, zeros, zeros)[2, 2, 2] == SKIN


def test_box_morphology_matches_matlab_for_even_sizes():
    from roast_py.segmentation.touchup import _box_dilate, _box_erode

    point = np.zeros((20, 20, 20), dtype=bool)
    point[10, 10, 10] = True
    dilated = np.argwhere(_box_dilate(point, 10))
    # MATLAB: centre floor((10+1)/2) = 5 -> offsets -4..+5 around the point
    assert dilated.min(axis=0).tolist() == [6, 6, 6] and dilated.max(axis=0).tolist() == [15, 15, 15]
    hole = ~point
    eroded_out = np.argwhere(~_box_erode(hole, 10))
    assert eroded_out.min(axis=0).tolist() == [5, 5, 5] and eroded_out.max(axis=0).tolist() == [14, 14, 14]
    assert _box_erode(np.ones((4, 4, 4), bool), 12).all()  # outside counts as foreground


def test_integer_smoothing_rounds_like_imfilter():
    from roast_py.segmentation.touchup import _smooth_slices

    volume = np.zeros((7, 7, 1), dtype=np.uint8)
    volume[3, 3, 0] = 255
    out = _smooth_slices(volume, 1.0)
    assert out.dtype == np.uint8
    as_float = _smooth_slices(volume.astype(float), 1.0)
    assert np.array_equal(out, np.floor(as_float + 0.5).astype(np.uint8))


# --------------------------------------------------------------------------
# runtime: location, license, installer
# --------------------------------------------------------------------------


def test_default_prefix_is_per_user_or_matlab_runtime_path(monkeypatch, tmp_path):
    monkeypatch.delenv(runtime.PREFIX_ENV_VAR, raising=False)
    assert str(runtime.default_prefix()).startswith(str(Path.home())) or "LOCALAPPDATA" in os.environ
    monkeypatch.setenv(runtime.PREFIX_ENV_VAR, str(tmp_path))
    assert runtime.default_prefix() == tmp_path


def test_pins_are_consistent():
    assert runtime.SPM_PACKAGES["spm-runtime-r2025b"].startswith("25.1.2")
    assert runtime.MATLAB_RELEASE == "R2025b"
    assert all("==" in spec for spec in runtime.install_command()[4:])


def _no_runtime(monkeypatch, tmp_path):
    monkeypatch.setenv(runtime.PREFIX_ENV_VAR, str(tmp_path / "rt"))
    monkeypatch.delenv(runtime.LICENSE_ENV_VAR, raising=False)
    monkeypatch.setattr(runtime, "find_runtime", lambda prefix=None: None)
    monkeypatch.setattr(runtime, "ensure_packages", lambda install_missing=True: None)


def test_license_is_never_accepted_silently(monkeypatch, tmp_path):
    _no_runtime(monkeypatch, tmp_path)
    monkeypatch.setattr(sys, "stdin", None)
    with pytest.raises(runtime.MatlabRuntimeLicenseNotAccepted, match=runtime.LICENSE_ENV_VAR):
        runtime.ensure_spm()
    with pytest.raises(runtime.MatlabRuntimeLicenseNotAccepted, match="declined"):
        runtime.ensure_spm(accept_license=False)


class _FakeTTY:
    def isatty(self):
        return True


def test_interactive_prompt_must_say_yes(monkeypatch, tmp_path):
    _no_runtime(monkeypatch, tmp_path)
    monkeypatch.setattr(sys, "stdin", _FakeTTY())
    monkeypatch.setattr("builtins.input", lambda prompt: "no")
    with pytest.raises(runtime.MatlabRuntimeLicenseNotAccepted, match="not accepted"):
        runtime.ensure_spm()


@pytest.mark.parametrize("how", ["argument", "env", "prompt"])
def test_accepted_license_installs_into_the_per_user_prefix(monkeypatch, tmp_path, how):
    _no_runtime(monkeypatch, tmp_path)
    calls, found = [], {}

    def fake_install(version, prefix=None, auto_answer=False):
        calls.append((version, prefix, auto_answer))
        Path(prefix, version).mkdir(parents=True)
        found["path"] = Path(prefix, version)

    import matlab_runtime

    monkeypatch.setattr(matlab_runtime, "install", fake_install)
    monkeypatch.setattr(runtime, "find_runtime", lambda prefix=None: found.get("path"))
    kwargs = {}
    if how == "argument":
        kwargs["accept_license"] = True
    elif how == "env":
        monkeypatch.setenv(runtime.LICENSE_ENV_VAR, "yes")
    else:
        monkeypatch.setattr(sys, "stdin", _FakeTTY())
        monkeypatch.setattr("builtins.input", lambda prompt: "yes")

    path = runtime.ensure_spm(**kwargs)
    assert calls == [("R2025b", str((tmp_path / "rt").resolve()), True)]
    assert path == (tmp_path / "rt").resolve() / "R2025b"
    assert not (tmp_path / "rt" / ".download").exists()


def test_a_failed_download_explains_what_to_do(monkeypatch, tmp_path):
    import urllib.error

    import matlab_runtime

    _no_runtime(monkeypatch, tmp_path)

    def unreachable(*args, **kwargs):
        raise urllib.error.URLError("Tunnel connection failed: 403 Forbidden")

    monkeypatch.setattr(matlab_runtime, "install", unreachable)
    with pytest.raises(RuntimeError, match=r"(?s)ssd\.mathworks\.com.*MATLAB_RUNTIME_PATH"):
        runtime.ensure_spm(accept_license=True)
    assert not (tmp_path / "rt" / ".download").exists()


def test_an_existing_runtime_is_used_without_asking(monkeypatch, tmp_path):
    existing = tmp_path / "MATLAB" / "R2025b"
    monkeypatch.setattr(runtime, "ensure_packages", lambda install_missing=True: None)
    monkeypatch.setattr(runtime, "find_runtime", lambda prefix=None: existing)
    monkeypatch.setattr(runtime, "install_runtime", lambda **k: pytest.fail("installed"))
    assert runtime.ensure_spm() == existing
    env = runtime.worker_env(existing)
    assert env[runtime.PREFIX_ENV_VAR] == str(tmp_path / "MATLAB")


def test_install_can_be_switched_off(monkeypatch, tmp_path):
    _no_runtime(monkeypatch, tmp_path)
    with pytest.raises(RuntimeError, match="automatic installation is off"):
        runtime.ensure_spm(install_missing=False)


def test_worker_reports_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(runtime.subprocess, "run", lambda cmd, env: subprocess.CompletedProcess(cmd, 3))
    with pytest.raises(RuntimeError, match="worker exit status 3"):
        runtime.run_job({"check": True}, runtime=tmp_path / "R2025b")


def test_worker_command_uses_mwpython_on_macos(monkeypatch):
    monkeypatch.setattr(runtime.platform, "system", lambda: "Darwin")
    cmd = runtime.worker_command("/j.json")
    assert "mwpython2" in cmd[2] and cmd[3:7] == ["-variant", "R2025b", "-m", "roast_py.spm.worker"]
    monkeypatch.setattr(runtime.platform, "system", lambda: "Linux")
    assert runtime.worker_command("/j.json") == [sys.executable, "-m", "roast_py.spm.worker", "/j.json"]


# --------------------------------------------------------------------------
# roast(): SPM options ask for the runtime before any work
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs",
    [
        {"multiaxial": "off"},
        {"t2": str(SUBJECT1).replace("subject1", "subject1_T2")},
        {"resampling": "on"},
    ],
)
def test_roast_sets_up_spm_first(monkeypatch, tmp_path, kwargs):
    import roast_py.segmentation.multiaxial as multiaxial

    _no_runtime(monkeypatch, tmp_path)
    monkeypatch.setattr(sys, "stdin", None)
    monkeypatch.setattr(multiaxial, "segment", lambda *a, **k: pytest.fail("started segmenting"))
    from roast_py import roast

    with pytest.raises(runtime.MatlabRuntimeLicenseNotAccepted):
        roast(str(SUBJECT1), **kwargs, install_missing=True, visualize=False)
