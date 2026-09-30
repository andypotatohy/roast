"""Ports runNiftyReg.m: affine registration of the head to MNI space.

Runs NiftyReg's `reg_aladin` (bundled under lib/NiftyReg/) to register the
subject's MRI to the MNI152 template (example/MNI152_T1_1mm.nii), and
stores the result the way SPM's `_seg8.mat` does, which is what the rest of
ROAST expects:

* ``affine`` -- SPM's ``Affine``: subject world (mm) -> MNI world (mm).
  reg_aladin's own matrix maps the other way (reference/MNI world ->
  floating/subject world, the direction resampling needs), so it is
  inverted, exactly as runNiftyReg.m does.
* ``image_affine`` / ``tpm_affine`` -- the voxel-to-world matrices of the
  subject MRI and of eTPM.nii (SPM's ``image(1).mat`` / ``tpm(1).mat``).

One convention differs from MATLAB: those voxel-to-world matrices are
nibabel affines, taking **0-based** voxel indices like everything else in
roast_py, where SPM's ``.mat`` takes 1-based ones. Composed matrices such
as ``mri2mni`` and ``tpm2mri`` therefore also work on 0-based indices.
"""

from __future__ import annotations

import json
import os
import platform
import stat
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _find_bundled(relative: str, what: str) -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / relative
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate the bundled {what} ({relative}). Pass its path explicitly.")


def find_reg_aladin(bin_path: str | os.PathLike | None = None) -> Path:
    """The bundled reg_aladin for this platform (runNiftyReg.m's switch on `computer('arch')`)."""
    if bin_path is not None:
        return Path(bin_path)
    system = platform.system()
    relative = {
        "Linux": "lib/NiftyReg/linux/reg_aladin",
        "Darwin": "lib/NiftyReg/mac/reg_aladin",  # x86_64: runs under Rosetta on Apple silicon
        "Windows": "lib/NiftyReg/win/reg_aladin.exe",
    }.get(system)
    if relative is None:
        raise RuntimeError(f"Unsupported operating system: {system!r}")
    binary = _find_bundled(relative, "NiftyReg reg_aladin binary")
    if system != "Windows" and not os.access(binary, os.X_OK):
        # The repo doesn't carry the executable bit for every platform;
        # runNiftyReg.m chmods it too.
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return binary


def find_mni_template() -> Path:
    """example/MNI152_T1_1mm.nii, runNiftyReg.m's registration reference."""
    return _find_bundled("example/MNI152_T1_1mm.nii", "MNI152 template")


def find_tpm() -> Path:
    """eTPM.nii, the tissue probability atlas the TPM landmarks are defined in."""
    return _find_bundled("eTPM.nii", "eTPM.nii tissue probability atlas")


@dataclass(frozen=True)
class Registration:
    """A subject-to-MNI affine registration, in _seg8.mat/_niftyReg.mat terms."""

    affine: np.ndarray  # SPM's Affine: subject world -> MNI world (mm)
    image_affine: np.ndarray  # subject voxel (0-based) -> subject world
    tpm_affine: np.ndarray  # eTPM voxel (0-based) -> MNI world

    @property
    def mri2mni(self) -> np.ndarray:
        """Subject voxel (0-based) -> MNI mm: roast.m's `Affine*image(1).mat`."""
        return self.affine @ self.image_affine

    @property
    def tpm2mri(self) -> np.ndarray:
        """eTPM voxel -> subject voxel (both 0-based):
        roast.m's `inv(image(1).mat)*inv(Affine)*tpm(1).mat`."""
        return np.linalg.inv(self.image_affine) @ np.linalg.inv(self.affine) @ self.tpm_affine

    def save(self, path: str | os.PathLike) -> None:
        with open(path, "w") as f:
            json.dump(
                {
                    "Affine": self.affine.tolist(),
                    "image_affine": self.image_affine.tolist(),
                    "tpm_affine": self.tpm_affine.tolist(),
                    "note": "Affine maps subject world -> MNI world (SPM convention); "
                    "the *_affine matrices are nibabel (0-based voxel) affines.",
                },
                f,
                indent=2,
            )

    @classmethod
    def load(cls, path: str | os.PathLike) -> Registration:
        with open(path) as f:
            data = json.load(f)
        return cls(
            np.asarray(data["Affine"], dtype=float),
            np.asarray(data["image_affine"], dtype=float),
            np.asarray(data["tpm_affine"], dtype=float),
        )


def registration_path(input_mri: str | os.PathLike, out_dir: str | os.PathLike | None = None) -> Path:
    """Where run_niftyreg() saves its result: `<mri>_niftyReg.json`
    (MATLAB: `<mri>_niftyReg.mat`)."""
    input_mri = Path(input_mri)
    base = input_mri.name.removesuffix(".gz").removesuffix(".nii")
    return Path(out_dir or input_mri.parent) / f"{base}_niftyReg.json"


def run_niftyreg(
    input_mri: str | os.PathLike,
    out_dir: str | os.PathLike | None = None,
    reference: str | os.PathLike | None = None,
    tpm: str | os.PathLike | None = None,
    bin_path: str | os.PathLike | None = None,
) -> Registration:
    """Ports runNiftyReg.m: registers `input_mri` to the MNI152 template.

    Runs `reg_aladin -ref MNI152 -flo input_mri -aff ... -res ... -voff`
    (a couple of minutes on CPU), inverts its matrix into SPM's `Affine`
    convention, and saves the result to `<mri>_niftyReg.json` in `out_dir`
    (default: next to the MRI). The resampled image and the raw matrix
    file are temporary, as in MATLAB.
    """
    import nibabel as nib

    input_mri = Path(input_mri).resolve()
    reference = Path(reference) if reference else find_mni_template()
    tpm = Path(tpm) if tpm else find_tpm()
    binary = find_reg_aladin(bin_path)
    out_path = registration_path(input_mri, out_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=out_path.parent, prefix=".niftyreg_") as tmp:
        forward = Path(tmp) / "tmp_forward.txt"
        resampled = Path(tmp) / "resampled.nii"
        cmd = [
            str(binary), "-ref", str(reference), "-flo", str(input_mri),
            "-aff", str(forward), "-res", str(resampled), "-voff",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not forward.exists():
            raise RuntimeError(
                f"niftyReg failed (exit status {result.returncode}).\ncommand: {' '.join(cmd)}\n"
                f"stdout:\n{result.stdout[-2000:]}\nstderr:\n{result.stderr[-2000:]}"
            )
        ref_to_flo = np.loadtxt(forward)

    registration = Registration(
        affine=np.linalg.inv(ref_to_flo),  # "to be consistent with SPM format"
        image_affine=np.asarray(nib.load(str(input_mri)).affine, dtype=float),
        tpm_affine=np.asarray(nib.load(str(tpm)).affine, dtype=float),
    )
    registration.save(out_path)
    return registration
