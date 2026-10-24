"""ROAST's SPM steps, run through spm-python (see runtime.py):

* resample_to_1mm()   -- resampToOneMM.m (the 'resampling' option)
* realign_t2()        -- realignT2.m (the 'T2' option)
* spm_segment()       -- start_seg.m + renameSPMres.m + segTouchup.m

Output names follow ROAST's MATLAB code, so a folder processed by either
looks the same: `<mri>_1mm.nii`, `<T1>_T2_aligned.nii`,
`<mri>_T1orT2_seg8.mat` (`_T1andT2` with a T2) and
`<mri>_T1orT2_SPM_masks.nii`. As in MATLAB, a step whose output already
exists is skipped.

ROAST's modified SPM12 also wrote `_rmask.mat` (skull holes, extra-cerebral
CSF/eyes and a marked white-matter surface, warped from atlas space) for
segTouchup's patching passes. Official SPM doesn't, so spm_segment() has
SPM write its forward deformation and samples the same atlas masks through
it the same way (rmask_volumes()).
"""

from __future__ import annotations

import os
import shutil
import warnings
from pathlib import Path

import nibabel as nib
import numpy as np

from .batches import coregister_batch, reslice_batch, segment_batch
from .runtime import run_job


def _split(path: str | os.PathLike) -> tuple[Path, str, str]:
    path = Path(path).resolve()
    name = path.name
    for ext in (".nii.gz", ".nii", ".hdr", ".img"):
        if name.endswith(ext):
            return path.parent, name[: -len(ext)], ext
    return path.parent, path.stem, path.suffix


def _nii_out(directory: Path, name: str) -> Path:
    """SPM writes NIfTI-1 .nii for any input format."""
    return directory / f"{name}.nii"


def _find_bundled(relative: str) -> Path:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / relative
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate {relative} (part of the ROAST repository).")


def find_tpm() -> Path:
    return _find_bundled("eTPM.nii")


# --------------------------------------------------------------------------
# resampling to 1 mm (resampToOneMM.m)
# --------------------------------------------------------------------------


def one_mm_grid(img: nib.Nifti1Image) -> tuple[np.ndarray, tuple[int, int, int]]:
    """The 1 mm grid resampToOneMM.m reslices onto: (0-based affine, shape).

    MATLAB: bb = spm_get_bbox(V) -- the world bounding box of the 8 corner
    voxel centers -- then VV.mat = spm_matrix([bb(1,:) 0 0 0 1 1 1]) *
    spm_matrix([-1 -1 -1]) and VV.dim = ceil(VV.mat \\ [bb(2,:) 1]' - 0.1).
    """
    shape = np.array(img.shape[:3])
    corners = np.array([[i, j, k] for i in (0, shape[0] - 1) for j in (0, shape[1] - 1) for k in (0, shape[2] - 1)])
    world = corners @ img.affine[:3, :3].T + img.affine[:3, 3]
    bb_min, bb_max = world.min(axis=0), world.max(axis=0)
    affine = np.eye(4)
    affine[:3, 3] = bb_min  # 0-based voxel 0 sits at bb_min (MATLAB's voxel 1)
    dim = np.ceil(bb_max - bb_min + 1 - 0.1).astype(int)
    return affine, tuple(int(d) for d in dim)


def resample_to_1mm(mri: str | os.PathLike, runtime=None) -> str:
    """Ports resampToOneMM.m with resampling on: reslices `mri` to 1 mm
    isotropic voxels (7th-degree B-spline) as `<mri>_1mm.nii`, and returns
    that path. An MRI that already has 1 mm voxels is returned unchanged.
    """
    directory, base, _ = _split(mri)
    img = nib.load(str(mri))
    if np.all(np.diag(img.affine)[:3] == 1):
        warnings.warn(f"The MRI {mri} already has a 1 mm isotropic resolution. No need to resample.", stacklevel=2)
        return str(mri)
    out = _nii_out(directory, f"{base}_1mm")
    if out.exists():
        warnings.warn(
            f"{mri} has already been resampled to 1mm resolution and saved as {out}. ROAST will use that file "
            "as the input.",
            stacklevel=2,
        )
        return str(out)

    print(f"Resampling {mri} to 1 mm isotropic resolution...")
    affine, shape = one_mm_grid(img)
    grid = directory / f"{base}_1mm_grid.nii"
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.uint8), affine), str(grid))
    try:
        run_job({"batches": [reslice_batch(str(grid), str(mri))]}, runtime=runtime)
        _nii_out(directory, f"r{base}").replace(out)
    finally:
        grid.unlink(missing_ok=True)
    print(f"{mri} has been resampled to 1 mm isotropic resolution, and is saved as:\n{out}")
    return str(out)


# --------------------------------------------------------------------------
# T2 alignment (realignT2.m)
# --------------------------------------------------------------------------


def _same_space(a: nib.Nifti1Image, b: nib.Nifti1Image) -> bool:
    """realignT2.m's test: identical dim(1:5), pixdim(1:4) and sform rows."""
    ha, hb = a.header, b.header
    return (
        np.array_equal(ha["dim"][:5], hb["dim"][:5])
        and np.array_equal(ha["pixdim"][:4], hb["pixdim"][:4])
        and all(np.array_equal(ha[f"srow_{r}"], hb[f"srow_{r}"]) for r in "xyz")
    )


def realign_t2(t2: str | os.PathLike, t1: str | os.PathLike, runtime=None) -> str:
    """Ports realignT2.m: if the T2 isn't in the T1's voxel space,
    coregisters and reslices it onto the T1 (SPM Coregister: Estimate &
    Reslice, NMI, 7th-degree B-spline) as `<T1 name>_T2_aligned.nii` next
    to the T2. Returns the T2 to use."""
    t1_img, t2_img = nib.load(str(t1)), nib.load(str(t2))
    if _same_space(t1_img, t2_img):
        return str(t2)

    _, t1_name, _ = _split(t1)
    t2_dir, _, t2_ext = _split(t2)
    aligned = _nii_out(t2_dir, f"{t1_name}_T2_aligned")
    if aligned.exists() and _same_space(t1_img, nib.load(str(aligned))):
        warnings.warn(
            f"Original T2 {t2} is not aligned to T1 image space, but an aligned T2 has been found ({aligned}). "
            "ROAST will use that file as the input.",
            stacklevel=2,
        )
        return str(aligned)

    warnings.warn(f"T2 image {t2} is not aligned to T1 image space. ROAST will align it now.", stacklevel=2)
    print("Aligning T2 to T1...")
    # Coregistration rewrites the source's header, so work on a copy.
    source = t2_dir / f"{t1_name}_T2{t2_ext}"
    shutil.copyfile(t2, source)
    if t2_ext == ".hdr":
        shutil.copyfile(Path(t2).with_suffix(".img"), source.with_suffix(".img"))
    run_job({"batches": [coregister_batch(str(t1), str(source))]}, runtime=runtime)
    _nii_out(t2_dir, f"r{t1_name}_T2").replace(aligned)
    print(f"Aligned T2 is saved as {aligned}, and will be used by ROAST.")
    return str(aligned)


# --------------------------------------------------------------------------
# segmentation (start_seg.m + renameSPMres.m + segTouchup.m)
# --------------------------------------------------------------------------

# ROAST's extra atlas masks (lib/spm12), in the eTPM.nii grid, and the
# names segTouchup.m knows them by.
RMASK_FILES = {
    "holes_vol": "lib/spm12/eTPM_Skull Holes.nii",
    "eyes_vol": "lib/spm12/eTPM_ExtraCerebral_CSF.nii",
    "WMexclude_vol": "lib/spm12/eTPM_Marked_WM_Surface.nii",
}


def sample_atlas_mask(mask: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Ports spm_load_mask.m + spm_sample_mask.m (ROAST's additions to SPM12).

    `mask` is an atlas volume; `coords` (..., 3) are 0-based voxel
    coordinates in it. As in MATLAB: the mask is clipped to [0, 1] and
    log-transformed (log(x + 1e-5)); its degree-1 B-spline coefficients
    (the samples themselves) are then evaluated as a degree-2 B-spline
    (`mask.deg = deg + 1` -- the original notes the mismatch too); the
    result is exponentiated, and points outside the atlas are 0.
    """
    from scipy.ndimage import map_coordinates

    coefficients = np.log(np.clip(np.asarray(mask, dtype=float), 0, 1) + 1e-5)
    shape = np.array(coefficients.shape)
    flat = coords.reshape(-1, 3)
    inside = np.all((flat >= 0) & (flat <= shape - 1), axis=1)
    out = np.zeros(flat.shape[0])
    values = map_coordinates(coefficients, flat[inside].T, order=2, prefilter=False, mode="mirror")
    out[inside] = np.exp(values)
    return out.reshape(coords.shape[:-1])


def rmask_volumes(deformation: str | os.PathLike) -> dict[str, np.ndarray]:
    """The three `_rmask.mat` volumes ROAST's modified SPM12 saved, from
    official SPM's forward deformation `y_<name>.nii`.

    That field holds, for every voxel of the segmented MRI, its position in
    the atlas (mm) -- the same mapping spm_preproc_write8.m computes
    internally (`[t1,t2,t3] = defs(...)`) and sampled the masks with.
    """
    y = nib.load(str(deformation))
    mm = np.asarray(y.dataobj, dtype=float).reshape(*y.shape[:3], 3)
    volumes = {}
    for name, relative in RMASK_FILES.items():
        atlas = nib.load(str(_find_bundled(relative)))
        to_voxel = np.linalg.inv(atlas.affine)
        coords = mm @ to_voxel[:3, :3].T + to_voxel[:3, 3]
        volumes[name] = sample_atlas_mask(np.asarray(atlas.dataobj, dtype=float).reshape(atlas.shape[:3]), coords)
    return volumes


def spm_output_name(t1: str | os.PathLike, with_t2: bool) -> str:
    """`<mri>_T1orT2` or `<mri>_T1andT2`: roast.m's subjRasRSPDspm name."""
    _, base, _ = _split(t1)
    return f"{base}_T1andT2" if with_t2 else f"{base}_T1orT2"


def seg8_path(t1: str | os.PathLike, with_t2: bool) -> Path:
    directory, _, _ = _split(t1)
    return directory / f"{spm_output_name(t1, with_t2)}_seg8.mat"


def spm_segment(t1: str | os.PathLike, t2: str | os.PathLike | None = None, runtime=None) -> str:
    """ROAST's SPM segmentation: returns the path of the 6-tissue mask
    `<mri>_T1orT2_SPM_masks.nii` (`_T1andT2` with a T2), writing
    `<...>_seg8.mat` (the registration roast() reads landmarks from) too.

    `t2` must already be in the T1's space (see realign_t2()).
    """
    from ..segmentation.touchup import touchup

    directory, base, _ = _split(t1)
    target = spm_output_name(t1, t2 is not None)
    masks_path = directory / f"{target}_SPM_masks.nii"
    seg8 = directory / f"{target}_seg8.mat"
    deformation = directory / f"y_{target}.nii"
    tissue_paths = [directory / f"c{i}{target}.nii" for i in range(1, 7)]

    if masks_path.exists() and seg8.exists():
        print(f"MRI segmented by SPM; using {masks_path}")
        return str(masks_path)

    if not (seg8.exists() and deformation.exists() and all(p.exists() for p in tissue_paths)):
        print(f"Segmenting {t1} with SPM" + (f" (with T2 {t2})" if t2 else "") + " ...")
        run_job({"batches": [segment_batch(str(t1), None if t2 is None else str(t2), tpm=str(find_tpm()))]}, runtime=runtime)
        # renameSPMres.m: c1<mri>.nii ... -> c1<mri>_T1orT2.nii ...
        for i, path in enumerate(tissue_paths, start=1):
            (directory / f"c{i}{base}.nii").replace(path)
        (directory / f"{base}_seg8.mat").replace(seg8)
        (directory / f"y_{base}.nii").replace(deformation)

    print("SPM segmentation touchup ...")
    # load_untouch_nii reads SPM's uint8 maps raw (0-255, not rescaled),
    # and segTouchup.m's smoothing works on those integers.
    tissues = [np.asanyarray(nib.load(str(p)).dataobj.get_unscaled()) for p in tissue_paths]
    gray, white, csf, bone, skin, air = tissues
    patch = rmask_volumes(deformation)
    labels = touchup(
        gray, white, csf, bone, skin, air,
        holes=patch["holes_vol"], eyes=patch["eyes_vol"], wm_exclude=patch["WMexclude_vol"],
    )
    reference = nib.load(str(t1))
    out = nib.Nifti1Image(labels.astype(np.uint8), reference.affine)
    out.header["descrip"] = b"tissue masks"
    nib.save(out, str(masks_path))
    # segTouchup.m deletes the per-tissue maps (and _rmask.mat) once done.
    for path in [*tissue_paths, deformation]:
        path.unlink(missing_ok=True)
    return str(masks_path)
