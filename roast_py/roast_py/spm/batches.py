"""SPM batch jobs, as ROAST's MATLAB code writes them.

Each builder returns one `matlabbatch{1}` entry as plain Python data, so
jobs can be written to JSON, sent to the worker process and checked in
tests without SPM installed:

* dict            -> MATLAB struct
* list of dicts   -> struct array (`tissue(1)`, `tissue(2)`, ...)
* list of str     -> cell array of strings (`{'file.nii,1'}`)
* list of numbers -> numeric row vector (`[1 0]`)
* str / number    -> char / scalar

to_matlabbatch() turns a list of them into spm-python's objects.

One deliberate difference from the MATLAB code: ROAST ships a modified
SPM12 whose `prefix` option appends instead of prepends (`_1mm`,
`_aligned`), so its outputs are `subject1_1mm.nii`. Official SPM prepends,
so these jobs use SPM's default prefix 'r' and roast_py.spm.routines renames
the outputs to ROAST's names.
"""

from __future__ import annotations

import os

# start_seg.m's Gaussians per tissue class (gray, white, CSF, bone, skin, air).
TISSUE_NGAUS = (1, 1, 2, 3, 4, 2)


def _vol(path: str) -> str:
    """SPM's `file,1` volume reference. Analyze pairs are given by their
    .img, as start_seg.m does."""
    path = os.path.abspath(path)
    if path.endswith(".hdr"):
        path = path[: -len(".hdr")] + ".img"
    return f"{path},1"


def segment_batch(t1: str, t2: str | None = None, tpm: str = "", write_deformation: bool = True) -> dict:
    """Ports start_seg.m (native space only, ROAST's default `norm=false`).

    `write_deformation` additionally writes SPM's inverse deformation
    `iy_<t1>.nii` (warp.write = [1 0]; start_seg.m has [0 0]): on the
    MRI's own grid, the atlas position (mm) of every voxel. It doesn't
    change the segmentation. roast_py samples ROAST's extra atlas masks
    through it, which ROAST's modified SPM12 did inside spm_preproc_write8.m
    (see routines.rmask_volumes).
    """
    channel = {"vols": [_vol(t1)], "biasreg": 0.001, "biasfwhm": 60, "write": [0, 0]}
    channels = [channel]
    if t2 is not None:
        channels.append({"vols": [_vol(t2)], "biasreg": 0.001, "biasfwhm": 60, "write": [0, 0]})
    tissues = [
        {"tpm": [f"{os.path.abspath(tpm)},{i}"], "ngaus": ngaus, "native": [1, 0], "warped": [0, 0]}
        for i, ngaus in enumerate(TISSUE_NGAUS, start=1)
    ]
    warp = {
        "reg": [0, 0.001, 0.5, 0.05, 0.2],
        "affreg": "mni",
        "samp": 3,
        "write": [1, 0] if write_deformation else [0, 0],
        "mrf": 0,
        "cleanup": 0,
        "fwhm": 0,
    }
    return {"spm": {"spatial": {"preproc": {"channel": channels, "tissue": tissues, "warp": warp}}}}


def coregister_batch(ref: str, source: str) -> dict:
    """Ports realignT2.m's coregister-and-reslice (estimate + write)."""
    return {
        "spm": {
            "spatial": {
                "coreg": {
                    "estwrite": {
                        "ref": [_vol(ref)],
                        "source": [_vol(source)],
                        "other": [""],
                        "eoptions": {
                            "cost_fun": "nmi",
                            "sep": [4, 2],
                            "tol": [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001],
                            "fwhm": [7, 7],
                        },
                        "roptions": {"interp": 7, "wrap": [0, 0, 0], "mask": 0, "prefix": "r"},
                    }
                }
            }
        }
    }


def reslice_batch(ref: str, source: str) -> dict:
    """Reslices `source` onto `ref`'s grid. With the 1 mm reference grid
    from routines.one_mm_grid(), this is resampToOneMM.m's
    `spm_reslice(VV, struct('mean',false,'which',1,'interp',7))` --
    Coregister: Reslice runs the same spm_reslice call."""
    return {
        "spm": {
            "spatial": {
                "coreg": {
                    "write": {
                        "ref": [_vol(ref)],
                        "source": [_vol(source)],
                        "roptions": {"interp": 7, "wrap": [0, 0, 0], "mask": 0, "prefix": "r"},
                    }
                }
            }
        }
    }


def to_matlabbatch(entries: list[dict]):
    """spm-python's `matlabbatch` (an spm.Cell) for a list of batch entries.

    Builds it attribute by attribute -- `batch(0).spm.spatial.preproc
    .tissue[2].ngaus = 2`.
    """
    import numpy as np
    from spm import Array, Cell

    def is_struct_array(value):
        return isinstance(value, list) and value and all(isinstance(v, dict) for v in value)

    def is_number(v):
        return isinstance(v, (int, float)) and not isinstance(v, bool)

    def convert(value):
        # Numbers go over as doubles, like MATLAB literals: a Python int
        # would arrive as int64, and SPM fails mixing that with doubles
        # (spm_reslice: "Integers can only be combined with integers of
        # the same class, or scalar doubles").
        if isinstance(value, list) and value and all(is_number(v) for v in value):
            return Array.from_any(np.asarray(value, dtype=float))
        if is_number(value):
            return float(value)
        if isinstance(value, list):  # cell array of strings
            return list(value)
        return value

    def leaves(spec: dict, path=()):
        """(path, value) for every leaf; path items are field names (str)
        or struct-array indices (int)."""
        for key, value in spec.items():
            if isinstance(value, dict):
                yield from leaves(value, (*path, key))
            elif is_struct_array(value):
                for i, element in enumerate(value):
                    yield from leaves(element, (*path, key, i))
            else:
                yield (*path, key), value

    batch = Cell()
    for i, entry in enumerate(entries):
        # Each assignment walks the whole path from the batch again --
        # `batch(i).spm.spatial.preproc.tissue[2].ngaus = 2` -- because
        # spm-python's placeholder objects can't be kept and reused.
        # matlabbatch itself is a cell (`batch(i)`, MATLAB's `{i}`); struct
        # arrays inside it use brackets (`tissue[2]`, MATLAB's `(i)`): in
        # spm-python, `tissue(2).f` would make a *cell* of structs, which
        # SPM rejects ("Can not initialise tissue value(s): job is not a
        # struct").
        for path, value in leaves(entry):
            target = batch(i)
            for step in path[:-1]:
                target = target[step] if isinstance(step, int) else getattr(target, step)
            setattr(target, path[-1], convert(value))
    return batch
