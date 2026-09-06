# roast_py — Python port of ROAST's core pipeline (work in progress)

This is an in-progress translation of ROAST's core simulation pipeline
(`roast()`, `roast_target()`, `reviewRes()`) from MATLAB to Python, with the
goal of removing the MATLAB license requirement. Scope, phasing, and the
two open technical risks (the CGAL mesher and SPM segmentation are not
plain standalone binaries the way getDP and NiftyReg are) are written up
in full in the project's plan; the short version is below.

**Status: Phase 0 (I/O) started. Everything else is not yet implemented.**
This is not a working `roast()` replacement yet — do not use it for actual
simulations.

## What's here so far

`roast_py/io/nifti.py` ports the pure-numerical MRI preprocessing helpers:

| MATLAB | Python | Notes |
|---|---|---|
| `convertToRAS.m` | `convert_to_ras()` | Uses nibabel's `as_closest_canonical` rather than reimplementing the sform/qform flip math by hand. Verified against `example/MNI152_T1_1mm.nii` (the LAS-oriented file the top-level README calls out) — the resulting affine matches ROAST's translation-update formula exactly. |
| `zeroPadding.m` | `zero_pad()` | Verified world-coordinates are preserved across padding on the same example file. |
| `alignHeader2mni.m` | `align_header_to_mni()` | Ports the `update_affine` helper (forces sform, writes `_MNI` suffix). |

`resampToOneMM.m` and `realignT2.m` are deliberately **not** ported yet —
both call into SPM (`spm_reslice`, `spm_jobman` coreg) to do real
resampling/registration, not just header math, so they're deferred to the
segmentation phase where the SPM dependency gets resolved as a whole
(likely via SPM12's free MATLAB-Runtime-based standalone build, the same
approach other Python neuroimaging pipelines use to call SPM without a
MATLAB license).

Tests in `tests/test_nifti.py` run against the real example MRIs bundled
in the repo's `example/` directory (not synthetic fixtures), so they
double as a numerical-parity check against ROAST's own documented example
data.

## Remaining phases (see the full plan for detail)

0. ~~I/O & preprocessing~~ (in progress)
1. Segmentation — port `lib/multiaxial` (already pure Python/TensorFlow) as
   the default path; SPM12-standalone as an optional secondary path
2. Electrode placement (`fitCap2individual.m`, `electrodePlacement.m`)
3. Meshing — resolve the CGAL-mesher-is-a-MEX-file-not-a-binary risk
4. FEM solve — `.pro` generation + `getdp` subprocess + `.pos` parsing
5. End-to-end numerical validation against MATLAB ROAST (gate before continuing)
6. Targeting (`roast_target()`, CVX → cvxpy)
7. Visualization (`reviewRes()`) & packaging

## Running the tests

```
cd roast_py
pip install -e ".[dev]"
pytest tests/ -v
```
