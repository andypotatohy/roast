# roast_py — Python port of ROAST's core pipeline (work in progress)

This is an in-progress translation of ROAST's core simulation pipeline
(`roast()`, `roast_target()`, `reviewRes()`) from MATLAB to Python, with the
goal of removing the MATLAB license requirement. Scope, phasing, and the
two open technical risks (the CGAL mesher and SPM segmentation are not
plain standalone binaries the way getDP and NiftyReg are) are written up
in full in the project's plan; the short version is below.

**Status: Phases 0-2 (I/O, segmentation, electrode placement) done. Phases
3+ not yet implemented.** This is not a working `roast()` replacement yet —
do not use it for actual simulations.

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

## Segmentation (Phase 1)

Two paths, both producing ROAST's 6-tissue label scheme (1=white, 2=gray,
3=csf, 4=bone, 5=skin, 6=air, 0=background):

- **`roast_py/segmentation/multiaxial.py` (default).** Ports
  `lib/multiaxial/{SEGMENT.py,preprocessing_lib.py,utils.py}` — already
  pure Python/TensorFlow — into an in-process function, replacing MATLAB's
  `runMultiaxial.m` + subprocess-into-a-separate-conda-env with a plain
  function call. **Verified end-to-end** against the real bundled model
  weights and `example/subject1.nii` (see `tests/test_multiaxial.py`,
  `@pytest.mark.slow`, ~2-3 min on CPU) — produces an anatomically sane
  6-tissue segmentation. Needs the legacy Keras 2 runtime to load the
  bundled `.h5` models under Keras 3 (`pip install roast_py[multiaxial]`;
  see `roast_py/segmentation/_keras_compat.py` for why).

  A Python-native alternative to SPM was considered for the *default* path
  too: nothing widely used does ROAST's exact 6-class full-head
  segmentation (ANTsPyNet/FSL only segment brain tissue, not skull/scalp/
  air, which matter for TES modeling; SimNIBS's `charm` is the closest real
  alternative but would mean importing a second large third-party
  simulation codebase). The already-bundled multiaxial CNN turns out to
  already be exactly that Python-native alternative, which is why it's the
  default rather than a fallback.

- **`roast_py/segmentation/spm_standalone.py` (optional/parity path).**
  Drives SPM12's official standalone build (compiled against the free
  MATLAB Runtime, no MATLAB license) to reproduce `start_seg.m`'s New
  Segment step, then ports `segTouchup.m`'s cleanup pipeline
  (`roast_py/segmentation/touchup.py`, unit-tested with synthetic data) to
  turn SPM's raw tissue-probability maps into the same 6-label scheme.
  **Not runtime-tested** — SPM standalone + the MATLAB Runtime aren't
  installable in this environment. **Known gap:** `segTouchup.m`'s three
  final patching passes (gray-matter/CSF/bone) depend on `eyes_vol`,
  `holes_vol`, and `WMexclude_vol` — extra volumes computed by ROAST's own
  patched copy of SPM's `lib/spm12/spm_preproc_write8.m` from extra classes
  in the extended `eTPM.nii` atlas. `touchup.py` implements everything else
  in the pipeline (smoothing, binarization, CSF-continuity fix,
  disconnected-voxel pruning, empty-voxel relabeling) but not those three
  passes yet — tracked as a follow-up task.

## Electrode placement & cap fitting (Phase 2)

`roast_py/geometry/` ports the whole electrode-placement pipeline:

| MATLAB | Python |
|---|---|
| `mask2EdgePointCloud.m`, `project2ClosestSurfacePoints.m`, `map2Points.m`, `convertToRASpointCloud.m` | `geometry/point_cloud.py` |
| `cylinder2P.m`, `drawCylinder.m`, `drawCuboid.m`, `drawLine.m` | `geometry/shapes.py` |
| `lib/ncs2daprox/ncs2dapprox.m` (+ its Merge/MaxSqDist helpers) | `geometry/spline.py` |
| `capInfo.xlsx` reader (`readtable(...)` calls) | `geometry/cap_info.py` |
| `fitCap2individual.m` | `geometry/cap_fitting.py` |
| `cleanScalp.m` | `geometry/scalp.py` |
| `placeNeckElec.m`, `generateElecMask.m`, `placeAndModelElectrodes.m`, `electrodePlacement.m`, `elecPreproc.m`'s classification | `geometry/placement.py` |

One deliberate API change from the MATLAB original: `electrode_placement()`
takes `elec_names`/`elec_paras` in whatever order the caller wants (it
internally reorders into MATLAB's predefined/neck/custom pool-sorted order
the way `elecPreproc.m`'s `ind2UI` does, then relabels the result back) —
a caller never has to pre-sort anything, unlike MATLAB's version, which
expects its caller (`roast.m`) to have already applied that permutation.
This is exercised in `tests/test_placement_integration.py` by
deliberately requesting electrodes out of pool order.

**Verified two ways:**
- `tests/test_cap_fitting.py`: a synthetic ellipsoid "head" with hand-placed
  landmarks — checks fitted electrodes land on the surface and that
  anatomically-named points end up anatomically placed (Fpz frontal, Oz
  occipital, Cz at the vertex, Fp1 left of Fp2).
- `tests/test_placement_integration.py` (`@pytest.mark.slow`, ~2.5 min):
  full pipeline (multiaxial segmentation → electrode placement) on the real
  `example/subject1.nii`, with landmarks derived heuristically from the
  segmentation's bounding box (real landmark detection isn't ported yet —
  see Remaining phases). Confirms Cz is the highest point, Fpz is more
  anterior than Oz, and gel never overlaps electrodes or other tissue.

One correctness fix worth calling out: `cleanScalp.m`'s morphological
close/open operations use structuring elements that grow up to `ones(30,30,30)`
or larger, applied directly to a MATLAB image — at full head resolution
(e.g. 192×256×256) a literal cube that size is computationally intractable
in scipy. `geometry/scalp.py` decomposes each into repeated
dilation/erosion with a 3×3×3 cube (`scipy.ndimage`'s optimized `iterations`
path), which is mathematically identical for an all-ones cube structuring
element but tractable — confirmed by the real-data test above completing in
seconds rather than hanging.

## Remaining phases (see the full plan for detail)

0. ~~I/O & preprocessing~~ done
1. ~~Segmentation~~ done (see above)
2. ~~Electrode placement & cap fitting~~ done (see above)
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
