# roast_py — Python port of ROAST's core pipeline (work in progress)

This is an in-progress translation of ROAST's core simulation pipeline
(`roast()`, `roast_target()`, `reviewRes()`) from MATLAB to Python, with the
goal of removing the MATLAB license requirement. Scope, phasing, and the
two open technical risks (the CGAL mesher and SPM segmentation are not
plain standalone binaries the way getDP and NiftyReg are) are written up
in full in the project's plan; the short version is below.

**Status: Phases 0-4 (I/O, segmentation, electrode placement, meshing, FEM
solve) done, wired into a top-level `roast()` (see Quickstart below), plus
the visualization half of Phase 7 (`reviewRes`/`visualizeRes` and the views
they call — see Visualization below). Phases 5 and 6 not yet implemented.** Runs end to end and produces physically
sane output, but **has not yet been numerically validated against MATLAB
ROAST** (that's Phase 5 — the actual accuracy gate). Landmarks now come
from registration to MNI space, as in MATLAB (see
[Landmarks](#landmarks-registration-to-mni)). Don't use it for real
simulations until Phase 5 passes.

## Install

roast_py runs on **one exact, tested environment**: the package versions
that were installed when `roast("../example/subject1.nii")` last ran end to
end, segmentation through FEM solve. They are listed in
`requirements-lock.txt` (identical to `TESTED_ENVIRONMENT` in
`roast_py/dependencies.py`):

| | |
|---|---|
| Python | **3.11.15** (3.11–3.13 supported; 3.10 and 3.14 are not) |
| Installer | **pip only** — no conda packages; conda, if used, only provides Python |
| Key packages | tensorflow 2.21.0, tf-keras 2.21.0, keras 3.15.1, numpy 2.4.6, scipy 1.17.1, nibabel 5.4.2, pandas 3.0.5, openpyxl 3.1.5, scikit-image 0.26.0, h5py 3.14.0, matplotlib 3.11.2, pyvista 0.49.0, vtk 9.7.1 (+ their dependencies, 58 pins total) |
| Platforms with wheels | Linux x86_64/aarch64 (glibc ≥ 2.27), macOS on Apple silicon, Windows x86_64 — **not** Intel Macs (no TensorFlow 2.21 build) |

**Recommended setup** — a fresh environment on the tested Python:

```
conda create -n roast_py python=3.11 -y
conda activate roast_py
```

then just call `roast()`:

```python
from roast_py import roast      # works with nothing installed
result = roast("../example/subject1.nii")
```

Before doing any work, `roast()` compares the running interpreter against
the tested environment and installs anything missing or at a different
version, with a single `python -m pip install name==version ...` into that
interpreter. It refuses up front (with the commands above) on an
unsupported Python. After installing it re-checks every version and
test-imports TensorFlow + tf-keras in a fresh process, so a broken install
fails right away with an explanation, not minutes into segmentation.

Importing `roast_py` deliberately needs nothing but the standard library,
which is what makes that possible: if getting hold of `roast` required the
packages it installs, the auto-install could never run.

To install up front instead (all equivalent):

```
cd roast_py
python -m roast_py.dependencies --install   # into the current interpreter
pip install -r requirements-lock.txt        # same pins, by hand
conda env create -f environment.yml         # new env: python=3.11 + the lock via pip
```

`python -m roast_py.dependencies` on its own reports how the current
environment differs from the tested one.

**Opting out of auto-install** (CI, an environment you manage yourself):
pass `roast(..., install_missing=False)` or set
`ROAST_PY_NO_AUTO_INSTALL=1`. Missing packages and a mismatched
TensorFlow/tf-keras pair then raise an actionable error; other version
differences only warn.

### Why exact pins, and why pip rather than conda

The bundled segmentation models are Keras 2 `.h5` files, loadable only
through the `tf-keras` compatibility package, and `tf-keras` X.Y only works
with `tensorflow` X.Y.\*. A mismatched pair fails at import, classically
with:

```
AttributeError: module 'tensorflow._api.v2.compat.v2.__internal__'
has no attribute 'register_load_context_function'
```

(tf-keras < 2.16 calls that function; TensorFlow removed it in 2.16.)
Version ranges kept producing such pairs in practice — e.g. conda-forge
shipping TensorFlow 2.22 while the newest tf-keras on PyPI is 2.21, so no
matching tf-keras exists at all. Pinning every package to the combination
that actually ran removes that whole class of problem, and pip is where
that combination came from (TensorFlow's official channel, and the only
one that publishes tf-keras reliably).

If your existing environment already has a conda-installed TensorFlow,
`roast()` will try to replace it with the pip build and tell you it's doing
so; if that can't be done cleanly it stops and points you at the fresh
environment above, which is the reliable route.

## Quickstart

Python equivalent of MATLAB's `roast('example/subject1.nii')` (default
recipe: anode Fp1 1 mA, cathode P4 -1 mA):

```python
from roast_py import roast

result = roast("/path/to/subject1.nii")  # copy it out of example/ first -- see examples/quickstart.py
```

That's the whole call — `roast()` chains segmentation → electrode
placement → meshing → FEM solve and saves `<subj>_v.nii` (voltage),
`<subj>_e.nii` (E-field, 3 components) and `<subj>_emag.nii` (E-field
magnitude) next to the input, exactly like MATLAB's `postGetDP.m`. It also
returns a `RoastResult` with those same volumes as numpy arrays
(`result.vol_v`, `result.vol_e`, `result.ef_mag`) plus the intermediate
tissue/electrode/gel masks, for inspection without re-reading the NIfTI
files. Like MATLAB's `roast()`, it then shows the results — see
[Visualization](#visualization-phase-7) — and returns once you close the
windows (`visualize=False` skips them).

A custom montage is just a dict of electrode name → current in mA (must
sum to ~0):

```python
result = roast("/path/to/subject1.nii", {"F3": 1.0, "F4": -1.0})
```

Runnable end-to-end example: `examples/quickstart.py`. It works straight
from a git clone with no install step — it puts the package root on
`sys.path` itself, and `roast()` installs the tested environment if needed:

```
python examples/quickstart.py                    # install the tested env if needed, then run
python examples/quickstart.py --check-deps       # just report what differs from it
python examples/quickstart.py --no-install-deps  # don't install anything
python examples/quickstart.py --no-visualize     # skip the result figures
```

Verified in `tests/test_roast.py` (`@pytest.mark.slow`) and by running
exactly `roast("../example/subject1.nii")` from `roast_py/` with the default
recipe: **~5.5–6.5 min total** on CPU (331 s in the tested environment; 377 s
in a fresh venv that started out with a broken TensorFlow 2.20 + tf-keras
2.15 pair, which `roast()` replaced with the pinned versions before running),
producing a voltage volume of 0 to 362 mV and an E-field magnitude up to
31 V/m across the full 192×256×256 grid — physically sane magnitudes for a
1 mA montage (not the 10^10-scale nonsense a mis-set-up model produces,
see the FEM solve section below), with `subject1_v.nii`/`_e.nii`/`_emag.nii`
written out alongside the intermediate `.msh`/`.pro`/`.pos` files.

Those subject1 numbers predate two later changes. Landmarks now come
from NiftyReg registration, which adds about 2 minutes and moves the
electrodes to their correct positions. Meshing now uses `roast.m`'s mesh
settings, which roughly triples subject1's mesh (see
[Meshing](#meshing-phase-3)).

With both in place, the full pipeline was run on MATLAB's own default
head, `example/MNI152_T1_1mm.nii` (converted to RAS first). It took
**938 s (~15.6 min)** on 4 CPU cores, with a 384,549-node mesh, voltage
from 0 to 280 mV and E-field up to 30 V/m. The electrode centers landed
at MNI (−36, 82, −4) for Fp1 and (59, −85, 52) for P4. The slow test
`tests/test_roast.py` now runs on this head, because subject1 at these
mesh settings needs more memory than some test machines have.

Current limitations of `roast()` itself, beyond Phase 5 validation:
input must already be RAS-oriented (run `roast_py.io.nifti.convert_to_ras`
first if not — `example/subject1.nii` already is); only disc electrodes
via the 10-05 cap are wired into `roast()`'s own keyword arguments so far
(pad/ring electrodes work at the `geometry.placement.ElectrodeParams`
level but aren't exposed as `roast()` options yet); no neck/custom
electrodes, no T2-assisted segmentation, no zero-padding option.

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
  bundled `.h5` models under Keras 3 (tf-keras, part of the tested environment;
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
- `tests/test_placement_integration.py` (`@pytest.mark.slow`, ~5 min):
  full pipeline (multiaxial segmentation → NiftyReg registration →
  landmarks → electrode placement) on the real `example/subject1.nii`.
  Confirms Cz, Fpz and Oz land where the 10-10 system puts them in MNI
  space, and gel never overlaps electrodes or other tissue.

## Landmarks (registration to MNI)

Electrode placement needs the nasion, inion, both ear points and two neck
points on the individual head. As in MATLAB with the Multiaxial
segmentation, `roast()` gets them by registration rather than by
detecting them in the image:

1. **`registration/niftyreg.py`** ports `runNiftyReg.m`. It runs the
   bundled `lib/NiftyReg/<os>/reg_aladin` to affinely register the head to
   `example/MNI152_T1_1mm.nii` (about 2 minutes on CPU). It inverts
   reg_aladin's matrix into SPM's `Affine` convention (subject world → MNI
   world), and saves it with the MRI and eTPM voxel-to-world matrices as
   `<subj>_niftyReg.json`, the counterpart of `_niftyReg.mat`.
2. **`geometry/landmarks.py`** holds `roast.m`'s `landmarksInTPM`
   verbatim. These are 16 points in `eTPM.nii` voxels: the six head
   landmarks, the scalp center, and nine 10-10 positions on the midline.
   It maps them onto the head with `roast.m`'s
   `tpm2mri = inv(image.mat)·inv(Affine)·tpm.mat`, rounding the way MATLAB
   does.

The same registration also gives `mri2mni` (`Affine·image.mat`). `roast()`
saves it in `_roastOptions.json` and returns it on `RoastResult`, and the
slice viewers use it to show and accept MNI coordinates, as MATLAB's do.

One convention differs from MATLAB. `roast_py`'s matrices are nibabel
affines on **0-based** voxel indices, where SPM's `.mat` uses 1-based
ones, so every landmark here is exactly MATLAB's minus one.
`tests/test_landmarks.py` checks this against MATLAB's 1-based formulas.
The same test file runs the real reg_aladin on subject1 (slow) and
checks the landmarks land on the scalp.

This replaced a bounding-box heuristic that had been a stand-in, and the
difference is large. On subject1 the heuristic put Cz about 5 cm too far
forward (MNI y = +36) and Fp1 at nose level (MNI z = −42). With the
registered landmarks Cz lands at MNI (−1, −15, 99) and Fp1 at
(−28, 83, −2), where the 10-10 system puts them.

**Not ported:**

- `checkLandmarks.m`, the manual landmark GUI that can refit the
  registration from clicked points.
- The SPM path's registration, which reads `Affine` from `_seg8.mat`.
  `Registration` takes that format too, so it plugs in once SPM
  segmentation runs.
- `alignHeader2mni` on the outputs (the `_MNI` header copies).
- Caching: unlike MATLAB, `roast()` re-runs the registration every time
  rather than reusing an existing `_niftyReg` file.

As in MATLAB, the back-neck landmark can fall below the image on scans
that stop at the skull base. It does for subject1. It only matters if
neck electrodes are requested, and placement then refuses with MATLAB's
error.

One correctness fix worth calling out: `cleanScalp.m`'s morphological
close/open operations use structuring elements that grow up to `ones(30,30,30)`
or larger, applied directly to a MATLAB image — at full head resolution
(e.g. 192×256×256) a literal cube that size is computationally intractable
in scipy. `geometry/scalp.py` decomposes each into repeated
dilation/erosion with a 3×3×3 cube (`scipy.ndimage`'s optimized `iterations`
path), which is mathematically identical for an all-ones cube structuring
element but tractable — confirmed by the real-data test above completing in
seconds rather than hanging.

## Meshing (Phase 3)

`roast_py/meshing/` ports the volumetric mesh generation step:

| MATLAB | Python |
|---|---|
| `saveinr.m`, `readmedit.m`, `sortmesh.m`, `savemsh.m` | `meshing/mesh_io.py` |
| `cgalv2m.m`, `meshByIso2mesh.m` | `meshing/cgal_mesher.py` |

**The CGAL-mesher risk flagged in the original plan turned out not to be
real.** `lib/iso2mesh/bin/cgalmesh.mexa64` *looks* like a MATLAB MEX shared
library from its extension (iso2mesh's own naming convention), but `file`/
`readelf` show it's actually a plain, statically-linked, directly
executable ELF binary — confirmed by running it directly (`./cgalmesh.mexa64`
prints a usage message and exits 0). MATLAB's own `cgalv2m.m` already
invokes it via `system()`, exactly like the getdp/NiftyReg subprocess
wrappers elsewhere in this port, so `cgal_mesher.py` does the same — no
fallback mesher (pygalmesh/TetGen) was needed after all.

One real behavior worth documenting (found empirically, not a bug):
cgalmesh assigns output region ids by *sorted distinct nonzero input
label*, not by literal label value — an input using only labels `{1, 7}`
comes back with elements labeled `{1, 2}`. `mesh_by_iso2mesh()`'s region
numbering (1-6 tissue, then gel, then electrodes) only lines up correctly
because a real head's 6 tissue labels are essentially always all present
with no gaps; this is an existing latent assumption in MATLAB ROAST too
(same underlying binary), not something the port introduced or needs to
work around.

**Verified two ways:**
- `tests/test_cgal_mesher.py` (`@pytest.mark.slow`, ~1s): runs the real
  bundled binary on small synthetic multi-region volumes, checking the
  produced tetrahedra are valid (positive volume, in-range node
  references) and that region numbering behaves as documented above.
- Manually verified against the real pipeline output: `example/subject1.nii`
  segmented (Phase 1) → 3 electrodes placed (Phase 2) → meshed at
  `maxvol=10` produced 223,926 nodes / ~1.33M tetrahedra / 897,814
  triangles across all 12 expected regions (6 tissue + 3 gel + 3
  electrode), in about a minute. That run used cgalv2m's own sizing
  defaults; see the next paragraph.

**Mesh options now match `roast.m`.** MATLAB's `roast()` always passes
its `meshOpt` to the mesher (`radbound` 5, `angbound` 30, `distbound` 0.3,
`reratio` 3, `maxvol` 10). The Python `roast()` passed nothing, so it
meshed with iso2mesh's generic defaults (`radbound` 6, `distbound` 0.5),
which are coarser. That went unnoticed until the registered landmarks
moved the electrodes. For subject1's default montage, CGAL then missed
both thin (~3-voxel) electrode regions entirely, and the solve stopped
with "Electrode was not meshed properly". With `roast.m`'s settings both
mesh, and the head mesh is about 2.6× finer (~584k nodes), as in MATLAB.
`roast(..., mesh_options={...})` ports MATLAB's `'meshOptions'`, with the
same defaults and validation. It replaces the old `maxvol` argument.

**Memory.** At these settings getDP's direct solver (MUMPS LU) needs
several GB: about **7.6 GB** for MNI152 (384,549 nodes), and more than
8.6 GB for subject1 (583,756 nodes). MATLAB ROAST needs the same for the
same mesh. The node count is set by the surface settings (`radbound`,
`distbound`), not by `maxvol`: `maxvol=20` gave subject1 583,626 nodes.

Two things keep this manageable:

- **Solve and post-processing run as two getDP calls.** MATLAB uses one
  (`-solve EleSta_v -pos Map`), which keeps the LU factorization in memory
  while it writes the E-field. That was enough to push MNI152 over this
  sandbox's 8.6 GB limit after the solve had already finished. `run_getdp`
  now runs `-solve`, then `-pos Map` reading the saved `.res`. On a test
  mesh the voltage output is byte-identical, and the E-field agrees to
  1.4×10⁻¹¹ V/m on fields up to 4,287 V/m (round-off from reloading the
  solution).
- **A clear error when memory runs out.** If getDP is killed anyway,
  `roast()` says it probably ran out of memory, rather than MATLAB's
  generic "cannot work properly on your system".

That real run also caught a genuine performance bug before it shipped:
`read_medit()`'s and `save_msh()`'s initial implementations converted
node/element data token-by-token / row-by-row in plain Python loops, which
would be a real bottleneck at the ~11 million numbers a real head mesh
involves. Both now parse/write in bulk via numpy (`str.split()` + single
`np.array(..., dtype=...)` calls for reading, `np.savetxt` for writing).

## FEM solve (Phase 4)

`roast_py/fem/` ports boundary-condition setup, the `.pro` (getDP problem
script) generation, the getdp subprocess, and `.pos` post-processing:

| MATLAB | Python |
|---|---|
| `prepareForGetDP.m` (boundary extraction: `TriRep`/`freeBoundary`) | `fem/boundary.py`, `fem/prepare.py` |
| `solveByGetDP.m`'s `.pro` generation | `fem/pro_writer.py` |
| `solveByGetDP.m`'s `system(cmd)` call | `fem/getdp_runner.py` |
| `postGetDP.m` (`.pos` parsing + `TriScatteredInterp`) | `fem/pos_parser.py` |
| the roast() (non-lead-field) path end to end | `fem/solve.py` |

This was the plan's flagged highest-precision-risk phase (PDE assembly,
boundary tags, element ordering all have to actually match, not just
"look right"), so it got the most direct kind of check available without
a MATLAB reference to diff against: **does the real getdp binary actually
solve the problem this code hands it, and is the answer physically
correct?**

`fem/boundary.py` replaces MATLAB's `TriRep(tets, node)` + `freeBoundary`
(which finds a tet subset's free surface, returning it in a locally
renumbered point list) with a from-scratch implementation that skips that
renumbering — every boundary face is built directly from the tet list's
own global node references, so there's no reason to translate to a local
numbering only to translate back for the `.msh` boundary-element section,
which needs global references anyway. Verified against a small hand-built
two-cube conforming mesh (`tests/test_fem_boundary.py`): a unit cube's
free boundary is exactly its 12 surface triangles (area 6), and excluding
an adjacent cube's shared interface face leaves exactly 5 faces (area 5).

**Verified two ways:**
- `tests/test_fem_solve.py` (`@pytest.mark.slow`, ~3s): a small synthetic
  two-electrode "head" (6 concentric tissue shells so all tissue labels
  are present with no gaps, plus disc-shaped gel/electrode pads poking
  outward on opposite sides — closer to `placeAndModelElectrodes.m`'s real
  electrode geometry than a fully gel-enclosed electrode, which turned out
  to produce zero usable outer surface area and caught a real geometry bug
  in an early draft of this test). Runs the **real getdp binary** through
  the full chain — mesh → boundary extraction → `.pro` → solve → `.pos`
  parsing → voxel-grid interpolation — and asserts the solved potential is
  actually higher near the current-injecting anode than the current-
  sinking cathode, which only holds if the boundary condition sign, region
  wiring, and area computation all have the right sign and magnitude
  together, not just individually.
- Manually run against the real pipeline output: Phase 1-3's segmented,
  electrode-placed, meshed `example/subject1.nii` (223,926 nodes, ~1.33M
  elements), solved with a real 3-electrode montage (a 1mA/-0.5mA/-0.5mA
  Oz/Fpz/Cz montage). **Completed successfully end to end** — meshed in
  ~80s, real getdp solve in ~101s, voltage interpolated onto the full
  192×256×256 grid in ~20s, landing in a physically sane range (0 to
  ~268 mV after re-referencing) — not the 10^10-scale nonsense the
  synthetic test's bad tissue geometry produced below. A full Python
  `roast()`-equivalent run (segmentation + placement + meshing + solve) on
  this subject now takes roughly 6 minutes total, no MATLAB involved
  anywhere.

  An early version of the small synthetic test above also caught a real
  modeling mistake worth noting: the first synthetic sphere's tissue-label
  ordering put "air" (conductivity 2.5e-14 S/m, about 10 orders of
  magnitude lower than any other tissue) as the *outermost* shell directly
  touching the current-injecting gel, which is anatomically backwards (air
  pockets are internal, e.g. sinuses) and produced absurd voltages
  (10^10-scale) as the solver forced current through a nearly-insulating
  boundary layer — not a bug in the Python code, but a reminder that this
  pipeline's numerical results are only as physically sane as the tissue
  geometry handed to it.

## Visualization (Phase 7)

`roast_py/viz/` ports ROAST's visualization. `roast()` calls it at the end
of every simulation, and `review_res()` redraws a finished one from disk:

```python
from roast_py import review_res

review_res("../example/subject1.nii")                    # what roast() shows
review_res("../example/subject1.nii", tissue="skin")     # results on another tissue
review_res("../example/subject1.nii", show=False, save_dir="figs")  # PNGs only
```

| MATLAB | Python | What it draws |
|---|---|---|
| `reviewRes.m` | `viz.review_res()` | Everything below, from a finished simulation's saved outputs. Same `tissue` (`'white'`, `'gray'`, `'csf'`, `'bone'`, `'skin'`, `'air'`, `'brain'`, `'all'`) and `fastRender` options. |
| `visualizeRes.m` | `viz.visualize_res()` | Voltage and E-field on the gray matter in 3D (electrodes colored by injected current, two color bars), and voltage/E-field slice views of the brain, the E-field with direction arrows. |
| `viewMRI.m` | `viz.view_mri()` | T1 (and T2) slice viewer. |
| `viewSeg.m` | `viz.view_seg()` | Segmentation slice viewer, with viewSeg's anatomical colormap. |
| `viewElectrodes.m` | `viz.view_electrodes()` | 3D scalp (translucent), gray matter, electrodes, gel and the four head landmarks. |
| `sliceshow.m` | `viz.sliceshow()` | The interactive viewer behind all the slice views: coronal/sagittal/axial panels; click to navigate or type voxel (and MNI) coordinates; value readout; optional vector-field arrows and brain crop. |
| `brainCrop.m` | `viz.brain_crop()` | White-matter bounding box used to crop the slice views. |

**How it's displayed.** Slice viewers are matplotlib windows; the 3D
views are PyVista (VTK). All 3D views share one window, a panel each, with
linked cameras; drag to rotate. While it's open the slice viewers stay
clickable, and they stay open after it's closed. `roast()` and
`review_res()` return once every window is closed. That's needed in a
script, where exiting would close them, and it differs from MATLAB, whose
figures outlive the call. In IPython with `%matplotlib` on, the slice
viewers don't block.

**No display** (a remote server, CI): the figures are saved as PNGs in
`<work_dir>/<subj>_figures/` instead and the path is printed. The 3D views
also need OpenGL, which headless machines often lack; if VTK can't get a
context they are skipped with a message rather than crashing, since
VTK's failure there can take the whole process down. The simulation
outputs are already on disk before any of this runs, so a display problem
never costs you the simulation.

**Differences from MATLAB:**

- Voxel coordinates in the slice viewers are **0-based**, like the rest of
  roast_py: MATLAB's voxel (129, 129, 129) is (128, 128, 128) here.
- MNI coordinates come from the NiftyReg registration (see
  [Landmarks](#landmarks-registration-to-mni)). Simulations run with an
  older roast_py have none saved, so their viewers show voxel coordinates
  only.
- No simulation tags, so `review_res(subj)` takes no `simTag`. A new
  `roast()` run on the same subject and work directory replaces the
  previous one.
- `roast_target()` isn't ported, so neither is reviewRes's targeting
  branch (`tarTag`, the montage topoplot).
- `fast_render=False` smooths the displayed surface with VTK's Laplacian
  smoother, in place of iso2mesh's `sms`.
- The 3D views are panels of one window rather than separate figures. The
  color bars sit under each panel, labeled at their ends, rather than
  beside it.

To make `review_res()` possible, `roast()` now also saves
`<subj>_mask_elec.nii`, `<subj>_mask_gel.nii`, `<subj>_mesh.npz` and
`<subj>_roastOptions.json`. These are the counterparts of MATLAB's
`_mask_elec.nii`, `_mask_gel.nii`, `<subj>_<tag>.mat` and
`_roastOptions.mat`.

**Verified** on `subject1.nii` with a full `roast()` run:

- **Rendering, under Xvfb.** The 3D panels show the voltage gradient
  between the Fp1 anode (left frontal) and the P4 cathode (right
  parietal), and E-field hot spots under both. The slice views put those
  hot spots on the correct sides.
- **Interactive session, with a real Qt backend under Xvfb.** A queued
  mouse click on a slice viewer was handled while the 3D window's event
  loop was running, and moved the crosshair to the clicked voxel. Closing
  the 3D window handed over to matplotlib's loop, and closing the viewers
  returned from `review_res()`.
- **Headless.** With no display the figures were saved as PNGs, and the
  3D views were skipped cleanly on a machine without OpenGL.
- **Unit tests.** `tests/test_viz.py` covers the viewer logic (clicks,
  typed coordinates, cropping, arrows), brainCrop's thresholds, the
  mesh-to-world mapping, color ranges and color bars (drawn into a
  recording fake plotter), and the save-instead-of-show fallback.

Needs a GUI toolkit for interactive slice viewers: tkinter (included with
conda's and python.org's Python) or Qt. Without one, matplotlib can only
save files, and the slice views are saved as PNGs instead (the 3D window
still opens).

## Remaining phases (see the full plan for detail)

0. ~~I/O & preprocessing~~ done
1. ~~Segmentation~~ done (see above)
2. ~~Electrode placement & cap fitting~~ done (see above)
3. ~~Meshing~~ done (see above)
4. ~~FEM solve~~ done (see above)
5. End-to-end numerical validation against MATLAB ROAST (gate before continuing)
6. Targeting (`roast_target()`, CVX → cvxpy)
7. ~~Visualization (`reviewRes()`)~~ done (see above); packaging remains

## Running the tests

```
cd roast_py
pip install -e ".[dev]"
pytest tests/ -v
```
