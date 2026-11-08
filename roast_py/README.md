# roast_py 0.1: ROAST in Python, no MATLAB needed

roast_py is a Python port of ROAST's
`roast()`: it builds a head model from a T1 MRI, places the electrodes
you ask for, and computes the voltage and electric field of a
transcranial electrical stimulation (TES) montage. It does what MATLAB
ROAST does, with the same options and output files, and needs no MATLAB
license.

**Status of 0.1.** The full `roast()` pipeline works: segmentation
(Multiaxial or SPM), electrode placement, meshing, the FEM solve, and
ROAST's result viewers (`reviewRes`). Electrode placement has been checked
against MATLAB ROAST's masks: pads match voxel for voxel, and discs and
rings match within the spread of MATLAB's own random rim orientation (see
[Electrode placement](#electrode-placement-and-cap-fitting)). The
E-field has been checked against MATLAB ROAST end to end on 8 of ROAST's
example simulations, and differs from MATLAB's by 4–9% on average (see
[Validation against MATLAB ROAST](#validation-against-matlab-roast)).
Not ported yet: `roast_target()`, the New York head, lead fields, and the
manual landmark GUI.

- [Quickstart](#quickstart)
- [Using roast_py](#using-roast_py)
- [Developer guide](#developer-guide)

## Quickstart

### What you need

- **The whole ROAST repository**, not just this folder. roast_py uses
  files that sit at the repository root: the getDP, CGAL and NiftyReg
  binaries and the Multiaxial models under `lib/`, `capInfo.xlsx`,
  `eTPM.nii` and `example/MNI152_T1_1mm.nii`
  ([full list](#files-roast_py-uses-from-the-rest-of-the-repository)).
- **Python 3.11, 3.12 or 3.13.** 3.11 is what roast_py is tested on.
- **Memory.** The FEM solve needs about **8 GB** for ROAST's default head
  (MNI152) and about **20 GB** for a head like `example/subject1.nii`
  at default mesh settings. A coarser mesh needs less (see
  [Memory](#memory)).
- **A computer roast_py has been tested on: Linux.** Bundled binaries also
  exist for Windows and macOS (Apple silicon, through Rosetta), but those
  haven't been tested. Intel Macs can't run the Multiaxial segmentation
  (no TensorFlow build).

### Install

```
git clone https://github.com/andypotatohy/roast.git
cd roast/roast_py
conda create -n roast_py python=3.11 -y      # or any Python 3.11-3.13 environment
conda activate roast_py
pip install -e .
```

Install with `-e` (editable). A regular `pip install .` copies only the
Python code, and roast_py then can't find the files it needs from the
rest of the repository. (Installing is optional too: from inside
`roast/roast_py`, `import roast_py` works without it.)

You don't need to install the other packages yourself. The first time you
call `roast()`, it checks your environment against the exact package
versions roast_py was tested with and installs whatever is missing with
pip (see [Dependencies](#dependencies)).

### Your first simulation

```python
from roast_py import roast

result = roast()
```

This is MATLAB's `roast()`: ROAST's default head
(`example/MNI152_T1_1mm.nii`) with its default montage, anode Fp1 at
1 mA and cathode P4 at −1 mA. It takes 10–15 minutes on 4 CPU cores.
When it's done, figures open showing the head, the electrodes, and the
voltage and E-field; `roast()` returns when you close them.

Your own head and montage:

```python
result = roast("/data/sub01/T1.nii", {"F3": 1.0, "F4": -1.0})
```

The recipe is a dict of electrode name to current in mA, and the currents
must add up to zero. Electrode names come from the 10-05 system by default
(`cap_type` for BioSemi or EGI); `Nk1`–`Nk4` are neck electrodes, and
`custom*` electrodes are read from `<subj>_customLocations`, as in
MATLAB.

### What you get

`roast()` saves its results next to the MRI (or in `work_dir=`), named as
MATLAB names them:

| File | Contents |
|---|---|
| `<subj>_<tag>_v.nii` | voltage (mV) |
| `<subj>_<tag>_e.nii` | E-field (V/m), 3 components |
| `<subj>_<tag>_emag.nii` | E-field magnitude (V/m) |
| `<subj>_<tag>_mask_elec.nii`, `_mask_gel.nii` | electrodes and gel, labeled in your recipe's order |
| `<mri>_multiaxial_masks.nii` or `<mri>_T1orT2_SPM_masks.nii` | the segmentation: 1 white matter, 2 gray matter, 3 CSF, 4 bone, 5 skin, 6 air |
| `<subj>_<tag>_roastOptions.json`, `<subj>_roastLog` | the options used, and a log of every simulation |

`<tag>` is `simulation_tag=` or, by default, the date and time
(`20261107T162433`). The results are in the space of the MRI that was
modeled, i.e. after RAS conversion, resampling and zero-padding, as in
MATLAB.

The returned `RoastResult` holds the same volumes as numpy arrays
(`result.vol_v`, `result.vol_e`, `result.ef_mag`), along with the
segmentation (`result.tissue_labels`), the electrode and gel masks, the
mesh, the landmarks, the MRI-to-MNI transform (`result.mri2mni`) and the
tag (`result.sim_tag`).

To look at a finished simulation again:

```python
from roast_py import review_res

review_res("/data/sub01/T1.nii")                          # the latest simulation
review_res("/data/sub01/T1.nii", sim_tag="mySim")          # a given one
review_res("/data/sub01/T1.nii", tissue="skin")            # results on another tissue
review_res("/data/sub01/T1.nii", show=False, save_dir="figs")  # save PNGs only
```

## Using roast_py

### Options

`roast()` takes MATLAB's options as keyword arguments, with the same
defaults and checks. MATLAB's `roast(subj, recipe, 'elecType', 'pad',
'zeroPadding', 20)` is `roast(subj, recipe, elec_type="pad",
zero_padding=20)`. Cell arrays become lists, `[]` becomes `None`, and
structs become dicts:

```python
roast(None, {"Fp1": 1, "FC4": 1, "POz": -2},
      elec_type=["pad", "disc", "pad"], elec_ori=["ap", None, [0, 0.71, 0.71]],
      elec_size=[[45, 25, 4], None, [50, 30, 3]],
      conductivities={"csf": 0.6}, mesh_options={"maxvol": 8},
      zero_padding=10, simulation_tag="mySim")
```

| MATLAB option | Python | Notes |
|---|---|---|
| `subj`, `recipe` | positional | `None` = default; the recipe is a dict (MATLAB's flat list also works). A non-RAS MRI is converted first (`<subj>_ras.nii`). |
| `'capType'` | `cap_type` | `"1020"`, `"1010"` (default), `"1005"`, `"biosemi"`, `"egi"` |
| `'elecType'`, `'elecSize'`, `'elecOri'` | `elec_type`, `elec_size`, `elec_ori` | disc, pad and ring; one value for all electrodes, or a list with one per electrode |
| `'T2'` | `t2` | a T2 image to segment with (SPM); aligned to the T1 if needed |
| `'multiaxial'` | `multiaxial` | `"on"`/`True` or `"off"`/`False`; see [Segmentation](#segmentation-multiaxial-or-spm) |
| `'meshOptions'` | `mesh_options` | dict of `radbound`, `angbound`, `distbound`, `reratio`, `maxvol` |
| `'conductivities'` | `conductivities` | dict, including per-electrode gel and electrode values |
| `'simulationTag'` | `simulation_tag` | see [Tags and resuming](#tags-and-resuming) |
| `'resampling'` | `resampling` | reslice to 1 mm voxels first (SPM) |
| `'zeroPadding'` | `zero_padding` | voxels of padding on every side |
| | `work_dir` | where to save the outputs (default: next to the MRI) |
| | `visualize` | `False` skips the figures |
| | `install_missing` | `False` installs no packages (see [Dependencies](#dependencies)) |
| | `accept_matlab_runtime_license` | for SPM in scripts; see below |
| `'manualGui'`, `'nyhead'`, `'leadfield'` | | not ported: `NotImplementedError` |

### Segmentation: Multiaxial or SPM

roast_py can segment the head two ways, as MATLAB ROAST can:

- **Multiaxial** (roast_py's default): a neural network that runs inside
  Python with TensorFlow. Nothing extra to install.
- **SPM** (MATLAB ROAST's default): used with `multiaxial="off"`, or
  whenever a `t2` is given. It runs a compiled SPM on the free MATLAB
  Runtime, which roast_py installs the first time you need it (about
  11 GB; keep about 20 GB free while it installs). Installing it means
  accepting MathWorks' license: roast_py asks once in a terminal, or pass
  `accept_matlab_runtime_license=True` (or set
  `ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE=yes`) in scripts. To install it
  ahead of time:

  ```
  python -m roast_py.spm install --accept-license
  python -m roast_py.spm check     # start SPM once and print its version
  ```

The default differs from MATLAB's so that the common case needs no MATLAB
Runtime. To get MATLAB's default, pass `multiaxial="off"`:

```python
roast("/data/sub01/T1.nii", {"F3": 1.0, "F4": -1.0}, multiaxial="off")
```

### Tags and resuming

Every simulation is saved under a tag, as in MATLAB. Running the same call
again finds the saved simulation and just loads its results. If a run is
interrupted, running it again resumes where it stopped: each of the six
steps (segmentation, registration, electrode placement, meshing, solving,
saving results) is skipped when its output already exists, and a new
montage on the same head starts at electrode placement. To redo a step,
delete its output; the steps after it are redone too. Reusing a tag for
different options is an error, as in MATLAB. Details:
[Simulation tags](#simulation-tags-and-resuming-details).

### Memory

Most of the memory goes to the FEM solve (getDP's direct solver), and it
grows with the mesh: about 8 GB for MNI152 (385k mesh nodes) and about
20 GB for `example/subject1.nii` (600–800k nodes). MATLAB ROAST needs the
same, since it runs the same solver. If getDP runs out of memory,
`roast()` says so. On a smaller machine, use a coarser mesh; the node
count is set mostly by `radbound` and `distbound`:

```python
roast(subj, recipe, mesh_options={"radbound": 8, "distbound": 0.4, "maxvol": 30})  # ~470k nodes for subject1
```

Don't go much coarser: the mesher can then miss small electrodes, and
`roast()` stops with "Electrode was not meshed properly".

### Figures

Slice viewers are matplotlib windows; the 3D views are one PyVista (VTK)
window with a panel each, and you can drag to rotate. `roast()` and
`review_res()` return once you close every window. Without a display
(a server, CI), the figures are saved as PNGs in
`<work_dir>/<subj>_figures/` instead. Interactive slice viewers need
tkinter or Qt.

### Validation against MATLAB ROAST

The E-field was compared end to end with MATLAB ROAST's on 8 simulations
from `runExamples.m` (see [Running ROAST's examples](#running-roasts-examples)):

| cmd | head | segmentation |
|---|---|---|
| 1 | MNI152 | SPM |
| 2 | MNI152 | Multiaxial |
| 5 | MNI152, zero-padded, 10-05 electrode | SPM |
| 6 | MNI152, zero-padded, 10-05 electrode | Multiaxial |
| 7 | subject1, four electrodes | SPM |
| 8 | subject1, four electrodes | Multiaxial |
| 14 | MNI152, disc, pad and ring | SPM |
| 25 | subject1 with T2, five electrode kinds | SPM |

roast_py's E-field differs from MATLAB's by **4–9% on average**, relative
to MATLAB's results. Two known differences could explain this:

- **MATLAB's `convertToRAS.m` has a small glitch** in the MRI header it
  writes, which roast_py doesn't reproduce.
- **Different SPM versions** (the SPM runs only): roast_py runs SPM 25.01,
  MATLAB ROAST its bundled, modified SPM12.

### Differences from MATLAB ROAST

- **The MRI header after RAS conversion.** roast_py reorients with
  nibabel and doesn't reproduce the glitch in MATLAB's `convertToRAS.m`
  (see [Validation](#validation-against-matlab-roast)).
- **Multiaxial is the default segmenter** (see above). With SPM, roast_py
  uses SPM 25.01, where MATLAB ROAST bundles a modified SPM12, so
  segmentations can differ slightly.
- **Voxel coordinates are 0-based** in the viewers and in `RoastResult`:
  MATLAB's voxel (129, 129, 129) is (128, 128, 128) here.
- **Electrode masks are labeled in your recipe's order.** MATLAB labels
  them in its own sorted order.
- **Discs and rings are reproducible.** MATLAB orients their rims at
  random (`rand` in `cylinder2P.m`), so every MATLAB run differs a little;
  roast_py uses a fixed orientation.
- **Loose bits of mesh are left out of the solve.** The mesher sometimes
  keeps a stray sliver of "skin" that doesn't touch the head, which can
  make the solve fail in MATLAB ROAST. roast_py leaves such pieces out and
  says so.
- **Resuming redoes the steps after one that ran;** MATLAB checks each
  file on its own and could reuse a stale mesh.
- **`review_res(subj)` shows the latest simulation;** MATLAB asks for a
  tag.
- **Figures block until closed** in a script, since the script would
  otherwise exit and close them.

### Running ROAST's examples

`examples/run_examples.py` is MATLAB's `runExamples.m`: the same calls,
in the same order. Like MATLAB, it logs failures to `errLog.txt` and
carries on.

```
python examples/run_examples.py            # all of them
python examples/run_examples.py 14 20 24   # just these
python examples/run_examples.py --list     # print the calls
```

Expect some to fail on purpose, with a message saying why: missing input
files (3, 4), a T2 with Multiaxial (10), the New York head (26–31, 33)
and lead fields (32). The `roast_target` calls and the reviews of their
results (34–37, 39–41) are left out. The calls on `subject1.nii` need about 20 GB of memory (see
[Memory](#memory)).

`examples/quickstart.py` runs `roast()` on ROAST's default head (MNI152,
about 8 GB of memory), working in a temporary folder:

```
python examples/quickstart.py                 # install the tested environment if needed, then run
python examples/quickstart.py --check-deps    # just report what differs from it
python examples/quickstart.py --no-visualize  # skip the figures
```

(c) [Yu (Andy) Huang](https://andypotatohy.github.io/), [University of Colorado at Colorado Springs](https://eas.uccs.edu/staff/dr-yu-andy-huang)

andypotatohy@gmail.com

Co-authored by Claude Code powered by Opus 5.5

November 2026

## Developer guide

Everything below is for working on roast_py: how it's laid out, how each
part was ported from MATLAB, and how it was checked.

### Layout

| Module | Ports |
|---|---|
| `roast_py/roast.py` | `roast.m`: the six steps, tags and resuming |
| `roast_py/options.py` | `roast.m`'s option handling, error messages included |
| `roast_py/simulation_log.py` | `writeRoastLog.m`, `isNewOptions.m` |
| `roast_py/io/` | `convertToRAS.m`, `zeroPadding.m`, `alignHeader2mni.m` |
| `roast_py/segmentation/` | Multiaxial (`lib/multiaxial`), `segTouchup.m` |
| `roast_py/spm/` | `start_seg.m`, `realignT2.m`, `resampToOneMM.m` through spm-python |
| `roast_py/registration/` | `runNiftyReg.m` |
| `roast_py/geometry/` | landmarks, cap fitting, electrode placement |
| `roast_py/morphology.py` | MATLAB's box `imdilate`/`imerode`/`imclose`/`imopen` |
| `roast_py/meshing/` | `meshByIso2mesh.m`, iso2mesh's `cgalv2m.m` and mesh I/O |
| `roast_py/fem/` | `prepareForGetDP.m`, `solveByGetDP.m`, `postGetDP.m` |
| `roast_py/viz/` | `reviewRes.m`, `visualizeRes.m`, `sliceshow.m` and the viewers |
| `roast_py/dependencies.py` | the tested environment and its auto-install |

### Files roast_py uses from the rest of the repository

roast_py finds these by searching the folders above its own source
files, which is why it has to run from a ROAST checkout (or an editable
install of one):

| File | Used for | Size |
|---|---|---|
| `lib/getdp-3.2.0/bin/` (`getdp`, `getdp.exe`, `getdpMac`) | FEM solve | 81 MB folder |
| `lib/iso2mesh/bin/` (`cgalmesh.mexa64`, `.mexmaci64`, `.exe` + DLLs) | meshing | 7 MB |
| `lib/NiftyReg/{linux,mac,win}/reg_aladin` | registration to MNI (Multiaxial path) | 5 MB |
| `lib/multiaxial/*.h5` | Multiaxial model weights | 13 MB |
| `lib/spm12/eTPM_*.nii` (3 files) | segTouchup's atlas masks (SPM path) | 6.5 MB |
| `eTPM.nii` | landmark atlas; SPM's tissue priors | 70 MB |
| `capInfo.xlsx` | EEG cap layouts | 44 KB |
| `example/MNI152_T1_1mm.nii` | default head; NiftyReg's template | 14 MB |

The examples and tests also use `example/subject1.nii` and friends. Each
binary can be overridden (`roast(..., getdp_bin=, cgalmesh_bin=,
niftyreg_bin=, model_dir=)`).

### Dependencies

roast_py runs on **one exact, tested environment**: the package versions
that ran the pipeline end to end, listed in `requirements-lock.txt`
(identical to `TESTED_ENVIRONMENT` in `roast_py/dependencies.py`):

| | |
|---|---|
| Python | **3.11.15** (3.11–3.13 supported; 3.10 and 3.14 are not) |
| Installer | **pip only**; conda, if used, only provides Python |
| Key packages | tensorflow 2.21.0, tf-keras 2.21.0, keras 3.15.1, numpy 2.4.6, scipy 1.17.1, nibabel 5.4.2, pandas 3.0.5, openpyxl 3.1.5, scikit-image 0.26.0, h5py 3.14.0, matplotlib 3.11.2, pyvista 0.49.0, vtk 9.7.1 (+ their dependencies, 58 pins total) |
| Platforms with wheels | Linux x86_64/aarch64 (glibc ≥ 2.27), macOS on Apple silicon, Windows x86_64; **not** Intel Macs (no TensorFlow 2.21 build) |

Before doing any work, `roast()` compares the running interpreter with
the tested environment and installs anything missing or at a different
version, with a single `python -m pip install name==version ...` into that
interpreter. It refuses up front on an unsupported Python. After
installing it re-checks every version and test-imports TensorFlow and
tf-keras in a fresh process, so a broken install fails right away with an
explanation, not minutes into segmentation.

Importing `roast_py` needs nothing but the standard library, which is
what makes that possible: if getting hold of `roast` required the
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

SPM's packages (spm-python and the compiled SPM, an 80 MB wheel) aren't
part of this environment; they're installed the first time a call needs
SPM (see [SPM through spm-python](#spm-through-spm-python)).

#### Why exact pins, and why pip rather than conda

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

### Simulation tags and resuming: details

`roast_py/simulation_log.py` ports writeRoastLog.m and isNewOptions.m. As
in MATLAB, every simulation is saved under a tag:

- **The tag** is `simulation_tag`, or by default the date and time,
  `20261025T143012` (MATLAB's `datestr(now,30)`). `RoastResult.sim_tag`
  says which one was used.
- **Options.** Once the MRI is segmented and registered, the options are
  compared with every simulation of that MRI already saved in the work
  directory. The comparison uses isNewOptions.m's list: the recipe, cap,
  each electrode's type and size (and orientation, for pads), T2,
  segmenter, the MNI registration, mesh options, conductivities,
  resampling and zero-padding.
  - If one matches, roast() carries on under its tag. If you passed a
    different tag, it warns that yours is ignored.
  - Otherwise the options are saved as `<subj>_<tag>_roastOptions.json`
    (MATLAB: `.mat`), and an entry in MATLAB's format is appended to
    `<subj>_roastLog`.
  - Reusing a tag for different options is an error, as in MATLAB.
- **Each step is skipped if its output exists:**

  | step | skipped when this exists |
  |---|---|
  | 1. segmentation | `<mri>_T1orT2_SPM_masks.nii` + `_seg8.mat` (SPM), `<mri>_multiaxial_masks.nii` (Multiaxial) |
  | 2. registration | `_seg8.mat` (SPM), `<mri>_niftyReg.json` (Multiaxial) |
  | 3. electrode placement | `<subj>_<tag>_mask_elec.nii` + `_mask_gel.nii` |
  | 4. meshing | `<subj>_<tag>_mesh.npz` + `<subj>_<tag>.msh` |
  | 5. solving | `<subj>_<tag>_v.pos` + `_e.pos`, both complete |
  | 6. saving the results | `<subj>_<tag>_v.nii`, `_e.nii`, `_emag.nii` |

  Steps 1–2 belong to the modeled MRI (`<mri>` is after RAS, resampling
  and padding), so every simulation of it shares them. Steps 3–6 belong
  to one tag. Running the same call twice just loads the results, and a
  new recipe on the same head starts at step 3.

Three things go beyond MATLAB, so that an interrupted run resumes
correctly:

- **Atomic writes.** Step outputs are written under a temporary name and
  renamed into place, so a file that exists is complete.
- **Complete solutions only.** The `.pos` files count only if getDP
  finished writing them. They are deleted when getDP fails or its
  solution is empty or NaN, so the next run solves again.
- **Redoing cascades.** When a step runs, every step after it runs too.
  Delete `<subj>_<tag>_mask_elec.nii` and the electrodes, mesh, solution
  and results are all redone. MATLAB checks each file on its own and
  would keep the old mesh.

### runExamples.m in Python

`examples/run_examples.py` translates the master branch's `runExamples.m`:
the same calls, in the same order and with the same numbers. Like MATLAB,
it appends each failure to `errLog.txt` ("error at running cmd N:" plus
the message) and moves on.

```
python examples/run_examples.py            # everything, unattended
python examples/run_examples.py 14 20 24   # just these
python examples/run_examples.py --list     # print the calls
python examples/run_examples.py --show 14  # figure windows (waits for you to close them)
```

The `roast_target` calls (34–37) and the `reviewRes` calls on their
results (39–41) are left out. After each successful `roast()` the script
also saves that run's figures in `example/run_examples_figures/cmdNN/`.
As in MATLAB, each call's outputs in `example/` carry its own tag, and
running a call again reuses them: the second run takes seconds.

Expect some calls to fail, quickly and with a message saying what's
missing. These stop in about a second with the expected error:

- missing input files: `example/bikson.hdr`/`.img` (3, 4)
- T2 plus Multiaxial (10): MATLAB's own error
- the New York head (26–31, 33)
- lead fields (32)

**The calls on subject1 need a machine with ~20 GB of memory.** At
roast.m's mesh settings, subject1 is too big to solve on the 14 GB
machine they were run on. Cmd 9 gives 801,767 nodes, and MUMPS asks
for ~19 GB. Cmd 11 gives 666,927 nodes and was killed at ~13 GB. So they
were also run there with a coarser mesh (`mesh_options={'radbound': 8,
'distbound': 0.4, 'maxvol': 30}`, about 470k nodes), passed by a separate
driver script; the examples themselves are unchanged:

- 9 (T1 + T2): 544 s.
- 11 as it then stood (SPM resampling, padding 20, and Multiaxial, since it
  left `multiaxial` out): 671 s. That run is now cmd 12. cmd 11 now passes
  `multiaxial='off'` and segments with SPM, as in MATLAB; it hasn't been run
  here since. cmd 12 at roast.m's mesh settings later stopped with MUMPS
  INFOG(1) = −10 (a singular matrix) because of a loose sliver of mesh.
  Those are now left out of the solve (see [Meshing](#meshing)), and it
  gets as far as running out of memory on the 14 GB machine.
- 22 (resampling, padding 10, T2, pads): 534 s.
- 25, as written (its own `radbound` 4 mesh): every step up to the solve
  runs, with all five electrodes meshed (823,387 nodes). It needs the
  ~20 GB machine to solve. With the T2, SPM labels part of the padding
  under the neck as skin (a slab); see the SPM section for why. roast_py
  leaves the segmentation as it is, so, as in MATLAB, Nk1 lands on the
  slab in front of the chin. It stays part of the head, so it meshes and
  solves.
- 38 (15 s): `review_res` of an `awesomeSimulation` made by cmd 9's call with
  that tag, since cmd 25 couldn't finish here.

The calls that run the pipeline were checked on the MNI152 head, which
fits this sandbox's ~8.6 GB memory limit; subject1's solve needs more.
Each took 9–10 minutes and saved figures:

- 5: zero-padding, plus 10-05 electrode Exx19
- 14: disc, pad and ring together
- 20: per-pad orientations given as MATLAB's cell array
- 24: per-electrode gel and electrode conductivities, which reach getDP's
  `.pro` file

13 (BioSemi plus neck and custom electrodes on subject1) stops at
placement with "MRI does not cover the neck". subject1's scan ends above
the back-of-neck landmark, which is the case MATLAB's own check guards
against. Its BioSemi and custom electrodes were checked separately: the
custom ones land on the scalp next to their `subject1_customLocations`
coordinates.

Each call runs in its own Python process, so memory is released between
calls and a crash or out-of-memory kill is logged as that call's error
rather than ending the batch. (`--in-process` runs them all in one.)

### MRI I/O

`roast_py/io/nifti.py` ports the pure-numerical MRI preprocessing helpers:

| MATLAB | Python | Notes |
|---|---|---|
| `convertToRAS.m` | `convert_to_ras()` | Uses nibabel's `as_closest_canonical` rather than reimplementing the sform/qform flip math by hand. Checked on `example/MNI152_T1_1mm.nii`, the LAS-oriented file the top-level README calls out. MATLAB's `convertToRAS.m` has a small glitch in the header it writes, which this doesn't reproduce; it is one likely source of the 4–9% E-field difference (see [Validation](#validation-against-matlab-roast)). |
| `zeroPadding.m` | `zero_pad()` | Verified world-coordinates are preserved across padding on the same example file. |
| `alignHeader2mni.m` | `align_header_to_mni()` | Ports the `update_affine` helper (forces sform, writes `_MNI` suffix). |

`resampToOneMM.m` and `realignT2.m` run SPM (`spm_reslice`, coregistration),
so they live in `roast_py/spm/` (see
[SPM through spm-python](#spm-through-spm-python)).

Tests in `tests/test_nifti.py` run against the real example MRIs bundled
in the repo's `example/` directory (not synthetic fixtures), so they
double as a numerical-parity check against ROAST's own documented example
data.

### Segmentation

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

- **SPM (`roast_py/spm/`), MATLAB ROAST's default path.** Runs ROAST's
  SPM steps through SPM's official Python interface,
  [spm-python](https://github.com/spm/spm-python) — see
  [SPM through spm-python](#spm-through-spm-python) below — and ports
  `segTouchup.m` in full (`roast_py/segmentation/touchup.py`) to turn SPM's
  tissue maps into the 6-label scheme.

### SPM through spm-python

`roast()` runs ROAST's SPM steps with SPM's official Python interface,
[spm-python](https://github.com/spm/spm-python). That package runs a
compiled SPM (25.01) on the free MATLAB Runtime, so no MATLAB license is
needed:

| MATLAB | Python | Used for |
|---|---|---|
| `start_seg.m`, `renameSPMres.m` | `spm.routines.spm_segment()` | segmentation with `multiaxial='off'` or a T2; writes `<mri>_T1orT2_seg8.mat` (`_T1andT2` with a T2), whose registration gives the landmarks, as in MATLAB |
| `segTouchup.m` | `segmentation.touchup.touchup()` | SPM's tissue maps → `<mri>_T1orT2_SPM_masks.nii` |
| `realignT2.m` | `spm.routines.realign_t2()` | `'T2'`: coregister and reslice the T2 onto the T1 → `<T1>_T2_aligned.nii` |
| `resampToOneMM.m` | `spm.routines.resample_to_1mm()` | `'resampling'`: reslice to 1 mm voxels → `<mri>_1mm.nii` |

Output names match MATLAB's, and a step whose output exists is skipped,
as in MATLAB.

**Installation is automatic.** The first time a call needs SPM, roast_py
does two things before any other work:

1. **Python packages.** It installs spm-python and the compiled SPM with
   pip, at pinned versions (`spm.runtime.SPM_PACKAGES`). They aren't part of
   the always-installed environment because only SPM calls need them, and
   the compiled SPM is an 80 MB wheel.
2. **MATLAB Runtime R2025b.** It finds an existing installation (MathWorks'
   default locations, `$MATLAB_RUNTIME_PATH`, or a MATLAB R2025b), or
   installs one. Installing uses the
   [matlab-runtime](https://github.com/balbasty/matlab-runtime) package
   from spm-python's authors, which downloads MathWorks' official installer
   and runs it silently. It goes into a per-user folder by default
   (`~/.local/share/roast_py/MATLAB_Runtime` on Linux,
   `~/Library/Application Support/...` on macOS, `%LOCALAPPDATA%\...` on
   Windows), so no administrator rights are needed on Linux or macOS. Set
   `$MATLAB_RUNTIME_PATH` to put it elsewhere. The installed runtime takes
   11 GB on Linux; keep about 20 GB free while installing, since MathWorks'
   installer archive and its unpacked copy need room too.

**The MATLAB Runtime license.** Installing the runtime means accepting
MathWorks' license, and roast_py never accepts it on its own. In a
terminal it asks once; for scripts and batch jobs, pass
`roast(..., accept_matlab_runtime_license=True)` or set
`ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE=yes`. Without either it stops
before doing any work and says how. To install up front:

```
python -m roast_py.spm install --accept-license   # packages + runtime
python -m roast_py.spm status                     # what's installed where
python -m roast_py.spm check                      # start SPM once, print its version
```

**How it runs.** SPM always runs in a separate Python process
(`roast_py.spm.worker`), never in the caller's. The MATLAB Runtime is
large and slow to start, and it ships its own copies of common libraries,
while roast_py's process usually has TensorFlow and VTK loaded already. On
macOS that process is started through MathWorks' `mwpython`, as MathWorks
requires, via matlab-runtime's `mwpython2`.

**Differences from MATLAB ROAST's SPM:**

- **SPM version.** spm-python ships SPM 25.01, while ROAST bundles a
  modified SPM12, so segmentations can differ slightly.
- **The `_rmask.mat` volumes.** ROAST's SPM12 also warped three atlas
  masks (skull holes, extra-cerebral CSF/eyes, a marked white-matter
  surface; `lib/spm12/eTPM_*.nii`) into the subject's space for
  `segTouchup`'s patching passes. Official SPM doesn't. Instead, roast_py
  has SPM write its inverse deformation (`iy_*.nii`) and samples the same
  masks through it, porting `spm_load_mask.m`/`spm_sample_mask.m`
  exactly, down to their degree-2 B-spline over the log mask.
- **Output names.** ROAST's SPM12 turns the `prefix` option into a
  suffix. Official SPM prepends it, so roast_py lets SPM write `r*.nii` and
  renames the result to ROAST's names.
- **With a T2, SPM can label empty space as skin.** The segmentation is
  left as SPM and segTouchup give it, as in MATLAB. SPM labels a voxel from
  the atlas prior alone when any channel is exactly zero. subject1's T2 is
  zero across much of its background (16% of voxels, against 0.1% in the
  T1), and dropping those air voxels from the fit lowers SPM's background
  weight (0.127 → 0.080 on cmd 25's image). Below the atlas, where the
  prior is about 40% soft tissue, prior-only voxels then come out as soft
  tissue (0.52 vs 0.46): part of the zero-padding under the neck becomes a
  slab of "skin". With T1 alone it stays background. ROAST's SPM12 has the
  same masking, and MATLAB's result for cmd 25 shows the slab too, with Nk1
  placed on it.

**`segTouchup.m` is now ported in full.** This includes the three patching
passes (gray matter, CSF, bone, then skull holes to skin) and the "remove
outside air" step, which the earlier port lacked. Three details were also
corrected to match MATLAB:

- SPM's maps are smoothed as uint8 integers, with MATLAB's rounding.
- Tie-breaking follows each call's argument order; air comes first when
  relabeling empty voxels.
- The even-sized 10³ and 12³ boxes are centered where MATLAB centers
  them, and erosion treats the outside of the volume as foreground.

**Verified here:**

- **Batches.** The jobs match `start_seg.m`, `realignT2.m` and
  `resampToOneMM.m` field for field, and they convert into spm-python's
  `matlabbatch` objects.
- **`_seg8.mat` reading.** It was checked on a real SPM output, the New
  York head's. The converted matrices equal the masks' and `eTPM.nii`'s
  headers, and the TPM landmarks it gives land within 3 mm of that head's
  scalp.
- **Atlas-mask sampling.** It matches the B-spline formula exactly. On the
  MNI152 head the sampled masks are left/right symmetric and anatomically
  placed.
- **Patching passes and morphology.** These are checked against MATLAB's
  definitions.
- **Installer logic.** License consent, install location, an existing
  runtime, a failed download and the macOS launcher are all covered by
  `tests/test_spm.py`.

**Verified by running SPM** (Linux, 4 CPUs, 14 GB memory):

- **Installation.** `python -m roast_py.spm install --accept-license`
  downloaded and installed MATLAB Runtime R2025b Update 7 (11 GB) in under
  a minute. `check` then started SPM25 in about 20 s.
- **Each routine on subject1.**
  - `resample_to_1mm` took 31 s; the mean intensity was unchanged.
  - `realign_t2` took 42 s. The aligned T2 correlates 0.986 with the
    original at the same world points.
  - Segment took 2–3 minutes, with T1 alone or T1 plus T2.
  - The touchup gives sensible tissue volumes and an intensity order of
    WM > skin > GM > CSF > bone > air. The eyes come out as CSF, so the
    atlas masks sampled through `iy_` work.
- **Bugs that only a real run could show,** all fixed:
  - Numbers reached MATLAB as int64, which `spm_reslice` rejects.
  - Channels and tissues were built as cells, where SPM needs struct
    arrays.
  - The masks were first sampled through `y_`, the forward deformation,
    which is on the atlas grid. They need `iy_`.
- **The full pipeline with SPM segmentation**
  (`roast('example/MNI152_T1_1mm.nii', multiaxial='off')`, all defaults)
  finished in 558 s. The brain E-field median is 0.11 V/m in WM and
  0.09 V/m in GM, peaking under Fp1 and P4.
- **runExamples' SPM calls on subject1:** see
  [runExamples.m in Python](#runexamplesm-in-python).

### Electrode placement and cap fitting

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

### Landmarks (registration to MNI)

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
- `alignHeader2mni` on the outputs (the `_MNI` header copies).

With SPM segmentation the registration comes from SPM instead: `Affine`
in `_seg8.mat`, read by the same `Registration` class. Either way it is
saved and reused, like the other steps (see
[Simulation tags](#simulation-tags-and-resuming-details)).

As in MATLAB, the back-neck landmark can fall below the image on scans
that stop at the skull base. It does for subject1. It only matters if
neck electrodes are requested, and placement then refuses with MATLAB's
error.

`cleanScalp.m` closes and opens the scalp with boxes `ones(se,se,se)`,
`se` = 10, 20, 30, … until the result is closed and smooth. A literal box
that size is slow in scipy at head resolution. `geometry/scalp.py` uses
`roast_py/morphology.py` instead, which applies each box as three 1-D
min/max filters (a box is separable). That is fast at any size and exact,
including MATLAB's centring of even-sized boxes and its border handling.

This replaced a bug. The first port iterated a 3×3×3 cube `se` times,
which is a (2·se+1)³ box: 21, 41, 61, … voxels where MATLAB uses 10, 20,
30, …, so roughly twice as big. On ordinary heads the difference was
small. On cmd 25 (T2, `zeroPadding` 30), SPM labels part of the padding
under the neck as skin, the loop runs to large `se`, and the doubled box
shaved the scalp into a block with flat faces. POz's pad then sat on a
face whose normal was parallel to its `'ap'` orientation, giving NaN and
"Electrode POz goes out of image boundary". With MATLAB's box sizes POz is
placed as in MATLAB, and Nk1 lands on the slab in front of the chin, as
MATLAB places it.

**Where the cap sits.** `fitCap2individual.m` tries 11 cap heights
(`factor` 1, 0.95, …, 0.5). It keeps the one whose central electrodes
(Oz…Fpz) best split the nasion–inion arc into ten equal parts, measuring
that arc along the scalp on the central sagittal slice. The choice can be
close (on MNI152, factors 1.0 and 0.95 differ by about 10% in fit error),
so that arc is now traced exactly as MATLAB traces it:
- `imclose`/`imopen` with MATLAB's box sizes and padding (`morphology.py`);
- MATLAB's `edge()`: Sobel, automatic threshold and thinning;
- `find()`'s column-major order, stable sorts, and the `linspace` sampling
  of the spline;
- MATLAB's rounding for the slice index.

The first port took the contour as an inner boundary in row-major order
with an unstable sort. That zig-zagged the path wherever it runs level, so
across the top of the head. On MNI152 the arc came out 500–560 mm instead
of about 417, the fit chose factor 1.0 where MATLAB chooses 0.95, and every
electrode sat 2–4 mm higher than MATLAB puts it.

**MATLAB's `imclose` pads the image.** That is, it pads, dilates, erodes,
then crops, so closing never grows a mask onto the image border. A plain
`imerode(imdilate(...))` does grow it there, because `imerode` counts the
outside as foreground. `morphology.imclose_box` now pads. On MNI152 the
unpadded closing grew the filled head onto the left and right faces, so
`cleanScalp.m`'s opening loop ran to `se` = 30 where MATLAB stops at 10.
That changed the cleaned scalp around the forehead and tilted Fp1's disc
5.8° from MATLAB's. The same closing is used in `fitCap2individual.m` and
segTouchup's outside-air step.

**Checked against MATLAB** on MATLAB's own segmentation (SPM12) and
landmarks for `roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'})`.
Each cell is the electrode-mask Dice with MATLAB's mask:

| | factor | Fp1 (disc) | FC4 (pad) | POz (ring) |
|---|---|---|---|---|
| first port | 1.0 | 0.78 | 0.87 | 0.23 |
| cap fit traced as MATLAB | 0.95, as MATLAB | 0.86 | **1.00** | 0.97 |
| + padded `imclose` | 0.95 | **0.98** | **1.00** | **0.98** |

The pad matches voxel for voxel. Discs and rings can't, even in MATLAB:
`cylinder2P.m` orients each circle with `rand(1,3)`, so every MATLAB run
voxelizes their rims a little differently. Twenty roast_py runs with that
random orientation agree with each other at Dice 0.97–0.99, the same
range as with MATLAB's masks, and every axis agrees with MATLAB's to 0.2°.
roast_py keeps a fixed orientation (`shapes.cylinder2p`), so its results
are reproducible. On roast_py's own SPM25 segmentation the electrodes
differ from MATLAB's by at most a voxel, which is the landmark difference
between SPM12 and SPM25.

**subject1, runExamples cmd 25** (all five electrode kinds: disc, rings,
a 4 mm pad, a neck ring and a custom pad), on MATLAB's segmentation. Two
more details had to follow MATLAB for the pads to match:
- **Pad layers.** A pad's gel and metal are where its cuboid meets shells
  made by dilating the scalp with `ones(h,h,h)`, `h` its thickness.
  `mask2edge_point_cloud` used scipy for that, which centres an even box
  one voxel off MATLAB, so POz's 4 mm pad sat a voxel too far out
  (Dice 0.40). It now uses `morphology.py`'s MATLAB-exact boxes.
- **Ties.** Points come in `find()`'s column-major order, directions are
  normalised in double and then single precision as in MATLAB, and every
  ranking is a stable sort. A custom electrode's local scalp patch is its
  100 nearest points, and on a voxel grid many distances tie. Ordering
  ties differently changed which points made the cut (custom1: Dice 0.996
  before, 1.000 now).

| | Fp1 (disc) | F8 (ring) | POz (4 mm pad) | Nk1 (ring) | custom1 (pad) |
|---|---|---|---|---|---|
| Dice with MATLAB | 0.98 | 0.97 | **1.00** | 0.94 | **1.00** |
| Dice between two random-orientation roast_py runs | 0.96–0.995 | 0.96–0.996 | 1.00 | 0.91–0.995 | 1.00 |

Both pads are identical to MATLAB's. The disc and rings fall within the
spread of MATLAB's own random orientation, with every axis within 0.6° of
MATLAB's. These runs used roast_py's SPM25 landmarks; MATLAB's `seg8.mat`
for this run wasn't available.

### Meshing

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
  segmented → 3 electrodes placed → meshed at
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
8.6 GB for subject1 (583,756 nodes). SPM's segmentation of subject1
meshes finer still (801,767 nodes), and for that MUMPS asks for
**~19 GB**. Its 1-mm resampled, 20-voxel padded version (cmd 11, 666,927
nodes) was killed at ~13 GB on a 14 GB machine. MATLAB ROAST needs the
same for the same mesh: it ships the same getDP, built on complex-valued
PETSc, so every factor entry takes 16 bytes. The node count is set by the surface settings (`radbound`,
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
  generic "cannot work properly on your system". Memory can also run out
  without a kill. When MUMPS can't allocate the factors up front
  (INFOG(1) = −13, which is what subject1's 19 GB request hits), PETSc
  fills the solution with NaN and getDP exits 0, and MATLAB then reports
  "getDP did not converge". roast_py passes PETSc's
  `-ksp_error_if_not_converged`, so getDP exits with MUMPS's error, and
  `roast()` reports it as running out of memory.
- **Loose pieces of the mesh are left out of the solve.** Not in MATLAB
  ROAST. Multiaxial's segmentation leaves stray bits of "skin" just
  outside the scalp, and the mesher sometimes meshes one as a sliver
  that doesn't touch the head. No current reaches it and no voltage is
  fixed on it, so it makes the matrix singular. Usually round-off hides
  that, but runExamples cmd 12 on subject1 got a 14-node sliver 0.7 mm
  off the scalp, and MUMPS stopped with INFOG(1) = −10, "numerically
  singular matrix" (MATLAB ROAST is exposed to the same luck: moving the
  electrodes 2–4 mm changed which slivers the mesher kept).
  `prepare_for_getdp` now leaves every piece with no gel or electrode in
  it out of `_ready.msh` and says how many tetrahedra it left out. The
  saved `.msh` is unchanged, and these pieces can't affect the field in
  the head.

That real run also caught a genuine performance bug before it shipped:
`read_medit()`'s and `save_msh()`'s initial implementations converted
node/element data token-by-token / row-by-row in plain Python loops, which
would be a real bottleneck at the ~11 million numbers a real head mesh
involves. Both now parse/write in bulk via numpy (`str.split()` + single
`np.array(..., dtype=...)` calls for reading, `np.savetxt` for writing).

### FEM solve

`roast_py/fem/` ports boundary-condition setup, the `.pro` (getDP problem
script) generation, the getdp subprocess, and `.pos` post-processing:

| MATLAB | Python |
|---|---|
| `prepareForGetDP.m` (boundary extraction: `TriRep`/`freeBoundary`) | `fem/boundary.py`, `fem/prepare.py` |
| `solveByGetDP.m`'s `.pro` generation | `fem/pro_writer.py` |
| `solveByGetDP.m`'s `system(cmd)` call | `fem/getdp_runner.py` |
| `postGetDP.m` (`.pos` parsing + `TriScatteredInterp`) | `fem/pos_parser.py` |
| the roast() (non-lead-field) path end to end | `fem/solve.py` |

This was the step with the highest precision risk (PDE assembly,
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
- Manually run against the real pipeline output (an early run, before the
  mesh settings matched roast.m): the segmented,
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

### Visualization

`roast_py/viz/` ports ROAST's visualization. `roast()` calls it at the end
of every simulation, and `review_res()` redraws a finished one from disk:

```python
from roast_py import review_res

review_res("../example/subject1.nii")                    # the latest simulation, as roast() shows it
review_res("../example/subject1.nii", sim_tag="awesomeSimulation")  # a given one
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
- `sim_tag` may be left out: `review_res(subj)` shows the subject's most
  recent simulation, where MATLAB asks for a tag. An unknown tag's error
  lists the saved ones.
- `roast_target()` isn't ported, so neither is reviewRes's targeting
  branch (`tarTag`, the montage topoplot).
- `fast_render=False` smooths the displayed surface with VTK's Laplacian
  smoother, in place of iso2mesh's `sms`.
- The 3D views are panels of one window rather than separate figures. The
  color bars sit under each panel, labeled at their ends, rather than
  beside it.

To make `review_res()` possible, `roast()` also saves
`<subj>_<tag>_mask_elec.nii`, `_mask_gel.nii`, `_mesh.npz` and
`_roastOptions.json`. These are the counterparts of MATLAB's
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

### Roadmap

Done in 0.1: MRI I/O, segmentation (Multiaxial and SPM), landmarks,
electrode placement, meshing, the FEM solve, tags and resuming, and
visualization, with the E-field checked against MATLAB ROAST end to end
on 8 examples.

Next:

1. **Targeting:** `roast_target()` (CVX → cvxpy) and reviewRes's
   targeting views.
2. **Packaging:** bundle the files roast_py uses from the rest of the
   repository (see [above](#files-roast_py-uses-from-the-rest-of-the-repository))
   so that roast_py installs on its own with a regular `pip install`.
3. The New York head, lead fields and the manual landmark GUI.

### Running the tests

```
cd roast_py
pip install -e ".[dev]"
pytest tests/ -m "not slow"   # about 30 s
pytest tests/                 # also the slow end-to-end tests (real models, binaries and getDP)
```

(c) [Yu (Andy) Huang](https://andypotatohy.github.io/), [University of Colorado at Colorado Springs](https://eas.uccs.edu/staff/dr-yu-andy-huang)

andypotatohy@gmail.com

Co-authored by Claude Code powered by Opus 5.5

November 2026
