"""Optional secondary segmentation path: drives SPM12's official standalone
build (compiled with the MATLAB Compiler SDK against the free MATLAB
Runtime -- no MATLAB license needed) to reproduce ROAST's default
New-Segment pipeline (start_seg.m), then runs the same tissue-cleanup
pipeline as segTouchup.m (roast_py.segmentation.touchup) to get a 6-label
mask matching ROAST's scheme.

Python-native (no MATLAB Runtime at all) alternative considered and ruled
out: nothing widely used does ROAST's exact 6-class *full-head*
segmentation (white/gray/csf/bone/skin/air). ANTsPyNet's `deep_atropos` and
FSL's `FAST` both only segment brain tissue -- no skull/scalp/air, which
matter a lot for TES current-flow modeling, since skull's conductivity is
~2 orders of magnitude lower than scalp's. SimNIBS's `charm` pipeline (a
MATLAB-free, Python + modified-SAMSEG full-head segmenter) is the closest
real alternative, but absorbing it would mean importing a second large
third-party simulation codebase with its own compiled dependencies --
out of scope here. The bundled `multiaxial` CNN
(roast_py.segmentation.multiaxial) already *is* a Python-native, MATLAB-free
full-head segmenter producing this exact label scheme, which is why it --
not this module -- is roast_py's default segmentation path; this module
exists for parity validation against MATLAB ROAST's own default path, and
for users who specifically want SPM-equivalent segmentation.

Status: NOT runtime-tested in this environment. SPM12's standalone build
and the MATLAB Runtime (multiple GB) aren't installable here, and require a
platform-specific download from https://www.fil.ion.ucl.ac.uk/spm/ that
this port can't fetch or validate offline. The batch-job text generation
and the post-segmentation cleanup (touchup.py) ARE tested; only the
`subprocess.run()` call into the real `run_spm12.sh`/`spm12.exe` binary is
unverified end-to-end -- validate that step against a real install before
relying on this path.

Why a subprocess call to the CLI rather than SPM's Python object API: SPM's
Python packaging (built via MATLAB Compiler SDK) is a moving target
version-to-version, and guessing at its object API without being able to
install and check it here would be worse than targeting the standalone
binary's long-stable, documented command-line interface (`run_spm12.sh
<mcr_root> batch <job.m>` on Linux/Mac, `spm12.exe batch <job.m>` on
Windows) -- the same style of interface getdp/reg_aladin already use
elsewhere in this port. Swap this subprocess call for the official Python
bindings if you'd rather do that, once you can test against them.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from .touchup import touchup

_BATCH_TEMPLATE = """\
matlabbatch{{1}}.spm.spatial.preproc.channel.vols = {{'{t1_path},1'}};
matlabbatch{{1}}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{{1}}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{{1}}.spm.spatial.preproc.channel.write = [0 0];
{tissue_lines}
matlabbatch{{1}}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{{1}}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{{1}}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{{1}}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{{1}}.spm.spatial.preproc.warp.mrf = 0;
matlabbatch{{1}}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{{1}}.spm.spatial.preproc.warp.fwhm = 0;
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
"""

# Ports start_seg.m's per-tissue ngaus (number of Gaussians in SPM's
# mixture model for that tissue class) -- not a free parameter, must match
# what the bundled eTPM.nii atlas was designed for.
_TISSUE_NGAUS = [1, 1, 2, 3, 4, 2]


def _find_default_tpm() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "eTPM.nii"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the bundled eTPM.nii tissue probability atlas. "
        "Pass template_path explicitly."
    )


def write_batch_job(t1_path: str, job_path: str, template_path: str | None = None, native: bool = True) -> str:
    """Writes an SPM batch script equivalent to start_seg.m's matlabbatch
    for a single T1 image (no T2 channel; native-space output only, i.e.
    ROAST's default `norm=false`)."""
    template = str(template_path) if template_path else str(_find_default_tpm())
    n = "[1 0]" if native else "[1 1]"
    w = "[0 0]"

    tissue_lines = []
    for i, ngaus in enumerate(_TISSUE_NGAUS, start=1):
        tissue_lines.append(
            f"matlabbatch{{1}}.spm.spatial.preproc.tissue({i}).tpm = {{'{template},{i}'}};\n"
            f"matlabbatch{{1}}.spm.spatial.preproc.tissue({i}).ngaus = {ngaus};\n"
            f"matlabbatch{{1}}.spm.spatial.preproc.tissue({i}).native = {n};\n"
            f"matlabbatch{{1}}.spm.spatial.preproc.tissue({i}).warped = {w};"
        )

    script = _BATCH_TEMPLATE.format(t1_path=t1_path, tissue_lines="\n".join(tissue_lines))
    Path(job_path).write_text(script)
    return job_path


def run_spm_standalone(job_path: str, spm_standalone_path: str, mcr_root: str | None = None) -> None:
    """Calls SPM12's standalone executable to run the batch job. On
    Linux/Mac `spm_standalone_path` is the `run_spm12.sh` script shipped
    with the standalone build and `mcr_root` (required there) is the path
    to the installed MATLAB Runtime; on Windows `spm_standalone_path` is
    `spm12.exe` directly and `mcr_root` is not used.
    """
    if mcr_root is not None:
        cmd = [spm_standalone_path, mcr_root, "batch", job_path]
    else:
        cmd = [spm_standalone_path, "batch", job_path]
    subprocess.run(cmd, check=True)


def segment(
    t1_path: str,
    out_path: str | None = None,
    spm_standalone_path: str | None = None,
    mcr_root: str | None = None,
    template_path: str | None = None,
    work_dir: str | None = None,
) -> str:
    """Segments a T1 head MRI via SPM12 standalone + touchup, producing a
    6-label mask matching ROAST's scheme (same output contract as
    roast_py.segmentation.multiaxial.segment). Requires
    ``spm_standalone_path`` (and, on Linux/Mac, ``mcr_root``) pointing at a
    real SPM12 standalone + MATLAB Runtime install -- see module docstring.
    """
    if spm_standalone_path is None:
        raise ValueError(
            "spm_standalone_path is required: path to SPM12's standalone "
            "run_spm12.sh (Linux/Mac, plus mcr_root) or spm12.exe (Windows). "
            "See module docstring for background."
        )

    t1_path = os.path.abspath(t1_path)
    dirname = os.path.dirname(t1_path)
    base = os.path.splitext(os.path.basename(t1_path))[0]
    work_dir = work_dir or dirname

    job_path = os.path.join(work_dir, f"{base}_spm_job.m")
    write_batch_job(t1_path, job_path, template_path=template_path)
    run_spm_standalone(job_path, spm_standalone_path, mcr_root=mcr_root)

    # SPM's New Segment writes c1<base>.nii .. c6<base>.nii (gray, white,
    # csf, bone, skin, air, in that order) next to the T1.
    tissue_imgs = [nib.load(os.path.join(dirname, f"c{i}{base}.nii")) for i in range(1, 7)]
    gray, white, csf, bone, skin, air = (np.asarray(img.dataobj) for img in tissue_imgs)

    mask = touchup(gray, white, csf, bone, skin, air)

    t1_img = nib.load(t1_path)
    out_img = nib.Nifti1Image(mask, t1_img.affine)
    if out_path is None:
        out_path = os.path.join(dirname, f"{base}_masks.nii")
    nib.save(out_img, out_path)
    return out_path
