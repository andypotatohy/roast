"""Multiaxial segmentation: default, MATLAB-free segmentation path for roast_py.

Ports lib/multiaxial/SEGMENT.py's ``main()`` (previously invoked by
runMultiaxial.m via a `system()` call into a separate conda environment) and
the inference-only pieces of lib/multiaxial/utils.py into an in-process
Python function -- there's no more MATLAB, and no more subprocess, needed
to run this segmentation path at all: it becomes roast_py's default rather
than the "for abnormal anatomy" opt-in it is in MATLAB ROAST, since it's
already free of the MATLAB/SPM dependency the rest of this port exists to
remove. The SPM12 path (roast_py.segmentation.spm_standalone) exists
alongside it for parity validation against MATLAB ROAST's default path.

Output tissue labels match ROAST's own scheme (see the top-level README's
`'conductivities'` option, listed in that order): 1=white, 2=gray, 3=csf,
4=bone, 5=skin, 6=air, 0=background. The models themselves predict a
different label order (see lib/multiaxial/README.md); _remap_to_roast_labels
translates between the two, ported unchanged from SEGMENT.py's `mapping`.
"""

from __future__ import annotations

import os
from pathlib import Path

import nibabel as nib
import numpy as np

from . import _keras_compat

_keras_compat.ensure_legacy_keras()

import tensorflow as tf  # noqa: E402  (must follow ensure_legacy_keras())

from ._losses import CUSTOM_OBJECTS  # noqa: E402
from ._preprocessing import preprocess_head_MRI, reshape_back_to_original  # noqa: E402

_MODEL_FILES = (
    "sagittal_model.h5",
    "axial_model.h5",
    "coronal_model.h5",
    "consensus_layer.h5",
)

# multiaxial labels: 0 background, 1 air, 2 gray, 3 white, 4 csf, 5 bone, 6 skin
# ROAST labels:      0 background, 1 white, 2 gray, 3 csf, 4 bone, 5 skin, 6 air
_MULTIAXIAL_TO_ROAST = {1: 6, 3: 1, 4: 3, 5: 4, 6: 5}

_models_cache: dict[str, tuple] = {}


def find_model_dir(model_dir: str | os.PathLike | None = None) -> Path:
    if model_dir is not None:
        return Path(model_dir)
    env = os.environ.get("ROAST_MULTIAXIAL_MODELS")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "lib" / "multiaxial"
        if (candidate / "sagittal_model.h5").exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the multiaxial model weights "
        f"({', '.join(_MODEL_FILES)}). Pass model_dir explicitly, or set "
        "the ROAST_MULTIAXIAL_MODELS environment variable."
    )


def _load_models(model_dir: Path) -> tuple:
    key = str(model_dir)
    if key in _models_cache:
        return _models_cache[key]

    def load(name, custom_objects=None):
        return tf.keras.models.load_model(
            model_dir / name, custom_objects=custom_objects, compile=False
        )

    models = (
        load("sagittal_model.h5", CUSTOM_OBJECTS),
        load("axial_model.h5", CUSTOM_OBJECTS),
        load("coronal_model.h5", CUSTOM_OBJECTS),
        load("consensus_layer.h5"),
    )
    _models_cache[key] = models
    return models


def segment_MRI(img, coords, model_sagittal, model_axial, model_coronal, model_layer):
    """Ported from lib/multiaxial/utils.py. Runs the three per-view 2D
    segmenters (each treating one axis of the 3D volume as a batch of 2D
    slices) and combines their per-voxel class probabilities with the
    consensus model.

    Unlike the original, exceptions from the consensus step are not
    silently swallowed (the original wraps this in a bare
    ``try/except: print('failed')``, which would otherwise return `None`
    with no indication anything went wrong).
    """
    yhat_sagittal = model_sagittal.predict(
        [np.expand_dims(img, -1), coords], batch_size=1, verbose=2
    )

    img_coronal = np.swapaxes(img, 0, 1)
    coords_coronal = np.swapaxes(coords, 0, 1)
    yhat_coronal = model_coronal.predict(
        [np.expand_dims(img_coronal, -1), coords_coronal], batch_size=1, verbose=2
    )
    yhat_coronal = np.swapaxes(yhat_coronal, 0, 1)

    img_axial = np.swapaxes(np.swapaxes(img, 1, 2), 0, 1)
    coords_axial = np.swapaxes(np.swapaxes(coords, 1, 2), 0, 1)
    yhat_axial = model_axial.predict(
        [np.expand_dims(img_axial, -1), coords_axial], batch_size=1, verbose=2
    )
    yhat_axial = np.swapaxes(np.swapaxes(yhat_axial, 0, 1), 1, 2)

    X_test = np.concatenate(
        [np.expand_dims(img / np.percentile(img, 95), -1), yhat_sagittal, yhat_coronal, yhat_axial],
        axis=-1,
    )
    yhat_new_model = model_layer.predict(np.expand_dims(X_test, 0), verbose=2)
    return np.argmax(yhat_new_model[0], -1).astype(np.int32)


def _remap_to_roast_labels(seg: np.ndarray) -> np.ndarray:
    out = seg.copy()
    for mp_label, roast_label in _MULTIAXIAL_TO_ROAST.items():
        out[seg == mp_label] = roast_label
    return out


def segment(
    t1_path: str,
    out_path: str | None = None,
    model_dir: str | os.PathLike | None = None,
    anterior_commissure: tuple | None = None,
) -> str:
    """Segments a T1 head MRI into ROAST's 6-tissue scheme, entirely in
    Python/TensorFlow (no MATLAB, no SPM). Mirrors runMultiaxial.m + SEGMENT.py's
    ``main()``, but runs in-process instead of shelling out to a separate
    conda environment.

    Returns the path to the saved `<subject>_multiaxial_masks.nii`.
    """
    resolved_model_dir = find_model_dir(model_dir)
    model_sagittal, model_axial, model_coronal, consensus_model = _load_models(resolved_model_dir)

    nii = nib.load(t1_path)
    nii_out, _, coords, _, reconstruction_parms = preprocess_head_MRI(
        nii, anterior_commissure=anterior_commissure, keep_parameters_for_reconstruction=True
    )

    seg = segment_MRI(
        nii_out.get_fdata(), coords, model_sagittal, model_axial, model_coronal, consensus_model
    )
    reconstructed = reshape_back_to_original(seg, nii, reconstruction_parms, resample_order=0).get_fdata()
    roast_masks = _remap_to_roast_labels(np.asarray(reconstructed, dtype=np.uint8))

    out_img = nib.Nifti1Image(roast_masks, nii.affine)
    if out_path is None:
        dirname = os.path.dirname(os.path.abspath(t1_path))
        base = os.path.splitext(os.path.basename(t1_path))[0]
        out_path = os.path.join(dirname, f"{base}_multiaxial_masks.nii")
    nib.save(out_img, out_path)
    return out_path
