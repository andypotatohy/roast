"""End-to-end electrode-placement test against a real segmented head.

Chains roast_py.segmentation.multiaxial.segment() (~2-3 min on CPU) with
roast_py.geometry.placement.electrode_placement() on the real
example/subject1.nii. Landmarks are derived heuristically from the
segmentation's bounding box (front/back/left/right-most scalp points near
the midline) rather than a real landmark-detection step (not yet ported --
see roast_py/README.md), so this checks the placement pipeline runs
correctly end-to-end and produces anatomically sane output on real data,
not MATLAB bit-parity.
"""

import shutil

import numpy as np
import pytest

from roast_py.geometry.cap_info import load_cap_info
from roast_py.geometry.placement import ElectrodeParams, electrode_placement

from .test_nifti import SUBJECT1


def _heuristic_landmarks(labels: np.ndarray):
    scalp_idx = np.argwhere(labels > 0)
    mid_x = int(round((scalp_idx[:, 0].min() + scalp_idx[:, 0].max()) / 2))

    slab = scalp_idx[np.abs(scalp_idx[:, 0] - mid_x) <= 3]
    nasion = slab[np.argmax(slab[:, 1])].astype(float)
    inion = slab[np.argmin(slab[:, 1])].astype(float)

    mid_y = int(round((scalp_idx[:, 1].min() + scalp_idx[:, 1].max()) / 2))
    mid_z = int(round((scalp_idx[:, 2].min() + scalp_idx[:, 2].max()) / 2))
    band = scalp_idx[(np.abs(scalp_idx[:, 1] - mid_y) <= 5) & (np.abs(scalp_idx[:, 2] - mid_z) <= 5)]
    right = band[np.argmax(band[:, 0])].astype(float)
    left = band[np.argmin(band[:, 0])].astype(float)

    return np.array([nasion, inion, right, left, nasion, inion])


@pytest.mark.slow
def test_electrode_placement_end_to_end_on_subject1(tmp_path):
    import nibabel as nib

    from roast_py.segmentation.multiaxial import segment

    t1_copy = tmp_path / "subject1.nii"
    shutil.copy(SUBJECT1, t1_copy)
    mask_path = segment(str(t1_copy))

    img = nib.load(mask_path)
    labels = np.asarray(img.dataobj, dtype=np.uint8)
    voxel_size = np.abs(np.diag(img.affine)[:3])

    landmarks = _heuristic_landmarks(labels)

    names, template = load_cap_info("1010")
    # Deliberately not in cap-sheet order, to exercise classify_electrodes'
    # pool-sort + relabel-back-to-request-order path.
    elec_names = ["Oz", "Fpz", "Cz"]
    elec_paras = [ElectrodeParams(elec_type="disc", elec_size=np.array([6.0, 2.0])) for _ in elec_names]

    elec_mask, gel_mask = electrode_placement(
        labels, landmarks, elec_names, elec_paras, names, template, voxel_size
    )

    assert elec_mask.shape == labels.shape
    for i in range(1, len(elec_names) + 1):
        assert np.sum(elec_mask == i) > 0
        assert np.sum(gel_mask == i) > 0

    # Anatomical sanity: Cz (vertex) should be the highest of the three;
    # Fpz (frontal) more anterior than Oz (occipital).
    centroids = {name: np.argwhere(elec_mask == i + 1).mean(axis=0) for i, name in enumerate(elec_names)}
    assert centroids["Cz"][2] > centroids["Fpz"][2]
    assert centroids["Cz"][2] > centroids["Oz"][2]
    assert centroids["Fpz"][1] > centroids["Oz"][1]

    # Gel must never overlap electrode voxels or any tissue voxel.
    assert not np.any((elec_mask > 0) & (gel_mask > 0))
    assert not np.any((gel_mask > 0) & (labels > 0))
