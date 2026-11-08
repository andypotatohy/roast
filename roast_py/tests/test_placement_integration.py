"""End-to-end electrode-placement test against a real segmented head.

Chains roast_py.segmentation.multiaxial.segment() (~2-3 min on CPU) with
roast_py.geometry.placement.electrode_placement() on the real
example/subject1.nii, with landmarks from the TPM mapped through a real
NiftyReg registration (~2 min) -- the same chain roast() runs. Checks the
placement is anatomically right on real data (in MNI coordinates), not
MATLAB bit-parity.
"""

import shutil

import numpy as np
import pytest

from roast_py.geometry.cap_info import load_cap_info
from roast_py.geometry.placement import ElectrodeParams, electrode_placement

from .test_nifti import SUBJECT1


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

    from roast_py.geometry.landmarks import tpm_landmarks_to_subject
    from roast_py.registration import run_niftyreg

    registration = run_niftyreg(t1_copy, out_dir=tmp_path)
    landmarks = tpm_landmarks_to_subject(registration.tpm2mri)

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

    # And where the 10-10 system puts them in MNI space: Cz at the vertex
    # (x ~ 0, y ~ -15, z ~ 100), Fpz front and low, Oz back and low.
    mni = {n: (registration.mri2mni @ np.append(c, 1))[:3] for n, c in centroids.items()}
    assert abs(mni["Cz"][0]) < 10 and -30 < mni["Cz"][1] < 0 and mni["Cz"][2] > 85
    assert mni["Fpz"][1] > 70 and abs(mni["Fpz"][2]) < 25
    assert mni["Oz"][1] < -100 and abs(mni["Oz"][2]) < 25

    # Gel must never overlap electrode voxels or any tissue voxel.
    assert not np.any((elec_mask > 0) & (gel_mask > 0))
    assert not np.any((gel_mask > 0) & (labels > 0))
