import numpy as np
import pytest

from roast_py.geometry.placement import classify_electrodes, generate_elec_mask


def test_classify_electrodes_splits_by_pool_and_sorts_by_pool_index():
    cap_names = ["Fp1", "Fp2", "Cz", "Oz"]
    elec_names = ["Oz", "Fp1", "nk2", "custom1"]
    custom_names = ["custom1", "custom2"]

    result = classify_electrodes(elec_names, cap_names, custom_names)

    # ind_p sorted by position in cap_names: Fp1 (0) before Oz (3).
    np.testing.assert_array_equal(result["ind_p"], [0, 3])
    np.testing.assert_array_equal(result["ind_n"], [1])  # nk2 -> index 1 in NECK_POOL
    np.testing.assert_array_equal(result["ind_c"], [0])  # custom1 -> index 0

    # order maps back to positions in the original elec_names list, in
    # predefined/neck/custom concatenation order.
    np.testing.assert_array_equal(result["order"], [1, 0, 2, 3])


def test_classify_electrodes_rejects_unknown_name():
    with pytest.raises(ValueError):
        classify_electrodes(["NotARealElectrode"], ["Fp1", "Fp2"])


def test_classify_electrodes_rejects_unlisted_custom_name():
    with pytest.raises(ValueError):
        classify_electrodes(["custom99"], ["Fp1"], custom_names=["custom1"])


def test_generate_elec_mask_places_labeled_voxels():
    coords = [
        np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
        np.array([[5.0, 5.0, 5.0]]),
    ]
    vol = generate_elec_mask(coords, (10, 10, 10), ["e1", "e2"], do_warn=True)
    assert vol[1, 1, 1] == 1
    assert vol[2, 2, 2] == 1
    assert vol[5, 5, 5] == 2
    assert vol.sum() == 1 + 1 + 2


def test_generate_elec_mask_raises_on_out_of_bounds_electrode():
    coords = [np.array([[100.0, 100.0, 100.0]])]
    with pytest.raises(ValueError, match="out of image boundary"):
        generate_elec_mask(coords, (10, 10, 10), ["e1"], do_warn=True)


def test_generate_elec_mask_raises_on_overlap():
    coords = [
        np.array([[1.0, 1.0, 1.0]]),
        np.array([[1.0, 1.0, 1.0]]),  # same voxel as e1
    ]
    with pytest.raises(ValueError, match="overlaps"):
        generate_elec_mask(coords, (10, 10, 10), ["e1", "e2"], do_warn=True)


def test_generate_elec_mask_warns_on_partial_out_of_bounds():
    coords = [np.array([[1.0, 1.0, 1.0], [100.0, 100.0, 100.0]])]
    with pytest.warns(UserWarning, match="Part of the electrode"):
        vol = generate_elec_mask(coords, (10, 10, 10), ["e1"], do_warn=True)
    assert vol[1, 1, 1] == 1
