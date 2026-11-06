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


def test_clean_scalp_opens_with_matlabs_box_size():
    """cleanScalp.m opens the scalp with ones(se,se,se), se = 10, 20, ...
    The port once used a (2*se+1)^3 cube instead, twice as big, which
    shaved a head with a slab of "skin" under it (runExamples cmd 25) into
    a block with flat faces. A 12-voxel bump survives MATLAB's first
    10-voxel opening; a 21-voxel cube would erase it."""
    from roast_py.geometry.point_cloud import mask2edge_point_cloud
    from roast_py.geometry.scalp import clean_scalp

    shape = (90, 90, 90)
    xx, yy, zz = np.meshgrid(*[np.arange(n) for n in shape], indexing="ij")
    r = np.sqrt((xx - 45) ** 2 + (yy - 45) ** 2 + (zz - 45) ** 2)
    scalp = (r <= 30) & (r >= 26)
    bump = np.zeros(shape, bool)
    bump[39:51, 39:51, 74:80] = True  # 12 x 12 voxels across, on top of the head
    scalp |= bump
    surface, _ = mask2edge_point_cloud(scalp, "erode", np.ones((3, 3, 3), dtype=bool))

    scalp_out, scalp_filled = clean_scalp(scalp, surface)

    assert scalp_filled[45, 45, 45]  # closed and filled
    assert scalp_out[bump].mean() > 0.9


def _disc_scalp(radius=60):
    """A head-like slab: on every sagittal slice a disc of `radius` voxels
    centred at z = 70, with the inion and nasion at its widest points, so
    the arc between them over the top is half the circle."""
    shape = (5, 160, 140)
    yy, zz = np.meshgrid(np.arange(shape[1]), np.arange(shape[2]), indexing="ij")
    scalp = np.zeros(shape, bool)
    scalp[:] = (yy - 80) ** 2 + (zz - 70) ** 2 <= radius**2
    rows = np.where(scalp[2, :, 70])[0]
    return scalp, np.array([2, rows[0], 70]), np.array([2, rows[-1], 70])


def test_central_sagittal_arc_length_follows_the_scalp():
    """fitCap2individual.m picks the cap's height by comparing electrode
    spacing with this arc. The port once ordered the contour row-major
    with an unstable sort, which zig-zagged it wherever it runs level (the
    top of the head), and traced an inner boundary instead of MATLAB's
    edge(): on MNI152 the arc came out ~20% long and the cap one step too
    high, every electrode 2-4 mm above MATLAB's. The arc over half a disc
    of radius 60 is pi*60 ~= 188.5."""
    from roast_py.geometry.cap_fitting import _central_sagittal_arc_length

    scalp, inion, nasion = _disc_scalp(60)
    distance_all, yi, zi = _central_sagittal_arc_length(scalp, 2, inion, nasion)
    assert abs(distance_all - np.pi * 60) < 0.03 * np.pi * 60
    assert zi.max() == pytest.approx(130, abs=1.5)  # passes over the top, z = 70 + 60
    assert np.all(np.abs(np.diff(yi)) < 2)  # no zig-zag: consecutive samples stay close


def test_matlab_edge_keeps_the_pixel_after_each_step():
    """MATLAB's edge() thins a step edge to the pixel after it (scanning
    left-to-right and top-to-bottom): inside the shape on its left and top
    sides, outside on its right and bottom sides."""
    from roast_py.geometry.cap_fitting import _matlab_edge

    a = np.zeros((20, 20), bool)
    a[5:15, 5:15] = True
    e = _matlab_edge(a)
    assert e[8:12, 5].all() and not e[8:12, 4].any()  # left side: inside
    assert e[8:12, 15].all() and not e[8:12, 14].any()  # right side: outside
    assert e[5, 8:12].all() and not e[4, 8:12].any()  # top: inside
    assert e[15, 8:12].all() and not e[14, 8:12].any()  # bottom: outside


def test_mask2edge_point_cloud_matches_matlab_for_even_boxes_and_order():
    """A 4 mm pad (runExamples cmd 25, POz) builds its gel/electrode layers
    by dilating the scalp with ones(4,4,4). MATLAB centres an even box at
    floor((n+1)/2); scipy centred it one voxel over, which moved that pad a
    voxel outward. Points also come in find()'s column-major order, which
    MATLAB's stable sorts use to break ties."""
    from roast_py.geometry.point_cloud import mask2edge_point_cloud
    from roast_py.morphology import imdilate_box

    mask = np.zeros((12, 12, 12), bool)
    mask[5:7, 5:7, 5:7] = True
    points, dilated = mask2edge_point_cloud(mask, "dilate", np.ones((4, 4, 4), bool))
    assert np.array_equal(dilated, imdilate_box(mask, 4))
    assert points[:, 0].min() == 4 and points[:, 0].max() == 8  # MATLAB: 1 voxel before, 2 after (scipy: 2 and 1)
    keys = points[:, 2] * 144 + points[:, 1] * 12 + points[:, 0]
    assert np.all(np.diff(keys) > 0)  # column-major: x varies fastest


def test_map2points_breaks_ties_like_matlab():
    """MATLAB's sort is stable, so among equally distant points the first
    in the cloud wins for 'closest' and 'closer', and the last for
    'farthest'. Custom electrodes take their local scalp patch from the
    100 'closer' points, where voxel-grid distances tie all the time."""
    from roast_py.geometry.point_cloud import map2points

    goal = np.array([[1.0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 2]])
    _, closest = map2points(np.zeros(3), goal, "closest")
    _, closer = map2points(np.zeros(3), goal, "closer", 3)
    _, farthest = map2points(np.zeros(3), goal[:4], "farthest")
    assert closest[0] == 0 and closer[:, 0].tolist() == [0, 1, 2] and farthest[0] == 3
