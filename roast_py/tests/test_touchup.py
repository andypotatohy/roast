import numpy as np

from roast_py.segmentation.touchup import (
    AIR,
    BONE,
    CSF,
    GRAY,
    SKIN,
    WHITE,
    binary_mask_generate,
    bwareaopen,
    size_of_object,
    touchup,
)


def test_binary_mask_generate_picks_argmax_and_flags_empty():
    a = np.array([1.0, 0.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0, 0.0])
    c = np.array([0.0, 0.0, 1.0, 0.0])  # last voxel: nothing claims it -> empty

    empty, mask_a, mask_b, mask_c = binary_mask_generate(a, b, c)
    np.testing.assert_array_equal(empty, [False, False, False, True])
    np.testing.assert_array_equal(mask_a, [True, False, False, False])
    np.testing.assert_array_equal(mask_b, [False, True, False, False])
    np.testing.assert_array_equal(mask_c, [False, False, True, False])


def test_size_of_object_sorts_components_descending():
    mask = np.zeros((10, 10, 10), dtype=bool)
    mask[0:3, 0:3, 0:3] = True  # 27-voxel cube
    mask[5, 5, 5] = True  # 1-voxel speck
    mask[8:10, 8:10, 8] = True  # 4-voxel block

    sizes, _ = size_of_object(mask, conn=26)
    np.testing.assert_array_equal(sizes, [27, 4, 1])


def test_bwareaopen_removes_only_small_components():
    mask = np.zeros((10, 10, 10), dtype=bool)
    mask[0:3, 0:3, 0:3] = True  # 27 voxels, kept
    mask[5, 5, 5] = True  # 1 voxel, removed

    cleaned = bwareaopen(mask, thres=2, conn=26)
    assert cleaned[1, 1, 1]  # inside the big cube
    assert not cleaned[5, 5, 5]
    assert cleaned.sum() == 27


def test_touchup_assigns_every_voxel_and_removes_specks():
    shape = (20, 20, 20)
    gray = np.zeros(shape)
    white = np.zeros(shape)
    csf = np.zeros(shape)
    bone = np.zeros(shape)
    skin = np.zeros(shape)
    air = np.zeros(shape)

    # A tiny uniform baseline everywhere (including the z=0/z=19 boundary
    # slices) stands in for the residual background probability real SPM
    # tissue-probability maps always carry -- touchup()'s empty-voxel
    # relabeling smooths each z-slice independently (matching
    # segTouchup.m's per-slice imfilter), so a slice that's exactly zero
    # for every tissue can never be filled in; that's a property of the
    # ported algorithm, not something to work around here.
    skin[:] = 0.01

    # Concentric-ish blobs: white core, gray shell, csf ring, bone ring, skin shell.
    white[8:12, 8:12, 8:12] = 1.0
    gray[6:14, 6:14, 6:14] = 0.8
    gray[8:12, 8:12, 8:12] = 0.0  # don't overlap white's core
    csf[5:15, 5:15, 5:15] = 0.3
    bone[3:17, 3:17, 3:17] = 0.2
    skin[1:19, 1:19, 1:19] = 0.1

    # A single-voxel speck of "bone" far away that should get pruned as a
    # disconnected component.
    bone[0, 0, 0] = 0.9

    result = touchup(gray, white, csf, bone, skin, air, is_smooth=False, conn=18)

    assert result.shape == shape
    assert result.dtype == np.uint8
    assert set(np.unique(result).tolist()) <= {0, WHITE, GRAY, CSF, BONE, SKIN, AIR}

    # The white matter core should indeed end up labeled white.
    assert result[10, 10, 10] == WHITE

    # The disconnected bone speck should have been pruned (reassigned to the
    # nearest real tissue during empty-voxel relabeling, not left as an
    # isolated bone voxel at the volume's corner).
    assert result[0, 0, 0] != BONE
