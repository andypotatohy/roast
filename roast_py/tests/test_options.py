"""Tests for roast_py.options (roast.m's option checking), using the
option values runExamples.m passes, and for roast()'s refusals of
options that aren't ported (which must fail fast, before any work)."""

import numpy as np
import pytest

from roast_py.options import (
    DEFAULT_ELEC_SIZE,
    on_off,
    read_custom_locations,
    resolve_conductivities,
    resolve_electrodes,
    resolve_recipe,
)

from .test_nifti import MNI152, SUBJECT1


def test_recipe_dict_and_matlab_flat_form():
    assert resolve_recipe({"Fp1": 1, "P4": -1}) == (["Fp1", "P4"], [1.0, -1.0])
    assert resolve_recipe(["F1", 0.3, "P2", 0.7, "C5", -0.6, "O2", -0.4])[0] == ["F1", "P2", "C5", "O2"]
    with pytest.raises(ValueError, match="not balanced"):
        resolve_recipe({"Fp1": 1, "P4": -0.5})
    with pytest.raises(ValueError, match="electrodeName-injectedCurrent"):
        resolve_recipe(["Fp1", 1, "P4"])
    with pytest.raises(NotImplementedError, match="leadfield"):
        resolve_recipe("leadfield")


def test_on_off():
    assert on_off("on", "x") is True and on_off("OFF", "x") is False and on_off(True, "x") is True
    with pytest.raises(ValueError, match="'on' or 'off' for option 'multiaxial'"):
        on_off("yes", "multiaxial")


def _types_sizes_oris(paras):
    return (
        [p.elec_type for p in paras],
        [p.elec_size.tolist() for p in paras],
        [p.elec_ori if isinstance(p.elec_ori, (str, type(None))) else p.elec_ori.tolist() for p in paras],
    )


def test_defaults_are_roast_ms():
    types, sizes, oris = _types_sizes_oris(resolve_electrodes(2))
    assert types == ["disc", "disc"] and sizes == [[6, 2], [6, 2]] and oris == [None, None]
    types, sizes, oris = _types_sizes_oris(resolve_electrodes(3, elec_type=["disc", "pad", "ring"]))  # cmd 14
    assert sizes == [list(DEFAULT_ELEC_SIZE[t]) for t in types]
    assert oris == [None, "lr", None]  # pads default to 'lr'


def test_run_examples_size_forms():
    # cmd 15: mixed types, one size per electrode
    _, sizes, _ = _types_sizes_oris(
        resolve_electrodes(3, elec_type=["disc", "pad", "ring"], elec_size=[[8, 2], [45, 25, 4], [5, 8, 2]])
    )
    assert sizes == [[8, 2], [45, 25, 4], [5, 8, 2]]
    # cmd 25: None entries take that type's default
    _, sizes, _ = _types_sizes_oris(
        resolve_electrodes(
            5, elec_type=["disc", "ring", "pad", "ring", "pad"], elec_size=[None, [7, 9, 3], [40, 20, 4], None, None]
        )
    )
    assert sizes == [[6, 2], [7, 9, 3], [40, 20, 4], [4, 6, 2], [50, 30, 3]]
    # cmd 31: single type, an N-row matrix
    _, sizes, _ = _types_sizes_oris(resolve_electrodes(3, elec_type="ring", elec_size=[[7, 10, 3], [6, 8, 3], [4, 6, 2]]))
    assert sizes == [[7, 10, 3], [6, 8, 3], [4, 6, 2]]
    # cmd 30: single type, one size for all
    _, sizes, _ = _types_sizes_oris(resolve_electrodes(2, elec_type="ring", elec_size=[7, 10, 3]))
    assert sizes == [[7, 10, 3], [7, 10, 3]]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        (dict(elec_type="disc", elec_size=[6, 2, 1]), "Redundant size info"),
        (dict(elec_type="pad", elec_size=[50, 30]), "Insufficient size info"),
        (dict(elec_type="pad", elec_size=[30, 50, 3]), "width of the pad should not be bigger"),
        (dict(elec_type="pad", elec_size=[50, 30, 2]), "at least be 3 mm"),
        (dict(elec_type="ring", elec_size=[6, 4, 2]), "inner radius should be smaller"),
        (dict(elec_type="disc", elec_size=[-1, 2]), "non-negative"),
        (dict(elec_type="disc", elec_size=[[6, 2], [6, 2], [6, 2]]), "different sizes for each electrode"),
        (dict(elec_type=["disc", "pad", "ring"]), "which type for each electrode"),
        (dict(elec_type="cube"), "Supported electrodes"),
        (dict(cap_type="1030"), "Supported cap types"),
        (dict(elec_type="pad", elec_ori="up"), "Unrecognized pad orientation"),
    ],
)
def test_bad_electrode_options_are_rejected_like_matlab(kwargs, message):
    with pytest.raises(ValueError, match=message):
        resolve_electrodes(2, **kwargs)


def test_run_examples_orientation_forms():
    pad3 = dict(elec_type="pad")
    # cmd 16 / 17: one keyword or one vector for every pad
    assert _types_sizes_oris(resolve_electrodes(2, **pad3, elec_ori="ap"))[2] == ["ap", "ap"]
    assert _types_sizes_oris(resolve_electrodes(2, **pad3, elec_ori=[0.71, 0.71, 0]))[2] == [[0.71, 0.71, 0]] * 2
    # cmd 18: an N-by-3 matrix, one row per pad
    rows = [[0.71, 0.71, 0], [-0.71, 0.71, 0], [0, 0.71, 0.71]]
    assert _types_sizes_oris(resolve_electrodes(3, **pad3, elec_ori=rows))[2] == rows
    mixed = dict(elec_type=["pad", "disc", "pad"])
    # cmd 19: mixed types, matrix rows go to the pads in order
    oris = _types_sizes_oris(resolve_electrodes(3, **mixed, elec_ori=[[0.71, 0.71, 0], [0, 0.71, 0.71]]))[2]
    assert oris == [[0.71, 0.71, 0], None, [0, 0.71, 0.71]]
    # cmd 20: MATLAB's cell array, one entry per electrode
    oris = _types_sizes_oris(resolve_electrodes(3, **mixed, elec_ori=["ap", None, [0, 0.71, 0.71]]))[2]
    assert oris == ["ap", None, [0, 0.71, 0.71]]
    # cmd 25: one keyword for the pads among other types
    oris = _types_sizes_oris(
        resolve_electrodes(5, elec_type=["disc", "ring", "pad", "ring", "pad"], elec_ori="ap")
    )[2]
    assert oris == [None, None, "ap", None, "ap"]


def test_orientation_for_non_pads_is_ignored_with_a_warning():
    with pytest.warns(UserWarning, match="not placing pad electrodes"):
        paras = resolve_electrodes(2, elec_type="disc", elec_ori="ap")
    assert [p.elec_ori for p in paras] == [None, None]


def test_conductivities():
    c = resolve_conductivities(None, 3)
    assert c.csf == 1.65 and c.gel == [0.3] * 3 and c.electrode == [5.9e7] * 3
    with pytest.warns(UserWarning, match="advanced options"):
        c = resolve_conductivities({"csf": 0.6, "electrode": 0.1}, 3)  # cmd 23
    assert c.csf == 0.6 and c.white == 0.126 and c.electrode == [0.1] * 3
    with pytest.warns(UserWarning):
        c = resolve_conductivities({"gel": [1, 0.3, 1], "electrode": [0.1, 5.9e7, 0.1]}, 3)  # cmd 24
    assert c.gel == [1, 0.3, 1] and c.electrode == [0.1, 5.9e7, 0.1]
    with pytest.raises(ValueError, match="Unrecognized tissue names"):
        resolve_conductivities({"brain": 0.3}, 3)
    with pytest.raises(ValueError, match="positive number for the csf"):
        resolve_conductivities({"csf": 0}, 3)
    with pytest.raises(ValueError, match="different gel conductivities"):
        resolve_conductivities({"gel": [1, 2]}, 3)


def test_custom_locations_are_read_next_to_the_mri():
    names, coords = read_custom_locations(str(SUBJECT1))
    assert names[:3] == ["custom1", "custom2", "custom3"]
    assert coords[0].tolist() == [147, 101, 222]
    with pytest.raises(FileNotFoundError, match="customLocations"):
        read_custom_locations(str(MNI152))


@pytest.mark.parametrize(
    "args, kwargs, error, message",
    [
        (("nyhead",), {}, NotImplementedError, "New York head"),
        (("missing.nii",), {}, FileNotFoundError, "does not exist"),
        ((None, "leadfield"), {}, NotImplementedError, "leadfield"),
        (
            (str(SUBJECT1),),
            {"t2": str(SUBJECT1).replace("subject1", "subject1_T2"), "multiaxial": "on"},
            ValueError,
            "Multiaxial cannot be run with both T1 and T2",
        ),
        ((), {"manual_gui": "on"}, NotImplementedError, "manual landmark GUI"),
        ((), {"zero_padding": 2.5}, ValueError, "positive integer"),
    ],
)
def test_unported_or_invalid_options_fail_before_any_work(args, kwargs, error, message, monkeypatch):
    import roast_py.segmentation.multiaxial as multiaxial

    monkeypatch.setattr(multiaxial, "segment", lambda *a, **k: pytest.fail("started segmenting"))
    from roast_py import roast

    with pytest.raises(error, match=message):
        roast(*args, **kwargs, install_missing=False, visualize=False)
