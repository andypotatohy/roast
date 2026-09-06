from roast_py.geometry.cap_info import find_default_cap_info, load_cap_info


def test_find_default_cap_info_locates_bundled_file():
    path = find_default_cap_info()
    assert path.name == "capInfo.xlsx"
    assert path.exists()


def test_load_cap_info_10_05_sheet_has_expected_landmarks():
    names, template = load_cap_info("1010")
    assert template.shape == (len(names), 3)
    assert "LPA" in names
    assert "RPA" in names
    assert "Nz" in names
    assert "Fp1" in names


def test_load_cap_info_biosemi_and_egi():
    names_bs, template_bs = load_cap_info("biosemi")
    assert "A1" in names_bs
    assert template_bs.shape[0] == len(names_bs)

    names_egi, template_egi = load_cap_info("egi")
    assert "FidNz" in names_egi
    assert template_egi.shape[0] == len(names_egi)


def test_load_cap_info_rejects_unknown_cap_type():
    import pytest

    with pytest.raises(ValueError):
        load_cap_info("not-a-real-cap")
