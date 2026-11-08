from roast_py.fem.pro_writer import Conductivities, write_pro_file


def test_write_pro_file_structure(tmp_path):
    sigma = Conductivities(gel=[0.3, 0.3], electrode=[5.9e7, 5.9e7])
    pro_path = tmp_path / "subject1_test.pro"
    write_pro_file(
        str(pro_path),
        current=[1.0, -1.0],
        sigma=sigma,
        area_elec_needed=[10.0, 12.0],
        ind_use=[1, 2],
    )
    text = pro_path.read_text()

    assert "white = Region[1];" in text
    assert "air = Region[6];" in text
    assert "gel1 = Region[7];" in text
    assert "gel2 = Region[8];" in text
    assert "elec1 = Region[9];" in text
    assert "elec2 = Region[10];" in text
    assert "usedElec1 = Region[11];" in text
    assert "usedElec2 = Region[12];" in text

    assert "sigma[gel1] = 0.3;" in text
    assert "sigma[elec1] = 5.9e+07;" in text

    # du_dn = 1000*current/area
    assert "du_dn1[] = 100;" in text  # 1000*1.0/10.0
    assert "du_dn2[] = -83.3333;" in text  # 1000*-1.0/12.0

    assert 'Print [ v, OnElementsOf DomainC, File "subject1_test_v.pos", Format NodeTable ];' in text
    assert 'Print [ e, OnElementsOf DomainC, Smoothing, File "subject1_test_e.pos", Format NodeTable ];' in text

    assert "DomainC = Region[{white, gray, csf, bone, skin, air, gel1, gel2, elec1, elec2}];" in text


def test_write_pro_file_lead_field_tag_omits_voltage_print(tmp_path):
    sigma = Conductivities(gel=[0.3], electrode=[5.9e7])
    pro_path = tmp_path / "subject1_test.pro"
    write_pro_file(
        str(pro_path),
        current=[1.0],
        sigma=sigma,
        area_elec_needed=[10.0],
        ind_use=[1],
        lf_tag="3",
    )
    text = pro_path.read_text()
    assert "File \"subject1_test_v.pos\"" not in text
    assert 'File "subject1_test_e3.pos"' in text
