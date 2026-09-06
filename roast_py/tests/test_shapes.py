import numpy as np

from roast_py.geometry.shapes import cylinder2p, draw_cuboid, draw_cylinder, draw_line


def test_cylinder2p_ring_radius_and_axis():
    r1 = np.array([0.0, 0.0, 0.0])
    r2 = np.array([0.0, 0.0, 10.0])
    X, Y, Z = cylinder2p(np.array([2.0, 2.0]), 100, r1, r2)

    # Every ring point should be exactly radius 2 from the z-axis.
    radii = np.sqrt(X**2 + Y**2)
    np.testing.assert_allclose(radii, 2.0, atol=1e-9)
    # z spans from 0 to 10 across the two rings.
    assert np.isclose(Z.min(), 0.0)
    assert np.isclose(Z.max(), 10.0)


def test_draw_line_endpoints_and_spacing():
    coords = draw_line([0, 0, 0], [1, 0, 0], length=10, density=2)
    assert np.allclose(coords[0], [0, 0, 0])
    assert np.allclose(coords[-1], [10, 0, 0])
    assert coords.shape[0] == 21  # 0, 0.5, 1.0, ..., 10.0


def test_draw_cylinder_disc_radius_and_height():
    top = np.array([0.0, 0.0, 0.0])
    bottom = np.array([0.0, 0.0, 5.0])
    coords = draw_cylinder(0, 3, top, bottom, density=2)

    assert coords.shape[0] > 0
    radii = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    assert radii.max() <= 3.0 + 1e-6
    assert coords[:, 2].min() >= -1e-9
    assert coords[:, 2].max() <= 5.0 + 1e-9


def test_draw_cylinder_ring_has_hollow_center():
    top = np.array([0.0, 0.0, 0.0])
    bottom = np.array([0.0, 0.0, 2.0])
    coords = draw_cylinder(2, 4, top, bottom, density=2)
    radii = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    assert radii.min() >= 2.0 - 1e-6
    assert radii.max() <= 4.0 + 1e-6


def test_draw_cuboid_bounding_box():
    center = np.array([0.0, 0.0, 0.0])
    long_axis = np.array([1.0, 0.0, 0.0])
    short_axis = np.array([0.0, 1.0, 0.0])
    normal_axis = np.array([0.0, 0.0, 1.0])
    dim = [10.0, 4.0, 2.0]  # length, width, thickness

    coords = draw_cuboid(center, dim, long_axis, short_axis, normal_axis, density=2)

    assert coords.shape[0] > 0
    assert coords[:, 0].min() >= -5.0 - 1e-6
    assert coords[:, 0].max() <= 5.0 + 1e-6
    assert coords[:, 1].min() >= -2.0 - 1e-6
    assert coords[:, 1].max() <= 2.0 + 1e-6
    assert coords[:, 2].min() >= -1.0 - 1e-6
    assert coords[:, 2].max() <= 1.0 + 1e-6
