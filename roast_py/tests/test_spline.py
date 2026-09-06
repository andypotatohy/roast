import numpy as np

from roast_py.geometry.spline import ncs2d_approx


def test_ncs2d_approx_fits_straight_line_with_two_breakpoints():
    x = np.linspace(0, 10, 50)
    y = np.zeros_like(x)  # a straight line needs no interior knots
    bx, by, breaks = ncs2d_approx(x, y)
    np.testing.assert_array_equal(breaks, [1, 50])


def test_ncs2d_approx_adds_breakpoints_for_curved_data_and_stays_within_tolerance():
    t = np.linspace(0, 2 * np.pi, 60)
    x = np.cos(t)
    y = np.sin(t)  # a circle -- can't be fit by a 2-knot spline

    bx, by, breaks = ncs2d_approx(x, y, max_allowed_sq_dist=0.01)
    assert len(breaks) > 2

    from scipy.interpolate import CubicSpline

    cs = CubicSpline(breaks, np.stack([bx, by], axis=1), bc_type="not-a-knot")
    fitted = cs(np.arange(1, len(x) + 1))
    sq_dist = (x - fitted[:, 0]) ** 2 + (y - fitted[:, 1]) ** 2
    assert sq_dist.max() <= 0.01 + 1e-9


def test_ncs2d_approx_breaks_are_sorted_and_bounded():
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.normal(size=40))
    y = np.cumsum(rng.normal(size=40))
    bx, by, breaks = ncs2d_approx(x, y, max_allowed_sq_dist=0.5)
    assert np.all(np.diff(breaks) > 0)
    assert breaks[0] == 1
    assert breaks[-1] == 40
