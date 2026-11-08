"""Ports fitCap2individual.m: places a standard EEG cap layout (10-05,
BioSemi, or EGI template coordinates) onto an individual's scalp surface.

Landmarks (nasion, inion, right ear, left ear) are 0-based (x, y, z) voxel
coordinates, matching this package's convention everywhere else.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.interpolate import CubicSpline

from ..morphology import imclose_box, imopen_box
from .landmarks import _matlab_round
from .point_cloud import map2points, project2closest_surface_points
from .spline import ncs2d_approx


def _round_1based(x) -> int:
    """MATLAB's round() of a 0-based coordinate, taken on the 1-based value
    as fitCap2individual.m does, returned 0-based. (Python's round() sends
    halves to even: round(90.5) is 90 where MATLAB's round(91.5) is 92.)"""
    return int(_matlab_round(np.asarray(x, dtype=float) + 1.0)) - 1


def _matlab_edge(a: np.ndarray) -> np.ndarray:
    """MATLAB's edge(a) with its defaults: Sobel gradients (imfilter with
    'replicate' borders), the automatic threshold 4*mean(gradient^2), and
    thinning. Thinning keeps a pixel whose gradient magnitude is >= the
    previous neighbour's and > the next one's, horizontally or vertically,
    so on a step edge it keeps the pixel after the step in each direction:
    outside the silhouette on some sides and inside on others, not always
    the inner boundary."""
    a = np.asarray(a, dtype=float)
    op = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]]) / 8.0  # fspecial('sobel')/8
    bx = ndimage.correlate(a, op.T, mode="nearest")
    by = ndimage.correlate(a, op, mode="nearest")
    b = bx * bx + by * by
    cutoff = 4 * b.mean()
    bx, by = np.abs(bx), np.abs(by)
    m, n = b.shape
    left, right = b[:, np.r_[0, 0 : n - 1]], b[:, np.r_[1:n, n - 1]]
    up, down = b[np.r_[0, 0 : m - 1], :], b[np.r_[1:m, m - 1], :]
    tol = 100 * np.finfo(float).eps
    return (b > cutoff) & (
        ((bx >= by - tol) & (left <= b) & (right < b)) | ((by >= bx - tol) & (up <= b) & (down < b))
    )


def _central_sagittal_arc_length(scalp: np.ndarray, central_sag: int, inion, nasion):
    """Ports the non-EGI branch of fitCap2individual.m: measures the
    nasion-inion distance along the scalp surface (not straight-line) on
    the central sagittal slice, by tracing that slice's boundary and
    fitting a compact smoothing spline to it.

    Returns (distance_all, yi, zi): the total arc length, and the fitted
    center-line points' (row, col) = (y, z) coordinates, used later to
    measure inter-electrode spacing along the same line.

    This arc length is what fitCap2individual.m compares electrode
    spacings against to choose the cap's height (`factor`), so every step
    follows MATLAB exactly: its box morphology, its edge(), find()'s
    column-major order, and stable sorts. Approximations here (a plain
    inner boundary, row-major order and an unstable sort, which zig-zags
    the path wherever the contour runs level) made the arc ~20% too long,
    which on MNI152 tipped the fit from MATLAB's factor 0.95 to 1.0 and
    put every electrode 2-4 mm higher than MATLAB does.
    """
    img_c = scalp[central_sag, :, :] > 0

    # Force-close the bottom of the slice: some MRIs cut off mid-face/neck,
    # which would otherwise leave the silhouette open.
    cols_with_content = np.where(img_c.any(axis=0))[0]
    indc = cols_with_content[0]
    rows_with_content = np.where(img_c[:, indc])[0]
    indr1, indr2 = rows_with_content[0], rows_with_content[-1]
    img_c = img_c.copy()
    img_c[indr1 : indr2 + 1, indc] = True

    ytemp, ztemp = np.nonzero(img_c)
    centroid = (_round_1based(ytemp.mean()), _round_1based(ztemp.mean()))

    se = 0
    is_filled = False
    im_test = img_c
    while not is_filled:
        se += 8
        im_test = ndimage.binary_fill_holes(imclose_box(img_c, se))  # imfill(imclose(img_c,ones(se,se)),'holes')
        is_filled = bool(im_test[centroid[0], centroid[1]])

    bw_c = _matlab_edge(imopen_box(im_test, 3))  # edge(imopen(im_test,ones(3,3)))

    c_c, r_c = np.nonzero(bw_c.T)  # [r_c,c_c] = find(bw_c): column-major, (row, col) = (y, z)

    matches_inion = np.where(c_c == inion[2])[0]
    matches_nasion = np.where(c_c == nasion[2])[0]
    if len(matches_inion) == 0 or len(matches_nasion) == 0:
        raise ValueError(
            "Could not find the inion/nasion z-coordinate on the central "
            "sagittal slice boundary -- landmarks may not lie exactly on "
            "the scalp surface."
        )
    idx_inion = matches_inion[0]
    idx_nasion = matches_nasion[-1]

    i_max = np.argmax(c_c)
    right_up = np.where((c_c >= c_c[idx_inion]) & (r_c < r_c[i_max]))[0]
    right_up = right_up[np.argsort(c_c[right_up], kind="stable")]  # MATLAB's sort() is stable
    right_down = np.where((c_c >= c_c[idx_nasion]) & (r_c >= r_c[i_max]))[0]
    right_down = right_down[np.argsort(-c_c[right_down], kind="stable")]  # sort(...,'descend'), stable
    index = np.concatenate([right_up, right_down])

    bx, by, breaks = ncs2d_approx(r_c[index].astype(float), c_c[index].astype(float))
    cs = CubicSpline(breaks, np.stack([bx, by], axis=1), bc_type="not-a-knot")
    n = breaks[-1]
    yizi = cs(np.linspace(1, n, int(np.floor(n))))  # linspace(1, finalbreaks(end), finalbreaks(end))
    yi, zi = yizi[:, 0], yizi[:, 1]

    distance_all = np.sum(np.sqrt(np.diff(yi) ** 2 + np.diff(zi) ** 2))
    return distance_all, yi, zi


def fit_cap_to_individual(
    scalp: np.ndarray,
    scalp_surface: np.ndarray,
    landmarks: np.ndarray,
    voxel_size: np.ndarray,
    cap_names: list[str],
    cap_template: np.ndarray,
    ind_need: np.ndarray,
    is_biosemi: bool = False,
    is_egi: bool = False,
):
    """Ports fitCap2individual.m.

    `landmarks` is (6, 3): nasion, inion, right, left, front_neck, back_neck
    (only the first 4 are used here). `voxel_size` is (3,) mm/voxel
    (diag of the MRI's affine, used to correct for non-1mm/anisotropic
    resolution the same way ROAST's `hdrInfo(1).mat` diagonal does).
    `cap_names`/`cap_template` come from geometry.cap_info.load_cap_info.
    `ind_need` is 0-based indices into `cap_names`/`cap_template` for the
    electrodes actually being placed.

    Returns (electrode_coord, center).
    """
    nasion, inion, right, left = landmarks[0], landmarks[1], landmarks[2], landmarks[3]

    L = np.linalg.norm(inion - nasion)
    line_center = (inion + nasion) / 2

    if not is_egi:
        central_sag = _round_1based(line_center[0])  # round(line_center(1)), 1-based in MATLAB
        distance_all, yi, zi = _central_sagittal_arc_length(scalp, central_sag, inion, nasion)

    if not is_egi:
        if not is_biosemi:
            central_elec = ["Oz", "POz", "Pz", "CPz", "Cz", "FCz", "Fz", "AFz", "Fpz"]
        else:
            central_elec = ["A19", "POzAid", "A6", "CPzAid", "A1", "FCzAid", "E17", "AFzAid", "E12"]
        ind_central_elec = np.array([cap_names.index(name) for name in central_elec])
    else:
        ind_central_elec = np.array([], dtype=int)

    ind_fit = np.concatenate([ind_central_elec, np.asarray(ind_need, dtype=int)]).astype(int)
    elec_template = cap_template[ind_fit] / voxel_size

    theta = 23
    alpha = ((360 - 10 * theta) / 2) * (np.pi / 180)
    h = (L / 2) * (1 / np.tan(alpha))

    s = right - left
    s = s / np.linalg.norm(s)
    c = nasion - inion
    c = c / np.linalg.norm(c)
    a = np.cross(s, c)
    a = a / np.linalg.norm(a)

    factor = np.arange(1, 0.5 - 1e-9, -0.05) if not is_egi else np.array([1.0])
    F = np.zeros(len(factor))
    centers = np.zeros((len(factor), 3))
    elec_coords = np.zeros((elec_template.shape[0], 3, len(factor)))

    scale = round(max(scalp.shape) / 2)

    for n_i, f in enumerate(factor):
        center = line_center + h * f * a
        centers[n_i] = center
        shift = center

        affine = scale * np.block(
            [
                [s.reshape(3, 1), c.reshape(3, 1), a.reshape(3, 1), (shift / scale).reshape(3, 1)],
                [np.zeros((1, 3)), np.array([[1.0 / scale]])],
            ]
        )

        elec_adjusted = np.vstack([elec_template.T, np.ones(elec_template.shape[0])])
        elec_adjusted[2, :] *= f
        elec_transformed = (affine @ elec_adjusted)[:3, :].T

        cosine_sorted, ind_on_scalp = project2closest_surface_points(elec_transformed, scalp_surface, center)
        idx = np.zeros(elec_transformed.shape[0], dtype=int)
        for i in range(len(idx)):
            best = cosine_sorted[:, i] == cosine_sorted[:, i].max()
            test_pts = scalp_surface[ind_on_scalp[best, i]]
            _, idx_farthest = map2points(center, test_pts, "farthest")
            idx[i] = ind_on_scalp[np.nonzero(best)[0][idx_farthest[0]], i]

        elec_interp = scalp_surface[idx]
        elec_coords[:, :, n_i] = elec_interp

        if not is_egi:
            n_central = len(ind_central_elec)
            center_points = np.vstack([inion, elec_interp[:n_central], nasion])
            center_fit = np.stack([np.full(len(yi), central_sag), yi, zi], axis=1)

            indxsave = 0
            distances = np.zeros(center_points.shape[0] - 1)
            for ii in range(1, center_points.shape[0]):
                d = np.sqrt(np.sum((center_points[ii] - center_fit) ** 2, axis=1))
                indx = int(np.argmin(d))
                distances[ii - 1] = np.sum(
                    np.sqrt(np.diff(yi[indxsave : indx + 1]) ** 2 + np.diff(zi[indxsave : indx + 1]) ** 2)
                )
                indxsave = indx
            F[n_i] = np.sum(np.abs(distances - distance_all / 10))

    best = int(np.argmin(F)) if not is_egi else 0
    n_central = len(ind_central_elec)
    if len(ind_need) > 0:
        electrode_coord = elec_coords[n_central:, :, best]
    else:
        electrode_coord = elec_coords[:n_central, :, best]
    center = centers[best]
    return electrode_coord, center
