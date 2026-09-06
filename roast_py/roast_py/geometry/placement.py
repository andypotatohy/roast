"""Ports placeNeckElec.m, generateElecMask.m, placeAndModelElectrodes.m,
and electrodePlacement.m -- placing and voxelizing electrodes/gel on the
scalp surface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from .cap_fitting import fit_cap_to_individual
from .point_cloud import map2points, mask2edge_point_cloud, project2closest_surface_points
from .scalp import clean_scalp
from .shapes import draw_cuboid, draw_cylinder

_ELEC_ORI_KEYWORDS = {"lr": [1, 0, 0], "ap": [0, 1, 0], "si": [0, 0, 1]}


@dataclass
class ElectrodeParams:
    """One electrode's placement options, mirroring a row of ROAST's
    `elecPara` struct array."""

    elec_type: str  # 'disc' | 'pad' | 'ring'
    elec_size: np.ndarray  # disc: [radius, height]; pad: [length, width, height]; ring: [r_in, r_out, height]
    elec_ori: np.ndarray | str = "lr"  # pad only
    cap_type: str = "1010"


def place_neck_electrodes(scalp_shape, scalp_surface: np.ndarray, landmarks: np.ndarray, ind_need: np.ndarray):
    """Ports placeNeckElec.m. `scalp_shape` is the scalp volume's shape,
    used (per the original) to offset the left/right neck electrodes from
    the neck center by half the head's x-extent."""
    front_neck, back_neck = landmarks[4], landmarks[5]
    neck_center = (front_neck + back_neck) / 2

    left_neck = neck_center.copy()
    left_neck[0] -= round(scalp_shape[0] / 2)
    right_neck = neck_center.copy()
    right_neck[0] += round(scalp_shape[0] / 2)

    neck_elec = np.stack([front_neck, back_neck, left_neck, right_neck])[np.asarray(ind_need, dtype=int)]

    idx = np.zeros(neck_elec.shape[0], dtype=int)
    cosine_sorted, ind_on_scalp = project2closest_surface_points(neck_elec, scalp_surface, neck_center)
    for i in range(len(idx)):
        thresh = np.percentile(cosine_sorted[:, i], 99.99)
        test_pts = scalp_surface[ind_on_scalp[cosine_sorted[:, i] > thresh, i]]
        _, idx_farthest = map2points(neck_center, test_pts, "farthest")
        idx[i] = ind_on_scalp[np.nonzero(cosine_sorted[:, i] > thresh)[0][idx_farthest[0]], i]

    neck_coord = scalp_surface[idx]
    return neck_coord, neck_center


def generate_elec_mask(elec_all_coord: list[np.ndarray], coord_range: tuple[int, int, int], elec_names: list[str], do_warn: bool):
    """Ports generateElecMask.m: rasterizes each electrode's point cloud
    into a single labeled volume (voxel value = 1-based electrode index),
    erroring the same way ROAST does on out-of-bounds or overlapping
    electrodes.
    """
    vol_elec = np.zeros(coord_range, dtype=np.int32)
    vol_is_elec = np.zeros(coord_range, dtype=bool)

    for i, coords in enumerate(elec_all_coord):
        name = elec_names[i]
        if coords is None or coords.shape[0] == 0:
            raise ValueError(
                f"Electrode {name} goes out of image boundary. ROAST cannot proceed "
                "without a properly placed electrode. Please expand the input MRI by "
                "specifying the 'zeroPadding' option."
            )
        idx = np.round(coords).astype(int)
        in_bounds = (
            (idx[:, 0] >= 0)
            & (idx[:, 0] < coord_range[0])
            & (idx[:, 1] >= 0)
            & (idx[:, 1] < coord_range[1])
            & (idx[:, 2] >= 0)
            & (idx[:, 2] < coord_range[2])
        )
        if not np.any(in_bounds):
            raise ValueError(
                f"Electrode {name} goes out of image boundary. ROAST cannot proceed "
                "without a properly placed electrode. Please expand the input MRI by "
                "specifying the 'zeroPadding' option."
            )
        if not np.all(in_bounds) and do_warn:
            import warnings

            warnings.warn(
                f"Part of the electrode {name} goes out of image boundary. ROAST can "
                "continue but results may not be accurate. It is recommended that you "
                "expand the input MRI by specifying the 'zeroPadding' option."
            )
        idx = idx[in_bounds]
        flat = np.ravel_multi_index((idx[:, 0], idx[:, 1], idx[:, 2]), coord_range)
        flat_unique, first_pos = np.unique(flat, return_index=True)
        if np.any(vol_is_elec.ravel()[flat_unique]):
            prev_name = elec_names[i - 1] if i > 0 else name
            raise ValueError(
                f"Electrode {name} overlaps with Electrode {prev_name}. ROAST cannot "
                "continue as overlapping electrodes will confuse ROAST when setting up "
                "the boundary conditions for the model. To avoid overlapping, please do "
                "not place two electrodes too close to each other, or reduce the size of "
                "any neighboring electrodes."
            )
        vol_elec.ravel()[flat_unique] = i + 1
        vol_is_elec.ravel()[flat_unique] = True

    return vol_elec


def place_and_model_electrodes(
    elec_loc: np.ndarray,
    elec_range: np.ndarray,
    scalp_clean_surface: np.ndarray,
    scalp_filled: np.ndarray,
    elec_names: list[str],
    elec_paras: list[ElectrodeParams],
    resolution: float,
):
    """Ports placeAndModelElectrodes.m (isDebug plotting omitted -- this is
    a headless port; see roast_py's later visualization phase for
    reviewRes.m instead)."""
    pad_h = np.array(
        [p.elec_size[2] / resolution if p.elec_type.lower() == "pad" else 0.0 for p in elec_paras]
    )
    ind_pad = np.nonzero(pad_h > 0)[0]
    gel_layer: dict[float, np.ndarray] = {}
    elec_layer: dict[float, np.ndarray] = {}
    ind2_layer = np.zeros(len(elec_paras), dtype=float)
    if len(ind_pad) > 0:
        all_ph = np.unique(pad_h[ind_pad])
        for ph in all_ph:
            se = np.ones((max(round(ph), 1),) * 3, dtype=bool)
            gel_pts, scalp_dilated = mask2edge_point_cloud(scalp_filled, "dilate", se)
            elec_pts, _ = mask2edge_point_cloud(scalp_dilated, "dilate", se)
            gel_layer[ph] = gel_pts
            elec_layer[ph] = elec_pts
        for i in ind_pad:
            ind2_layer[i] = pad_h[i]

    nx, ny, nz = scalp_filled.shape
    scalp_filled = scalp_filled.copy()
    scalp_filled[:, :, 0] = False
    scalp_filled[:, :, nz - 1] = False
    scalp_filled[:, 0, :] = False
    scalp_filled[:, ny - 1, :] = False
    scalp_filled[0, :, :] = False
    scalp_filled[nx - 1, :, :] = False

    elec_all_coord: list[np.ndarray] = [None] * len(elec_paras)
    gel_all_coord: list[np.ndarray] = [None] * len(elec_paras)

    for i, para in enumerate(elec_paras):
        lcl = scalp_clean_surface[elec_range[i]]
        cov = np.cov(lcl, rowvar=False)
        eigvals, eigvecs = np.linalg.eigh(cov)
        normal = eigvecs[:, np.argmin(eigvals)]
        normal = normal / np.linalg.norm(normal)

        len_try = 1
        while True:
            test_in = np.round(elec_loc[i] - len_try * normal).astype(int)
            test_out = np.round(elec_loc[i] + len_try * normal).astype(int)
            both = np.stack([test_in, test_out])
            if not (np.all(both.min(axis=0) >= 0) and np.all(both.max(axis=0) <= np.array(scalp_filled.shape) - 1)):
                break
            in_val = scalp_filled[test_in[0], test_in[1], test_in[2]]
            out_val = scalp_filled[test_out[0], test_out[1], test_out[2]]
            if bool(in_val) != bool(out_val):
                break
            len_try += 1
        both = np.stack(
            [
                np.round(elec_loc[i] - len_try * normal).astype(int),
                np.round(elec_loc[i] + len_try * normal).astype(int),
            ]
        )
        if np.all(both.min(axis=0) >= 0) and np.all(both.max(axis=0) <= np.array(scalp_filled.shape) - 1):
            test_in = both[0]
            if not scalp_filled[test_in[0], test_in[1], test_in[2]]:
                normal = -normal

        kind = para.elec_type.lower()
        den = 2
        if kind == "pad":
            pad_length = para.elec_size[0] / resolution
            pad_width = para.elec_size[1] / resolution
            ori = para.elec_ori
            if isinstance(ori, str):
                ori = _ELEC_ORI_KEYWORDS[ori.lower()]
            ori = np.asarray(ori, dtype=float)

            pad_ori_short = np.cross(normal, ori)
            pad_ori_short = pad_ori_short / np.linalg.norm(pad_ori_short)
            pad_ori_long = np.cross(normal, pad_ori_short)
            pad_ori_long = pad_ori_long / np.linalg.norm(pad_ori_long)

            dim_try = np.mean([pad_length, pad_width])
            pad_coor = draw_cuboid(elec_loc[i], [pad_length, pad_width, dim_try], pad_ori_long, pad_ori_short, normal, den)
            pad_coor = np.unique(np.round(pad_coor), axis=0)

            ph = ind2_layer[i]
            gel_coor = _intersect_rows(pad_coor, gel_layer[ph])
            elec_coor = _intersect_rows(pad_coor, elec_layer[ph])

        elif kind == "disc":
            disc_radius = para.elec_size[0] / resolution
            disc_height = para.elec_size[1] / resolution

            gel_out = elec_loc[i] + 2 * disc_height * normal
            electrode = gel_out + disc_height * normal
            gel_in = gel_out - disc_radius * normal

            gel_coor = np.unique(np.round(draw_cylinder(0, disc_radius, gel_in, gel_out, den)), axis=0)
            elec_coor = np.unique(np.round(draw_cylinder(0, disc_radius, gel_out, electrode, den)), axis=0)

        elif kind == "ring":
            ring_radius_in = para.elec_size[0] / resolution
            ring_radius_out = para.elec_size[1] / resolution
            ring_height = para.elec_size[2] / resolution

            dim_try = np.mean([ring_radius_out, ring_radius_in])
            gel_out = elec_loc[i] + 2 * ring_height * normal
            electrode = gel_out + ring_height * normal
            gel_in = gel_out - dim_try * normal

            gel_coor = np.unique(np.round(draw_cylinder(ring_radius_in, ring_radius_out, gel_in, gel_out, den)), axis=0)
            elec_coor = np.unique(np.round(draw_cylinder(ring_radius_in, ring_radius_out, gel_out, electrode, den)), axis=0)
        else:
            raise ValueError(f"Unknown electrode type {para.elec_type!r}")

        elec_all_coord[i] = elec_coor
        gel_all_coord[i] = gel_coor

    return elec_all_coord, gel_all_coord


NECK_POOL = ["nk1", "nk2", "nk3", "nk4"]


def classify_electrodes(elec_names: list[str], cap_names: list[str], custom_names: list[str] | None = None):
    """Ports elecPreproc.m's classification of requested electrode names
    into predefined-cap / neck / custom groups. Returns a dict with
    0-based index arrays `ind_p` (into cap_names), `ind_n` (into
    NECK_POOL), `ind_c` (into custom_names), each sorted the way ROAST
    sorts them (by position in their respective pool), plus `order`: the
    positions in `elec_names` each classified electrode came from, in the
    same concatenation order roast_py uses elsewhere (predefined, neck,
    custom) -- mirroring elecPreproc.m's `ind2UI`.
    """
    custom_names = custom_names or []
    lower_neck_pool = NECK_POOL
    ind_p_raw, ind_n_raw, ind_c_raw = [], [], []
    order_p, order_n, order_c = [], [], []

    for pos, name in enumerate(elec_names):
        if name in cap_names:
            ind_p_raw.append(cap_names.index(name))
            order_p.append(pos)
        elif name.lower() in lower_neck_pool:
            ind_n_raw.append(lower_neck_pool.index(name.lower()))
            order_n.append(pos)
        elif "custom" in name.lower():
            if name not in custom_names:
                raise ValueError(f"Unrecognized electrode {name}")
            ind_c_raw.append(custom_names.index(name))
            order_c.append(pos)
        else:
            raise ValueError(f"Unrecognized electrode {name}")

    def _sorted(ind_raw, order_raw):
        ind_raw = np.array(ind_raw, dtype=int)
        order_raw = np.array(order_raw, dtype=int)
        sort_i = np.argsort(ind_raw, kind="stable")
        return ind_raw[sort_i], order_raw[sort_i]

    ind_p, order_p = _sorted(ind_p_raw, order_p)
    ind_n, order_n = _sorted(ind_n_raw, order_n)
    ind_c, order_c = _sorted(ind_c_raw, order_c)

    return {
        "ind_p": ind_p,
        "ind_n": ind_n,
        "ind_c": ind_c,
        "order": np.concatenate([order_p, order_n, order_c]),
    }


def electrode_placement(
    tissue_labels: np.ndarray,
    landmarks: np.ndarray,
    elec_names: list[str],
    elec_paras: list[ElectrodeParams],
    cap_names: list[str],
    cap_template: np.ndarray,
    voxel_size: np.ndarray,
    is_biosemi: bool = False,
    is_egi: bool = False,
    custom_names: list[str] | None = None,
    custom_coords: np.ndarray | None = None,
):
    """Ports electrodePlacement.m: the full pipeline from tissue labels +
    landmarks + a requested electrode montage to labeled electrode/gel
    voxel masks.

    `tissue_labels` is ROAST's 6-tissue label volume (see
    roast_py.segmentation). `elec_names`/`elec_paras` are given in
    whatever order the caller wants electrodes in -- internally this
    reorders into predefined/neck/custom, pool-sorted order the way
    MATLAB's electrodePlacement.m expects its own (pre-sorted-by-the-caller)
    inputs, then relabels the result back to `elec_names`' original order
    before returning, so a caller never has to think about that detour.

    Returns (elec_mask, gel_mask): uint8 volumes, 0 = not electrode/gel,
    i+1 = the i-th electrode in `elec_names`.
    """
    scalp = tissue_labels > 0
    scalp_surface, _ = mask2edge_point_cloud(scalp, "erode", np.ones((3, 3, 3), dtype=bool))

    classes = classify_electrodes(elec_names, cap_names, custom_names)
    ind_p, ind_n, ind_c = classes["ind_p"], classes["ind_n"], classes["ind_c"]

    if len(ind_p) > 0:
        electrode_coord_p, center_p = fit_cap_to_individual(
            scalp, scalp_surface, landmarks, voxel_size, cap_names, cap_template, ind_p, is_biosemi, is_egi
        )
    else:
        electrode_coord_p, center_p = np.empty((0, 3)), None

    if len(ind_n) > 0:
        if np.any(landmarks[4:6, 2] <= 0):
            raise ValueError(
                "MRI does not cover the neck, so cannot place electrodes on the neck."
            )
        electrode_coord_n, center_n = place_neck_electrodes(scalp.shape, scalp_surface, landmarks, ind_n)
    else:
        electrode_coord_n, center_n = np.empty((0, 3)), None

    if len(ind_c) > 0:
        if custom_coords is None:
            raise ValueError("custom_coords is required when custom electrodes are requested")
        elec_loc_c = custom_coords[ind_c]
        _, ind_on_scalp = map2points(elec_loc_c, scalp_surface, "closest")
        electrode_coord_c = scalp_surface[ind_on_scalp]
    else:
        electrode_coord_c = np.empty((0, 3))

    scalp_clean, scalp_filled = clean_scalp(scalp, scalp_surface)
    scalp_clean_surface, _ = mask2edge_point_cloud(scalp_clean, "erode", np.ones((3, 3, 3), dtype=bool))

    elec_ranges = []
    if len(ind_p) > 0:
        _, ind_on_scalp = project2closest_surface_points(electrode_coord_p, scalp_clean_surface, center_p)
        elec_ranges.append(ind_on_scalp[:100].T)
    if len(ind_n) > 0:
        _, ind_on_scalp = project2closest_surface_points(electrode_coord_n, scalp_clean_surface, center_n)
        elec_ranges.append(ind_on_scalp[:100].T)
    if len(ind_c) > 0:
        _, ind_on_scalp = map2points(electrode_coord_c, scalp_clean_surface, "closer", 100)
        elec_ranges.append(ind_on_scalp.T)
    elec_range = np.concatenate(elec_ranges, axis=0) if elec_ranges else np.empty((0, 100), dtype=int)

    electrode_coord = np.concatenate([electrode_coord_p, electrode_coord_n, electrode_coord_c], axis=0)

    # electrode_coord (and elec_range) are in predefined/neck/custom,
    # pool-sorted order (matching ROAST's own internal ind2UI reordering
    # before calling electrodePlacement.m) -- not necessarily the order
    # `elec_names`/`elec_paras` were given in, so those need the same
    # permutation applied before they're zipped together below.
    order = classes["order"]
    elec_names_sorted = [elec_names[i] for i in order]
    elec_paras_sorted = [elec_paras[i] for i in order]

    resolution = float(np.mean(voxel_size))
    elec_all_coord, gel_all_coord = place_and_model_electrodes(
        electrode_coord, elec_range, scalp_clean_surface, scalp_filled, elec_names_sorted, elec_paras_sorted, resolution
    )

    coord_range = scalp.shape
    vol_elec_c = generate_elec_mask(elec_all_coord, coord_range, elec_names_sorted, True)
    vol_gel_c = generate_elec_mask(gel_all_coord, coord_range, elec_names_sorted, False)

    vol_elec = vol_elec_c > 0
    vol_gel = vol_gel_c > 0
    vol_gel = vol_gel & ~vol_elec  # remove gel overlapping the electrode
    for tissue in range(1, 7):
        vol_gel = vol_gel & ~(tissue_labels == tissue)  # remove gel that goes into other tissue

    elec_mask_sorted = (vol_elec_c * vol_elec).astype(np.uint8)
    gel_mask_sorted = (vol_gel_c * vol_gel).astype(np.uint8)

    # Relabel back to `elec_names`' original order, so a caller sees
    # elec_mask == i+1 for elec_names[i], regardless of the pool-sorting
    # detour above.
    sorted_to_original = np.zeros(len(order) + 1, dtype=np.uint8)  # index 0 stays background
    for sorted_pos, original_pos in enumerate(order):
        sorted_to_original[sorted_pos + 1] = original_pos + 1
    elec_mask = sorted_to_original[elec_mask_sorted]
    gel_mask = sorted_to_original[gel_mask_sorted]
    return elec_mask, gel_mask


def _intersect_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if a.size == 0 or b.size == 0:
        n_cols = a.shape[1] if a.ndim == 2 else 3
        return np.empty((0, n_cols), dtype=a.dtype if a.size else float)
    view_dtype = {"names": [f"f{i}" for i in range(a.shape[1])], "formats": a.shape[1] * [a.dtype]}
    a_view = np.ascontiguousarray(a).view(view_dtype)
    b_view = np.ascontiguousarray(b).view(view_dtype)
    inter = np.intersect1d(a_view, b_view)
    return inter.view(a.dtype).reshape(-1, a.shape[1])
