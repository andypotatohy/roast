"""Ports roast.m's option checking: turns roast()'s keyword arguments into
the per-electrode parameters, conductivities and switches the pipeline
uses, with roast.m's defaults, validation and (where it makes sense) its
error and warning messages.

Python spellings of MATLAB's values:

* recipe: a dict ``{"Fp1": 1, "P4": -1}``, or MATLAB's flat form
  ``["Fp1", 1, "P4", -1]``.
* elec_type: ``"disc"`` or a list, one entry per electrode.
* elec_size: one size ``[6, 2]``; per electrode for a single type, an
  N-row matrix ``[[6, 2], [8, 2]]``; for mixed types, a list with one
  entry per electrode where ``None`` means that type's default (MATLAB's
  cell array with ``[]`` entries).
* elec_ori (pads only): ``"lr"``/``"ap"``/``"si"``, a direction vector,
  an N-by-3 matrix, or a list mixing keywords, ``None`` and vectors
  (MATLAB's cell array).
* conductivities: a dict with any of roast.m's field names (missing ones
  take the defaults), or a ``Conductivities``.
* 'on'/'off' switches: ``True``/``False`` (or the strings).
"""

from __future__ import annotations

import warnings

import numpy as np

from .fem.pro_writer import Conductivities
from .geometry.placement import ElectrodeParams

CAP_TYPES = ("1020", "1010", "1005", "biosemi", "egi")
ELEC_TYPES = ("disc", "pad", "ring")
PAD_ORIENTATIONS = ("lr", "ap", "si")
DEFAULT_ELEC_SIZE = {"disc": (6.0, 2.0), "pad": (50.0, 30.0, 3.0), "ring": (4.0, 6.0, 2.0)}
DEFAULT_CONDUCTIVITIES = {
    "white": 0.126, "gray": 0.276, "csf": 1.65, "bone": 0.01,
    "skin": 0.465, "air": 2.5e-14, "gel": 0.3, "electrode": 5.9e7,
}

_BAD_ORI = (
    "Unrecognized pad orientation. Please enter 'lr', 'ap', or 'si' for pad orientation; "
    "or just enter the direction vector of the long axis of the pad"
)


def on_off(value, name: str) -> bool:
    """roast.m's 'on'/'off' option values; Python booleans also accepted."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in ("on", "off"):
        return value.lower() == "on"
    raise ValueError(f"Unrecognized option value. Please enter 'on' or 'off' for option '{name}'.")


def resolve_recipe(recipe) -> tuple[list[str], list[float]]:
    """(electrode names, injected currents in mA), in the order given."""
    if isinstance(recipe, str):
        if recipe.lower() == "leadfield":
            raise NotImplementedError(
                "recipe='leadfield' (lead field generation for roast_target) is not ported "
                "to roast_py yet; roast_target() isn't either."
            )
        raise ValueError("Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair.")
    if isinstance(recipe, dict):
        names, currents = list(recipe.keys()), list(recipe.values())
    else:
        recipe = list(recipe)
        if len(recipe) % 2:
            raise ValueError(
                "Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair."
            )
        names, currents = recipe[0::2], recipe[1::2]
    currents = [float(c) for c in currents]
    if abs(sum(currents)) > 1e-9:
        raise ValueError(
            "Electric currents going in and out of the head not balanced. Please make sure they sum to 0."
        )
    return [str(n) for n in names], currents


def _is_numeric(value) -> bool:
    if value is None or isinstance(value, str):
        return False
    try:
        np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return False
    return True


def _check_one_size(kind: str, size: np.ndarray) -> None:
    """roast.m's per-type size checks, for an (n, k) array of sizes."""
    if np.any(size <= 0):
        raise ValueError("Please enter non-negative values for electrode size.")
    if size.shape[1] not in (2, 3):
        raise ValueError(
            "Unrecognized electrode sizes. Please specify as [radius height] for disc, "
            "[length width height] for pad, and [innerRadius outterRadius height] for ring electrode."
        )
    if kind == "disc" and size.shape[1] == 3:
        raise ValueError("Redundant size info for Disc electrodes. Please enter as [radius height]")
    if kind in ("pad", "ring") and size.shape[1] == 2:
        raise ValueError(
            "Insufficient size info for Pad or Ring electrodes. Please specify as [length width height] "
            "for pad, and [innerRadius outterRadius height] for ring electrode."
        )
    if kind == "pad":
        if np.any(size[:, 0] < size[:, 1]):
            raise ValueError(
                "For Pad electrodes, the width of the pad should not be bigger than its length. "
                "Please enter as [length width height]"
            )
        if np.any(size[:, 2] < 3):
            raise ValueError("For Pad electrodes, the thickness should at least be 3 mm.")
        if np.any(size > 80):
            warnings.warn(
                "You're placing large pad electrodes (one of its dimensions is bigger than 8 cm). For "
                "large pads, the size will not be exact in the model because they will be bent to fit "
                "the scalp surface.",
                stacklevel=4,
            )
    if kind == "ring" and np.any(size[:, 0] >= size[:, 1]):
        raise ValueError(
            "For Ring electrodes, the inner radius should be smaller than outter radius. Please enter "
            "as [innerRadius outterRadius height]"
        )


def _check_ori_value(ori):
    if isinstance(ori, str):
        if ori.lower() not in PAD_ORIENTATIONS:
            raise ValueError(_BAD_ORI)
        return ori.lower()
    arr = np.atleast_2d(np.asarray(ori, dtype=float))
    if arr.shape[1] != 3:
        raise ValueError(_BAD_ORI)
    return arr


def resolve_electrodes(
    n: int, cap_type: str = "1010", elec_type="disc", elec_size=None, elec_ori=None
) -> list[ElectrodeParams]:
    """One ElectrodeParams per electrode (recipe order), porting roast.m's
    capType/elecType/elecSize/elecOri handling."""
    cap_type = str(cap_type).lower()
    if cap_type not in CAP_TYPES:
        raise ValueError("Supported cap types are: '1020', '1010', '1005', 'BioSemi' and 'EGI'.")

    # -- types ---------------------------------------------------------------
    mixed = not isinstance(elec_type, str)
    if mixed:
        types = [str(t).lower() for t in elec_type]
        if len(types) != n:
            raise ValueError(
                "You want to place more than 1 type of electrodes, but did not tell ROAST which type for "
                "each electrode. Please provide the type for each electrode respectively, as the value for "
                "option 'elecType', in a list of length equals to the number of electrodes to be placed."
            )
    else:
        types = [elec_type.lower()] * n
    if any(t not in ELEC_TYPES for t in types):
        raise ValueError("Supported electrodes are: 'disc', 'pad' and 'ring'.")

    # -- sizes ---------------------------------------------------------------
    if elec_size is None:
        sizes = [np.array(DEFAULT_ELEC_SIZE[t]) for t in types]
    elif not mixed:
        if not _is_numeric(elec_size):  # MATLAB: a cell array with one type
            warnings.warn(
                "Looks like you're placing only 1 type of electrodes. ROAST will only use the 1st entry of "
                "'elecSize'. If you meant different sizes for different electrodes of the same type, enter "
                "'elecSize' as an N-by-2 or N-by-3 matrix, where N is number of electrodes to be placed.",
                stacklevel=3,
            )
            elec_size = list(elec_size)[0]
        size = np.atleast_2d(np.asarray(elec_size, dtype=float))
        _check_one_size(types[0], size)
        if size.shape[0] > 1 and size.shape[0] != n:
            raise ValueError(
                "You want different sizes for each electrode. Please tell ROAST the size for each electrode "
                "respectively, in a N-row matrix, where N is the number of electrodes to be placed."
            )
        sizes = [size[i if size.shape[0] > 1 else 0] for i in range(n)]
    else:
        if isinstance(elec_size, np.ndarray) or (_is_numeric(elec_size) and np.ndim(elec_size) == 1):
            raise ValueError(
                "You want to place at least 2 types of electrodes, but only provided size info for 1 type. "
                "Please provide complete size info for all types of electrodes in a list as the value for "
                "option 'elecSize', or just use defaults by not specifying 'elecSize' option."
            )
        elec_size = list(elec_size)
        if len(elec_size) != n:
            raise ValueError(
                "You want to place more than 1 type of electrodes. Please tell ROAST the size for each "
                "electrode respectively, as the value for option 'elecSize', in a list of length equals "
                "to the number of electrodes to be placed."
            )
        sizes = []
        for kind, entry in zip(types, elec_size):
            if entry is None or (np.size(entry) == 0):
                sizes.append(np.array(DEFAULT_ELEC_SIZE[kind]))
                continue
            size = np.atleast_2d(np.asarray(entry, dtype=float))
            if size.shape[0] > 1:
                raise ValueError(
                    "You're placing more than 1 type of electrodes. Please put size info for each electrode "
                    "as a 1-row vector in the list for option 'elecSize'."
                )
            _check_one_size(kind, size)
            sizes.append(size[0])

    # -- pad orientations ------------------------------------------------------
    pads = [i for i, t in enumerate(types) if t == "pad"]
    oris: list = ["lr" if t == "pad" else None for t in types]
    if elec_ori is not None:
        if not mixed and types[0] != "pad":
            warnings.warn(
                "You're not placing pad electrodes; customized orientation options will be ignored.",
                stacklevel=3,
            )
        elif not mixed:
            if not isinstance(elec_ori, str) and not _is_numeric(elec_ori):  # MATLAB cell
                entries = list(elec_ori)
                if all(isinstance(e, str) for e in entries):
                    if len(entries) != n:
                        raise ValueError(
                            "You want different orientations for each pad electrode by using pre-defined "
                            "keywords in a list. Please make sure the list has a length equal to the number "
                            "of pad electrodes."
                        )
                    oris = [_check_ori_value(e) for e in entries]
                else:
                    warnings.warn(
                        "Looks like you're only placing pad electrodes. ROAST will only use the 1st entry "
                        "of 'elecOri'. If you meant different orientations for different pads, enter it as "
                        "an N-by-3 matrix, or as a list of 'lr', 'ap' or 'si', with N the number of pads.",
                        stacklevel=3,
                    )
                    elec_ori = entries[0]
            if isinstance(elec_ori, str) or _is_numeric(elec_ori):
                value = _check_ori_value(elec_ori)
                if isinstance(value, str) or value.shape[0] == 1:
                    oris = [value if isinstance(value, str) else value[0]] * n
                elif value.shape[0] == n:
                    oris = list(value)
                else:
                    raise ValueError(
                        "You want different orientations for each pad electrode. Please tell ROAST the "
                        "orientation for each pad respectively, in a N-by-3 matrix, where N is the number of "
                        "pads to be placed."
                    )
        elif isinstance(elec_ori, str) or _is_numeric(elec_ori):  # mixed types, one value or a matrix
            value = _check_ori_value(elec_ori)
            if isinstance(value, str):
                rows = [value] * len(pads)
            elif value.shape[0] == 1:
                rows = [value[0]] * len(pads)
            elif value.shape[0] == len(pads):
                rows = list(value)
            else:
                raise ValueError(
                    "You want different orientations for each pad electrode. Please tell ROAST the orientation "
                    "for each pad respectively, in a N-by-3 matrix, where N is the number of pads to be placed."
                )
            for i, row in zip(pads, rows):
                oris[i] = row
        else:  # mixed types, MATLAB cell: one entry per electrode
            entries = list(elec_ori)
            if len(entries) != n:
                raise ValueError(
                    "You want to place another type of electrodes aside from pad. Please tell ROAST the "
                    "orienation for each electrode respectively, as the value for option 'elecOri', in a list "
                    "of length equals to the number of electrodes to be placed (put None for non-pad "
                    "electrodes)."
                )
            for i in pads:
                entry = entries[i]
                if entry is None or (not isinstance(entry, str) and np.size(entry) == 0):
                    continue  # stays 'lr'
                value = _check_ori_value(entry)
                if not isinstance(value, str):
                    if value.shape[0] > 1:
                        raise ValueError(
                            "You're placing more than 1 type of electrodes. Please put orientation info for "
                            "each pad electrode as a 1-by-3 vector or one of the three keywords 'lr', 'ap', or "
                            "'si' in a list for option 'elecOri'."
                        )
                    value = value[0]
                oris[i] = value

    return [
        ElectrodeParams(elec_type=t, elec_size=np.asarray(s, dtype=float), elec_ori=o, cap_type=cap_type)
        for t, s, o in zip(types, sizes, oris)
    ]


def resolve_conductivities(conductivities, n: int) -> Conductivities:
    """roast.m's conductivity handling: defaults for missing tissues, and
    gel/electrode values expanded to one per electrode (recipe order)."""
    if conductivities is None:
        values = dict(DEFAULT_CONDUCTIVITIES)
    else:
        if isinstance(conductivities, Conductivities):
            given = {k: v for k, v in vars(conductivities).items() if not (isinstance(v, list) and not v)}
        elif isinstance(conductivities, dict):
            given = dict(conductivities)
        else:
            raise TypeError(
                "Unrecognized format of conductivity values. Please enter as a dict, with keys 'white', "
                "'gray', 'csf', 'bone', 'skin', 'air', 'gel' and 'electrode'."
            )
        unknown = set(given) - set(DEFAULT_CONDUCTIVITIES)
        if not given or unknown:
            raise ValueError(
                "Unrecognized tissue names detected. Supported tissue names in the conductivity option are "
                "'white', 'gray', 'csf', 'bone', 'skin', 'air', 'gel' and 'electrode'."
            )
        values = dict(DEFAULT_CONDUCTIVITIES)
        for name, value in given.items():
            arr = np.ravel(np.asarray(value, dtype=float))
            if arr.size == 0 or np.any(arr <= 0):
                raise ValueError(f"Please enter a positive number for the {name} conductivity.")
            if name in ("gel", "electrode"):
                if arr.size > 1 and arr.size != n:
                    raise ValueError(
                        f"You want to assign different {name} conductivities under different electrodes, "
                        "but didn't tell ROAST clearly which conductivity each electrode should use. Please "
                        f"follow the order of electrodes in 'recipe', one value per electrode, for '{name}'."
                    )
                values[name] = arr.tolist()
            else:
                if arr.size > 1:
                    raise ValueError(
                        "Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity."
                    )
                values[name] = float(arr[0])
        if conductivities is not None and not isinstance(conductivities, Conductivities):
            warnings.warn(
                "You're changing the advanced options of ROAST. Unless you know what you're doing, please "
                "keep conductivity values default.",
                stacklevel=3,
            )
    for name in ("gel", "electrode"):
        v = np.ravel(np.asarray(values[name], dtype=float))
        values[name] = (np.repeat(v, n) if v.size == 1 else v).tolist()
    return Conductivities(**values)


def read_custom_locations(subj: str) -> tuple[list[str], np.ndarray]:
    """`<subj>_customLocations` next to the MRI: lines of `name x y z`, in
    MATLAB's 1-based voxel coordinates of that MRI."""
    import os

    dirname, filename = os.path.split(os.path.abspath(subj))
    base = filename.removesuffix(".gz").removesuffix(".nii").removesuffix(".hdr").removesuffix(".img")
    path = os.path.join(dirname, f"{base}_customLocations")
    if not os.path.exists(path):
        raise FileNotFoundError(
            "You specified customized electrode locations but did not provide the location coordinates. "
            f"Please put together all the coordinates in a text file '{base}_customLocations' and store it "
            "under the subject folder."
        )
    names, coords = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                names.append(parts[0])
                coords.append([float(p) for p in parts[1:4]])
    return names, np.array(coords, dtype=float).reshape(-1, 3)
