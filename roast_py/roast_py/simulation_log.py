"""Simulation tags and the options log: ports writeRoastLog.m and
isNewOptions.m (their 'roast' branches).

Every roast() simulation is saved under a tag: the user's
`simulation_tag`, or by default a date-time stamp like MATLAB's
`datestr(now,30)` (`20261025T143012`). Its options go to
`<subj>_<tag>_roastOptions.json` (MATLAB: `_roastOptions.mat`), and a
readable entry is appended to `<subj>_roastLog`, both in the work
directory. Before electrode placement, roast() compares the new options
with every saved simulation of the same MRI. If one matches, roast()
carries on under that tag and so picks up its saved outputs; otherwise
the options are saved under a new tag.

What counts as "the same options" is isNewOptions.m's list: the recipe,
cap type, each electrode's type and size (and orientation, for pads), T2,
segmentation (multiaxial or SPM), the MNI registration (`Affine`), mesh
options, conductivities, resampling and zero-padding.
"""

from __future__ import annotations

import glob
import json
import os
import warnings
from datetime import datetime

import numpy as np

# capType values ROAST treats as one cap: 10-05, 10-10 and 10-20 come from
# the same electrode list (isNewOptions.m's capNew/capOld).
_CAP_GROUPS = {"1020": "10-x", "1010": "10-x", "1005": "10-x"}


def timestamp_tag(now: datetime | None = None) -> str:
    """MATLAB's `datestr(now,30)`: 'yyyymmddTHHMMSS'."""
    return (now or datetime.now()).strftime("%Y%m%dT%H%M%S")


def options_path(work_dir: str, base: str, tag: str) -> str:
    return os.path.join(work_dir, f"{base}_{tag}_roastOptions.json")


def log_path(work_dir: str, base: str) -> str:
    return os.path.join(work_dir, f"{base}_roastLog")


def check_tag(tag: str | None) -> None:
    if tag is None:
        return
    if not isinstance(tag, str) or not tag or any(c in tag for c in '/\\:*?"<>|'):
        raise ValueError(
            f"Invalid simulation tag {tag!r}: it becomes part of the output file names, so use a "
            "non-empty name without / \\ : * ? \" < > |."
        )


def _same_numbers(a, b) -> bool:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return a.shape == b.shape and bool(np.allclose(a, b, rtol=1e-9, atol=1e-12))


def _same_orientation(a, b) -> bool:
    if isinstance(a, str) and isinstance(b, str):
        return a.lower() == b.lower()
    if a is None or b is None or isinstance(a, str) or isinstance(b, str):
        return False  # isNewOptions.m: a keyword vs a vector counts as new
    return _same_numbers(a, b)


def is_new_options(new: dict, old: dict) -> bool:
    """isNewOptions.m: True unless `old` (a saved simulation's options) is
    the same simulation as `new`. Options saved by an older roast_py that
    lack a field count as new."""
    try:
        if list(new["recipe"]) != list(old["recipe"]):
            return True
        if not _same_numbers(list(new["recipe"].values()), list(old["recipe"].values())):
            return True
        if _CAP_GROUPS.get(new["cap_type"], new["cap_type"]) != _CAP_GROUPS.get(old["cap_type"], old["cap_type"]):
            return True
        if len(new["electrodes"]) != len(old["electrodes"]):
            return True
        for e_new, e_old in zip(new["electrodes"], old["electrodes"]):
            if e_new["elec_type"].lower() != e_old["elec_type"].lower():
                return True
            if not _same_numbers(e_new["elec_size"], e_old["elec_size"]):
                return True
            if e_new["elec_type"].lower() == "pad" and not _same_orientation(e_new["elec_ori"], e_old["elec_ori"]):
                return True
        if (new["t2"] or None) != (old["t2"] or None):
            return True
        if new["segmentation"] != old["segmentation"]:
            return True
        if not _same_numbers(new["Affine"], old["Affine"]):
            return True
        if set(new["mesh_options"]) != set(old["mesh_options"]):
            return True
        if any(not _same_numbers(v, old["mesh_options"][k]) for k, v in new["mesh_options"].items()):
            return True
        if set(new["conductivities"]) != set(old["conductivities"]):
            return True
        if any(not _same_numbers(v, old["conductivities"][k]) for k, v in new["conductivities"].items()):
            return True
        if bool(new["resampling"]) != bool(old["resampling"]):
            return True
        if int(new["zero_padding"]) != int(old["zero_padding"]):
            return True
    except (KeyError, TypeError, AttributeError, ValueError):
        return True
    return False


def saved_simulations(work_dir: str, base: str, subj: str) -> list[dict]:
    """The options of every simulation saved for MRI `subj` (named `base`)
    in `work_dir`, oldest first. MATLAB's `dir([subjName '_*_roastOptions.mat'])`;
    the stored `subj` keeps another MRI whose name starts with `base_` out."""
    found = []
    for path in glob.glob(os.path.join(glob.escape(work_dir), f"{glob.escape(base)}_*_roastOptions.json")):
        try:
            with open(path) as f:
                opt = json.load(f)
        except (OSError, ValueError):
            continue
        tag = opt.get("simulation_tag") if isinstance(opt, dict) else None
        if tag and opt.get("subj") == subj and os.path.basename(path) == f"{base}_{tag}_roastOptions.json":
            found.append((opt.get("created") or "", os.path.getmtime(path), opt))
    found.sort(key=lambda t: (t[0], t[1]))
    return [opt for _, _, opt in found]


def latest_simulation(work_dir: str, base: str, subj: str) -> dict | None:
    sims = saved_simulations(work_dir, base, subj)
    return sims[-1] if sims else None


def register_simulation(options: dict, work_dir: str, base: str) -> tuple[dict, bool]:
    """roast.m's log tracking. `options["simulation_tag"]` is the user's
    tag or None. Returns (options with the tag to use, reused), where
    `reused` means a saved simulation had the same options and its tag is
    carried on; otherwise the options were just saved under a new tag and
    logged."""
    tag = options.get("simulation_tag")
    check_tag(tag)
    for old in saved_simulations(work_dir, base, options["subj"]):
        if not is_new_options(options, old):
            if tag and tag != old["simulation_tag"]:
                warnings.warn(
                    f"The simulation with the same options has been run before under tag "
                    f"'{old['simulation_tag']}'. The new tag you specified '{tag}' will be ignored.",
                    stacklevel=3,
                )
            return dict(options, simulation_tag=old["simulation_tag"], created=old.get("created")), True

    if tag:
        if os.path.exists(options_path(work_dir, base, tag)):
            raise ValueError(
                "You're about to run a simulation using options that you never used before, but forgot "
                f"to use a new tag for it: tag '{tag}' is already used by another simulation of {base}. "
                "ROAST will get confused when managing different simulations with the same tag. "
                "Please use a new tag."
            )
    else:
        tag = timestamp_tag()
        stamp, n = tag, 2
        while os.path.exists(options_path(work_dir, base, tag)):  # two runs within one second
            tag, n = f"{stamp}-{n}", n + 1
    options = dict(options, simulation_tag=tag, created=datetime.now().isoformat(timespec="seconds"))

    path = options_path(work_dir, base, tag)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(options, f, indent=2)
    os.replace(tmp, path)
    with open(log_path(work_dir, base), "a") as f:
        f.write(log_entry(options))
    return options, False


def _vector(values, fmt: str) -> str:
    return "[" + ",".join(fmt % v for v in np.ravel(values)) + "]"


def log_entry(opt: dict) -> str:
    """writeRoastLog.m's entry for one simulation."""
    lines = [f"{opt['simulation_tag']}:", f"recipe:\t{opt['config_txt']}", f"capType:\t{opt['cap_type']}"]
    lines.append("elecType:\t" + "\t".join(e["elec_type"] for e in opt["electrodes"]))
    lines.append("elecSize:\t" + "; ".join(_vector(e["elec_size"], "%.1f") for e in opt["electrodes"]))
    oris = []
    for e in opt["electrodes"]:
        ori = e["elec_ori"]
        oris.append("N/A" if ori is None else ori if isinstance(ori, str) else _vector(ori, "%.4f"))
    lines.append("elecOri:\t" + "; ".join(oris))
    lines.append(f"T2:\t{opt['t2'] or 'none'}")
    lines.append(f"multiaxial:\t{'on' if opt['segmentation'] == 'multiaxial' else 'off'}")
    m = opt["mesh_options"]
    lines.append(
        f"meshOpt:\tradbound: {m['radbound']:g}; angbound: {m['angbound']:g}; distbound: {m['distbound']:.1f}; "
        f"reratio: {m['reratio']:g}; maxvol: {m['maxvol']:g}"
    )
    c = opt["conductivities"]
    lines.append(
        f"conductivities:\twhite: {c['white']:.3f}; gray: {c['gray']:.3f}; CSF: {c['csf']:.3f}; "
        f"bone: {c['bone']:.3f}; skin: {c['skin']:.3f}; air: {c['air']:.1e}; "
        + "gel: " + "".join(f"{g:.3f}; " for g in c["gel"])
        + "electrode: " + "".join(f"{g:.1e}; " for g in c["electrode"])
    )
    lines.append(f"reSampling:\t{'on' if opt['resampling'] else 'off'}")
    lines.append(f"zeroPadding:\t{opt['zero_padding'] if opt['zero_padding'] > 0 else 'none'}")
    lines.append(f"original MRI in RAS?:\t{'no, and is re-oriented into RAS' if opt.get('is_non_ras') else 'yes'}")
    modeled = f"MRI that was modeled by ROAST:\t{opt['model_mri']}"
    if opt["t2"]:
        modeled += f"\n                         \t{opt['t2']}"
    lines.append(modeled)
    lines.append(f"mri2mni matrix:\t {'[' + '; '.join(','.join(f'{v:.1f}' for v in row) for row in opt['mri2mni']) + ']'}")
    return "\n".join(lines) + "\n\n"
