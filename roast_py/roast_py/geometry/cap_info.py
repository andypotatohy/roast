"""Reads capInfo.xlsx -- the standard EEG electrode template layouts
(10-05, BioSemi, EGI) bundled at the ROAST repo root -- replacing
`readtable('capInfo.xlsx','Sheet',...)`.

The sheets have no header row (row 1 is real data -- the first electrode,
e.g. 'LPA' for 10-05), so this reads with header=None to match ROAST's
own usage (`capInfo = table2cell(readtable(...))`, using every row).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

_SHEET_BY_CAP_TYPE = {
    "1020": "10-05",
    "1010": "10-05",
    "1005": "10-05",
    "biosemi": "BioSemi",
    "egi": "EGI",
}


def find_default_cap_info() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "capInfo.xlsx"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the bundled capInfo.xlsx. Pass cap_info_path explicitly."
    )


def load_cap_info(cap_type: str, cap_info_path: str | Path | None = None):
    """Returns (electrode_names, template_xyz) for the given cap type
    ('1020'/'1010'/'1005'/'biosemi'/'egi'). `template_xyz` is an (N, 3)
    array of the normalized template coordinates from the spreadsheet.
    """
    cap_type = cap_type.lower()
    if cap_type not in _SHEET_BY_CAP_TYPE:
        raise ValueError(f"Unknown capType {cap_type!r}")

    path = Path(cap_info_path) if cap_info_path else find_default_cap_info()
    df = pd.read_excel(path, sheet_name=_SHEET_BY_CAP_TYPE[cap_type], header=None)
    names = df[0].tolist()
    template = df[[1, 2, 3]].to_numpy(dtype=float)
    return names, template
