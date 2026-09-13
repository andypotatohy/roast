"""Ports cgalv2m.m and meshByIso2mesh.m.

The one open technical risk flagged in the original translation plan --
that iso2mesh's `cgalmesh` mesher ships as a MATLAB MEX file on Linux/Mac,
not a standalone binary -- turned out not to be real: despite the
`.mexa64`/`.mexmaci64` filenames (iso2mesh's own naming convention,
apparently reused from its actual MEX tools), `lib/iso2mesh/bin/cgalmesh.mexa64`
is a plain, statically-linked, directly-executable ELF binary (verified
with `file`/`readelf`: `ELF 64-bit LSB EXECUTABLE ... statically linked`,
and it runs and prints a usage message when invoked directly). MATLAB's
own cgalv2m.m already calls it via `system()`, exactly the same style as
the getdp/NiftyReg subprocess wrappers elsewhere in this port -- so this
module does the same, with no fallback library needed.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np

from .mesh_io import save_inr, save_msh, sort_mesh, read_medit

# Ports cgalv2m.m's defaults.
_DEFAULT_ANG = 30
_DEFAULT_SSIZE = 6
_DEFAULT_APPROX = 0.5
_DEFAULT_RERATIO = 3


def find_cgalmesh_binary(bin_path: str | os.PathLike | None = None) -> Path:
    if bin_path is not None:
        return Path(bin_path)
    here = Path(__file__).resolve()
    for parent in here.parents:
        for candidate_name in ("cgalmesh.mexa64", "cgalmesh.mexmaci64", "cgalmesh.exe", "cgalmesh"):
            candidate = parent / "lib" / "iso2mesh" / "bin" / candidate_name
            if candidate.exists():
                return candidate
    raise FileNotFoundError(
        "Could not locate the bundled cgalmesh binary under lib/iso2mesh/bin/. "
        "Pass bin_path explicitly."
    )


def run_cgal_mesher(
    vol: np.ndarray,
    work_dir: str | os.PathLike,
    radbound: float = _DEFAULT_SSIZE,
    angbound: float = _DEFAULT_ANG,
    distbound: float = _DEFAULT_APPROX,
    reratio: float = _DEFAULT_RERATIO,
    maxvol: float = 10.0,
    randseed: int = 0x623F9A9E,
    bin_path: str | os.PathLike | None = None,
):
    """Ports cgalv2m.m: writes `vol` (a uint8 multi-domain label volume) to
    INR format, runs the cgalmesh binary, and reads back the resulting
    surface+volume mesh.

    Returns (node, elem, face) with node offset by +0.5 (ports cgalv2m.m's
    own `node=node+0.5` at the very end -- CGAL's mesher outputs voxel-
    corner-aligned coordinates; this recenters them to voxel centers) and
    sorted by roast_py.meshing.mesh_io.sort_mesh for cache locality, same
    as the original.
    """
    if vol.dtype != np.uint8 and vol.dtype != np.bool_:
        raise ValueError("cgalmesher can only handle uint8 volumes")
    if not np.any(vol):
        raise ValueError("no labeled regions found in the input volume")

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    inr_path = work_dir / "pre_cgalmesh.inr"
    mesh_path = work_dir / "post_cgalmesh.mesh"
    mesh_path.unlink(missing_ok=True)

    save_inr(vol, str(inr_path))

    binary = find_cgalmesh_binary(bin_path)
    cmd = [
        str(binary),
        str(inr_path),
        str(mesh_path),
        f"{angbound:f}",
        f"{radbound:f}",
        f"{distbound:f}",
        f"{reratio:f}",
        f"{maxvol:f}",
        str(randseed),
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    if not mesh_path.exists():
        raise RuntimeError(f"cgalmesh did not produce an output file; command was: {' '.join(cmd)}")

    node, elem, face = read_medit(str(mesh_path))

    if node.shape[0] > 0:
        node, elem, face = sort_mesh(node[0, :3], node, elem, [0, 1, 2, 3], face, [0, 1, 2])

    node = node.copy()
    node[:, :3] = node[:, :3] + 0.5
    return node, elem, face


def mesh_by_iso2mesh(
    tissue_labels: np.ndarray,
    elec_mask: np.ndarray,
    gel_mask: np.ndarray,
    voxel_size: np.ndarray,
    work_dir: str | os.PathLike,
    out_path: str | os.PathLike,
    radbound: float = _DEFAULT_SSIZE,
    angbound: float = _DEFAULT_ANG,
    distbound: float = _DEFAULT_APPROX,
    reratio: float = _DEFAULT_RERATIO,
    maxvol: float = 10.0,
    bin_path: str | os.PathLike | None = None,
):
    """Ports meshByIso2mesh.m: combines the 6-tissue label volume with the
    per-electrode gel/electrode masks (roast_py.geometry.placement's
    output) into one multi-domain volume, meshes it, converts node
    coordinates from voxel-corner to physical (mm) space, and writes a
    Gmsh .msh file.

    Region ids in the output mesh, matching meshByIso2mesh.m exactly:
    1-6 = white/gray/csf/bone/skin/air, 7..7+numGel-1 = per-electrode gel,
    7+numGel..7+numGel+numElec-1 = per-electrode metal.

    This numbering scheme relies on cgalmesh assigning consecutive region
    ids to *sorted distinct nonzero input label values* -- confirmed
    empirically, e.g. an input using only labels {1, 7} (a gap) comes back
    with elements labeled {1, 2}, not {1, 7}. That's harmless here only
    because a real head's 6 tissue labels are essentially always all
    present (no gaps 1-6), so gel starting at 7 and electrodes right after
    stay contiguous with them; MATLAB's own meshByIso2mesh.m makes the
    same assumption against the same underlying binary, so this isn't a
    behavior gap introduced by the port. A tissue label that happens to be
    entirely absent (e.g. zero air voxels) would silently shift every
    later region id down by one in both MATLAB and here.

    Returns (node, elem, face) -- node already in physical/mm space.
    """
    num_tissue = 6
    all_mask = tissue_labels.copy().astype(np.uint8)

    num_gel = int(gel_mask.max()) if gel_mask.size else 0
    for i in range(1, num_gel + 1):
        all_mask[gel_mask == i] = num_tissue + i

    num_elec = int(elec_mask.max()) if elec_mask.size else 0
    for i in range(1, num_elec + 1):
        all_mask[elec_mask == i] = num_tissue + num_gel + i

    node, elem, face = run_cgal_mesher(
        all_mask, work_dir, radbound, angbound, distbound, reratio, maxvol, bin_path=bin_path
    )

    # Voxel-corner mesh coordinates -> voxel-center, then scaled into
    # physical (mm) space by the MRI's voxel size -- ports
    # meshByIso2mesh.m's `node(:,1:3)=node(:,1:3)+0.5` (already applied by
    # run_cgal_mesher) followed by `node(:,i)=node(:,i)*imgHdr.mat(i,i)`.
    node = node.copy()
    for i in range(3):
        node[:, i] = node[:, i] * voxel_size[i]

    region_names = ["WHITE", "GRAY", "CSF", "BONE", "SKIN", "AIR"]
    region_names += [f"GEL{i}" for i in range(1, num_gel + 1)]
    region_names += [f"ELEC{i}" for i in range(1, num_elec + 1)]

    save_msh(node[:, :3], elem, str(out_path), region_names)
    return node, elem, face
