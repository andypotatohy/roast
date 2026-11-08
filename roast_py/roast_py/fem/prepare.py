"""Ports prepareForGetDP.m's file-rewriting step: appends 2D boundary
(triangle) elements for each electrode's outer surface to the .msh file,
producing the `_ready.msh` getDP actually solves on.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from .boundary import extract_electrode_outer_surface


def attached_to_electrodes(elem: np.ndarray, num_tissue: int = 6) -> np.ndarray:
    """Which tetrahedra belong to a piece of the mesh that contains a gel
    or electrode element (pieces = tetrahedra connected through shared
    nodes). Not in MATLAB ROAST.

    Anything else floats free: no current enters it and no voltage is set
    on it, so its block of getDP's matrix is singular on its own. Usually
    round-off hides that, but not always. runExamples cmd 12 (Multiaxial
    segmentation of subject1) meshes three loose specks of "skin" from the
    segmentation, and MUMPS stops with INFOG(1)=-10, "numerically singular
    matrix". Such a piece cannot affect the field in the head.
    """
    tets = np.asarray(elem[:, :4], dtype=np.int64) - 1
    n = int(tets.max()) + 1
    rows = np.repeat(tets[:, 0], 3)
    cols = tets[:, 1:].ravel()
    graph = sp.coo_matrix((np.ones(rows.size, dtype=np.int8), (rows, cols)), shape=(n, n))
    _, piece = connected_components(graph, directed=False)
    tet_piece = piece[tets[:, 0]]
    with_electrodes = np.unique(tet_piece[elem[:, 4] > num_tissue])
    return np.isin(tet_piece, with_electrodes)


def prepare_for_getdp(msh_path: str, ready_msh_path: str, node: np.ndarray, elem: np.ndarray, elec_names: list[str]):
    """Computes each electrode's outer-surface boundary elements and area,
    then rewrites `msh_path` (as produced by
    roast_py.meshing.cgal_mesher.mesh_by_iso2mesh) into `ready_msh_path`
    with those boundary elements appended before `$EndElements`, bumping
    the declared element count accordingly. Pieces of the mesh with no
    electrode in them are left out of `ready_msh_path`, as they would make
    the system singular (see attached_to_electrodes).

    Returns (element_elec_needed, area_elec_needed): per-electrode lists
    matching elec_names' order -- outer-surface faces (1-based global node
    index triples) and total surface area (mm^2, physical space).
    """
    num_tissue = 6
    num_elec = len(elec_names)

    element_elec_needed = []
    area_elec_needed = np.zeros(num_elec)

    for i in range(1, num_elec + 1):
        gel_tets = elem[elem[:, 4] == num_tissue + i, :4]
        elec_tets = elem[elem[:, 4] == num_tissue + num_elec + i, :4]
        try:
            outer_faces, area = extract_electrode_outer_surface(gel_tets, elec_tets, node[:, :3])
        except ValueError as e:
            raise ValueError(f"{e} (electrode {elec_names[i - 1]!r})") from e
        element_elec_needed.append(outer_faces)
        area_elec_needed[i - 1] = area

    num_of_part = len(np.unique(elem[:, 4]))
    total_extra = sum(f.shape[0] for f in element_elec_needed)
    keep = attached_to_electrodes(elem, num_tissue)
    if not keep.all():
        print(f"Leaving out {int((~keep).sum())} tetrahedra in pieces of the mesh not attached to the head")

    with open(msh_path) as fin, open(ready_msh_path, "w") as fout:
        lines = iter(fin)
        for raw_line in lines:
            s = raw_line.rstrip("\n")
            if s == "$Elements":
                fout.write(s + "\n")
                num_listed = int(next(lines).rstrip("\n"))
                if num_listed != len(elem):
                    raise ValueError(f"{msh_path} lists {num_listed} elements, the mesh has {len(elem)}")
                num_elem = int(keep.sum())
                fout.write(str(num_elem + total_extra) + "\n")
                eid = 0
                for k in range(num_listed):
                    line = next(lines)
                    if keep[k]:
                        eid += 1
                        fout.write(f"{eid} {line.split(None, 1)[1]}")
            elif s == "$EndElements":
                offset = 0
                for j in range(num_elec):
                    faces = element_elec_needed[j]
                    region = num_of_part + j + 1
                    for i in range(faces.shape[0]):
                        eid = num_elem + offset + i + 1
                        n1, n2, n3 = (int(x) for x in faces[i])
                        fout.write(f"{eid} 2 2 {region} {region} {n1} {n2} {n3} \n")
                    offset += faces.shape[0]
                fout.write(s + "\n")
            else:
                fout.write(s + "\n")

    return element_elec_needed, area_elec_needed
