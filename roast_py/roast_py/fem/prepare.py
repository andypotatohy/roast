"""Ports prepareForGetDP.m's file-rewriting step: appends 2D boundary
(triangle) elements for each electrode's outer surface to the .msh file,
producing the `_ready.msh` getDP actually solves on.
"""

from __future__ import annotations

import numpy as np

from .boundary import extract_electrode_outer_surface


def prepare_for_getdp(msh_path: str, ready_msh_path: str, node: np.ndarray, elem: np.ndarray, elec_names: list[str]):
    """Computes each electrode's outer-surface boundary elements and area,
    then rewrites `msh_path` (as produced by
    roast_py.meshing.cgal_mesher.mesh_by_iso2mesh) into `ready_msh_path`
    with those boundary elements appended before `$EndElements`, bumping
    the declared element count accordingly.

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

    with open(msh_path) as fin, open(ready_msh_path, "w") as fout:
        lines = iter(fin)
        for raw_line in lines:
            s = raw_line.rstrip("\n")
            if s == "$Elements":
                fout.write(s + "\n")
                num_elem = int(next(lines).rstrip("\n"))
                fout.write(str(num_elem + total_extra) + "\n")
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
