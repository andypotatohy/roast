"""Ports the boundary-extraction core of prepareForGetDP.m: for each
electrode, finds the part of the electrode-metal region's outer surface
that does NOT touch the gel region (the true "outward-facing" surface,
where ROAST applies the injected-current Neumann boundary condition), and
its area.

MATLAB does this with `TriRep(tets, node)` + `freeBoundary`, which
computes the free boundary of a tetrahedral sub-mesh and returns it with a
*compacted, locally renumbered* point list (a MATLAB TriRep/freeBoundary
API convention). This port skips that renumbering: since every boundary
face is built directly from the tet list's own (already-global) node
references, there's no reason to translate to a local numbering only to
translate back -- `free_boundary_faces` returns faces in terms of the
mesh's real global (1-based) node indices throughout, which is also
exactly the form the .msh boundary-element section needs them in.
"""

from __future__ import annotations

import numpy as np

# The 4 triangular faces of a tetrahedron (local vertex index triples),
# in the same order MATLAB's TriRep uses.
_TET_LOCAL_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def free_boundary_faces(tets: np.ndarray) -> np.ndarray:
    """Ports freeBoundary(TriRep(tets, node)): the triangular faces that
    belong to exactly one tetrahedron in `tets` (a face shared by two
    tets in the list is internal, not a boundary). `tets` is (n, 4),
    1-based global node indices. Returns (m, 3), same node-index
    convention, in each face's original tet-local vertex order.
    """
    tets = np.asarray(tets, dtype=np.int64)
    all_faces = np.concatenate([tets[:, list(f)] for f in _TET_LOCAL_FACES], axis=0)

    sorted_faces = np.sort(all_faces, axis=1)
    _, inverse, counts = np.unique(sorted_faces, axis=0, return_inverse=True, return_counts=True)
    inverse = inverse.ravel()
    is_boundary = counts[inverse] == 1
    return all_faces[is_boundary]


def triangle_areas(faces: np.ndarray, node_xyz: np.ndarray) -> np.ndarray:
    """Areas of triangles given as 1-based node-index triples into
    `node_xyz` (physical/mm coordinates)."""
    p0 = node_xyz[faces[:, 0] - 1]
    p1 = node_xyz[faces[:, 1] - 1]
    p2 = node_xyz[faces[:, 2] - 1]
    cross = np.cross(p1 - p0, p2 - p0)
    return 0.5 * np.sqrt(np.sum(cross**2, axis=1))


def extract_electrode_outer_surface(gel_tets: np.ndarray, elec_tets: np.ndarray, node_xyz: np.ndarray):
    """Ports the per-electrode body of prepareForGetDP.m's loop: given one
    electrode's gel-region and electrode-metal-region tetrahedra, returns
    (outer_faces, area) -- the electrode surface not in contact with the
    gel, and its total physical-space area.
    """
    if gel_tets.shape[0] == 0:
        raise ValueError(
            "Gel under this electrode was not meshed properly. Reasons may be: 1) "
            "electrode size is too small so the mesher cannot capture it; 2) mesh "
            "resolution is not high enough. Consider using bigger electrodes or "
            "increasing the mesh resolution."
        )
    if elec_tets.shape[0] == 0:
        raise ValueError(
            "Electrode was not meshed properly. Reasons may be: 1) electrode size is "
            "too small so the mesher cannot capture it; 2) mesh resolution is not high "
            "enough. Consider using bigger electrodes or increasing the mesh resolution."
        )

    gel_faces = free_boundary_faces(gel_tets)
    elec_faces = free_boundary_faces(elec_tets)

    gel_boundary_nodes = set(np.unique(gel_faces).tolist())
    is_gel_contact = np.array(
        [all(v in gel_boundary_nodes for v in face) for face in elec_faces]
    )
    outer_faces = elec_faces[~is_gel_contact]

    area = float(np.sum(triangle_areas(outer_faces, node_xyz)))
    return outer_faces, area
