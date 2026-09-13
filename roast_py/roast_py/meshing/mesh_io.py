"""Ports saveinr.m, readmedit.m, sortmesh.m, and savemsh.m -- the file
formats at the mesher's input/output boundary and the Gmsh output ROAST
hands to getDP.
"""

from __future__ import annotations

import numpy as np

_INR_HEADER_SIZE = 256

_INR_DTYPES = {
    "uint8": ("unsigned fixed", np.uint8, 8),
    "uint16": ("unsigned fixed", np.uint16, 16),
    "float32": ("float", np.float32, 32),
    "float64": ("float", np.float64, 64),
}


def save_inr(vol: np.ndarray, fname: str) -> None:
    """Ports saveinr.m: writes a volume in INR-4 format (a fixed 256-byte
    ASCII header, newline-padded, followed by raw voxel data), the format
    the cgalmesh binary reads as input.
    """
    if vol.dtype == np.bool_:
        key = "uint8"
        vol = vol.astype(np.uint8)
    elif str(vol.dtype) in _INR_DTYPES:
        key = str(vol.dtype)
    else:
        raise ValueError(f"volume format not supported: {vol.dtype}")

    btype, np_dtype, bitlen = _INR_DTYPES[key]
    nx, ny, nz = vol.shape
    header = (
        f"#INRIMAGE-4#{{\nXDIM={nx}\nYDIM={ny}\nZDIM={nz}\nVDIM=1\nTYPE={btype}\n"
        f"PIXSIZE={bitlen} bits\nCPU=decm\nVX=1\nVY=1\nVZ=1\n"
    )
    pad = "\n" * (_INR_HEADER_SIZE - 4 - len(header))
    header = header + pad + "##}\n"
    assert len(header) == _INR_HEADER_SIZE, len(header)

    with open(fname, "wb") as f:
        f.write(header.encode("ascii"))
        # INR stores the volume in Fortran (column-major) order, i.e. x
        # varying fastest -- matches MATLAB's own native array layout, so
        # this transpose is the actual translation work (not a no-op).
        f.write(np.asarray(vol, dtype=np_dtype).tobytes(order="F"))


def read_medit(filename: str):
    """Ports readmedit.m: reads MEDIT .mesh format (cgalmesh's output).

    Faithfully replicates the original's parsing quirk: for any section
    keyword it doesn't recognize (e.g. the "MeshVersionFormatted 1" /
    "Dimension 3" header lines cgalmesh always emits), it reads exactly
    one following integer and moves on, rather than skipping that
    section's actual data block -- which only works because those two
    header lines are single key-value pairs, not key+count+data blocks.

    Tokenizes the whole file up front with str.split() (a single C-level
    pass) and slices out each data block as one array via
    numpy.array(..., dtype=...) instead of converting token-by-token in a
    Python loop -- a real head mesh has hundreds of thousands to millions
    of nodes/elements, where the naive per-token approach would be far too
    slow.

    Returns (node, elem, face): node is (N,4) -- x,y,z,ref (1-based mesh
    convention, not roast_py's usual 0-based voxel convention -- see
    package docstring); elem is (N,5) -- 4 node indices + region; face is
    (N,4) -- 3 node indices + boundary id.
    """
    node = np.empty((0, 4))
    elem = np.empty((0, 5))
    face = np.empty((0, 4))

    with open(filename) as fid:
        tokens = fid.read().split()

    i = 0
    n = len(tokens)
    while i < n:
        key = tokens[i]
        i += 1
        if key == "End":
            break
        val = int(tokens[i])
        i += 1
        if key == "Vertices":
            block = tokens[i : i + 4 * val]
            node = np.array(block, dtype=float).reshape(val, 4)
            i += 4 * val
        elif key == "Triangles":
            block = tokens[i : i + 4 * val]
            face = np.array(block, dtype=int).reshape(val, 4)
            i += 4 * val
        elif key == "Tetrahedra":
            block = tokens[i : i + 5 * val]
            elem = np.array(block, dtype=int).reshape(val, 5)
            i += 5 * val
    return node, elem, face


def sort_mesh(origin, node: np.ndarray, elem: np.ndarray, ecol=None, face: np.ndarray | None = None, fcol=None):
    """Ports sortmesh.m: reorders nodes (and remaps element/face node
    references accordingly) by spherical coordinates around `origin`, to
    improve cache locality for the downstream FEM solve. Node references
    in `elem`/`face` are 1-based (see module docstring).

    Returns (node_sorted, elem_sorted, face_sorted_or_None).
    """
    node = np.asarray(node, dtype=float)
    if origin is None:
        origin = node[0, :3]
    origin = np.asarray(origin, dtype=float)

    d = node[:, :3] - origin
    r = np.sqrt(np.sum(d**2, axis=1))
    # MATLAB's cart2sph: theta = azimuth = atan2(y, x), phi = elevation = atan2(z, hypot(x, y)).
    theta = np.arctan2(d[:, 1], d[:, 0])
    phi = np.arctan2(d[:, 2], np.sqrt(d[:, 0] ** 2 + d[:, 1] ** 2))
    sort_key = np.stack([r, phi, theta], axis=1)

    nodemap = np.lexsort((sort_key[:, 2], sort_key[:, 1], sort_key[:, 0]))  # sortrows == lexsort on reversed keys
    no = node[nodemap]

    nidx = np.argsort(nodemap)  # inverse permutation: old 0-based index -> new 0-based index

    if ecol is None:
        ecol = list(range(elem.shape[1]))
    el = elem.copy()
    remapped = nidx[(el[:, ecol].astype(int) - 1)] + 1  # translate 1-based refs through the permutation
    el[:, ecol] = np.sort(remapped, axis=1)
    el = el[np.lexsort(tuple(el[:, c] for c in reversed(ecol)))]

    if face is not None:
        if fcol is None:
            fcol = list(range(face.shape[1]))
        fc = face.copy()
        remapped_f = nidx[(fc[:, fcol].astype(int) - 1)] + 1
        fc[:, fcol] = np.sort(remapped_f, axis=1)
        fc = fc[np.lexsort(tuple(fc[:, c] for c in reversed(fcol)))]
        return no, el, fc

    return no, el, None


def save_msh(node: np.ndarray, elem: np.ndarray, fname: str, region_names: list[str] | None = None) -> None:
    """Ports savemsh.m: writes a tetrahedral mesh in Gmsh MSH 2.2 ASCII
    format (element type 4 = 4-node tetrahedron; 2 tags per element, both
    set to the region id, matching the original's `buffer(4,:) =
    buffer(5,:) = M.Elements.region(et)`).

    Uses numpy.savetxt for the bulk node/element blocks rather than a
    Python-level per-row write loop -- a real head mesh has hundreds of
    thousands of nodes and over a million elements (e.g. ~224k nodes /
    1.3M elements for a real subject in this port's own testing), where a
    plain per-row f.write() loop is a real bottleneck.
    """
    node = np.asarray(node, dtype=float)
    elem = np.asarray(elem)
    if elem.shape[1] < 5:
        elem = np.column_stack([elem, np.ones(elem.shape[0], dtype=int)])

    n_nodes = node.shape[0]
    n_elem = elem.shape[0]

    with open(fname, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        f.write("$Nodes\n")
        f.write(f"{n_nodes}\n")
        node_rows = np.column_stack([np.arange(1, n_nodes + 1), node[:, :3]])
        np.savetxt(f, node_rows, fmt=["%d", "%.10f", "%.10f", "%.10f"])
        f.write("$EndNodes\n")

        f.write("$Elements\n")
        f.write(f"{n_elem}\n")
        region = elem[:, 4].astype(int)
        elem_rows = np.column_stack(
            [np.arange(1, n_elem + 1), np.full(n_elem, 4), np.full(n_elem, 2), region, region, elem[:, :4].astype(int)]
        )
        np.savetxt(f, elem_rows, fmt="%d")
        f.write("$EndElements\n")
