"""Tetrahedral meshing (ROAST Phase 3 port).

Note on indexing convention: unlike the rest of roast_py (0-based voxel
coordinates everywhere), mesh *node references* inside `elem`/`face`
arrays here are kept 1-based, matching the MEDIT and Gmsh MSH file formats
this package reads/writes -- there's no reason to renumber them 0-based
internally only to renumber back for every file boundary.
"""
