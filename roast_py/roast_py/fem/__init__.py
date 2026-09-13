"""FEM solve (ROAST Phase 4 port): boundary-condition setup, .pro
generation, the getdp subprocess, and .pos post-processing.

Like roast_py.meshing, mesh node references here are kept 1-based to
match the Gmsh/getDP file formats, not roast_py's usual 0-based voxel
convention.
"""
