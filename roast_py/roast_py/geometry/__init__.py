"""Electrode-placement geometry (ROAST Phase 2 port).

Convention used consistently throughout this package: voxel coordinates
are 0-based numpy indices (i, j, k), unlike MATLAB's 1-based
ind2sub/sub2ind -- there is no reason to carry MATLAB's 1-based offset
into a Python port, since nothing here writes back into a MATLAB-indexed
file format. Physical/world coordinates are only ever mentioned
explicitly as such.
"""
