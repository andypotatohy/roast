"""Subprocess wrapper around the getdp binary (ports the `system(cmd)`
call at the end of solveByGetDP.m). getdp is already a plain, per-platform
executable (lib/getdp-3.2.0/bin/{getdp,getdp.exe,getdpMac}) -- unlike the
CGAL mesher, there was never any MEX-file ambiguity here.
"""

from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path


def find_getdp_binary(bin_path: str | os.PathLike | None = None) -> Path:
    if bin_path is not None:
        return Path(bin_path)

    system = platform.system()
    name = {"Windows": "getdp.exe", "Linux": "getdp", "Darwin": "getdpMac"}.get(system)
    if name is None:
        raise RuntimeError(f"Unsupported operating system: {system!r}")

    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "lib" / "getdp-3.2.0" / "bin" / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Could not locate the bundled getdp binary ({name}) under lib/getdp-3.2.0/bin/. "
        "Pass bin_path explicitly."
    )


def run_getdp(pro_path: str | os.PathLike, ready_msh_path: str | os.PathLike, bin_path: str | os.PathLike | None = None) -> None:
    """Ports solveByGetDP.m's `system(cmd)` call: runs getdp's `EleSta_v`
    resolution + `Map` post-operation against the given .pro/_ready.msh.

    MATLAB does both in one getdp call. Here they are two: `-solve` (which
    saves the solution to the .res file), then `-pos Map` reading that
    .res back. The results match -- the voltage .pos byte for byte, the
    E-field to ~1e-11 V/m (round-off from reloading the solution) -- but
    the peak memory is lower, because the one-call form keeps the direct
    solver's LU factorization in memory through the post-processing. For
    MNI152 at roast.m's mesh settings that's the difference between
    fitting in ~7.6 GB and being killed above 8.6 GB.

    Runs with cwd set to the .pro file's own directory, since the .pro
    file's `Print[...File "bare_name.pos"...]` directives (see
    pro_writer.py) use bare filenames, not paths -- matching where
    postGetDP.m expects to find them afterward.
    """
    pro_path = Path(pro_path).resolve()
    ready_msh_path = Path(ready_msh_path).resolve()
    res_path = pro_path.with_suffix(".res")
    binary = find_getdp_binary(bin_path)

    steps = [
        [str(binary), str(pro_path), "-solve", "EleSta_v", "-msh", str(ready_msh_path)],
        [str(binary), str(pro_path), "-msh", str(ready_msh_path), "-res", str(res_path), "-pos", "Map"],
    ]
    try:
        for cmd in steps:
            _run_step(cmd, cwd=pro_path.parent)
    finally:
        # Ports the original's clean-up of getdp's intermediate .pre/.res files.
        for ext in (".pre", ".res"):
            pro_path.with_suffix(ext).unlink(missing_ok=True)


def _run_step(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode in (-9, 137):  # SIGKILL (137 = 128 + 9 through a shell)
        # Almost always the kernel's out-of-memory killer: getDP's direct
        # (MUMPS LU) solve needs several GB for a full head mesh, the same
        # as in MATLAB ROAST, and dies without printing anything.
        raise RuntimeError(
            "getDP was killed while solving, most likely for running out of memory "
            "(its direct solver needs several GB for a full-head mesh: ~7.6 GB for "
            "the MNI152 head at ROAST's default mesh settings, more for bigger "
            "heads -- the same as in MATLAB ROAST). Run it on a machine with more "
            "memory, or use a coarser mesh via roast(..., mesh_options=...). A head "
            "mesh's size is set mostly by the surface settings (radbound, "
            "distbound), not maxvol; too coarse and the mesher can miss small "
            "electrodes.\n"
            f"command: {' '.join(cmd)}\nlast output:\n{result.stdout[-1500:]}"
        )
    if result.returncode != 0:
        raise RuntimeError(
            "getDP solver cannot work properly on your system. Please check any error "
            f"message you got.\ncommand: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
