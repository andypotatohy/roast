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

    Runs with cwd set to the .pro file's own directory, since the .pro
    file's `Print[...File "bare_name.pos"...]` directives (see
    pro_writer.py) use bare filenames, not paths -- matching where
    postGetDP.m expects to find them afterward.
    """
    pro_path = Path(pro_path).resolve()
    ready_msh_path = Path(ready_msh_path).resolve()
    binary = find_getdp_binary(bin_path)

    cmd = [str(binary), str(pro_path), "-solve", "EleSta_v", "-msh", str(ready_msh_path), "-pos", "Map"]
    result = subprocess.run(cmd, cwd=pro_path.parent, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "getDP solver cannot work properly on your system. Please check any error "
            f"message you got.\ncommand: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    # Ports the original's clean-up of getdp's intermediate .pre/.res files.
    for ext in (".pre", ".res"):
        stale = pro_path.with_suffix(ext)
        stale.unlink(missing_ok=True)
