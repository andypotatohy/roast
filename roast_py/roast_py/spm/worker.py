"""Runs SPM jobs through spm-python, in a process of its own.

    python -m roast_py.spm.worker job.json

`job.json` is {"batches": [<entry>, ...]} with entries from
roast_py.spm.batches (run together, as one matlabbatch, the way ROAST's
MATLAB code runs them), or {"check": true} to only start SPM and print
its version. roast_py.spm.runtime.run_job() writes the file and starts
this, with $MATLAB_RUNTIME_PATH pointing at the runtime (on macOS through
mwpython).
"""

from __future__ import annotations

import json
import sys


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: python -m roast_py.spm.worker job.json", file=sys.stderr)
        return 2
    with open(argv[0]) as f:
        job = json.load(f)

    import spm  # starts the MATLAB Runtime on first call into SPM

    from .batches import to_matlabbatch

    print(f"SPM: {spm.spm('Ver')}", flush=True)
    if job.get("check"):
        return 0
    spm.spm_jobman("initcfg", nargout=0)
    spm.spm_jobman("run", to_matlabbatch(job["batches"]), nargout=0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
