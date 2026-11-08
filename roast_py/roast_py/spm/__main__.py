"""python -m roast_py.spm {status,install,check} -- see roast_py.spm."""

from __future__ import annotations

import argparse
import sys

from . import runtime


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m roast_py.spm", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status", help="report the SPM packages and MATLAB Runtime found")
    install = sub.add_parser("install", help="install the SPM packages and the MATLAB Runtime")
    install.add_argument("--accept-license", action="store_true", help="accept the MATLAB Runtime license")
    install.add_argument("--prefix", default=None, help="where to install the runtime (default: per-user folder)")
    sub.add_parser("check", help="start SPM once and print its version")
    args = parser.parse_args(argv)

    if args.command == "status":
        for name, wanted in runtime.SPM_PACKAGES.items():
            have = runtime.installed_version(name)
            print(f"{name:20s} {have or '(not installed)':16s} pinned {wanted}")
        found = runtime.find_runtime() if runtime.installed_version("matlab-runtime") else None
        print(f"MATLAB Runtime {runtime.MATLAB_RELEASE}: {found or 'not installed'} "
              f"(default location {runtime.default_prefix()})")
        return 0 if found and not runtime.package_drift() else 1

    if args.command == "install":
        runtime.ensure_packages()
        found = runtime.find_runtime(args.prefix)
        if found:
            print(f"MATLAB Runtime {runtime.MATLAB_RELEASE} already installed: {found}")
        else:
            found = runtime.install_runtime(args.prefix, accept_license=True if args.accept_license else None)
            print(f"MATLAB Runtime {runtime.MATLAB_RELEASE} installed: {found}")
        return 0

    found = runtime.ensure_spm()
    runtime.run_job({"check": True}, runtime=found)
    return 0


if __name__ == "__main__":
    sys.exit(main())
