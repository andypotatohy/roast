"""Python translation of ROAST's runExamples.m.

Runs the same calls, in the same order and with the same numbers, through
roast_py. As in MATLAB, each one is wrapped so a failure is appended to
errLog.txt ("error at running cmd N:" + the message) and the script moves
on to the next.

    python examples/run_examples.py              # all of them, unattended
    python examples/run_examples.py 1 14 21      # just these
    python examples/run_examples.py --list       # print the calls
    python examples/run_examples.py --show 14    # open the figure windows

Each call runs in its own Python process, so memory is released between
calls and a crash (a segmentation fault, or the out-of-memory killer) is
logged as that call's error instead of ending the whole run. Pass
--in-process to run them all in this one.

Inputs and outputs are under the repository's example/ folder, as in
MATLAB, so example/ fills up with results. Untagged runs on the same
subject overwrite each other there (roast_py doesn't give each run a
timestamp tag the way MATLAB does), so after each successful roast() call
the script also saves its figures as PNGs in
example/run_examples_figures/cmdNN/. Pass --no-figures to skip that, or
--show to open the windows instead, which waits for you to close them
before the next call (MATLAB's commented-out `pause`).

Left out: the roast_target() calls (cmds 34-37) and the reviewRes() calls
that show their results (cmds 39-41). The lead-field roast() calls that
only feed roast_target (cmds 32-33) are kept; they fail fast because lead
fields aren't ported.

The calls with T2 images, resampling or multiaxial 'off' use SPM, which
needs the MATLAB Runtime (installed on first use, after you accept its
license; set ROAST_PY_ACCEPT_MATLAB_RUNTIME_LICENSE=yes to run unattended
-- see roast_py.spm). Calls that leave 'multiaxial' out use roast_py's
default, Multiaxial, where MATLAB's default is SPM. The New York head and
lead fields aren't ported and fail fast with NotImplementedError. Each
roast() run that goes through takes ~10-15 minutes on CPU and several GB
of memory for the solve: ~7.6 GB for the MNI152 head, up to ~19 GB for
subject1 (with SPM's segmentation).
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]  # the roast_py/ project directory
REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = REPO_ROOT / "example"

# Work straight from a git clone, with no `pip install -e .` step.
if importlib.util.find_spec("roast_py") is None:
    sys.path.insert(0, str(PACKAGE_ROOT))

from roast_py import review_res, roast  # noqa: E402  (importable with nothing installed)

SUBJECT1 = str(EXAMPLE / "subject1.nii")
SUBJECT1_T2 = str(EXAMPLE / "subject1_T2.nii")

FOUR_ELECTRODES = {"F1": 0.3, "P2": 0.7, "C5": -0.6, "O2": -0.4}
THREE_ELECTRODES = {"Fp1": 1, "FC4": 1, "POz": -2}

# (MATLAB cmd number, function, positional args, keyword args), in
# runExamples.m's order. `None` stands for MATLAB's [] / '' (use the default).
EXAMPLES = [
    (1, roast, (), {}),
    (2, roast, (None, None), dict(multiaxial="on")),
    (3, roast, (str(EXAMPLE / "bikson.hdr"),), {}),
    (4, roast, (str(EXAMPLE / "bikson.img"),), {}),
    (5, roast, (None, {"Exx19": 1, "C4": -1}), dict(zero_padding=20)),
    (6, roast, (None, {"Exx19": 1, "C4": -1}), dict(zero_padding=20, multiaxial="on")),
    (7, roast, (SUBJECT1, FOUR_ELECTRODES), {}),
    (8, roast, (SUBJECT1, FOUR_ELECTRODES), dict(multiaxial="on")),
    (9, roast, (SUBJECT1, None), dict(t2=SUBJECT1_T2)),
    (10, roast, (SUBJECT1, None), dict(t2=SUBJECT1_T2, multiaxial="on")),
    (11, roast, (SUBJECT1, FOUR_ELECTRODES), dict(resampling="on", zero_padding=20)),
    (12, roast, (SUBJECT1, FOUR_ELECTRODES), dict(resampling="on", zero_padding=20, multiaxial="on")),
    # roast('example/subject1.nii',{'G12',1,'J7',-1},'captype','biosemi')  -- commented out in MATLAB too
    (
        13,
        roast,
        (SUBJECT1, {"G12": 0.25, "J7": -0.25, "Nk1": 0.5, "Nk3": -0.5, "custom1": 0.25, "custom3": -0.25}),
        dict(cap_type="biosemi"),
    ),
    (14, roast, (None, THREE_ELECTRODES), dict(elec_type=["disc", "pad", "ring"])),
    (
        15,
        roast,
        (None, THREE_ELECTRODES),
        dict(elec_type=["disc", "pad", "ring"], elec_size=[[8, 2], [45, 25, 4], [5, 8, 2]]),
    ),
    (16, roast, (None, None), dict(elec_type="pad", elec_ori="ap")),
    (17, roast, (None, None), dict(elec_type="pad", elec_ori=[0.71, 0.71, 0])),
    # roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',{'ap','lr','si'})
    #   -- commented out in MATLAB too
    (
        18,
        roast,
        (SUBJECT1, THREE_ELECTRODES),
        dict(elec_type="pad", elec_ori=[[0.71, 0.71, 0], [-0.71, 0.71, 0], [0, 0.71, 0.71]]),
    ),  # MATLAB's note: the direction vectors for pads aren't really intuitive
    (
        19,
        roast,
        (None, THREE_ELECTRODES),
        dict(elec_type=["pad", "disc", "pad"], elec_ori=[[0.71, 0.71, 0], [0, 0.71, 0.71]]),
    ),
    (
        20,
        roast,
        (None, THREE_ELECTRODES),
        dict(elec_type=["pad", "disc", "pad"], elec_ori=["ap", None, [0, 0.71, 0.71]]),
    ),
    (21, roast, (None, None), dict(mesh_options={"radbound": 4, "maxvol": 8})),
    (
        22,
        roast,
        (SUBJECT1, None),
        dict(resampling="on", zero_padding=10, t2=SUBJECT1_T2, elec_type="pad"),
    ),
    (23, roast, (None, THREE_ELECTRODES), dict(conductivities={"csf": 0.6, "electrode": 0.1})),
    (
        24,
        roast,
        (None, THREE_ELECTRODES),
        dict(
            elec_type=["pad", "disc", "pad"],
            conductivities={"gel": [1, 0.3, 1], "electrode": [0.1, 5.9e7, 0.1]},
        ),
    ),
    (
        25,
        roast,
        (SUBJECT1, {"Fp1": 0.3, "F8": 0.2, "POz": -0.4, "Nk1": 0.5, "custom1": -0.6}),
        dict(
            elec_type=["disc", "ring", "pad", "ring", "pad"],
            elec_size=[None, [7, 9, 3], [40, 20, 4], None, None],
            elec_ori="ap",
            t2=SUBJECT1_T2,
            mesh_options={"radbound": 4, "maxvol": 8},
            conductivities={"csf": 0.6, "skin": 1.0},
            resampling="on",
            zero_padding=30,
            simulation_tag="awesomeSimulation",
        ),
    ),
    (26, roast, ("nyhead",), {}),
    (27, roast, ("nyhead", None), dict(multiaxial="on")),
    (28, roast, ("nyhead", None), dict(resampling="on", zero_padding=25)),
    (29, roast, ("nyhead", None), dict(resampling="on", zero_padding=25, multiaxial="on")),
    (
        30,
        roast,
        ("nyhead", None),
        dict(zero_padding=5, elec_type="ring", elec_size=[7, 10, 3], multiaxial="on"),
    ),
    (
        31,
        roast,
        ("nyhead", THREE_ELECTRODES),
        dict(elec_type="ring", elec_size=[[7, 10, 3], [6, 8, 3], [4, 6, 2]]),
    ),
    (
        32,
        roast,
        (None, "leadfield"),
        dict(multiaxial="on", zero_padding=10, elec_type="ring", simulation_tag="LFwithMA"),
    ),
    (
        33,
        roast,
        ("nyhead", "leadfield"),
        dict(multiaxial="on", zero_padding=10, elec_type="ring", simulation_tag="nyLFwithMA"),
    ),
    # cmds 34-37: roast_target(...) -- left out (roast_target isn't ported)
    (38, review_res, (SUBJECT1,), dict(sim_tag="awesomeSimulation", tissue="all")),
    # cmds 39-41: reviewRes(..., tarTag) of the roast_target results -- left out
    # roast('','','manualgui','on')  -- commented out in MATLAB too
]


def describe(fn, args, kwargs) -> str:
    def short(value):
        text = repr(value)
        return text.replace(str(EXAMPLE), "example") if isinstance(value, str) else text

    parts = [short(a) for a in args] + [f"{k}={short(v)}" for k, v in kwargs.items()]
    return f"{fn.__name__}({', '.join(parts)})"


def log_error(log_path: Path, number: int, error: BaseException) -> None:
    """MATLAB's catch block: append the command number and message."""
    with open(log_path, "a") as f:
        f.write(f"error at running cmd {number}:\n")
        f.write(f"{type(error).__name__}: {error}\n")


def run_one(number, fn, args, kwargs, figures: str, log_path: Path) -> bool:
    print(f"\n{'=' * 70}\ncmd {number}: {describe(fn, args, kwargs)}\n{'=' * 70}", flush=True)
    start = time.time()
    try:
        if fn is roast:
            result = roast(*args, **kwargs, visualize=figures == "show")
            if figures == "save":
                save_dir = EXAMPLE / "run_examples_figures" / f"cmd{number:02d}"
                review_res(result.subj, sim_tag=result.sim_tag, show=False, save_dir=str(save_dir))
        else:
            fn(*args, **kwargs, show=figures == "show",
               save_dir=str(EXAMPLE / "run_examples_figures" / f"cmd{number:02d}") if figures == "save" else None)
    except Exception as error:  # noqa: BLE001 -- runExamples.m's try/catch
        traceback.print_exc()
        log_error(log_path, number, error)
        print(f"cmd {number}: FAILED after {time.time() - start:.0f} s ({type(error).__name__})", flush=True)
        return False
    finally:
        try:  # MATLAB's `close all`
            import matplotlib.pyplot as plt

            plt.close("all")
        except ImportError:
            pass
    print(f"cmd {number}: OK in {time.time() - start:.0f} s", flush=True)
    return True


def run_in_child(number: int, figures: str, log_path: Path) -> bool:
    """Runs one cmd in a fresh interpreter; logs it if the child dies."""
    cmd = [sys.executable, str(Path(__file__).resolve()), str(number), "--child", "--log", str(log_path)]
    cmd += {"show": ["--show"], "none": ["--no-figures"], "save": []}[figures]
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(PACKAGE_ROOT), env.get("PYTHONPATH")]))
    returncode = subprocess.run(cmd, env=env).returncode
    if returncode < 0:  # killed by a signal: the child couldn't log it itself
        name = signal.Signals(-returncode).name
        hint = " (most likely out of memory)" if -returncode == signal.SIGKILL else ""
        log_error(log_path, number, RuntimeError(f"Python process killed by {name}{hint}"))
        print(f"cmd {number}: FAILED, process killed by {name}{hint}", flush=True)
    return returncode == 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Python translation of ROAST's runExamples.m")
    parser.add_argument("numbers", nargs="*", type=int, help="cmd numbers to run (default: all)")
    parser.add_argument("--list", action="store_true", help="print the calls and exit")
    parser.add_argument("--log", default="errLog.txt", help="error log to append to (default: %(default)s)")
    view = parser.add_mutually_exclusive_group()
    view.add_argument("--show", action="store_true", help="open the figure windows (waits for you to close them)")
    view.add_argument("--no-figures", action="store_true", help="don't draw or save figures")
    parser.add_argument(
        "--in-process", action="store_true", help="run every cmd in this process instead of one process each"
    )
    parser.add_argument("--child", action="store_true", help=argparse.SUPPRESS)  # one cmd for run_in_child
    args = parser.parse_args(argv)

    known = {number for number, *_ in EXAMPLES}
    unknown = sorted(set(args.numbers) - known)
    if unknown:
        parser.error(f"no such cmd: {unknown} (available: {sorted(known)}; 34-37 and 39-41 are roast_target)")
    selected = [e for e in EXAMPLES if not args.numbers or e[0] in args.numbers]

    if args.list:
        for number, fn, call_args, kwargs in selected:
            print(f"{number:3d}  {describe(fn, call_args, kwargs)}")
        return 0

    figures = "show" if args.show else "none" if args.no_figures else "save"
    log_path = Path(args.log).resolve()
    if args.in_process or args.child:
        results = {number: run_one(number, fn, a, k, figures, log_path) for number, fn, a, k in selected}
    else:
        results = {number: run_in_child(number, figures, log_path) for number, *_ in selected}

    if args.child:
        return 0 if all(results.values()) else 1
    passed = [n for n, ok in results.items() if ok]
    failed = [n for n, ok in results.items() if not ok]
    print(f"\n{len(passed)} of {len(results)} succeeded: {passed}")
    if failed:
        print(f"failed: {failed}; see {log_path}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
