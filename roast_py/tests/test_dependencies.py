"""Tests the dependency preflight and installer.

The drift test matters: the bug this module exists to fix was
pyproject.toml declaring a set of dependencies that didn't match what the
code actually imports, so `pip install -e .` produced an install that
couldn't run roast().

The conda paths are exercised by simulating a conda environment
(a directory containing conda-meta/ plus a stubbed conda executable),
since the container these tests run in has no conda.
"""

import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

from roast_py.dependencies import (
    NO_AUTO_INSTALL_ENV_VAR,
    REQUIRED,
    Dependency,
    auto_install_disabled,
    check_dependencies,
    conda_install_command,
    ensure_dependencies,
    format_missing,
    in_conda_environment,
    install_commands,
    missing_dependencies,
    pip_install_command,
    use_conda_for,
)

PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"

FAKE = Dependency("definitely_not_a_real_module_xyz", "not-real-xyz", "testing", "not-real-xyz")
FAKE_PIP_ONLY = Dependency("also_not_real_xyz", "also-not-real", "testing")  # conda_name None


@pytest.fixture
def simulated_conda_env(tmp_path, monkeypatch):
    """Makes the running interpreter look like it lives in a conda env."""
    (tmp_path / "conda-meta").mkdir()
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    monkeypatch.setattr("roast_py.dependencies.find_conda", lambda: "/fake/bin/conda")
    return tmp_path


# --------------------------------------------------------------------------
# declaration drift
# --------------------------------------------------------------------------


def _declared_runtime_deps() -> set[str]:
    with open(PYPROJECT, "rb") as f:
        data = tomllib.load(f)
    names = set()
    for spec in data["project"]["dependencies"]:
        for sep in (">=", "==", "<=", "~=", ">", "<", "!="):
            spec = spec.split(sep)[0]
        names.add(spec.strip())
    return names


def test_required_list_matches_pyproject_runtime_dependencies():
    declared = _declared_runtime_deps()
    listed = {dep.pip_name for dep in REQUIRED}
    assert listed == declared, (
        "roast_py/dependencies.py REQUIRED and pyproject.toml [project] dependencies "
        f"have drifted apart.\n  only in dependencies.py: {sorted(listed - declared)}\n"
        f"  only in pyproject.toml: {sorted(declared - listed)}"
    )


def test_every_required_dependency_is_actually_imported_by_the_package():
    """Guards the other direction: nothing declared that the code doesn't use.

    (openpyxl is the documented exception -- pandas imports it internally
    as its .xlsx engine, so it never appears in an `import` statement here
    but is genuinely required to read capInfo.xlsx.)
    """
    package_root = Path(__file__).resolve().parents[1] / "roast_py"
    sources = "\n".join(
        p.read_text() for p in package_root.rglob("*.py") if "dependencies.py" not in str(p)
    )
    for dep in REQUIRED:
        if dep.import_name == "openpyxl":
            continue
        imported = f"import {dep.import_name}" in sources or f"from {dep.import_name}" in sources
        assert imported, f"{dep.pip_name} is declared required but never imported by the package"


# --------------------------------------------------------------------------
# detection
# --------------------------------------------------------------------------


def test_all_required_dependencies_are_importable_in_this_environment():
    assert missing_dependencies() == []
    assert check_dependencies() == []


def test_missing_dependencies_reports_unimportable_packages():
    assert missing_dependencies((FAKE,)) == [FAKE]


def test_in_conda_environment_detects_conda_meta_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    assert in_conda_environment() is False  # no conda-meta yet

    (tmp_path / "conda-meta").mkdir()
    assert in_conda_environment() is True


def test_in_conda_environment_ignores_activated_env_vars(tmp_path, monkeypatch):
    """CONDA_PREFIX describes the *shell's* env, not this interpreter's."""
    monkeypatch.setattr(sys, "prefix", str(tmp_path))  # no conda-meta
    monkeypatch.setenv("CONDA_PREFIX", "/somewhere/else")
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "someenv")
    assert in_conda_environment() is False


# --------------------------------------------------------------------------
# conda / pip split and command construction
# --------------------------------------------------------------------------


def test_outside_conda_everything_goes_to_pip():
    conda_deps, pip_deps = use_conda_for([FAKE, FAKE_PIP_ONLY])
    assert conda_deps == []
    assert pip_deps == [FAKE, FAKE_PIP_ONLY]


def test_inside_conda_only_packages_with_a_conda_name_go_to_conda(simulated_conda_env):
    conda_deps, pip_deps = use_conda_for([FAKE, FAKE_PIP_ONLY])
    assert conda_deps == [FAKE]
    assert pip_deps == [FAKE_PIP_ONLY]


def test_tensorflow_and_tf_keras_are_always_pip_installed(simulated_conda_env):
    """pip is TensorFlow's official channel; tf-keras conda availability is unverified."""
    conda_deps, pip_deps = use_conda_for(list(REQUIRED))
    pip_names = {d.pip_name for d in pip_deps}
    assert {"tensorflow", "tf-keras"} <= pip_names
    # ...while the scientific stack does go through conda.
    conda_names = {d.pip_name for d in conda_deps}
    assert {"numpy", "scipy", "pandas", "scikit-image", "nibabel"} <= conda_names


def test_pip_install_command_targets_the_running_interpreter():
    cmd = pip_install_command([FAKE])
    # Must be sys.executable -m pip, not a bare `pip`, so it installs into
    # the active environment rather than whatever `pip` resolves to first.
    assert cmd[:4] == [sys.executable, "-m", "pip", "install"]
    assert cmd[-1] == "not-real-xyz"


def test_conda_install_command_targets_this_interpreters_prefix(simulated_conda_env):
    cmd = conda_install_command([FAKE])
    assert cmd[0] == "/fake/bin/conda"
    assert cmd[1] == "install"
    # --prefix sys.prefix is the conda equivalent of `sys.executable -m pip`:
    # it writes into the env this interpreter belongs to, not the activated one.
    assert "--prefix" in cmd
    assert cmd[cmd.index("--prefix") + 1] == sys.prefix
    assert "conda-forge" in cmd
    assert "-y" in cmd
    assert cmd[-1] == "not-real-xyz"


def test_conda_install_command_uses_conda_name_when_it_differs():
    dep = Dependency("skimage", "scikit-image", "testing", "scikit-image")
    assert conda_install_command([dep], conda_exe="/fake/conda")[-1] == "scikit-image"


def test_install_commands_runs_conda_before_pip(simulated_conda_env):
    commands = install_commands([FAKE, FAKE_PIP_ONLY])
    assert len(commands) == 2
    # conda first, pip last -- so conda's solver sees the environment
    # before pip writes into it.
    assert commands[0][0] == "/fake/bin/conda"
    assert commands[1][:3] == [sys.executable, "-m", "pip"]


def test_install_commands_outside_conda_is_pip_only():
    commands = install_commands([FAKE, FAKE_PIP_ONLY])
    assert len(commands) == 1
    assert commands[0][:3] == [sys.executable, "-m", "pip"]


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def test_check_dependencies_raises_with_actionable_message(monkeypatch):
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (FAKE,))

    with pytest.raises(ImportError) as excinfo:
        check_dependencies()

    message = str(excinfo.value)
    assert "not-real-xyz" in message
    assert "pip install" in message
    assert "python -m roast_py.dependencies --install" in message


def test_format_missing_shows_the_conda_command_inside_conda(simulated_conda_env):
    text = format_missing([FAKE])
    assert "conda" in text
    assert "--prefix" in text


def test_format_missing_lists_reason_for_each_package():
    text = format_missing([Dependency("x", "x-pkg", "doing the thing")])
    assert "x-pkg" in text
    assert "doing the thing" in text


# --------------------------------------------------------------------------
# ensure_dependencies (what roast() calls)
# --------------------------------------------------------------------------


def test_ensure_dependencies_is_a_noop_when_nothing_is_missing():
    ensure_dependencies(install_missing=True)  # must not attempt any install


def test_ensure_dependencies_raises_instead_of_installing_when_opted_out(monkeypatch):
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (FAKE,))
    called = []
    monkeypatch.setattr(
        "roast_py.dependencies.install_dependencies", lambda *a, **k: called.append(a)
    )

    with pytest.raises(ImportError):
        ensure_dependencies(install_missing=False)
    assert called == [], "must not install when install_missing=False"


def test_ensure_dependencies_respects_the_environment_opt_out(monkeypatch):
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (FAKE,))
    monkeypatch.setenv(NO_AUTO_INSTALL_ENV_VAR, "1")
    called = []
    monkeypatch.setattr(
        "roast_py.dependencies.install_dependencies", lambda *a, **k: called.append(a)
    )

    with pytest.raises(ImportError):
        ensure_dependencies(install_missing=True)
    assert called == [], f"must not install when {NO_AUTO_INSTALL_ENV_VAR} is set"


def test_ensure_dependencies_installs_when_asked(monkeypatch):
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (FAKE,))
    monkeypatch.delenv(NO_AUTO_INSTALL_ENV_VAR, raising=False)
    called = []
    monkeypatch.setattr(
        "roast_py.dependencies.install_dependencies", lambda deps, **k: called.append(deps)
    )

    ensure_dependencies(install_missing=True)
    assert called == [[FAKE]]


@pytest.mark.parametrize(
    "value,disabled",
    [("", False), ("0", False), ("false", False), ("False", False), ("1", True), ("yes", True)],
)
def test_auto_install_disabled_parses_the_env_var(monkeypatch, value, disabled):
    monkeypatch.setenv(NO_AUTO_INSTALL_ENV_VAR, value)
    assert auto_install_disabled() is disabled


# --------------------------------------------------------------------------
# bootstrapping properties
# --------------------------------------------------------------------------


def test_dependencies_module_imports_without_any_third_party_packages():
    """The whole point of this module: it must import in a bare environment."""
    package_parent = str(Path(__file__).resolve().parents[1])
    code = (
        "import sys;"
        f"sys.path.insert(0, {package_parent!r});"
        "import roast_py.dependencies as d;"
        "print(len(d.REQUIRED))"
    )
    result = subprocess.run([sys.executable, "-S", "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(len(REQUIRED))


def test_roast_is_importable_with_no_dependencies_installed():
    """`from roast_py import roast` must work in a bare environment.

    roast() is what installs the missing dependencies, so if importing it
    required them first, the auto-install could never run. This is the
    property that keeps roast.py's imports inside the function body.
    """
    package_parent = str(Path(__file__).resolve().parents[1])
    code = f"""
import sys
BLOCKED = {{"numpy","scipy","nibabel","pandas","openpyxl","skimage","tensorflow","tf_keras"}}
class Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCKED:
            raise ImportError("blocked: " + name)
        return None
sys.meta_path.insert(0, Blocker())
sys.path.insert(0, {package_parent!r})
from roast_py import roast
print(type(roast).__name__)
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "function"


def test_lazy_attribute_access_returns_the_function_not_the_submodule():
    """roast_py.roast is both a submodule and a function name.

    Importing the submodule binds it as a package attribute, so without
    explicit caching the first `from roast_py import roast` would give the
    function and every later one the module.
    """
    package_parent = str(Path(__file__).resolve().parents[1])
    code = (
        "import sys;"
        f"sys.path.insert(0, {package_parent!r});"
        "from roast_py import roast as r1;"
        "from roast_py import roast as r2;"
        "import roast_py;"
        "print(type(r1).__name__, type(r2).__name__, r1 is r2, type(roast_py.roast).__name__)"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "function function True function"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def test_cli_reports_missing_and_exits_nonzero(monkeypatch, capsys):
    from roast_py.dependencies import _main

    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (FAKE,))
    assert _main([]) == 1
    assert "not-real-xyz" in capsys.readouterr().out


def test_cli_reports_success_and_exits_zero(capsys):
    from roast_py.dependencies import _main

    assert _main([]) == 0
    assert "All roast_py dependencies are installed" in capsys.readouterr().out


def test_conda_failure_falls_back_to_pip_for_those_packages(simulated_conda_env, monkeypatch):
    """A package that isn't on conda-forge must not dead-end the install.

    This is the safety net for the fact that conda channel contents can't
    be assumed: if conda can't install something, pip gets it instead.
    """
    monkeypatch.setattr("roast_py.dependencies.missing_dependencies", lambda *a, **k: [])
    runs = []

    class Result:
        def __init__(self, rc):
            self.returncode = rc

    def fake_run(cmd, **kwargs):
        runs.append(cmd)
        # conda fails, pip succeeds
        return Result(1 if "conda" in cmd[0] else 0)

    monkeypatch.setattr("roast_py.dependencies.subprocess.run", fake_run)

    from roast_py.dependencies import install_dependencies

    install_dependencies([FAKE, FAKE_PIP_ONLY], quiet=True)

    assert len(runs) == 2, "expected a conda attempt then a pip retry"
    assert "conda" in runs[0][0]
    pip_cmd = runs[1]
    assert pip_cmd[:3] == [sys.executable, "-m", "pip"]
    # The package conda failed on must appear in the pip retry.
    assert "not-real-xyz" in pip_cmd
    assert "also-not-real" in pip_cmd


# --------------------------------------------------------------------------
# TensorFlow / tf-keras version pairing
#
# tf-keras X.Y requires tensorflow X.Y.*, and tf-keras < 2.16 calls
# tf.compat.v2.__internal__.register_load_context_function, which
# TensorFlow removed in 2.16. A mismatched pair therefore fails with a bare
# AttributeError at import time.
# --------------------------------------------------------------------------


def _fake_versions(monkeypatch, tf_version, keras_version):
    def fake(dist):
        return {"tensorflow": tf_version, "tf-keras": keras_version}.get(dist)

    monkeypatch.setattr("roast_py.dependencies.installed_version", fake)


def test_tf_keras_spec_pins_to_the_installed_tensorflow(monkeypatch):
    from roast_py.dependencies import tf_keras_spec

    _fake_versions(monkeypatch, "2.16.1", None)
    assert tf_keras_spec() == "tf-keras>=2.16,<2.17"

    _fake_versions(monkeypatch, "2.19.0", None)
    assert tf_keras_spec() == "tf-keras>=2.19,<2.20"


def test_tf_keras_spec_falls_back_to_a_floor_when_tensorflow_is_absent(monkeypatch):
    from roast_py.dependencies import MIN_TF_KERAS, tf_keras_spec

    _fake_versions(monkeypatch, None, None)
    assert tf_keras_spec() == f"tf-keras>={MIN_TF_KERAS[0]}.{MIN_TF_KERAS[1]}"
    # Never below 2.16, whose predecessors call the removed TF internal.
    assert MIN_TF_KERAS >= (2, 16)


def test_matching_versions_are_not_reported_as_a_mismatch(monkeypatch):
    from roast_py.dependencies import tensorflow_keras_mismatch

    _fake_versions(monkeypatch, "2.21.0", "2.21.0")
    assert tensorflow_keras_mismatch() is None
    # Patch releases may differ; only major.minor has to agree.
    _fake_versions(monkeypatch, "2.16.2", "2.16.0")
    assert tensorflow_keras_mismatch() is None


def test_the_reported_error_case_is_detected(monkeypatch):
    """Old tf-keras beside a modern TensorFlow: the register_load_context_function case."""
    from roast_py.dependencies import tensorflow_keras_mismatch

    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    problem = tensorflow_keras_mismatch()

    assert problem is not None
    assert "2.19.0" in problem and "2.15.0" in problem
    assert "register_load_context_function" in problem  # names the symptom
    assert "tf-keras>=2.19,<2.20" in problem  # and the exact fix


def test_mismatch_is_not_reported_when_one_side_is_absent(monkeypatch):
    from roast_py.dependencies import tensorflow_keras_mismatch

    _fake_versions(monkeypatch, "2.19.0", None)
    assert tensorflow_keras_mismatch() is None
    _fake_versions(monkeypatch, None, "2.15.0")
    assert tensorflow_keras_mismatch() is None


def test_pip_install_command_pins_tf_keras(monkeypatch):
    from roast_py.dependencies import REQUIRED, pip_install_command

    _fake_versions(monkeypatch, "2.18.0", None)
    tf_keras_dep = next(d for d in REQUIRED if d.pip_name == "tf-keras")
    cmd = pip_install_command([tf_keras_dep])
    assert cmd[-1] == "tf-keras>=2.18,<2.19", "must pin, or pip can backtrack to a broken release"


def test_check_dependencies_raises_on_a_version_mismatch(monkeypatch):
    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    with pytest.raises(ImportError, match="register_load_context_function"):
        check_dependencies()


def test_ensure_dependencies_repairs_a_mismatch_when_allowed(monkeypatch):
    """The user's case: everything installed, but TF and tf-keras disagree."""
    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    monkeypatch.delenv(NO_AUTO_INSTALL_ENV_VAR, raising=False)

    runs = []

    class Result:
        returncode = 0

    def fake_run(cmd, **kwargs):
        runs.append(cmd)
        # After the repair, report matching versions.
        _fake_versions(monkeypatch, "2.19.0", "2.19.0")
        return Result()

    monkeypatch.setattr("roast_py.dependencies.subprocess.run", fake_run)

    ensure_dependencies(install_missing=True, quiet=True)

    assert len(runs) == 1, "expected exactly one repair install"
    assert "tf-keras>=2.19,<2.20" in runs[0]


def test_ensure_dependencies_reports_mismatch_instead_of_repairing_when_opted_out(monkeypatch):
    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    monkeypatch.setattr(
        "roast_py.dependencies.subprocess.run",
        lambda *a, **k: pytest.fail("must not install when opted out"),
    )

    with pytest.raises(ImportError, match="register_load_context_function"):
        ensure_dependencies(install_missing=False)


# --------------------------------------------------------------------------
# TensorFlow newer than any released tf-keras
#
# tf-keras trails TensorFlow, and conda-forge can ship a TensorFlow ahead of
# PyPI, so pinning tf-keras to the installed TensorFlow's minor can name a
# release that does not exist (e.g. TF 2.22 -> 'tf-keras>=2.22,<2.23',
# which pip fails on).
# --------------------------------------------------------------------------


def test_alignment_specs_pin_neither_side(monkeypatch):
    from roast_py.dependencies import MIN_TF_KERAS, tensorflow_alignment_specs

    specs = tensorflow_alignment_specs()
    floor = f"{MIN_TF_KERAS[0]}.{MIN_TF_KERAS[1]}"
    assert specs == [f"tensorflow>={floor}", f"tf-keras>={floor}"]
    # No hard-coded ceiling: tf-keras's own tensorflow<X.(Y+1) requirement
    # caps TensorFlow, so this can't go stale when a new tf-keras ships.
    assert not any("<" in s for s in specs)


def test_unsatisfiable_pin_falls_back_to_joint_resolution(monkeypatch):
    """The reported failure: pip install 'tf-keras>=2.22,<2.23' exits 1."""
    from roast_py.dependencies import install_matching_tf_keras

    _fake_versions(monkeypatch, "2.22.0", None)
    runs = []
    state = {"aligned": False}

    class Result:
        def __init__(self, rc):
            self.returncode = rc

    def fake_run(cmd, **kwargs):
        runs.append(cmd)
        if "tf-keras>=2.22,<2.23" in cmd:
            return Result(1)  # no such release
        # The joint resolve succeeds and lands on a matching pair.
        state["aligned"] = True
        _fake_versions(monkeypatch, "2.21.0", "2.21.0")
        return Result(0)

    monkeypatch.setattr("roast_py.dependencies.subprocess.run", fake_run)

    install_matching_tf_keras(quiet=True)

    assert len(runs) == 2, "expected the pin attempt, then the joint resolve"
    assert "tf-keras>=2.22,<2.23" in runs[0]
    assert "tensorflow>=2.16" in runs[1] and "tf-keras>=2.16" in runs[1]
    assert state["aligned"]


def test_both_approaches_failing_explains_instead_of_raising_calledprocesserror(monkeypatch):
    from roast_py.dependencies import install_matching_tf_keras

    _fake_versions(monkeypatch, "2.22.0", None)

    class Result:
        returncode = 1

    monkeypatch.setattr("roast_py.dependencies.subprocess.run", lambda *a, **k: Result())

    with pytest.raises(RuntimeError) as excinfo:
        install_matching_tf_keras(quiet=True)

    message = str(excinfo.value)
    assert "2.22.0" in message
    assert "conda install" in message  # conda-native way down
    assert "tensorflow>=2.16" in message  # and the pip way
    assert "pypi.org/project/tf-keras" in message
    # Must not be a bare subprocess error.
    assert "non-zero exit status" not in message


def test_a_pin_that_works_does_not_touch_tensorflow(monkeypatch):
    """Don't downgrade a (possibly conda-managed) TensorFlow unnecessarily."""
    from roast_py.dependencies import install_matching_tf_keras

    _fake_versions(monkeypatch, "2.19.0", None)
    runs = []

    class Result:
        returncode = 0

    def fake_run(cmd, **kwargs):
        runs.append(cmd)
        _fake_versions(monkeypatch, "2.19.0", "2.19.0")
        return Result()

    monkeypatch.setattr("roast_py.dependencies.subprocess.run", fake_run)
    install_matching_tf_keras(quiet=True)

    assert len(runs) == 1, "a working pin should not escalate"
    assert "tf-keras>=2.19,<2.20" in runs[0]
    assert not any("tensorflow" in part for part in runs[0][4:]), "must not reinstall tensorflow"
