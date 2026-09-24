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
