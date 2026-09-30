"""Tests the dependency preflight and installer.

roast_py pins its runtime to one exact environment (TESTED_ENVIRONMENT /
requirements-lock.txt) that ran the full pipeline end to end, and roast()
installs it with pip when the running interpreter differs. These tests
cover the declarations staying consistent with each other and with the
code, the version checks, and the install/opt-out decisions -- the latter
with pip stubbed out, since a real install downloads TensorFlow.
"""

import re
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

import roast_py.dependencies as deps
from roast_py.dependencies import (
    NO_AUTO_INSTALL_ENV_VAR,
    REQUIRED,
    SUPPORTED_PYTHON,
    TESTED_ENVIRONMENT,
    TESTED_PYTHON,
    Dependency,
    auto_install_disabled,
    check_dependencies,
    ensure_dependencies,
    format_missing,
    install_command,
    install_dependencies,
    lock_file,
    lock_specs,
    missing_dependencies,
    python_version_problem,
    tensorflow_keras_mismatch,
    version_drift,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = PROJECT_ROOT / "pyproject.toml"
ENVIRONMENT_YML = PROJECT_ROOT / "environment.yml"

FAKE = Dependency("definitely_not_a_real_module_xyz", "not-real-xyz", "testing")


def _normalize(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


# --------------------------------------------------------------------------
# declarations stay consistent
# --------------------------------------------------------------------------


def _declared_runtime_deps() -> list[str]:
    with open(PYPROJECT, "rb") as f:
        return tomllib.load(f)["project"]["dependencies"]


def _requirement_name(spec: str) -> str:
    return re.match(r"[A-Za-z0-9_.-]+", spec).group(0)


def test_required_list_matches_pyproject_runtime_dependencies():
    declared = {_requirement_name(s) for s in _declared_runtime_deps()}
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
    package_root = PROJECT_ROOT / "roast_py"
    sources = "\n".join(
        p.read_text() for p in package_root.rglob("*.py") if "dependencies.py" not in str(p)
    )
    for dep in REQUIRED:
        if dep.import_name == "openpyxl":
            continue
        imported = f"import {dep.import_name}" in sources or f"from {dep.import_name}" in sources
        assert imported, f"{dep.pip_name} is declared required but never imported by the package"


def test_lock_file_matches_tested_environment():
    lines = [
        line.strip()
        for line in lock_file().read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert lines == lock_specs(), (
        "requirements-lock.txt and TESTED_ENVIRONMENT in roast_py/dependencies.py "
        "have drifted apart"
    )


def test_every_direct_dependency_is_pinned():
    pinned = {_normalize(name) for name in TESTED_ENVIRONMENT}
    for dep in REQUIRED:
        assert _normalize(dep.pip_name) in pinned, f"{dep.pip_name} has no tested version"


def test_lock_satisfies_pyproject_ranges():
    requirements = pytest.importorskip("packaging.requirements")
    pinned = {_normalize(k): v for k, v in TESTED_ENVIRONMENT.items()}
    for spec in _declared_runtime_deps():
        req = requirements.Requirement(spec)
        version = pinned[_normalize(req.name)]
        assert req.specifier.contains(version), f"lock pins {req.name}=={version}, outside {spec}"


def test_lock_pins_tensorflow_and_tf_keras_to_the_same_minor():
    tf = TESTED_ENVIRONMENT["tensorflow"].split(".")[:2]
    keras = TESTED_ENVIRONMENT["tf-keras"].split(".")[:2]
    assert tf == keras


def test_tested_python_is_inside_the_supported_range():
    tested = tuple(int(p) for p in TESTED_PYTHON.split("."))
    assert python_version_problem(tested) is None


def test_pyproject_requires_python_matches_supported_range():
    specifiers = pytest.importorskip("packaging.specifiers")
    with open(PYPROJECT, "rb") as f:
        spec = specifiers.SpecifierSet(tomllib.load(f)["project"]["requires-python"])
    lo, hi = SUPPORTED_PYTHON
    assert f"{lo[0]}.{lo[1]}.0" in spec
    assert f"{hi[0]}.{hi[1]}.9" in spec
    assert f"{lo[0]}.{lo[1] - 1}.9" not in spec
    assert f"{hi[0]}.{hi[1] + 1}.0" not in spec


def test_environment_yml_uses_the_tested_python_and_the_lock_file():
    text = ENVIRONMENT_YML.read_text()
    tested_minor = ".".join(TESTED_PYTHON.split(".")[:2])
    assert f"- python={tested_minor}\n" in text
    assert "-r requirements-lock.txt" in text


# --------------------------------------------------------------------------
# detection
# --------------------------------------------------------------------------


def test_this_environment_matches_the_tested_one():
    """The environment these tests run in is the one the lock was taken from."""
    assert version_drift() == []
    assert missing_dependencies() == []
    assert check_dependencies() == []


def test_missing_dependencies_reports_unimportable_packages():
    assert missing_dependencies((FAKE,)) == [FAKE]


@pytest.mark.parametrize(
    "version, supported",
    [
        ((3, 10, 14), False),
        ((3, 11, 0), True),
        ((3, 12, 5), True),
        ((3, 13, 1), True),
        ((3, 14, 0), False),
    ],
)
def test_python_version_problem(version, supported):
    problem = python_version_problem(version)
    assert (problem is None) == supported
    if problem:
        assert "conda create -n roast_py python=3.11" in problem


def test_version_drift_reports_differences_and_absences(monkeypatch):
    scipy_version = deps.installed_version("scipy")
    monkeypatch.setattr(
        deps,
        "TESTED_ENVIRONMENT",
        {"numpy": "0.0.1", "not-real-xyz": "1.0", "scipy": scipy_version},
    )
    drift = version_drift()
    assert ("numpy", deps.installed_version("numpy"), "0.0.1") in drift
    assert ("not-real-xyz", None, "1.0") in drift
    assert all(name != "scipy" for name, _, _ in drift)


def _fake_versions(monkeypatch, tf_version, keras_version):
    def fake(dist):
        return {"tensorflow": tf_version, "tf-keras": keras_version}.get(dist)

    monkeypatch.setattr(deps, "installed_version", fake)


def test_matching_versions_are_not_reported_as_a_mismatch(monkeypatch):
    _fake_versions(monkeypatch, "2.21.0", "2.21.0")
    assert tensorflow_keras_mismatch() is None
    _fake_versions(monkeypatch, "2.16.2", "2.16.0")  # only major.minor must agree
    assert tensorflow_keras_mismatch() is None


def test_the_register_load_context_function_case_is_detected(monkeypatch):
    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    problem = tensorflow_keras_mismatch()
    assert problem is not None
    assert "2.19.0" in problem and "2.15.0" in problem
    assert "register_load_context_function" in problem
    assert "python -m roast_py.dependencies --install" in problem


def test_mismatch_is_not_reported_when_one_side_is_absent(monkeypatch):
    _fake_versions(monkeypatch, "2.19.0", None)
    assert tensorflow_keras_mismatch() is None


# --------------------------------------------------------------------------
# commands and messages
# --------------------------------------------------------------------------


def test_install_command_pins_everything_into_the_running_interpreter():
    cmd = install_command()
    assert cmd[:4] == [sys.executable, "-m", "pip", "install"]
    assert cmd[4:] == lock_specs()
    assert all("==" in spec for spec in cmd[4:])
    assert "tensorflow==2.21.0" in cmd and "tf-keras==2.21.0" in cmd


def test_format_missing_lists_reason_for_each_package():
    text = format_missing([FAKE])
    assert "not-real-xyz" in text and "testing" in text
    assert "python -m roast_py.dependencies --install" in text


@pytest.mark.parametrize(
    "value, disabled",
    [("", False), ("0", False), ("false", False), ("False", False), ("1", True), ("yes", True)],
)
def test_auto_install_disabled_parses_the_env_var(monkeypatch, value, disabled):
    monkeypatch.setenv(NO_AUTO_INSTALL_ENV_VAR, value)
    assert auto_install_disabled() is disabled


# --------------------------------------------------------------------------
# ensure_dependencies: what roast() calls first
# --------------------------------------------------------------------------

SOME_DRIFT = [("numpy", "1.26.4", TESTED_ENVIRONMENT["numpy"])]


@pytest.fixture
def no_pip(monkeypatch):
    """Fails the test if anything tries to run a subprocess."""

    def forbidden(*args, **kwargs):
        raise AssertionError(f"unexpected subprocess: {args}")

    monkeypatch.setattr(deps.subprocess, "run", forbidden)
    monkeypatch.delenv(NO_AUTO_INSTALL_ENV_VAR, raising=False)


def test_ensure_dependencies_is_a_noop_in_the_tested_environment(no_pip):
    ensure_dependencies()


def test_ensure_dependencies_installs_when_versions_differ(no_pip, monkeypatch):
    calls = []
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    monkeypatch.setattr(deps, "install_dependencies", lambda quiet=False: calls.append(quiet))
    ensure_dependencies()
    assert calls == [False]


def test_ensure_dependencies_refuses_unsupported_python_before_installing(no_pip, monkeypatch):
    monkeypatch.setattr(deps, "python_version_problem", lambda: "needs Python 3.11-3.13")
    monkeypatch.setattr(deps, "install_dependencies", lambda quiet=False: pytest.fail("installed"))
    with pytest.raises(RuntimeError, match="3.11-3.13"):
        ensure_dependencies()


def test_opted_out_missing_packages_raise(no_pip, monkeypatch):
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    monkeypatch.setattr(deps, "REQUIRED", (FAKE,))
    with pytest.raises(ImportError, match="not-real-xyz"):
        ensure_dependencies(install_missing=False)


def test_opted_out_tensorflow_mismatch_raises(no_pip, monkeypatch):
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    _fake_versions(monkeypatch, "2.19.0", "2.15.0")
    with pytest.raises(ImportError, match="register_load_context_function"):
        ensure_dependencies(install_missing=False)


def test_opted_out_other_drift_only_warns(no_pip, monkeypatch):
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    with pytest.warns(UserWarning, match="numpy"):
        ensure_dependencies(install_missing=False)


def test_env_var_opt_out_behaves_like_install_missing_false(no_pip, monkeypatch):
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    monkeypatch.setenv(NO_AUTO_INSTALL_ENV_VAR, "1")
    monkeypatch.setattr(deps, "install_dependencies", lambda quiet=False: pytest.fail("installed"))
    with pytest.warns(UserWarning):
        ensure_dependencies()


# --------------------------------------------------------------------------
# install_dependencies (pip stubbed)
# --------------------------------------------------------------------------


class _FakePip:
    def __init__(self, returncode=0):
        self.returncode = returncode
        self.commands = []

    def __call__(self, cmd, *args, **kwargs):
        self.commands.append(cmd)
        return subprocess.CompletedProcess(cmd, self.returncode)


def _drift_sequence(monkeypatch, *results):
    """version_drift() returns each of `results` in turn, then the last forever."""
    results = list(results)

    def fake():
        return results.pop(0) if len(results) > 1 else results[0]

    monkeypatch.setattr(deps, "version_drift", fake)


def test_install_runs_one_pip_command_with_every_pin(monkeypatch, capsys):
    pip = _FakePip()
    monkeypatch.setattr(deps.subprocess, "run", pip)
    monkeypatch.setattr(deps, "smoke_test", lambda: None)
    _drift_sequence(monkeypatch, SOME_DRIFT, [])

    install_dependencies()

    assert pip.commands == [install_command()]
    assert "tested environment is installed" in capsys.readouterr().out


def test_install_is_skipped_when_nothing_differs(monkeypatch):
    monkeypatch.setattr(deps.subprocess, "run", lambda *a, **k: pytest.fail("ran pip"))
    _drift_sequence(monkeypatch, [])
    install_dependencies(quiet=True)


def test_failed_pip_install_recommends_a_fresh_environment(monkeypatch):
    monkeypatch.setattr(deps.subprocess, "run", _FakePip(returncode=1))
    monkeypatch.setattr(deps, "smoke_test", lambda: None)
    _drift_sequence(monkeypatch, SOME_DRIFT, SOME_DRIFT)

    with pytest.raises(RuntimeError) as excinfo:
        install_dependencies(quiet=True)
    message = str(excinfo.value)
    assert "pip exited with status 1" in message
    assert "numpy" in message
    assert "conda create -n roast_py python=3.11" in message


def test_broken_imports_after_install_are_reported(monkeypatch):
    monkeypatch.setattr(deps.subprocess, "run", _FakePip())
    monkeypatch.setattr(
        deps, "smoke_test", lambda: "AttributeError: register_load_context_function"
    )
    _drift_sequence(monkeypatch, SOME_DRIFT, [])

    with pytest.raises(RuntimeError, match="register_load_context_function"):
        install_dependencies(quiet=True)


def test_install_notes_conda_owned_packages_it_replaces(monkeypatch, capsys):
    monkeypatch.setattr(deps.subprocess, "run", _FakePip())
    monkeypatch.setattr(deps, "smoke_test", lambda: None)
    monkeypatch.setattr(deps, "installed_by", lambda name: "conda")
    _drift_sequence(monkeypatch, [("tensorflow", "2.22.0", "2.21.0")], [])

    install_dependencies()
    assert "replacing conda-installed tensorflow" in capsys.readouterr().out


def test_install_refuses_unsupported_python(monkeypatch):
    monkeypatch.setattr(deps, "python_version_problem", lambda: "needs Python 3.11-3.13")
    monkeypatch.setattr(deps.subprocess, "run", lambda *a, **k: pytest.fail("ran pip"))
    with pytest.raises(RuntimeError, match="3.11-3.13"):
        install_dependencies()


@pytest.mark.slow
def test_smoke_test_passes_in_the_tested_environment():
    assert deps.smoke_test() is None


# --------------------------------------------------------------------------
# importable with nothing installed
# --------------------------------------------------------------------------


def test_dependencies_module_imports_without_any_third_party_packages():
    """The whole point of this module: it must import in a bare environment."""
    code = (
        "import sys;"
        f"sys.path.insert(0, {str(PROJECT_ROOT)!r});"
        "import roast_py.dependencies as d;"
        "print(len(d.TESTED_ENVIRONMENT))"
    )
    result = subprocess.run([sys.executable, "-S", "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(len(TESTED_ENVIRONMENT))


def test_roast_is_importable_with_no_dependencies_installed():
    """`from roast_py import roast` (and review_res) must work in a bare environment.

    roast() is what installs the dependencies, so if importing it required
    them first, the auto-install could never run. This is the property
    that keeps roast.py's imports inside the function body.
    """
    code = f"""
import sys
BLOCKED = {{"numpy","scipy","nibabel","pandas","openpyxl","skimage","tensorflow","tf_keras",
           "matplotlib","pyvista","vtk"}}
class Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCKED:
            raise ImportError("blocked: " + name)
        return None
sys.meta_path.insert(0, Blocker())
sys.path.insert(0, {str(PROJECT_ROOT)!r})
from roast_py import roast, review_res
import roast_py.viz
print(type(roast).__name__, type(review_res).__name__)
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "function function"


@pytest.mark.parametrize(
    "statements",
    [
        "from roast_py import DEFAULT_RECIPE, roast",  # another name resolved first
        "import roast_py.roast; from roast_py.viz import review; from roast_py import roast",
        "from roast_py.viz import SliceViewer, sliceshow",
    ],
)
def test_public_functions_are_never_shadowed_by_submodules(statements):
    code = (
        "import sys, types;"
        f"sys.path.insert(0, {str(PROJECT_ROOT)!r});"
        f"{statements};"
        "import roast_py, roast_py.viz;"
        "print(isinstance(roast_py.roast, types.FunctionType),"
        " isinstance(roast_py.viz.sliceshow, types.FunctionType))"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "True True"


def test_lazy_attribute_access_returns_the_function_not_the_submodule():
    """Repeated `from roast_py import roast` must keep returning the function."""
    code = (
        "import sys;"
        f"sys.path.insert(0, {str(PROJECT_ROOT)!r});"
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


def test_cli_reports_drift_and_exits_nonzero(monkeypatch, capsys):
    monkeypatch.setattr(deps, "version_drift", lambda: SOME_DRIFT)
    assert deps._main([]) == 1
    out = capsys.readouterr().out
    assert "numpy" in out and "--install" in out


def test_cli_reports_success_and_exits_zero(capsys):
    assert deps._main([]) == 0
    assert "tested environment is installed" in capsys.readouterr().out
