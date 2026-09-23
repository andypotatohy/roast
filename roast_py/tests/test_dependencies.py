"""Tests the dependency preflight.

The drift test matters: the bug this module exists to fix was
pyproject.toml declaring a set of dependencies that didn't match what the
code actually imports, so `pip install -e .` produced an install that
couldn't run roast().
"""

import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

from roast_py.dependencies import (
    REQUIRED,
    Dependency,
    check_dependencies,
    format_missing,
    install_command,
    missing_dependencies,
)

PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def _declared_runtime_deps() -> set[str]:
    with open(PYPROJECT, "rb") as f:
        data = tomllib.load(f)
    # Strip version specifiers: "numpy>=1.24" -> "numpy"
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


def test_all_required_dependencies_are_importable_in_this_environment():
    # The test environment installs the package, so nothing should be missing.
    assert missing_dependencies() == []
    assert check_dependencies() == []


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
        # Both spellings count: `import skimage` and `from skimage.transform import ...`
        imported = (
            f"import {dep.import_name}" in sources or f"from {dep.import_name}" in sources
        )
        assert imported, f"{dep.pip_name} is declared required but never imported by the package"


def test_missing_dependencies_reports_unimportable_packages():
    fake = Dependency("definitely_not_a_real_module_xyz", "not-real-xyz", "testing")
    missing = missing_dependencies((fake,))
    assert missing == [fake]


def test_check_dependencies_raises_with_actionable_message(monkeypatch):
    fake = Dependency("definitely_not_a_real_module_xyz", "not-real-xyz", "testing")
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (fake,))

    with pytest.raises(ImportError) as excinfo:
        check_dependencies()

    message = str(excinfo.value)
    assert "not-real-xyz" in message
    assert "pip install" in message
    assert "python -m roast_py.dependencies --install" in message


def test_install_command_targets_the_running_interpreter():
    fake = Dependency("x", "x-pkg", "testing")
    cmd = install_command([fake])
    # Must be sys.executable -m pip, not a bare `pip`, so it installs into
    # the active environment rather than whatever `pip` resolves to first.
    assert cmd[:4] == [sys.executable, "-m", "pip", "install"]
    assert cmd[-1] == "x-pkg"


def test_format_missing_lists_reason_for_each_package():
    fake = Dependency("x", "x-pkg", "doing the thing")
    text = format_missing([fake])
    assert "x-pkg" in text
    assert "doing the thing" in text


def test_dependencies_module_imports_without_any_third_party_packages():
    """The whole point of this module: it must import in a bare environment.

    Runs in a subprocess with a stubbed-out import path so that numpy,
    tensorflow et al. are unavailable, proving the module has no
    third-party imports of its own.
    """
    package_parent = str(Path(__file__).resolve().parents[1])
    code = (
        "import sys;"
        f"sys.path.insert(0, {package_parent!r});"
        "import roast_py.dependencies as d;"
        "print(len(d.REQUIRED))"
    )
    result = subprocess.run(
        [sys.executable, "-S", "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(len(REQUIRED))


def test_cli_reports_missing_and_exits_nonzero(monkeypatch, capsys):
    from roast_py.dependencies import _main

    fake = Dependency("definitely_not_a_real_module_xyz", "not-real-xyz", "testing")
    monkeypatch.setattr("roast_py.dependencies.REQUIRED", (fake,))

    assert _main([]) == 1
    assert "not-real-xyz" in capsys.readouterr().out


def test_cli_reports_success_and_exits_zero(capsys):
    from roast_py.dependencies import _main

    assert _main([]) == 0
    assert "All roast_py dependencies are installed" in capsys.readouterr().out
