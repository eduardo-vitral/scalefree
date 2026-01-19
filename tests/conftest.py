from __future__ import annotations

from pathlib import Path
import pytest


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def scalefree_exe() -> Path:
    exe = repo_root() / "fortran_src" / "scalefree.e"
    if not exe.exists():
        pytest.skip(
            f"ScaleFree Fortran executable not found at {exe}. "
            "Compile it first in CI or locally (gfortran ...)."
        )
    return exe


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    d = repo_root() / "tests" / "data"
    if not d.exists():
        pytest.skip("Missing tests/data directory with reference outputs.")
    return d
