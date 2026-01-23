from __future__ import annotations

import os
from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--scalefree-exe",
        action="store",
        default=None,
        help=(
            "Path to the compiled scalefree Fortran executable. "
            "Defaults to <repo>/fortran_src/scalefree.e if present."
        ),
    )


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def scalefree_exe(pytestconfig: pytest.Config, repo_root: Path) -> Path:
    """Return a path to the compiled Fortran executable.

    CI builds it via:
      gfortran -O2 -std=legacy -o fortran_src/scalefree.e fortran_src/scalefree.f

    Locally you can either compile the same way, set SCALEFREE_EXE, or pass
    --scalefree-exe.
    """

    # 1) CLI override
    cli = pytestconfig.getoption("--scalefree-exe")
    if cli:
        p = Path(cli).expanduser().resolve()
        if not p.exists():
            pytest.skip(f"scalefree executable not found at: {p}")
        return p

    # 2) Env var override
    env = os.environ.get("SCALEFREE_EXE")
    if env:
        p = Path(env).expanduser().resolve()
        if not p.exists():
            pytest.skip(f"SCALEFREE_EXE points to missing file: {p}")
        return p

    # 3) Default expected location
    p = (repo_root / "fortran_src" / "scalefree.e").resolve()
    if not p.exists():
        pytest.skip(
            "Compiled scalefree executable not found. "
            "Build it with: gfortran -O2 -std=legacy -o fortran_src/scalefree.e fortran_src/scalefree.f "
            "or pass --scalefree-exe / set SCALEFREE_EXE."
        )
    return p


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    return Path(__file__).resolve().parent / "data"
