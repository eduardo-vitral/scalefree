"""Pytest configuration and shared fixtures.

Notes
-----
- The Fortran executable is expected at: `fortran_src/scalefree.e`.
- CI compiles it in `.gitlab-ci.yml`.
- When running locally, tests can compile it if it is missing.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


def _repo_root() -> Path:
    # tests/ is at <repo>/tests
    return Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def build_fortran_executable() -> Path:
    """Build the Fortran executable if it is not already present."""

    root = _repo_root()
    exe = root / "fortran_src" / "scalefree.e"
    src = root / "fortran_src" / "scalefree.f"

    # If CI (or the user) has already compiled it, do nothing.
    if exe.exists():
        return exe

    # If the user has provided a custom path, do not attempt to compile.
    env_exe = os.environ.get("SCALEFREE_EXE")
    if env_exe:
        p = Path(env_exe).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"SCALEFREE_EXE is set but not found: {p}")
        return p

    if not src.exists():
        raise FileNotFoundError(f"Missing Fortran source file: {src}")

    exe.parent.mkdir(parents=True, exist_ok=True)

    # Keep flags conservative and consistent with CI.
    cmd = [
        "gfortran",
        "-O2",
        "-std=legacy",
        "-o",
        str(exe),
        str(src),
    ]

    try:
        subprocess.run(cmd, check=True, cwd=str(root))
    except FileNotFoundError as e:
        raise RuntimeError(
            "gfortran not found. Install gfortran (or rely on CI) to run tests locally."
        ) from e

    if not exe.exists():
        raise RuntimeError(f"Fortran build reported success but executable not found: {exe}")

    return exe


@pytest.fixture(scope="session")
def scalefree_exe(build_fortran_executable: Path) -> Path:
    """Return the path to the compiled Fortran executable."""

    env_exe = os.environ.get("SCALEFREE_EXE")
    if env_exe:
        p = Path(env_exe).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"SCALEFREE_EXE is set but not found: {p}")
        return p
    return build_fortran_executable


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    """Directory containing stored vprofile reference outputs."""

    return Path(__file__).resolve().parent / "data"


@pytest.fixture(scope="session")
def _seed_rngs() -> None:
    """Keep RNGs deterministic in tests that use stochastic sampling.

    vprofile regression tests are deterministic, but other tests may depend on RNG state.
    """

    os.environ.setdefault("PYTHONHASHSEED", "0")
