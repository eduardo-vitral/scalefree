from __future__ import annotations

from pathlib import Path
import os
import random
import shutil
import subprocess

import numpy as np
import pytest


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def resources() -> Path:
    """Path to on-disk resources used by regression tests.

    Reference outputs are stored under tests/data/.
    """
    return repo_root() / "tests" / "data"


@pytest.fixture(scope="session")
def ref_dir(resources: Path) -> Path:
    """Alias for ``resources`` used by the regression tests."""
    return resources


@pytest.fixture(scope="session", autouse=True)
def _seed_rngs() -> None:
    """Seed RNGs to reduce test flakiness.

    Tests are intended to be deterministic. If any underlying code path uses
    randomness (directly or indirectly), we seed both Python's `random` and
    NumPy's global RNG.

    Set SCALEFREE_TEST_SEED to override the default.
    """

    seed_env = os.environ.get("SCALEFREE_TEST_SEED", "0")
    try:
        seed = int(seed_env)
    except ValueError:
        seed = 0

    random.seed(seed)
    np.random.seed(seed)


def pytest_addoption(parser):
    parser.addoption(
        "--build",
        action="store_true",
        default=False,
        help="Build Fortran executable before running tests",
    )


@pytest.fixture(scope="session")
def build_fortran_executable(pytestconfig):
    """Optionally build the Fortran executable before tests run."""
    if not pytestconfig.getoption("--build"):
        return

    root = repo_root()
    src = root / "fortran_src" / "scalefree.f"
    # Keep in sync with CI (GitLab) build step.
    exe = root / "fortran_src" / "scalefree.e"

    if not src.exists():
        raise FileNotFoundError(f"Missing Fortran source: {src}")

    # Use gfortran if available
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        raise RuntimeError(
            "gfortran is required to build the Fortran executable. "
            "Install gfortran or run tests without --build."
        )

    cmd = [gfortran, "-O2", "-std=legacy", "-o", str(exe), str(src)]
    subprocess.check_call(cmd, cwd=str(root))


@pytest.fixture(scope="session")
def scalefree_exe(build_fortran_executable):
    """Return the path to the Fortran executable, if present.

    Historically the executable was named 'fitvp'. The CI pipeline currently
    builds 'scalefree.e'. We accept either to keep the test structure stable.
    """

    root = repo_root() / "fortran_src"
    for candidate in (root / "scalefree.e", root / "fitvp"):
        if candidate.exists():
            return candidate

    pytest.skip(
        "Fortran executable not found (expected 'fortran_src/scalefree.e' or "
        "'fortran_src/fitvp'). Build it locally or run tests with --build."
    )
