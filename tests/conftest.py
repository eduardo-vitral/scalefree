from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess

import pytest


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def scalefree_exe() -> Path:
    """Provide path to the ScaleFree Fortran executable.

    Behaviour:
      - If fortran_src/scalefree.e exists -> use it.
      - Else, compile fortran_src/scalefree.f into scalefree.e (requires gfortran).
      - If gfortran is not available -> skip tests that require the executable.

    Notes:
      - CI runners are typically clean, so compilation will happen there.
      - Locally, you may already have a built executable.
    """

    root = repo_root()
    exe = root / "fortran_src" / "scalefree.e"
    src = root / "fortran_src" / "scalefree.f"

    if exe.exists():
        return exe

    if not src.exists():
        pytest.skip(f"Missing Fortran source at {src}")

    gfortran = shutil.which("gfortran")
    if not gfortran:
        pytest.skip(
            "gfortran not available, cannot build ScaleFree executable. "
            "Install gfortran (e.g. apt-get install gfortran) or provide a prebuilt fortran_src/scalefree.e."
        )

    cmd = [
        gfortran,
        "-O2",
        "-std=legacy",
        "-ffixed-line-length-none",
        "-o",
        str(exe),
        str(src),
    ]

    try:
        subprocess.run(
            cmd,
            cwd=str(src.parent),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "Failed to compile ScaleFree Fortran backend.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Output:\n{e.stdout}"
        ) from e

    if not exe.exists():
        raise RuntimeError(f"Compilation reported success, but executable not found at {exe}")

    return exe


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    """Directory holding golden reference outputs."""
    d = repo_root() / "tests" / "data"
    if not d.exists():
        pytest.skip("Missing tests/data directory. Generate refs with: python tests/make_vprofile_refs.py")
    return d


@pytest.fixture(scope="session")
def include_both() -> bool:
    """Whether to include the optional kinematics='both' regression cases.

    Enabled by setting environment variable:
      SCALEFREE_TEST_INCLUDE_BOTH=1
    """
    return os.environ.get("SCALEFREE_TEST_INCLUDE_BOTH", "0") == "1"
