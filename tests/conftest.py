from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def scalefree_exe() -> Path:
    """
    Provide the path to the ScaleFree Fortran executable.

    Behaviour:
      - If fortran_src/scalefree.e exists -> use it.
      - Else, if gfortran is available,
      compile fortran_src/scalefree.f into scalefree.e.
      - Else, skip tests that depend on the executable.
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
            "Install gfortran (e.g. apt-get install gfortran) "
            "or commit a prebuilt fortran_src/scalefree.e."
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
        raise RuntimeError(
            f"Compilation reported success, but executable not found at {exe}"
        )

    return exe


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    """
    Directory holding golden reference outputs.
    We do not auto-generate refs here; tests will skip if refs are missing and
    tell you how to generate them using tests/make_vprofile_refs.py.
    """
    d = repo_root() / "tests" / "data"
    if not d.exists():
        pytest.skip(
            "Missing tests/data directory. "
            + "Create it and "
            + "add reference outputs."
        )
    return d
