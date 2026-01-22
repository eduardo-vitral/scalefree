from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def scalefree_exe(tmp_path_factory) -> Path:
    """
    Provide the path to the ScaleFree Fortran executable.

    Behaviour:
      - If fortran_src/scalefree.e exists -> use it.
      - Else, if gfortran is available, compile scalefree.f into a session tmp dir.
      - Else, skip tests that depend on the executable.
    """
    root = repo_root()
    prebuilt = root / "fortran_src" / "scalefree.e"
    src = root / "fortran_src" / "scalefree.f"

    if prebuilt.exists():
        return prebuilt

    if not src.exists():
        pytest.skip(f"Missing Fortran source at {src}")

    gfortran = shutil.which("gfortran")
    if not gfortran:
        pytest.skip(
            "gfortran not available, cannot build ScaleFree executable. "
            "Install gfortran (e.g. apt-get install gfortran) "
            "or commit a prebuilt fortran_src/scalefree.e."
        )

    build_dir = tmp_path_factory.mktemp("scalefree_build")
    exe = build_dir / "scalefree.e"

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
    """
    Directory holding golden reference outputs.

    Used only when SCALEFREE_STRICT_TESTS=1.
    """
    d = repo_root() / "tests" / "data"
    if not d.exists():
        # Do not skip: strict-mode tests will skip per-file if missing.
        # Keep fixture available.
        return d
    return d
