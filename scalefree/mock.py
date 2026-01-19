from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Tuple, Union

import os
import sys
import subprocess
import tempfile

import numpy as np


@dataclass(frozen=True)
class MockStats:
    binney_beta: float
    rms_method1: float
    rms_method2: float
    mean_vphi: float


def _fortran_src_dir() -> Path:
    """
    Locate fortran_src relative to the installed package.
    Expects repository layout:
      scalefree/
        fortran_src/
        scalefree/
          __init__.py
          ...
    i.e. this file is scalefree/scalefree/mock.py
    """
    pkg_dir = Path(__file__).resolve().parent
    return (pkg_dir.parent / "fortran_src").resolve()


def _parse_kind_block(stdout: str, kind: str) -> np.ndarray:
    """
    Parse a structured stdout block:

      # kind=<kind>
      # columns: ...
      <rows...>
      # kind=...

    Returns numeric array from rows.
    """
    lines = stdout.splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.strip() == f"# kind={kind}":
            start = i + 1
            break
    if start is None:
        raise RuntimeError(f"Could not find '# kind={kind}' in MCdraw output.")

    # skip header lines beginning with '#'
    j = start
    while j < len(lines) and lines[j].lstrip().startswith("#"):
        j += 1

    data_lines = []
    while j < len(lines):
        ln = lines[j].strip()
        if not ln:
            j += 1
            continue
        if ln.startswith("# kind="):
            break
        if ln.startswith("#"):
            j += 1
            continue
        data_lines.append(ln)
        j += 1

    if not data_lines:
        return np.empty((0, 0), dtype=float)

    return np.loadtxt(StringIO("\n".join(data_lines)))


def _build_executable(build_dir: Path) -> Path:
    """
    Build MCdraw_gh.e in an isolated temp directory.

    This avoids writing compiled artifacts into the repo/package directory.
    """
    src = _fortran_src_dir()
    needed = ["MCdraw_gh.f90", "mcdraw_legacy.f", "gh_velocity_sampler.py"]

    for name in needed:
        p = src / name
        if not p.exists():
            raise FileNotFoundError(
                f"Missing required Fortran source/helper: {p}",
            )

    # Copy into build dir
    for name in needed:
        (build_dir / name).write_bytes((src / name).read_bytes())

    # Compile legacy object
    subprocess.run(
        [
            "gfortran",
            "-O2",
            "-std=legacy",
            "-c",
            "mcdraw_legacy.f",
            "-o",
            "mcdraw_legacy.o",
        ],
        cwd=build_dir,
        check=True,
        capture_output=True,
        text=True,
    )

    # Compile+link main
    subprocess.run(
        [
            "gfortran",
            "-O2",
            "-std=f2008",
            "MCdraw_gh.f90",
            "mcdraw_legacy.o",
            "-o",
            "MCdraw_gh.e",
        ],
        cwd=build_dir,
        check=True,
        capture_output=True,
        text=True,
    )

    exe = build_dir / "MCdraw_gh.e"
    if not exe.exists():
        raise RuntimeError("Build succeeded but MCdraw_gh.e not found.")
    return exe


def mock(
    *,
    potential: int = 1,
    gamma: float = 4.0,
    q: float = 1.0,
    df: int = 1,
    beta: float = 0.0,
    s: float = 0.5,
    t: float = 0.0,
    nsamples: int = 50,
    seed: int = -101,
    rin: float = 1.0,
    rout: float = 1000.0,
    debug: bool = False,
    return_stats: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, MockStats]]:
    """
    Generate a mock 6D phase-space sample (x,y,z,vx,vy,vz)
    from the ScaleFree model.

    Notes
    -----
    - This currently uses GH sampling via the Python helper, but with h3=h4=0
      unless you later implement gh_moments_hook() to populate those.
    - The Fortran executable is built in a temporary directory per call.
      (We can later add a cache to avoid recompilation.)
    """
    # Prepare stdin in the interactive order expected by MCdraw_gh.e
    input_lines = [
        str(int(potential)),
        f"{float(gamma)}",
        f"{float(q)}",
        str(int(df)),
        f"{float(beta)}",
        f"{float(s)}",
        f"{float(t)}",
        str(int(nsamples)),
        str(int(seed)),
        f"{float(rin)}",
        f"{float(rout)}",
    ]
    stdin = "\n".join(input_lines) + "\n"

    with tempfile.TemporaryDirectory(prefix="scalefree_mcdraw_") as td:
        build_dir = Path(td)
        exe = _build_executable(build_dir)

        env = os.environ.copy()
        env["PYTHON"] = sys.executable  # ensure the same venv runs balrogo

        p = subprocess.run(
            [str(exe)],
            cwd=build_dir,
            input=stdin,
            capture_output=True,
            text=True,
            env=env,
            check=True,
        )

        if debug:
            print(p.stdout)
            if p.stderr.strip():
                print(p.stderr, file=sys.stderr)

        data = _parse_kind_block(p.stdout, "mock")
        if data.ndim == 1 and data.size == 0:
            data = data.reshape((0, 6))

        if data.shape[1] != 6:
            raise RuntimeError(
                f"Expected 6 columns in mock output; got shape {data.shape}"
            )

        if not return_stats:
            return data

        stats_arr = _parse_kind_block(p.stdout, "mock_stats").reshape(-1)
        if stats_arr.size < 4:
            raise RuntimeError("mock_stats block did not contain 4 values.")

        stats = MockStats(
            binney_beta=float(stats_arr[0]),
            rms_method1=float(stats_arr[1]),
            rms_method2=float(stats_arr[2]),
            mean_vphi=float(stats_arr[3]),
        )
        return data, stats
