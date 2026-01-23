from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess

import scalefree


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def build_exe(*, force: bool = False) -> Path:
    """Build (or reuse) the Fortran executable.

    By default we rebuild if:
      - executable is missing, OR
      - source is newer than executable.

    Use --force-rebuild to always rebuild.
    """
    root = repo_root()
    exe = root / "fortran_src" / "scalefree.e"
    src = root / "fortran_src" / "scalefree.f"

    if not src.exists():
        raise FileNotFoundError(f"Missing Fortran source at {src}")

    if exe.exists() and not force:
        try:
            if exe.stat().st_mtime >= src.stat().st_mtime:
                return exe
        except OSError:
            pass

    gfortran = shutil.which("gfortran")
    if not gfortran:
        raise RuntimeError(
            "gfortran not found. Install it, then re-run:\n"
            "  python tests/make_vprofile_refs.py"
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
    subprocess.run(cmd, cwd=str(src.parent), check=True)
    return exe


def case_cfg(algorithm: int) -> dict:
    # Keep aligned with tests.
    if algorithm == 1:
        return dict(maxmom=4)
    if algorithm == 2:
        return dict(maxmom=8, vp_reg_param=1.0)
    if algorithm == 3:
        return dict(maxmom=8, vp_smooth_eps=0.0)
    raise ValueError(f"Unsupported algorithm={algorithm}")


def ref_name(*, kinematics: str, average: bool, algorithm: int) -> str:
    stem = "avg" if average else "point"
    return f"{kinematics}_{stem}_alg{algorithm}_ref.txt"


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Generate/update golden reference outputs under tests/data/.\n\n"
            "Run this locally after changing Fortran output or the Python parser, "
            "then commit the updated *_ref.txt files."
        )
    )
    p.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Always rebuild fortran_src/scalefree.e from fortran_src/scalefree.f.",
    )
    p.add_argument(
        "--include-both",
        action="store_true",
        help="Also generate references for kinematics='both' (optional).",
    )

    args = p.parse_args()

    exe = build_exe(force=args.force_rebuild)
    root = repo_root()
    data_dir = root / "tests" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    runner = scalefree.ScaleFreeRunner(
        exe_path=exe,
        workdir=(root / "tests" / "_work"),
    )
    runner.workdir.mkdir(parents=True, exist_ok=True)

    kinematics_list = ["intrinsic", "projected"]
    if args.include_both:
        kinematics_list.append("both")

    for algorithm in (1, 2, 3):
        cfg = case_cfg(algorithm)
        for kinematics in kinematics_list:
            for average in (False, True):
                res = runner.vprofile(
                    potential="logarithmic",
                    gamma=2.0,
                    q=0.608,
                    df=1,
                    beta=0.189,
                    s=0.5,
                    t=0.0,
                    inclination=57.1,
                    xi=0.0,
                    theta=0.0,
                    integration=1,
                    ngl_or_eps=0,
                    algorithm=algorithm,
                    maxmom=cfg["maxmom"],
                    average=average,
                    kinematics=kinematics,
                    usevp=True,
                    verbose_vp=0,
                    output_path=None,  # keep file-free; res.raw_text is the reference
                    debug_prompts=False,
                    parse_stdout_fallback=False,
                    vp_reg_param=cfg.get("vp_reg_param", 1.0),
                    vp_smooth_eps=cfg.get("vp_smooth_eps", 0.0),
                )

                pth = data_dir / ref_name(
                    kinematics=kinematics,
                    average=average,
                    algorithm=algorithm,
                )
                pth.write_text(res.raw_text, encoding="utf-8")
                print(f"Wrote {pth.relative_to(root)}")

    print("Done. Commit the updated tests/data/*_ref.txt files.")
    if args.include_both:
        print("Note: 'both' references are optional; enable tests with SCALEFREE_TEST_INCLUDE_BOTH=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
