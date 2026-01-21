from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import scalefree


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def build_exe() -> Path:
    root = repo_root()
    exe = root / "fortran_src" / "scalefree.e"
    src = root / "fortran_src" / "scalefree.f"

    if exe.exists():
        return exe

    if not src.exists():
        raise FileNotFoundError(f"Missing Fortran source at {src}")

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
    if algorithm == 1:
        return dict(maxmom=4)
    if algorithm == 2:
        return dict(maxmom=8, vp_reg_param=1.0)
    if algorithm == 3:
        return dict(maxmom=8, vp_smooth_eps=0.0)
    raise ValueError(f"Unsupported algorithm={algorithm}")


def out_path(*, average: bool, algorithm: int) -> Path:
    stem = "out_avg" if average else "out_point"
    return repo_root() / "tests" / "data" / f"{stem}_alg{algorithm}_ref.txt"


def main() -> int:
    exe = build_exe()
    root = repo_root()
    data_dir = root / "tests" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    runner = scalefree.ScaleFreeRunner(
        exe_path=exe,
        workdir=(root / "tests" / "_work"),
    )
    runner.workdir.mkdir(parents=True, exist_ok=True)

    for algorithm in (1, 2, 3):
        cfg = case_cfg(algorithm)
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
                usevp=True,
                verbose_vp=0,
                output_path=None,  # STDOUT
                debug_prompts=False,
                parse_stdout_fallback=False,
                vp_reg_param=cfg.get("vp_reg_param", 1.0),
                vp_smooth_eps=cfg.get("vp_smooth_eps", 0.0),
            )

            p = out_path(average=average, algorithm=algorithm)
            p.write_text(res.raw_text, encoding="utf-8")
            print(f"Wrote {p.relative_to(root)}")

    print("Done. Commit the updated tests/data/*_ref.txt files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
