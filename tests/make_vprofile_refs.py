from __future__ import annotations

import argparse
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np

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


def _extract_numeric_tokens(text: str) -> np.ndarray:
    nums: list[float] = []
    for tok in text.replace(",", " ").split():
        try:
            nums.append(float(tok))
        except ValueError:
            continue
    return np.asarray(nums, dtype=float)


def _round_sig(x: float, sig: int = 5) -> float:
    if math.isnan(x) or math.isinf(x) or x == 0.0:
        return x
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)


def _compare_texts(ref_text: str, new_text: str, sig: int = 5) -> tuple[bool, str]:
    a = _extract_numeric_tokens(ref_text)
    b = _extract_numeric_tokens(new_text)

    if a.shape != b.shape:
        return False, f"shape mismatch: ref has {a.size} nums, new has {b.size} nums"

    for i, (ra, rb) in enumerate(zip(a, b)):
        if math.isnan(ra) and math.isnan(rb):
            continue
        if math.isinf(ra) and math.isinf(rb) and (ra > 0) == (rb > 0):
            continue
        if _round_sig(float(ra), sig) != _round_sig(float(rb), sig):
            return False, f"diff at index {i}: ref={ra} new={rb} (sig={sig})"

    return True, "match"


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Generate/update golden reference outputs under tests/data/.\n\n"
            "Default behaviour is to compare the newly-generated outputs against the\n"
            "existing *_ref.txt files up to 5 significant digits. Use --overwrite\n"
            "to replace the reference files."
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
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing *_ref.txt files with newly-generated outputs.",
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

    any_mismatch = False

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
                    # IMPORTANT: provide a short output filename.
                    # Some Fortran builds use fixed-length CHARACTER buffers for
                    # filenames; long autogenerated names may be truncated,
                    # causing the Python side to look for the wrong path.
                    output_path=runner.workdir / "vp.txt",
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

                if pth.exists() and not args.overwrite:
                    ok, msg = _compare_texts(pth.read_text(encoding="utf-8"), res.raw_text, sig=5)
                    if ok:
                        print(f"OK   {pth.relative_to(root)}")
                    else:
                        any_mismatch = True
                        cand = pth.with_name(pth.stem + "_candidate" + pth.suffix)
                        cand.write_text(res.raw_text, encoding="utf-8")
                        print(f"DIFF {pth.relative_to(root)} -> wrote {cand.relative_to(root)} ({msg})")
                else:
                    pth.write_text(res.raw_text, encoding="utf-8")
                    print(f"Wrote {pth.relative_to(root)}")

    if args.overwrite:
        print("Done. Commit the updated tests/data/*_ref.txt files.")
    else:
        if any_mismatch:
            print("\nOne or more references differed at 5 significant digits.")
            print("Review the *_candidate.txt files; if changes are intended, re-run with --overwrite and commit.")
            return 2
        print("\nAll references match at 5 significant digits.")

    if args.include_both:
        print("Note: 'both' references are optional; enable tests with SCALEFREE_TEST_INCLUDE_BOTH=1")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
