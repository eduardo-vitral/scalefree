#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _force_repo_import(repo_root: Path) -> None:
    """
    Force Python to import 'scalefree' from the working tree (repo checkout),
    not from a pip/conda installed distribution.
    """
    repo_root = repo_root.resolve()
    sys.path.insert(0, str(repo_root))

    # If scalefree was already imported (common in notebooks), drop it.
    for k in list(sys.modules.keys()):
        if k == "scalefree" or k.startswith("scalefree."):
            del sys.modules[k]


def potential_logarithmic() -> int:
    return 2  # Fortran: Kepler (1) or Logarithmic (2)


def print_table(title: str, columns, data):
    print(f"\n--- {title} ---")
    if columns:
        print("Columns:", " ".join(columns))
    if data is None or getattr(data, "size", 0) == 0:
        print("(no data)")
        return
    for row in data:
        try:
            print(" ".join(f"{x:.16g}" for x in row))
        except TypeError:
            # fallback if row has non-numerics
            print(" ".join(str(x) for x in row))


def _print_vp_summary_block(block_name: str, blk: dict):
    # Expected: columns + data; may also have by_iproj/by_icomp
    cols = blk.get("columns", [])
    data = blk.get("data")
    print_table(block_name, cols, data)

    # Convenience: print a short per-component summary if present
    by_key = None
    if "by_iproj" in blk:
        by_key = "by_iproj"
    if "by_icomp" in blk:
        by_key = "by_icomp"

    if by_key:
        by = blk.get(by_key, {})
        if by:
            print(f"\n{block_name}: parsed index -> row mapping ({by_key})")
            for kk in sorted(by.keys()):
                row = by[kk]

                def _fmt(x):
                    # Accept numbers or numeric strings; fall back to str
                    try:
                        xf = float(x)
                        return f"{xf:.16g}"
                    except Exception:
                        return str(x)

                print(f"  {kk}: " + " ".join(_fmt(x) for x in row))


def print_results(tag: str, res):
    print("\n==============================")
    print(f"CASE: {tag}")
    print(f"Blocks found: {sorted(res.blocks.keys())}")
    print("==============================")

    # Core kinematics blocks
    for k in (
        "intrinsic_point",
        "intrinsic_shell_average",
        "projected_point",
        "projected_circle_average",
    ):
        if k in res.blocks:
            blk = res.blocks[k]
            print_table(k, blk.get("columns", []), blk.get("data"))

    # VP summary blocks (projected + intrinsic)
    for k in ("vp", "vp_intrinsic"):
        if k in res.blocks:
            _print_vp_summary_block(k, res.blocks[k])

    # VP tables (projected)
    if "vp_table" in res.blocks:
        vpt = res.blocks["vp_table"]
        for iproj in sorted(vpt.keys()):
            tbl = vpt[iproj]
            print_table(
                f"vp_table iproj={iproj}",
                tbl.get("columns", []),
                tbl.get("data"),
            )

    # VP tables (intrinsic)
    if "vp_table_intrinsic" in res.blocks:
        vpt = res.blocks["vp_table_intrinsic"]
        for icomp in sorted(vpt.keys()):
            tbl = vpt[icomp]
            print_table(
                f"vp_table_intrinsic icomp={icomp}",
                tbl.get("columns", []),
                tbl.get("data"),
            )


def run_case(
    runner,
    *,
    average: bool,
    kinematics,
    usevp: bool,
    verbose: int,
    debug: bool,
    maxmom: int,
):
    return runner.vprofile(
        potential=potential_logarithmic,
        gamma=2.0,
        q=0.608,
        df=1,
        beta=0.189,
        s=0.5,
        t=0.0,
        inclination=57.1,
        xi=0.0,
        theta=0.0,
        integration=1,  # Gauss-Legendre
        ngl_or_eps=0,  # 0 => Fortran default
        algorithm=3,
        maxmom=maxmom,
        average=average,
        kinematics=kinematics,
        usevp=usevp,
        verbose_vp=verbose,
        output_path=None,  # no persistent file
        debug_prompts=debug,  # prints Fortran conversation
    )


def main():
    p = argparse.ArgumentParser(
        description=(
            "Local smoke tests for scalefree.vmoments kinematics modes,"
            " including VP diagnostics.\n\n"
            "This script exercises: intrinsic/projected/both × point/average,"
            " and prints VP blocks if enabled."
        ),
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Print full Fortran prompt/response flow.",
    )
    p.add_argument(
        "--no-vp",
        action="store_true",
        help="Disable VP diagnostics (usevp=False for all cases).",
    )
    p.add_argument(
        "--maxmom",
        type=int,
        default=10,
        help="Maximum moment order passed to Fortran (default: 10).",
    )
    p.add_argument(
        "--verbose-vp",
        type=int,
        default=0,
        help="Fortran VP verbosity flag (default: 0).",
    )

    p.add_argument(
        "--exe",
        type=str,
        default=None,
        help=(
            "Path to Fortran backend executable. "
            "If it does not exist, scalefree will compile it there. "
            "If omitted, this script forces a rebuild into "
            "<repo>/fortran_src/scalefree_intrvp.e."
        ),
    )
    args, _unknown = p.parse_known_args()
    # ignore ipykernel args if run in notebooks

    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner  # noqa: E402

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    # Force use of a fresh backend build so Fortran changes
    # (e.g., vp_intrinsic)
    # are definitely picked up even if a cached executable exists.
    if args.exe is None:
        exe_path = repo_root / "fortran_src" / "scalefree_intrvp.e"
    else:
        exe_path = Path(args.exe).expanduser().resolve()

    runner = ScaleFreeRunner(exe_path=exe_path)

    usevp = not args.no_vp

    # 1) Intrinsic point (iwhat=0 when average=False)
    res = run_case(
        runner,
        average=False,
        kinematics="intrinsic",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("intrinsic / point (average=False)", res)

    # 2) Intrinsic mass-weighted shell average (iwhat=2 when average=True)
    res = run_case(
        runner,
        average=True,
        kinematics="intrinsic",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("intrinsic / shell-average (average=True)", res)

    # 3) Projected point (iwhat=1 when average=False)
    res = run_case(
        runner,
        average=False,
        kinematics="projected",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("projected / point (average=False)", res)

    # 4) Projected circle-average (iwhat=3 when average=True)
    res = run_case(
        runner,
        average=True,
        kinematics="projected",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("projected / circle-average (average=True)", res)

    # 5) Both (intrinsic + projected) — point mode
    res = run_case(
        runner,
        average=False,
        kinematics="both",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("both (intrinsic + projected), average=False", res)

    # 6) Both (intrinsic + projected) — average mode
    res = run_case(
        runner,
        average=True,
        kinematics="both",
        usevp=usevp,
        verbose=args.verbose_vp,
        debug=args.debug,
        maxmom=args.maxmom,
    )
    print_results("both (intrinsic + projected), average=True", res)

    print("\nFinished.")


if __name__ == "__main__":
    main()
