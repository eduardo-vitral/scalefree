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
        print(" ".join(f"{x:.16g}" for x in row))


def print_results(tag: str, res):
    print("\n==============================")
    print(f"CASE: {tag}")
    print(f"Blocks found: {sorted(res.blocks.keys())}")
    print("==============================")

    for k in (
        "intrinsic_point",
        "intrinsic_shell_average",
        "projected_point",
        "projected_circle_average",
        "vp",
    ):
        if k in res.blocks:
            blk = res.blocks[k]
            print_table(k, blk.get("columns", []), blk.get("data"))

    if "vp_table" in res.blocks:
        vpt = res.blocks["vp_table"]
        for iproj in sorted(vpt.keys()):
            tbl = vpt[iproj]
            print_table(
                f"vp_table iproj={iproj}",
                tbl.get(
                    "columns",
                    [],
                ),
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
        maxmom=10,
        average=average,
        kinematics=kinematics,
        usevp=usevp,
        verbose_vp=verbose,
        output_path=None,  # no persistent file
        debug_prompts=debug,  # prints Fortran conversation
    )


def main():
    p = argparse.ArgumentParser(
        description="Local smoke tests for scalefree.vmoments kinematics modes.",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Print full Fortran prompt/response flow.",
    )
    args, _unknown = p.parse_known_args()  
    # ignore ipykernel args if run in notebooks

    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner  # noqa: E402

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    runner = ScaleFreeRunner()

    # 1) Intrinsic point (iwhat=0 when average=False)
    res = run_case(
        runner,
        average=False,
        kinematics="intrinsic",
        usevp=False,
        verbose=0,
        debug=args.debug,
    )
    print_results("intrinsic / point (average=False)", res)

    # 2) Intrinsic mass-weighted shell average (iwhat=2 when average=True)
    res = run_case(
        runner,
        average=True,
        kinematics="intrinsic",
        usevp=False,
        verbose=0,
        debug=args.debug,
    )
    print_results("intrinsic / shell-average (average=True)", res)

    # 3) Projected point (iwhat=1 when average=False)
    res = run_case(
        runner,
        average=False,
        kinematics="projected",
        usevp=True,
        verbose=0,
        debug=args.debug,
    )
    print_results("projected / point (average=False)", res)

    # 4) Projected circle-average (iwhat=3 when average=True)
    res = run_case(
        runner,
        average=True,
        kinematics="projected",
        usevp=True,
        verbose=0,
        debug=args.debug,
    )
    print_results("projected / circle-average (average=True)", res)

    # 5) Legacy: both (intrinsic + projected)
    res = run_case(
        runner,
        average=False,
        kinematics="both",
        usevp=True,
        verbose=0,
        debug=args.debug,
    )
    print_results("both (intrinsic + projected), average=False", res)

    print("\nFinished.")


if __name__ == "__main__":
    main()
