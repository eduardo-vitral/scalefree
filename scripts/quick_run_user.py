#!/usr/bin/env python3
from __future__ import annotations

from scalefree import ScaleFreeRunner


def potential_logarithmic():
    return 2  # Fortran: Kepler (1) or Logarithmic (2)


def print_table(title: str, columns, data):
    print(f"\n--- {title} ---")
    if columns:
        print("Columns:", " ".join(columns))
    if data is None:
        print("(no data)")
        return
    if getattr(data, "size", 0) == 0:
        print("(empty)")
        return
    for row in data:
        print(" ".join(f"{x:.16g}" for x in row))


def print_results(tag: str, res):
    print("\n==============================")
    print(f"CASE: {tag}")
    print(f"Output file: {res.output_path}")
    print(f"Blocks found: {sorted(res.blocks.keys())}")
    print("==============================")

    if "projected_point" in res.blocks:
        blk = res.blocks["projected_point"]
        print_table("projected_point", blk.get("columns", []), blk.get("data"))

    if "projected_circle_average" in res.blocks:
        blk = res.blocks["projected_circle_average"]
        print_table("projected_circle_average", blk.get("columns", []), blk.get("data"))

    if "vp" in res.blocks:
        blk = res.blocks["vp"]
        print_table("vp summary", blk.get("columns", []), blk.get("data"))
    else:
        print("\n(No 'vp' summary block found.)")

    if "vp_table" in res.blocks:
        vpt = res.blocks["vp_table"]
        for iproj in sorted(vpt.keys()):
            tbl = vpt[iproj]
            print_table(
                f"vp_table iproj={iproj}", tbl.get("columns", []), tbl.get("data")
            )
    else:
        print("\n(No 'vp_table' blocks found.)")


def run_case(runner: ScaleFreeRunner, *, average: bool, outname: str):
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
        algorithm=3,  # default algorithm
        maxmom=4,
        average=average,
        usevp=True,
        verbose_vp=0,
        output_path=outname,
        debug_prompts=False,
    )


def main():
    # Let scalefree resolve/build the backend automatically.
    # If gfortran is missing, scalefree will raise a clear error message.
    runner = ScaleFreeRunner()

    res_point = run_case(runner, average=False, outname="out_point.txt")
    print_results("average=False (point)", res_point)

    res_avg = run_case(runner, average=True, outname="out_avg.txt")
    print_results("average=True (circle-average)", res_avg)

    print("\nFinished.")


if __name__ == "__main__":
    main()
