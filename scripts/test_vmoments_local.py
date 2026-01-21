#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _force_repo_import(repo_root: Path) -> None:
    """
    Force Python to import 'scalefree' from the working tree (repo checkout),
    not from a pip/conda installed distribution.
    """
    repo_root = repo_root.resolve()
    sys.path.insert(0, str(repo_root))

    # If scalefree was already imported (very common in notebooks), drop it.
    for k in list(sys.modules.keys()):
        if k == "scalefree" or k.startswith("scalefree."):
            del sys.modules[k]


def potential_logarithmic() -> int:
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

    if "_fortran" in res.blocks:
        meta = res.blocks["_fortran"]
        print("\n--- fortran meta ---")
        print("parsed_from:", meta.get("parsed_from"))
        print("exe_path:", meta.get("exe_path"))
        # Print full captured stdout/stderr only if debug was enabled in vprofile
        # (vprofile will already print them when debug_prompts=True)
        if meta.get("stderr", "").strip():
            print("\n--- fortran stderr (captured; non-empty) ---")
            print(meta["stderr"])

    if "projected_point" in res.blocks:
        blk = res.blocks["projected_point"]
        print_table("projected_point", blk.get("columns", []), blk.get("data"))

    if "projected_circle_average" in res.blocks:
        blk = res.blocks["projected_circle_average"]
        print_table("projected_circle_average", blk.get("columns", []), blk.get("data"))

    if "vp" in res.blocks:
        blk = res.blocks["vp"]
        print_table("vp summary", blk.get("columns", []), blk.get("data"))

        # Hard check: header/data width match
        data = blk.get("data")
        cols = blk.get("columns", [])
        if getattr(data, "ndim", 0) == 2 and data.size:
            assert data.shape[1] == len(cols), (
                f"VP header/data mismatch: {len(cols)} columns vs {data.shape[1]} data width.\n"
                f"columns={cols}"
            )
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


def run_case(runner, *, average: bool, debug: bool):
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
        ngl_or_eps=0,  # 0 => Fortran default (and avoids "eps" prompt paths)
        algorithm=3,  # recommend 3 for physical VPs; change to 1 if you explicitly want it
        maxmom=10,
        average=average,
        usevp=True,
        verbose_vp=1,
        output_path=None,
        debug_prompts=True,
    )


def main():
    p = argparse.ArgumentParser(
        description="Local test of working-tree scalefree.vmoments parsing."
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug_prompts in vmoments and print full captured output.",
    )
    args, _unknown = p.parse_known_args()  # IMPORTANT: ignore Jupyter/ipykernel flags

    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner  # noqa: E402
    from scalefree import mock

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    runner = ScaleFreeRunner()

    res_point = run_case(runner, average=False, debug=False)
    print_results("average=False (point)", res_point)

    # xyzv = mock(
    #     potential=2,
    #     gamma=2.0,
    #     q=0.608,
    #     df=1,
    #     beta=0.189,
    #     s=0.5,
    #     t=0.0,
    #     nsamples=10,
    #     seed=101,  # use non-negative (or keep -101 if you applied the seed fix)
    #     rin=1.0,
    #     rout=10.0,
    #     inclination=57.1,
    #     xi=0.0,
    #     algorithm=1,
    #     nbins=6,  # key: few bins => few vprofile calls
    #     nsig=6,  # smaller PDF grid support is fine for a smoke check
    #     grid_n=801,  # smaller grid speeds up BALRoGO sampling
    #     debug=True,  # prints bin occupancy summary
    # )
    # print(xyzv.shape)
    # print(xyzv[:3])

    print("\nFinished.")


if __name__ == "__main__":
    main()
