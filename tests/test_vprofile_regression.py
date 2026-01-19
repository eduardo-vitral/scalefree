from __future__ import annotations

from pathlib import Path
import numpy as np

import scalefree
from scalefree.vmoments import parse_scalefree_output


def _run_case(exe_path: Path, *, average: bool, outname: str, workdir: Path):
    runner = scalefree.ScaleFreeRunner(exe_path=exe_path, workdir=workdir)

    # fixed test case (your quick_run numbers)
    return runner.vprofile(
        potential="logarithmic",  # maps to 2 in vmoments.py
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
        algorithm=1,
        maxmom=4,
        average=average,
        usevp=True,
        verbose_vp=0,
        output_path=outname,
        debug_prompts=False,
        parse_stdout_fallback=False,
    )


def _assert_block_close(new_blk, ref_blk, *, rtol=1e-8, atol=1e-10):
    # Compare columns exactly
    assert new_blk.get("columns", []) == ref_blk.get("columns", [])

    new_data = np.asarray(new_blk.get("data"))
    ref_data = np.asarray(ref_blk.get("data"))

    assert new_data is not None and ref_data is not None
    assert new_data.shape == ref_data.shape

    # Treat subnormals/denormals as zero (compiler/CPU dependent)
    tiny = 1e-300
    new_data = np.where(np.abs(new_data) < tiny, 0.0, new_data)
    ref_data = np.where(np.abs(ref_data) < tiny, 0.0, ref_data)

    # Use slightly looser tolerance (still very strict for scientific code)
    assert np.allclose(
        new_data,
        ref_data,
        rtol=rtol,
        atol=atol,
        equal_nan=True,
    )


def _compare_outputs(new_blocks, ref_blocks):
    # Only enforce a minimal set of blocks that should exist for regression.
    # Adjust these keys if your Fortran output names differ.
    must_have = []
    if "projected_point" in ref_blocks:
        must_have.append("projected_point")
    if "projected_circle_average" in ref_blocks:
        must_have.append("projected_circle_average")
    if "vp" in ref_blocks:
        must_have.append("vp")

    for k in must_have:
        assert k in new_blocks, f"Missing block '{k}' in new output"
        _assert_block_close(new_blocks[k], ref_blocks[k])

    # VP tables: compare only if present in reference
    if "vp_table" in ref_blocks:
        assert "vp_table" in new_blocks
        for iproj, ref_tbl in ref_blocks["vp_table"].items():
            assert iproj in new_blocks["vp_table"]
            _assert_block_close(new_blocks["vp_table"][iproj], ref_tbl)


def test_vprofile_point_regression(scalefree_exe, ref_dir, tmp_path):
    # Run new output
    res = _run_case(
        scalefree_exe,
        average=False,
        outname="out_point_new.txt",
        workdir=tmp_path,
    )
    new_blocks = res.blocks

    # Load reference output text and parse
    ref_path = ref_dir / "out_point_ref.txt"
    assert ref_path.exists(), f"Missing reference file: {ref_path}"
    ref_blocks = parse_scalefree_output(
        ref_path.read_text(encoding="utf-8", errors="replace")
    )

    _compare_outputs(new_blocks, ref_blocks)


def test_vprofile_average_regression(scalefree_exe, ref_dir, tmp_path):
    # Run new output
    res = _run_case(
        scalefree_exe,
        average=True,
        outname="out_avg_new.txt",
        workdir=tmp_path,
    )
    new_blocks = res.blocks

    # Load reference output text and parse
    ref_path = ref_dir / "out_avg_ref.txt"
    assert ref_path.exists(), f"Missing reference file: {ref_path}"
    ref_blocks = parse_scalefree_output(
        ref_path.read_text(encoding="utf-8", errors="replace")
    )

    _compare_outputs(new_blocks, ref_blocks)
