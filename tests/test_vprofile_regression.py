from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import pytest

import scalefree
from scalefree.vmoments import parse_scalefree_output


def _case_cfg(algorithm: int) -> dict:
    """Algorithm-dependent settings aligned with make_vprofile_refs.py."""
    if algorithm == 1:
        return dict(maxmom=4)
    if algorithm == 2:
        return dict(maxmom=8, vp_reg_param=1.0)
    if algorithm == 3:
        return dict(maxmom=8, vp_smooth_eps=0.0)
    raise ValueError(f"Unsupported algorithm={algorithm}")


def _run_case(exe_path: Path, *, average: bool, workdir: Path, algorithm: int, kinematics: str):
    runner = scalefree.ScaleFreeRunner(exe_path=exe_path, workdir=workdir)
    cfg = _case_cfg(algorithm)

    return runner.vprofile(
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
        output_path=None,
        debug_prompts=False,
        parse_stdout_fallback=False,
        vp_reg_param=cfg.get("vp_reg_param", 1.0),
        vp_smooth_eps=cfg.get("vp_smooth_eps", 0.0),
    )


def _assert_block_close(new_blk, ref_blk, *, rtol=1e-6, atol=5e-8):
    assert new_blk.get("columns", []) == ref_blk.get("columns", [])

    new_data = np.asarray(new_blk.get("data"))
    ref_data = np.asarray(ref_blk.get("data"))

    assert new_data is not None and ref_data is not None
    assert new_data.shape == ref_data.shape

    tiny = 1e-300
    new_data = np.where(np.abs(new_data) < tiny, 0.0, new_data)
    ref_data = np.where(np.abs(ref_data) < tiny, 0.0, ref_data)

    assert np.allclose(
        new_data,
        ref_data,
        rtol=rtol,
        atol=atol,
        equal_nan=True,
    )


def _compare_outputs(new_blocks, ref_blocks):
    # Compare every block that exists in refs.
    for k, ref_blk in ref_blocks.items():
        if k in ("vp_table", "vp_table_intrinsic"):
            continue

        assert k in new_blocks, f"Missing block '{k}' in new output"
        _assert_block_close(new_blocks[k], ref_blk)

    # Compare nested VP tables (projected)
    if "vp_table" in ref_blocks:
        assert "vp_table" in new_blocks
        for iproj, ref_tbl in ref_blocks["vp_table"].items():
            assert iproj in new_blocks["vp_table"]
            _assert_block_close(new_blocks["vp_table"][iproj], ref_tbl)

    # Compare nested VP tables (intrinsic)
    if "vp_table_intrinsic" in ref_blocks:
        assert "vp_table_intrinsic" in new_blocks
        for icomp, ref_tbl in ref_blocks["vp_table_intrinsic"].items():
            assert icomp in new_blocks["vp_table_intrinsic"]
            _assert_block_close(new_blocks["vp_table_intrinsic"][icomp], ref_tbl)


def _ref_path(ref_dir: Path, *, average: bool, algorithm: int, kinematics: str) -> Path:
    stem = "avg" if average else "point"
    return ref_dir / f"{kinematics}_{stem}_alg{algorithm}_ref.txt"


def _require_ref(ref_path: Path):
    if ref_path.exists():
        return

    msg = (
        f"Missing reference file: {ref_path}\n"
        "Generate/refresh refs by running:\n"
        "  python tests/make_vprofile_refs.py\n"
        "and commit the resulting tests/data/*_ref.txt files."
    )

    # Fail on CI, skip locally.
    if os.environ.get("CI"):
        raise AssertionError(msg)
    pytest.skip(msg)


@pytest.mark.parametrize("algorithm", [1, 2, 3])
@pytest.mark.parametrize("average", [False, True])
@pytest.mark.parametrize("kinematics", ["intrinsic", "projected"])
def test_vprofile_regression(
    scalefree_exe,
    ref_dir,
    tmp_path,
    algorithm,
    average,
    kinematics,
):
    res = _run_case(
        scalefree_exe,
        average=average,
        workdir=tmp_path,
        algorithm=algorithm,
        kinematics=kinematics,
    )
    new_blocks = res.blocks

    ref_path = _ref_path(ref_dir, average=average, algorithm=algorithm, kinematics=kinematics)
    _require_ref(ref_path)

    ref_blocks = parse_scalefree_output(ref_path.read_text(encoding="utf-8", errors="replace"))
    _compare_outputs(new_blocks, ref_blocks)

    # Guardrail: VP blocks should not be empty when present
    for vp_key in ("vp", "vp_intrinsic"):
        if vp_key in ref_blocks:
            assert new_blocks[vp_key]["data"].size > 0, f"{vp_key} block is empty"


@pytest.mark.parametrize("algorithm", [1, 2, 3])
@pytest.mark.parametrize("average", [False, True])
def test_vprofile_regression_both_optional(
    scalefree_exe,
    ref_dir,
    tmp_path,
    algorithm,
    average,
    include_both,
):
    if not include_both:
        pytest.skip("Optional: set SCALEFREE_TEST_INCLUDE_BOTH=1 to enable kinematics='both' regression")

    res = _run_case(
        scalefree_exe,
        average=average,
        workdir=tmp_path,
        algorithm=algorithm,
        kinematics="both",
    )
    new_blocks = res.blocks

    ref_path = _ref_path(ref_dir, average=average, algorithm=algorithm, kinematics="both")
    _require_ref(ref_path)

    ref_blocks = parse_scalefree_output(ref_path.read_text(encoding="utf-8", errors="replace"))
    _compare_outputs(new_blocks, ref_blocks)
