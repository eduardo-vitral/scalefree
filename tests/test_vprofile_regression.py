from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

import scalefree
from scalefree.vmoments import parse_scalefree_output


def _case_cfg(algorithm: int) -> dict:
    """
    Minimal algorithm-dependent settings.

    - alg=1: keep original behaviour (maxmom=4)
    - alg=2/3: use higher maxmom to avoid REGMAT2
    "Matrix too small" for some models
    """
    if algorithm == 1:
        return dict(maxmom=4)
    if algorithm == 2:
        return dict(maxmom=8, vp_reg_param=1.0)
    if algorithm == 3:
        return dict(maxmom=8, vp_smooth_eps=0.0)
    raise ValueError(f"Unsupported algorithm={algorithm}")


def _run_case(exe_path: Path, *, average: bool, workdir: Path, algorithm: int):
    runner = scalefree.ScaleFreeRunner(exe_path=exe_path, workdir=workdir)

    cfg = _case_cfg(algorithm)

    # fixed test case (your quick_run numbers)
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
        usevp=True,
        verbose_vp=0,  # keep output small/stable
        output_path=None,  # file-free strategy (STDOUT)
        debug_prompts=False,
        parse_stdout_fallback=False,
        # Only used when algorithm != 1; harmless otherwise
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

    if "vp_table" in ref_blocks:
        assert "vp_table" in new_blocks
        for iproj, ref_tbl in ref_blocks["vp_table"].items():
            assert iproj in new_blocks["vp_table"]
            _assert_block_close(new_blocks["vp_table"][iproj], ref_tbl)


def _ref_path(ref_dir: Path, *, average: bool, algorithm: int) -> Path:
    stem = "out_avg" if average else "out_point"
    return ref_dir / f"{stem}_alg{algorithm}_ref.txt"


def _skip_if_missing_ref(ref_path: Path):
    if not ref_path.exists():
        pytest.skip(
            f"Missing reference file: {ref_path}\n"
            "Generate/refresh refs by running:\n"
            "  python tests/make_vprofile_refs.py\n"
            "and commit the resulting tests/data/*_ref.txt files."
        )


@pytest.mark.parametrize("algorithm", [1, 2, 3])
def test_vprofile_point_regression(
    scalefree_exe,
    ref_dir,
    tmp_path,
    algorithm,
):
    res = _run_case(
        scalefree_exe,
        average=False,
        workdir=tmp_path,
        algorithm=algorithm,
    )
    new_blocks = res.blocks

    ref_path = _ref_path(ref_dir, average=False, algorithm=algorithm)
    _skip_if_missing_ref(ref_path)

    ref_blocks = parse_scalefree_output(
        ref_path.read_text(encoding="utf-8", errors="replace")
    )

    _compare_outputs(new_blocks, ref_blocks)

    # Guardrail: VP block should not be empty when present in refs
    if "vp" in ref_blocks:
        assert (
            new_blocks["vp"]["data"].size > 0
        ), "VP block is empty (backend likely STOPped mid-run)"


@pytest.mark.parametrize("algorithm", [1, 2, 3])
def test_vprofile_average_regression(
    scalefree_exe,
    ref_dir,
    tmp_path,
    algorithm,
):
    res = _run_case(
        scalefree_exe,
        average=True,
        workdir=tmp_path,
        algorithm=algorithm,
    )
    new_blocks = res.blocks

    ref_path = _ref_path(ref_dir, average=True, algorithm=algorithm)
    _skip_if_missing_ref(ref_path)

    ref_blocks = parse_scalefree_output(
        ref_path.read_text(encoding="utf-8", errors="replace")
    )

    _compare_outputs(new_blocks, ref_blocks)

    if "vp" in ref_blocks:
        assert (
            new_blocks["vp"]["data"].size > 0
        ), "VP block is empty (backend likely STOPped mid-run)"
