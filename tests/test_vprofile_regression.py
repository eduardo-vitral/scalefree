from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import scalefree
from scalefree.vmoments import parse_scalefree_output


STRICT = os.getenv("SCALEFREE_STRICT_TESTS", "0").strip() == "1"


def _case_cfg(algorithm: int) -> dict:
    """
    Minimal algorithm-dependent settings.

    Notes:
      - alg=1: original behaviour (few moments)
      - alg=2/3: more moments to avoid REGMAT2 issues in some regimes
    """
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
        kinematics=kinematics,   # <-- explicit: avoids ambiguity
        usevp=(kinematics == "projected"),
        verbose_vp=0,
        output_path=None,
        debug_prompts=False,
        parse_stdout_fallback=False,
        vp_reg_param=cfg.get("vp_reg_param", 1.0),
        vp_smooth_eps=cfg.get("vp_smooth_eps", 0.0),
    )


# -----------------------------
# Smoke / contract assertions
# -----------------------------

def _assert_block_has_columns_and_data(block: dict, *, min_rows: int = 1, min_cols: int = 1):
    assert isinstance(block, dict)
    cols = block.get("columns", [])
    data = block.get("data", None)

    assert isinstance(cols, list) and len(cols) >= min_cols
    assert data is not None

    arr = np.asarray(data)
    assert arr.ndim == 2
    assert arr.shape[0] >= min_rows
    assert arr.shape[1] >= min_cols


def _assert_projected_contract(blocks: dict, *, average: bool):
    kind = "projected_circle_average" if average else "projected_point"
    assert kind in blocks, f"Missing '{kind}' block"

    blk = blocks[kind]
    _assert_block_has_columns_and_data(blk, min_rows=3, min_cols=6)

    # Columns in Fortran structured output are expected to be:
    # iproj rho_p v1 v2 v3 v4
    expected = ["iproj", "rho_p", "v1", "v2", "v3", "v4"]
    assert blk["columns"] == expected

    arr = np.asarray(blk["data"], dtype=float)

    # iproj should be 1..3
    iproj = arr[:, 0].astype(int)
    assert set(iproj.tolist()) == {1, 2, 3}

    # rho_p should be positive
    rho_p = arr[:, 1]
    assert np.all(rho_p > 0)

    # v2 (second moment) should be >= 0
    v2 = arr[:, 3]
    assert np.all(v2 >= 0)


def _assert_vp_contract(blocks: dict):
    # VP blocks only exist when projected mode & usevp=True
    assert "vp" in blocks, "Missing 'vp' block (backend may have stopped mid-run)"
    _assert_block_has_columns_and_data(blocks["vp"], min_rows=3, min_cols=7)

    # vp_table should exist and have iproj 1..3 tables
    assert "vp_table" in blocks, "Missing 'vp_table' block"
    vpt = blocks["vp_table"]
    assert isinstance(vpt, dict) and vpt, "vp_table is empty"

    for ip in (1, 2, 3):
        assert ip in vpt, f"vp_table missing iproj={ip}"
        tbl = vpt[ip]
        _assert_block_has_columns_and_data(tbl, min_rows=5, min_cols=2)
        assert tbl["columns"] == ["v", "vp"]


def _assert_intrinsic_contract(blocks: dict, *, average: bool):
    kind = "intrinsic_shell_average" if average else "intrinsic_point"
    assert kind in blocks, f"Missing '{kind}' block"

    blk = blocks[kind]
    if average:
        assert blk["columns"] == ["rho", "vphi", "vr2", "vth2", "vphi2", "beta"]
        _assert_block_has_columns_and_data(blk, min_rows=1, min_cols=6)
    else:
        assert blk["columns"] == ["rho", "vphi", "vr2", "vth2", "vphi2"]
        _assert_block_has_columns_and_data(blk, min_rows=1, min_cols=5)

    arr = np.asarray(blk["data"], dtype=float)
    assert np.all(np.isfinite(arr[:, :5])), "Intrinsic moments contain NaN/Inf unexpectedly"
    assert np.all(arr[:, 0] > 0), "Intrinsic rho must be positive"


# -----------------------------
# Optional strict regression
# -----------------------------

def _ref_path(ref_dir: Path, *, average: bool, algorithm: int) -> Path:
    stem = "out_avg" if average else "out_point"
    return ref_dir / f"{stem}_alg{algorithm}_ref.txt"


def _skip_if_missing_ref(ref_path: Path):
    if not ref_path.exists():
        pytest.skip(
            f"Missing reference file: {ref_path}\n"
            "Generate/refresh refs by running:\n"
            "  python tests/make_vprofile_refs.py\n"
            "and commit tests/data/*_ref.txt."
        )


def _assert_block_close(new_blk, ref_blk, *, rtol=5e-5, atol=1e-7):
    """
    Relaxed tolerances to accommodate platform/compiler differences.
    """
    assert new_blk.get("columns", []) == ref_blk.get("columns", [])

    new_data = np.asarray(new_blk.get("data"), dtype=float)
    ref_data = np.asarray(ref_blk.get("data"), dtype=float)

    assert new_data.shape == ref_data.shape

    # Handle extremely tiny underflows and denormals consistently.
    tiny = 1e-300
    new_data = np.where(np.abs(new_data) < tiny, 0.0, new_data)
    ref_data = np.where(np.abs(ref_data) < tiny, 0.0, ref_data)

    assert np.allclose(new_data, ref_data, rtol=rtol, atol=atol, equal_nan=True)


def _compare_outputs(new_blocks, ref_blocks):
    # Only compare projected outputs (refs are projected)
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


# -----------------------------
# Tests
# -----------------------------

@pytest.mark.parametrize("algorithm", [1, 2, 3])
def test_projected_point_smoke(scalefree_exe, tmp_path, algorithm, ref_dir):
    res = _run_case(
        scalefree_exe,
        average=False,
        workdir=tmp_path,
        algorithm=algorithm,
        kinematics="projected",
    )
    blocks = res.blocks

    _assert_projected_contract(blocks, average=False)
    _assert_vp_contract(blocks)

    if STRICT:
        ref_path = _ref_path(ref_dir, average=False, algorithm=algorithm)
        _skip_if_missing_ref(ref_path)
        ref_blocks = parse_scalefree_output(ref_path.read_text(encoding="utf-8", errors="replace"))
        _compare_outputs(blocks, ref_blocks)


@pytest.mark.parametrize("algorithm", [1, 2, 3])
def test_projected_average_smoke(scalefree_exe, tmp_path, algorithm, ref_dir):
    res = _run_case(
        scalefree_exe,
        average=True,
        workdir=tmp_path,
        algorithm=algorithm,
        kinematics="projected",
    )
    blocks = res.blocks

    _assert_projected_contract(blocks, average=True)
    _assert_vp_contract(blocks)

    if STRICT:
        ref_path = _ref_path(ref_dir, average=True, algorithm=algorithm)
        _skip_if_missing_ref(ref_path)
        ref_blocks = parse_scalefree_output(ref_path.read_text(encoding="utf-8", errors="replace"))
        _compare_outputs(blocks, ref_blocks)


@pytest.mark.parametrize("average", [False, True])
def test_intrinsic_smoke(scalefree_exe, tmp_path, average):
    res = _run_case(
        scalefree_exe,
        average=average,
        workdir=tmp_path,
        algorithm=3,            # algorithm irrelevant for intrinsic moments; choose a stable default
        kinematics="intrinsic",
    )
    blocks = res.blocks
    _assert_intrinsic_contract(blocks, average=average)

    # Ensure projected-only blocks are not leaking into intrinsic-only mode
    assert not any(k.startswith("projected") for k in blocks.keys())
    assert "vp" not in blocks
    assert "vp_table" not in blocks
