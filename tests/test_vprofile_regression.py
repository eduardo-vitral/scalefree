"""Regression tests for vprofile.

This test compares the current vprofile output against stored reference files.
Because small numerical differences can arise across compilers/architectures,
comparisons are performed at 5 *significant digits*.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from scalefree import ScaleFreeRunner


_NUMERIC_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _extract_numeric_tokens(text: str) -> np.ndarray:
    """Return all numeric tokens found in a text blob as a float array."""

    toks = _NUMERIC_RE.findall(text)
    if not toks:
        return np.array([], dtype=float)
    return np.array([float(t) for t in toks], dtype=float)


def _round_sig(x: np.ndarray, sig: int = 5) -> np.ndarray:
    """Round to `sig` significant digits (element-wise)."""

    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)

    finite = np.isfinite(x)
    out[~finite] = x[~finite]

    xf = x[finite]
    if xf.size == 0:
        return out

    # For each value, decimals = sig - 1 - floor(log10(|x|))
    mags = np.floor(np.log10(np.abs(xf)))
    mags[np.isneginf(mags)] = 0.0  # handle zeros
    decimals = (sig - 1 - mags).astype(int)

    for i, (v, d) in enumerate(zip(xf, decimals)):
        # numpy.round accepts negative decimals as well
        out[np.where(finite)[0][i]] = np.round(v, d)

    return out


def _assert_text_close(ref_text: str, cur_text: str, sig: int = 5) -> None:
    """Compare reference vs current text outputs up to `sig` significant digits."""

    ref_nums = _extract_numeric_tokens(ref_text)
    cur_nums = _extract_numeric_tokens(cur_text)

    # Basic sanity: ensure we did not change the structure drastically.
    assert ref_nums.size == cur_nums.size, (
        f"Different number of numeric tokens: ref={ref_nums.size}, cur={cur_nums.size}.\n"
        "This usually indicates a formatting/structure change in the Fortran output."
    )

    ref_r = _round_sig(ref_nums, sig=sig)
    cur_r = _round_sig(cur_nums, sig=sig)

    # Compare token-by-token; NaNs must align.
    both_nan = np.isnan(ref_r) & np.isnan(cur_r)
    neq = ~both_nan & (ref_r != cur_r)

    if np.any(neq):
        idx = np.flatnonzero(neq)[:10]
        diffs = [
            f"[{i}] ref={ref_nums[i]} -> {ref_r[i]} | cur={cur_nums[i]} -> {cur_r[i]}"
            for i in idx
        ]
        raise AssertionError(
            "Numeric mismatch after rounding to 5 significant digits. First differences:\n"
            + "\n".join(diffs)
        )


def test_vprofile_regression(scalefree_exe: Path, ref_dir: Path, tmp_path: Path) -> None:
    """Test vprofile output against stored reference files (5 significant digits)."""

    tests = {
        "projected_point": {
            "project": "projected_point",
            "algorithm": 3,
            "rmax": 50,
            "rmin": 0.01,
            "vmax": 50,
            "vmin": -50,
            "order": 4,
            "grid": "1d",
        },
        "projected_line": {
            "project": "projected_line",
            "algorithm": 3,
            "rmax": 50,
            "rmin": 0.01,
            "vmax": 50,
            "vmin": -50,
            "order": 4,
            "grid": "1d",
        },
        "projected_plane": {
            "project": "projected_plane",
            "algorithm": 3,
            "rmax": 50,
            "rmin": 0.01,
            "vmax": 50,
            "vmin": -50,
            "order": 4,
            "grid": "1d",
        },
    }

    runner = ScaleFreeRunner(scalefree_exe, workdir=tmp_path)

    for name, args in tests.items():
        ref_path = ref_dir / f"{name}_alg3_ref.txt"
        ref_text = ref_path.read_text(encoding="utf-8")

        res = runner.vprofile(**args)
        cur_text = res.raw_text

        # Helpful debugging artifact: write current output into tmp_path
        # so it can be inspected in CI artifacts if needed.
        (tmp_path / f"{name}_alg3_cur.txt").write_text(cur_text, encoding="utf-8")

        _assert_text_close(ref_text, cur_text, sig=5)
