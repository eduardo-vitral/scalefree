from __future__ import annotations

from pathlib import Path
import math

import numpy as np

from scalefree.vmoments import ScaleFreeRunner


def _extract_numeric_tokens(text: str) -> np.ndarray:
    """Extract numeric tokens from runner text output.

    Reference files and current outputs are plain text. For a stable,
    low-maintenance regression check, we compare the sequence of numeric
    tokens (float-parsable values) extracted in order, ignoring labels.
    """
    nums: list[float] = []
    for tok in text.replace(",", " ").split():
        try:
            nums.append(float(tok))
        except ValueError:
            # Non-numeric token (e.g., labels like 'PROJ', 'mean').
            continue
    return np.asarray(nums, dtype=float)


def _round_sig(x: float, sig: int = 5) -> float:
    """Round a float to a given number of significant digits."""
    if math.isnan(x) or math.isinf(x) or x == 0.0:
        return x
    # round(x, ndigits) uses decimal digits; convert from significant digits.
    ndigits = sig - int(math.floor(math.log10(abs(x)))) - 1
    return round(x, ndigits)


def _round_sig_array(arr: np.ndarray, sig: int = 5) -> np.ndarray:
    return np.vectorize(lambda v: _round_sig(float(v), sig=sig), otypes=[float])(arr)


def _assert_equal_to_5sig(ref_text: str, cur_text: str, *, context: str) -> None:
    ref = _extract_numeric_tokens(ref_text)
    cur = _extract_numeric_tokens(cur_text)

    assert ref.size == cur.size, (
        f"{context}: token count differs (ref={ref.size}, cur={cur.size}). "
        "If the output format changed intentionally, regenerate reference files."
    )

    ref_r = _round_sig_array(ref, sig=5)
    cur_r = _round_sig_array(cur, sig=5)

    # Compare with NaN/Inf semantics.
    both_nan = np.isnan(ref_r) & np.isnan(cur_r)
    both_posinf = np.isposinf(ref_r) & np.isposinf(cur_r)
    both_neginf = np.isneginf(ref_r) & np.isneginf(cur_r)
    ok_special = both_nan | both_posinf | both_neginf

    eq = (ref_r == cur_r) | ok_special
    if not np.all(eq):
        bad = np.where(~eq)[0]
        # Provide a compact debugging message with the first few mismatches.
        k = min(10, bad.size)
        idxs = bad[:k]
        details = ", ".join(
            f"i={i}: ref={ref[i]:.16g} cur={cur[i]:.16g} (ref5={ref_r[i]:.16g} cur5={cur_r[i]:.16g})"
            for i in idxs
        )
        raise AssertionError(
            f"{context}: numeric tokens differ beyond 5 significant digits. "
            f"First {k} mismatches: {details}"
        )


def test_vprofile_regression(resources: Path) -> None:
    """Test vprofile output against stored reference files."""

    root = Path(__file__).resolve().parent
    data_dir = root / "data"

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

    runner = ScaleFreeRunner(resources)

    for name, args in tests.items():
        ref_path = data_dir / f"{name}_alg3_ref.txt"
        with open(ref_path, "r", encoding="utf-8") as f:
            ref_text = f.read()

        cur_text = runner.vprofile_text(**args)

        _assert_equal_to_5sig(ref_text, cur_text, context=name)
