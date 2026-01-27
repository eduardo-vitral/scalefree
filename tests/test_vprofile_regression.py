from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pytest

from scalefree.vmoments import ScaleFreeRunner


# -----------------------------------------------------------------------------
# Numerical comparison helpers (cross-platform / compiler tolerant)
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class SigSpec:
    """Specification for ~N significant-digit comparisons."""

    sig: int = 5
    tiny: float = 1e-300
    factor: float = 2.0


_NUM_RE = re.compile(
    r"""^
    [+-]?
    (?:
        (?:\d+\.?\d*)
        |
        (?:\d*\.\d+)
    )
    (?:
        [eEdD][+-]?\d+
        |
        [+-]\d+   # Fortran sometimes prints like 0.1416-319 (no 'E')
    )?
    $""",
    re.VERBOSE,
)


def _parse_fortran_float(tok: str) -> Optional[float]:
    """Parse a float token robustly across Fortran formatting variants."""
    t = tok.strip()
    if not t or not _NUM_RE.match(t):
        return None

    # Handle Fortran 'D' exponent
    t = t.replace("D", "E").replace("d", "e")

    # Handle missing 'E': e.g. 0.1416-319 or 1.23+05
    m = re.match(r"^([+-]?(?:\d+(?:\.\d*)?|\d*\.\d+))([+-]\d+)$", t)
    if m and ("e" not in t.lower()):
        t = f"{m.group(1)}e{m.group(2)}"

    try:
        return float(t)
    except ValueError:
        return None


# -----------------------------------------------------------------------------
# Tokenization helpers
# -----------------------------------------------------------------------------
# NOTE:
# The Fortran backend prints some numeric tables in fixed-width formats like:
#   (I3,1X,6E24.16,1X,7F10.5)
# Those formats do NOT guarantee whitespace between adjacent fields. On some
# compilers/platforms, two consecutive numbers can become "glued" (e.g.
# '0.1234567890123456E+00-0.1111111111111111E+00'), which breaks naive
# `str.split()` tokenization. The regex below extracts float-like substrings
# even when fields are adjacent.

_FLOAT_FIND_RE = re.compile(
    r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eEdD][+-]?\d+|[+-]\d+)?"
)


def _tokenize_relaxed(line: str) -> list[str]:
    """Tokenize a vprofile numeric row robustly.

    Keeps the leading integer index (if present) and then extracts all
    subsequent numeric fields using a regex, tolerating adjacent fields.
    """
    s = line.rstrip()
    if not s:
        return []

    # Preserve leading integer row id when present
    m = re.match(r"^\s*([+-]?\d+)\s+(.*)$", s)
    if m:
        head = [m.group(1)]
        rest = m.group(2)
        nums = _FLOAT_FIND_RE.findall(rest)
        return head + nums

    return _FLOAT_FIND_RE.findall(s)


def _sig_tol(x: float, sig: int, factor: float) -> float:
    """Absolute tolerance corresponding to ~sig significant digits."""
    ax = abs(x)
    if ax == 0.0 or not math.isfinite(ax):
        return 0.0
    exp10 = math.floor(math.log10(ax))
    ulp = 10 ** (exp10 - sig + 1)
    return factor * 0.5 * ulp


def _close_sig(a: float, b: float, spec: SigSpec) -> bool:
    """Return True if a and b agree to ~spec.sig significant digits."""
    if not (math.isfinite(a) and math.isfinite(b)):
        return a == b

    if abs(a) < spec.tiny and abs(b) < spec.tiny:
        return True

    scale = max(abs(a), abs(b))
    tol = _sig_tol(scale, spec.sig, spec.factor)
    return abs(a - b) <= max(tol, spec.tiny)


def _text_equal_sig(
    ref_text: str,
    got_text: str,
    spec: SigSpec,
) -> tuple[bool, str]:
    """Compare two vprofile raw_text blocks token-by-token.

    - Non-numeric tokens must match exactly.
    - Numeric tokens are compared to ~spec.sig significant digits.
    """
    ref_lines = [
        ln.rstrip()
        for ln in ref_text.replace(
            "\r\n",
            "\n",
        ).split("\n")
    ]
    got_lines = [
        ln.rstrip()
        for ln in got_text.replace(
            "\r\n",
            "\n",
        ).split("\n")
    ]

    # Drop trailing blank lines
    while ref_lines and ref_lines[-1] == "":
        ref_lines.pop()
    while got_lines and got_lines[-1] == "":
        got_lines.pop()

    if len(ref_lines) != len(got_lines):
        return (
            False,
            "Line"
            + " count"
            + " differs: "
            + f"ref={len(ref_lines)} "
            + f"got={len(got_lines)}",
        )

    for i, (rln, gln) in enumerate(zip(ref_lines, got_lines), start=1):
        rtoks = rln.split()
        gtoks = gln.split()

        # If naive whitespace tokenization disagrees (common with fixed-width
        # Fortran formats), fall back to a relaxed numeric tokenizer.
        if len(rtoks) != len(gtoks):
            rt_rel = _tokenize_relaxed(rln)
            gt_rel = _tokenize_relaxed(gln)
            if rt_rel and gt_rel and (len(rt_rel) == len(gt_rel)):
                rtoks, gtoks = rt_rel, gt_rel
            else:
                return (
                    False,
                    "Token count differs at line "
                    f"{i}: ref={len(rtoks)} got={len(gtoks)}\n"
                    f"ref: {rln}\n"
                    f"got: {gln}",
                )
        for j, (rt, gt) in enumerate(zip(rtoks, gtoks), start=1):
            ra = _parse_fortran_float(rt)
            ga = _parse_fortran_float(gt)
            if ra is not None and ga is not None:
                if not _close_sig(ra, ga, spec):
                    return (
                        False,
                        "Numeric mismatch at line "
                        f"{i}, "
                        f"token {j}: "
                        f"ref={rt} "
                        f"got={gt}",
                    )
            else:
                if rt != gt:
                    return (
                        False,
                        "Token mismatch at line "
                        f"{i}, token "
                        f"{j}: "
                        f"ref={rt!r} "
                        f"got={gt!r}",
                    )

    return True, ""


# -----------------------------------------------------------------------------
# Regression test
# -----------------------------------------------------------------------------


def test_vprofile_regression(
    scalefree_exe: Path, ref_dir: Path, tmp_path: Path
) -> None:
    """
    Compare vprofile raw output to stored references
    (~5 significant digits).
    """

    # IMPORTANT: keep these kwargs aligned with tests/make_vprofile_refs.py.
    # For algorithm=3 we exercise the vp-table path (usevp=True) and use a
    # smaller moment order (maxmom=8), which is both faster and more stable.
    base = dict(
        # model
        potential="logarithmic",
        gamma=2.0,
        q=0.608,
        beta=0.189,
        s=0.5,
        t=0.0,
        inclination=57.1,
        xi=0.0,
        theta=0.0,
        df=1,
        # run control
        maxmom=8,
        algorithm=3,
        kinematics="projected",
        integration=1,
        ngl_or_eps=0,
        usevp=True,
        verbose_vp=0,
        vp_smooth_eps=0.0,
        vp_reg_param=1.0,
        parse_stdout_fallback=False,
        debug_prompts=False,
        _skip_df1_beta_correction=True,
    )

    cases = {
        "projected_point": {
            "average": False,
            "ref": "projected_point_alg3_ref.txt",
        },
        "projected_avg": {
            "average": True,
            "ref": "projected_avg_alg3_ref.txt",
        },
    }

    runner = ScaleFreeRunner(scalefree_exe, workdir=tmp_path)
    # NOTE: in CI we sometimes see sub-e-8 variations
    # (and even underflow/formatting
    # differences like "0.1416-319" vs "0"). We treat |x|<1e-8 as zero.
    spec = SigSpec(sig=5, tiny=1e-8, factor=2.0)

    for name, cfg in cases.items():
        ref_path = ref_dir / cfg["ref"]
        if not ref_path.exists():
            pytest.fail(f"Missing reference file: {ref_path}")
        ref_text = ref_path.read_text(encoding="utf-8")

        # Provide a *short* output filename to avoid Fortran character-length
        # truncation issues on some platforms.
        out_path = tmp_path / f"vp_{name}.txt"

        res = runner.vprofile(
            **base,
            average=cfg["average"],
            output_path=out_path,
        )

        ok, msg = _text_equal_sig(ref_text, res.raw_text, spec)
        assert ok, msg
