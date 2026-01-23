"""Regression tests for the Fortran-backed vprofile interface.

The Fortran program prints floating-point values whose last digits can vary
across platforms/compilers and (in rare cases) across runs. To keep the test
stable and aligned with the project requirement, we normalise *every* numeric
token to **5 significant digits** before comparison.
"""

from __future__ import annotations

import re
from pathlib import Path

from scalefree.vmoments import ScaleFreeRunner


# Matches floats/ints in fixed or scientific notation.
_NUM_RE = re.compile(
    r"""(?x)
    (?P<num>
        [+-]?
        (?:
            (?:\d+\.\d*)|(?:\.\d+)|(?:\d+)
        )
        (?:[eE][+-]?\d+)?
    )
    """
)


def _norm_5sig(text: str) -> str:
    """Normalise numeric tokens to 5 significant digits."""

    def repl(match: re.Match[str]) -> str:
        tok = match.group("num")
        # Keep plain integers as integers (still stable and clearer in diffs)
        if re.fullmatch(r"[+-]?\d+", tok):
            return tok
        try:
            val = float(tok)
        except ValueError:
            return tok
        # Format with 5 significant digits; keep exponent where helpful.
        return f"{val:.5g}"

    return _NUM_RE.sub(repl, text).strip() + "\n"


def test_vprofile_regression(scalefree_exe: Path, ref_dir: Path, tmp_path: Path) -> None:
    """Compare vprofile raw output to stored references (5 significant digits)."""

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
    )

    cases = {
        "projected_point": {"average": False, "ref": "projected_point_alg3_ref.txt"},
        "projected_avg": {"average": True, "ref": "projected_avg_alg3_ref.txt"},
    }

    runner = ScaleFreeRunner(scalefree_exe, workdir=tmp_path)

    for name, cfg in cases.items():
        ref_text = (ref_dir / cfg["ref"]).read_text(encoding="utf-8")

        # IMPORTANT: provide a *short* output filename to avoid Fortran
        # character-length truncation issues on some platforms.
        out_path = tmp_path / f"vp_{name}.txt"

        res = runner.vprofile(**base, average=cfg["average"], output_path=out_path)

        got = _norm_5sig(res.raw_text)
        exp = _norm_5sig(ref_text)
        assert got == exp
