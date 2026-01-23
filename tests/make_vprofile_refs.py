#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate and (optionally) update vprofile reference outputs.

This script is meant to be run manually by maintainers before committing new
reference files under ``tests/data``.

Important: different Fortran compilers and platforms may format very small
numbers differently (e.g. ``0.1416-319`` vs ``0``). For that reason, reference
updates and comparisons are performed using *numerical* equality at roughly
5 significant digits, rather than strict text equality.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import math
import re
from typing import Tuple

from scalefree.vmoments import ScaleFreeRunner


# ------------------------------
# Numeric comparison utilities
# ------------------------------


_NUM_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")
_FORTRAN_NOE_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$")


def _parse_float_token(tok: str) -> float | None:
    """Parse a token that may be a Fortran/Python float.

    Supports:
      - standard E/e exponents (1.0e-3)
      - Fortran D exponent (1.0D-3)
      - Fortran output without 'E' (0.1416-319 meaning 0.1416e-319)
    """
    t = tok.strip()
    if not t:
        return None
    if _NUM_RE.match(t):
        try:
            return float(t.replace("D", "E").replace("d", "e"))
        except ValueError:
            return None
    m = _FORTRAN_NOE_RE.match(t)
    if m:
        base, exp = m.group(1), m.group(2)
        try:
            return float(f"{base}e{exp}")
        except ValueError:
            return None
    return None


@dataclass(frozen=True)
class SigSpec:
    sig: int = 5
    # Treat ultra-small values as zero for cross-platform stability.
    # (e.g. "0.1416-319" vs "0" in some Fortran/GFortran builds)
    tiny: float = 1e-8
    factor: float = 2.0  # tolerance multiplier (helps cross-platform)


def _sig_ulp_tol(x: float, sig: int, factor: float) -> float:
    """Half-ULP at the given significant digits (scaled by factor)."""
    ax = abs(x)
    if ax == 0.0:
        return 0.0
    exp10 = math.floor(math.log10(ax))
    # one unit in the last place for `sig` significant digits
    ulp = 10.0 ** (exp10 - sig + 1)
    return factor * 0.5 * ulp


def _close_sig(a: float, b: float, spec: SigSpec) -> bool:
    # Treat ultra-small values as zero for stability.
    if max(abs(a), abs(b)) < spec.tiny:
        return True
    # If one side is tiny, still compare against the other side's scale.
    scale = max(abs(a), abs(b))
    tol = _sig_ulp_tol(scale, spec.sig, spec.factor)
    return abs(a - b) <= max(tol, spec.tiny)


def _tokenise_lines(text: str) -> list[list[str]]:
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    # strip trailing whitespace; keep internal spacing irrelevant
    while lines and not lines[-1].strip():
        lines.pop()
    return [line.rstrip().split() for line in lines]


def texts_equivalent(
    ref_text: str,
    got_text: str,
    spec: SigSpec,
) -> Tuple[bool, str]:
    """Return (ok, message) comparing two outputs at ~sig digits."""
    ref_lines = _tokenise_lines(ref_text)
    got_lines = _tokenise_lines(got_text)
    if len(ref_lines) != len(got_lines):
        return (
            False,
            "line-count differs: "
            + f"expected {len(ref_lines)} "
            + f"got {len(got_lines)}",
        )

    for i, (rtoks, gtoks) in enumerate(zip(ref_lines, got_lines), start=1):
        if len(rtoks) != len(gtoks):
            return (
                False,
                "token-count differs on "
                + f"line {i}: "
                + f"expected {len(rtoks)} "
                + f"got {len(gtoks)}",
            )
        for j, (r, g) in enumerate(zip(rtoks, gtoks), start=1):
            rf = _parse_float_token(r)
            gf = _parse_float_token(g)
            if (rf is not None) and (gf is not None):
                if not _close_sig(rf, gf, spec):
                    return (
                        False,
                        "numeric mismatch on"
                        + f" line {i},"
                        + f" token {j}:"
                        + f" expected {r}"
                        + f" got {g}",
                    )
            else:
                if r != g:
                    return (
                        False,
                        "text mismatch on"
                        + f" line {i},"
                        + f" token {j}:"
                        + f" expected {r!r}"
                        + f" got {g!r}",
                    )
    return True, "ok"


# ------------------------------
# Reference generation
# ------------------------------


def _runner(exe: Path, workdir: Path) -> ScaleFreeRunner:
    return ScaleFreeRunner(exe, workdir=workdir)


def _common_kwargs(*, algorithm: int, kinematics: str) -> dict:
    # These are deliberately conservative parameters, selected to be stable.
    # Keep aligned with tests/test_vprofile_regression.py.
    base: dict = {
        # model
        "potential": "logarithmic",
        "gamma": 2.0,
        "q": 0.608,
        "beta": 0.189,
        "s": 0.5,
        "t": 0.0,
        "inclination": 57.1,
        "xi": 0.0,
        "theta": 0.0,
        "df": 1,
        # run control
        "algorithm": algorithm,
        "kinematics": kinematics,
        "integration": 1,
        "ngl_or_eps": 0,
        "debug_prompts": False,
    }

    if algorithm == 3:
        # algorithm=3 supports vp-table mode;
        # using fewer moments is both faster
        # and tends to be more stable across compilers.
        base.update(
            {
                "maxmom": 8,
                "usevp": True,
                "verbose_vp": 0,
                "vp_smooth_eps": 0.0,
                "vp_reg_param": 1.0,
                "parse_stdout_fallback": False,
            }
        )
    else:
        base.update({"maxmom": 20})

    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--exe",
        type=Path,
        default=Path(
            __file__,
        )
        .resolve()
        .parents[1]
        / "fortran_src"
        / "scalefree.e",
        help="Path to scalefree "
        + "Fortran executable "
        + "(default: repo fortran_src/scalefree.e)",
    )
    ap.add_argument(
        "--sig",
        type=int,
        default=5,
        help="Significant digits for comparison (default: 5)",
    )
    ap.add_argument(
        "--update",
        action="store_true",
        help="Overwrite reference files when"
        + " they differ at the "
        + "chosen precision.",
    )
    args = ap.parse_args()

    root = Path(__file__).resolve().parent
    data_dir = root / "data"
    workdir = root / "_work"
    data_dir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    runner = _runner(args.exe, workdir)
    spec = SigSpec(sig=args.sig)

    changed = False
    for algorithm in (1, 2, 3):
        for kinematics in ("intrinsic", "projected"):
            for average in (False, True):
                tag = "avg" if average else "point"
                out_name = f"{kinematics}_{tag}_alg{algorithm}_ref.txt"
                out_path = data_dir / out_name

                kwargs = _common_kwargs(
                    algorithm=algorithm,
                    kinematics=kinematics,
                )
                # Keep output paths short
                # (some Fortran builds truncate file names).
                tmp_out = workdir / "vp.txt"

                res = runner.vprofile(
                    average=average,
                    output_path=tmp_out,
                    **kwargs,
                )
                got_text = res.raw_text

                if out_path.exists():
                    ref_text = out_path.read_text(encoding="utf-8")
                    ok, msg = texts_equivalent(ref_text, got_text, spec)
                    if ok:
                        print(f"OK  {out_path.as_posix()}")
                        continue
                    if not args.update:
                        print(f"DIFF {out_path.name}: {msg}")
                        changed = True
                        continue

                # Write / update
                out_path.write_text(got_text, encoding="utf-8")
                print(f"Wrote {out_path.as_posix()}")
                changed = True

    if not changed:
        print(f"\nAll references match at ~{args.sig} significant digits.")
        return 0

    if args.update:
        print("\nReferences updated. Re-run without --update to confirm.")
        return 0

    print(f"\nSome references differ at ~{args.sig} significant digits.")
    print("Re-run with --update to overwrite reference files.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
