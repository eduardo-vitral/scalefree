#!/usr/bin/env python3
"""
gh_velocity_sampler.py

Reads per-sample Gauss-Hermite parameters and generates (vr, vtheta, vphi)
using balrogo.dynamics.mom_sample_generator.

Input format
------------
Ignore lines starting with '#'.
Each remaining line must have 18 columns:

  r theta phi x y z
  vr_mean vr_sig vr_h3 vr_h4
  vth_mean vth_sig vth_h3 vth_h4
  vph_mean vph_sig vph_h3 vph_h4

Output format
-------------
Header lines (# ...) and then N lines:
  vr  vtheta  vphi

Notes
-----
- In the scalefree context we set convolution errors eps=0.
- This implementation samples row-by-row (parameters vary with position).
  If you need speed for very large N, you can cache repeated parameter tuples.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def _read_table(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append([float(x) for x in s.split()])
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    arr = np.asarray(rows, dtype=float)
    if arr.shape[1] != 18:
        raise ValueError(f"Expected 18 columns, got {arr.shape[1]} in {path}")
    return arr


def _mom_sample(
    mean: float,
    sigma: float,
    h3: float,
    h4: float,
    n: int,
) -> np.ndarray:
    from balrogo import dynamics  # type: ignore

    mom_stats = np.zeros((4, 2), dtype=float)
    mom_stats[:, 0] = [mean, sigma, h3, h4]
    eps = np.zeros(n, dtype=float)

    samples = dynamics.mom_sample_generator(mom_stats, eps, nsig=10)

    if samples is None:
        # Fallback: Gaussian with mean/sigma
        return np.random.normal(loc=mean, scale=sigma, size=n).astype(float)

    samples = np.asarray(samples)
    if samples.ndim == 2 and samples.shape[1] == 1:
        samples = samples[:, 0]
    return samples.astype(float, copy=False)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(
            "Usage: gh_velocity_sampler.py <in_table.txt> <out_vel.txt>",
            file=sys.stderr,
        )
        return 2

    in_path = Path(argv[1])
    out_path = Path(argv[2])

    tab = _read_table(in_path)
    n = tab.shape[0]

    vr_params = tab[:, 6:10]
    vth_params = tab[:, 10:14]
    vph_params = tab[:, 14:18]

    vr = np.empty(n, dtype=float)
    vth = np.empty(n, dtype=float)
    vph = np.empty(n, dtype=float)

    for i in range(n):
        vr[i] = _mom_sample(*vr_params[i], n=1)[0]
        vth[i] = _mom_sample(*vth_params[i], n=1)[0]
        vph[i] = _mom_sample(*vph_params[i], n=1)[0]

    with out_path.open("w", encoding="utf-8") as f:
        f.write("# kind=mcdraw_gh_output\n")
        f.write("# columns: vr vtheta vphi\n")
        for i in range(n):
            f.write(f"{vr[i]: .16e} {vth[i]: .16e} {vph[i]: .16e}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
