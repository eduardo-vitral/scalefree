#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""scripts/check_mock_end2end.py

End-to-end check for the intrinsic-θ mock generator.

What this script validates
--------------------------
The current mock generator (``scalefree.mock.mock``) returns samples in the
model Cartesian frame. This check therefore:

1) Runs ``mock()`` to get (x,y,z,vx,vy,vz).
2) Converts the Cartesian velocities back to intrinsic spherical components
   (vr, vtheta, vphi) at each sampled position.
3) Compares *shell-averaged* intrinsic second moments against the analytic
   ScaleFree backend result (``vprofile(..., kinematics='intrinsic', average=True)``).

This intentionally avoids LOS / POSr / POSt geometry and focuses on the
intrinsic formalism consistency.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


def _force_repo_import(repo_root: Path) -> None:
    """Ensure we import *this repo* version of scalefree, not an installed wheel."""
    repo_root = repo_root.resolve()
    sys.path.insert(0, str(repo_root))
    for k in list(sys.modules.keys()):
        if k == "scalefree" or k.startswith("scalefree."):
            del sys.modules[k]


def _cart_to_sph_angles(
    x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.sqrt(x * x + y * y + z * z)
    # theta in [0, pi], phi in [0, 2pi)
    theta = np.arccos(np.clip(z / np.where(r > 0, r, 1.0), -1.0, 1.0))
    phi = (np.arctan2(y, x) + 2.0 * np.pi) % (2.0 * np.pi)
    return r, theta, phi


def _cart_to_sph_v(
    *,
    theta: np.ndarray,
    phi: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project Cartesian velocity onto (e_r, e_theta, e_phi)."""
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    # e_r     = ( st cp,  st sp,  ct)
    # e_theta = ( ct cp,  ct sp, -st)
    # e_phi   = (   -sp,    cp,   0)
    vr = vx * st * cp + vy * st * sp + vz * ct
    vtheta = vx * ct * cp + vy * ct * sp - vz * st
    vphi = -vx * sp + vy * cp
    return vr, vtheta, vphi


def _summarize_intrinsic(
    vr: np.ndarray, vtheta: np.ndarray, vphi: np.ndarray
) -> Dict[str, float]:
    vr = np.asarray(vr, dtype=float)
    vtheta = np.asarray(vtheta, dtype=float)
    vphi = np.asarray(vphi, dtype=float)

    out: Dict[str, float] = {}
    out["vr_mean"] = float(np.mean(vr))
    out["vtheta_mean"] = float(np.mean(vtheta))
    out["vphi_mean"] = float(np.mean(vphi))

    out["vr2"] = float(np.mean(vr * vr))
    out["vth2"] = float(np.mean(vtheta * vtheta))
    out["vphi2"] = float(np.mean(vphi * vphi))

    # velocity anisotropy beta_v = 1 - (vth^2 + vphi^2)/(2 vr^2)
    denom = max(out["vr2"], 1e-30)
    out["beta"] = float(1.0 - (out["vth2"] + out["vphi2"]) / (2.0 * denom))
    return out


def _read_intrinsic_shell_average(res) -> Dict[str, float]:
    blk = res.blocks.get("intrinsic_shell_average")
    if blk is None:
        raise RuntimeError("vprofile did not return 'intrinsic_shell_average'.")

    cols = list(blk.get("columns", []))
    data = np.asarray(blk.get("data"))
    if data.ndim != 2 or data.shape[0] < 1:
        raise RuntimeError("'intrinsic_shell_average' block is empty or malformed.")
    row = data[0]

    def _get(name: str) -> float:
        if name not in cols:
            raise KeyError(
                f"Missing column '{name}' in intrinsic_shell_average cols={cols}"
            )
        return float(row[cols.index(name)])

    return {
        "rho": _get("rho"),
        "vphi_mean": _get("vphi"),
        "vr2": _get("vr2"),
        "vth2": _get("vth2"),
        "vphi2": _get("vphi2"),
        "beta": _get("beta"),
    }


def _rel_err(a: float, b: float) -> float:
    denom = max(abs(b), 1e-30)
    return abs(a - b) / denom


def main() -> int:
    p = argparse.ArgumentParser(
        description="End-to-end check for scalefree.mock (intrinsic version)."
    )
    p.add_argument(
        "--repo",
        type=str,
        default=str(Path(__file__).resolve().parents[1]),
        help="Repo root (for local import).",
    )

    # Model
    p.add_argument(
        "--potential",
        type=int,
        default=1,
        choices=[1, 2],
        help="1=Kepler, 2=Logarithmic",
    )
    p.add_argument("--gamma", type=float, default=4.0)
    p.add_argument("--q", type=float, default=1.0)
    p.add_argument("--df", type=int, default=1, choices=[1, 2])
    p.add_argument("--beta", type=float, default=0.0)
    p.add_argument("--s", type=float, default=0.5)
    p.add_argument("--t", type=float, default=0.0)

    # Mock settings
    p.add_argument("--nsamples", type=int, default=50_000)
    p.add_argument("--theta-bins", type=int, default=45)
    p.add_argument("--rin", type=float, default=1.0)
    p.add_argument("--rout", type=float, default=1_000.0)
    p.add_argument("--seed", type=int, default=1234)
    p.add_argument("--maxmom", type=int, default=4)
    p.add_argument("--vp-smooth-eps", type=float, default=0.0)
    p.add_argument("--integration", type=int, default=1, choices=[0, 1])
    p.add_argument("--ngl-or-eps", type=float, default=0.0)
    p.add_argument("--nsig", type=int, default=10)

    # Validation tolerances
    p.add_argument(
        "--rtol",
        type=float,
        default=0.07,
        help="Relative tolerance on 2nd moments/beta (default 7%).",
    )
    p.add_argument(
        "--atol-mean",
        type=float,
        default=0.05,
        help="Absolute tolerance on means (default 0.05).",
    )

    p.add_argument(
        "--plot", action="store_true", help="Show histograms vs reference Gaussians."
    )

    # args = p.parse_args()
    # Jupyter injects extra args like --f=...json; ignore anything unknown
    args, _unknown = p.parse_known_args()

    _force_repo_import(Path(args.repo))

    from scalefree import ScaleFreeRunner
    from scalefree.mock import mock as generate_mock

    # -----------------
    # Generate mock
    # -----------------
    X = generate_mock(
        potential=lambda: int(args.potential),
        gamma=float(args.gamma),
        q=float(args.q),
        df=int(args.df),
        beta=float(args.beta),
        s=float(args.s),
        t=float(args.t),
        nsamples=int(args.nsamples),
        theta_bins=int(args.theta_bins),
        rin=float(args.rin),
        rout=float(args.rout),
        seed=int(args.seed),
        integration=int(args.integration),
        ngl_or_eps=float(args.ngl_or_eps),
        maxmom=int(args.maxmom),
        vp_smooth_eps=float(args.vp_smooth_eps),
        nsig=int(args.nsig),
        debug=False,
    )

    x, y, z, vx, vy, vz = X.T
    _r, theta, phi = _cart_to_sph_angles(x, y, z)
    vr, vtheta, vphi = _cart_to_sph_v(theta=theta, phi=phi, vx=vx, vy=vy, vz=vz)
    samp = _summarize_intrinsic(vr, vtheta, vphi)

    # -----------------
    # Analytic reference
    # -----------------
    runner = ScaleFreeRunner()
    res = runner.vprofile(
        potential=lambda: int(args.potential),
        gamma=float(args.gamma),
        q=float(args.q),
        df=int(args.df),
        beta=float(args.beta),
        s=float(args.s),
        t=float(args.t),
        inclination=90.0,
        xi=0.0,
        theta=0.0,
        integration=int(args.integration),
        ngl_or_eps=float(args.ngl_or_eps),
        algorithm=3,
        vp_smooth_eps=float(args.vp_smooth_eps),
        maxmom=int(args.maxmom),
        average=True,
        kinematics="intrinsic",
        usevp=False,
        verbose_vp=0,
        output_path=None,
        debug_prompts=False,
    )
    ref = _read_intrinsic_shell_average(res)

    # -----------------
    # Report + checks
    # -----------------
    print("\n=== Model ===")
    print(
        f"potential={args.potential} gamma={args.gamma} q={args.q} df={args.df} "
        f"beta={args.beta} s={args.s} t={args.t}"
    )
    print("\n=== Mock settings ===")
    print(
        f"nsamples={args.nsamples} theta_bins={args.theta_bins} rin={args.rin} rout={args.rout} "
        f"seed={args.seed} maxmom={args.maxmom}"
    )

    print("\n=== Intrinsic shell-average moments: sample vs reference ===")
    rows = [
        ("<vphi>", samp["vphi_mean"], ref["vphi_mean"]),
        ("<vr^2>", samp["vr2"], ref["vr2"]),
        ("<vtheta^2>", samp["vth2"], ref["vth2"]),
        ("<vphi^2>", samp["vphi2"], ref["vphi2"]),
        ("beta_v", samp["beta"], ref["beta"]),
    ]
    ok = True
    for name, v_s, v_r in rows:
        if name == "<vphi>":
            err = abs(v_s - v_r)
            passed = err <= float(args.atol_mean)
            print(
                f"{name:10s}  sample={v_s: .6e}  ref={v_r: .6e}  abs_err={err: .3e}  {'OK' if passed else 'FAIL'}"
            )
        else:
            err = _rel_err(v_s, v_r)
            passed = err <= float(args.rtol)
            print(
                f"{name:10s}  sample={v_s: .6e}  ref={v_r: .6e}  rel_err={err: .3e}  {'OK' if passed else 'FAIL'}"
            )
        ok = ok and passed

    if args.plot:
        import matplotlib.pyplot as plt

        def _plot_hist(vals: np.ndarray, mu: float, sig2: float, title: str) -> None:
            v = np.asarray(vals, dtype=float)
            v = v[np.isfinite(v)]
            sig = float(np.sqrt(max(sig2, 1e-30)))
            plt.figure()
            plt.hist(v, bins=100, density=True)
            grid = np.linspace(mu - 5 * sig, mu + 5 * sig, 400)
            g = (1.0 / (np.sqrt(2.0 * np.pi) * sig)) * np.exp(
                -0.5 * ((grid - mu) / sig) ** 2
            )
            plt.plot(grid, g)
            plt.title(title)
            plt.xlabel("v")
            plt.ylabel("pdf")

        _plot_hist(vr, 0.0, ref["vr2"], "vr: histogram vs N(0, <vr^2>)")
        _plot_hist(vtheta, 0.0, ref["vth2"], "vtheta: histogram vs N(0, <vtheta^2>)")
        _plot_hist(
            vphi,
            ref["vphi_mean"],
            max(ref["vphi2"] - ref["vphi_mean"] ** 2, 1e-30),
            "vphi: histogram vs Gaussian moments",
        )
        plt.show()

    if not ok:
        print("\nFAILED: at least one tolerance check did not pass.")
        return 1

    print(
        "\nPASSED: mock shell-averaged intrinsic moments are consistent with vprofile."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
