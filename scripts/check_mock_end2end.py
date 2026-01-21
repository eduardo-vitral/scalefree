#!/usr/bin/env python3
"""
Local end-to-end validation of scalefree mock generation (mock-focused).

Key features
-----------
- Forces import from the *repo working tree* (not pip install).
- Uses algorithm=1 for speed.
- Uses <=10k samples by default.
- Produces quick diagnostic plots into scripts/mock_check_outputs/.

What it checks
--------------
A) Runs ScaleFreeRunner.vprofile once (algorithm=1) and extracts:
   - vp summary rows for iproj=1,2,3 (gauss_V, gauss_sig, h3, h4 + true_* if present)
   - projected_circle_average moments (if present)
   - vp_table for iproj=1,2,3 (if present)

B) Density sanity check (q=1):
   - empirical shell density vs r^{-gamma} shape overlay

C) Flattening sanity check (q<1):
   - scatter in x-z plane with ellipse overlay

D) Velocity sanity check:
   - generate mock at rin=rout=1 (single radius)
   - compute (vlos, vposr, vpost) from returned (vx,vy,vz) and plot histograms
   - overlay ScaleFree vp_table if available
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Helpers: force repo import
# -----------------------------


def _force_repo_import(repo_root: Path) -> None:
    """
    Force Python to import 'scalefree' from the working tree (repo checkout),
    not from a pip/conda installed distribution.
    """
    repo_root = repo_root.resolve()
    sys.path.insert(0, str(repo_root))

    # If scalefree was already imported (common in notebooks), drop it.
    for k in list(sys.modules.keys()):
        if k == "scalefree" or k.startswith("scalefree."):
            del sys.modules[k]


# -----------------------------
# Density utilities
# -----------------------------


def volume_density_powerlaw_shape(gamma: float, r: np.ndarray) -> np.ndarray:
    """Shape-only 3D density for rho(r) ∝ r^{-gamma}."""
    return np.power(r, -gamma)


def shell_volume_density(
    r: np.ndarray, nbins: int = 25
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Empirical 3D density in spherical shells: rho_hat ~ counts / shell_volume.
    Uses log-spaced bins.
    Returns rmid (geometric midpoint), rho_hat for bins with counts>0.
    """
    r = np.asarray(r, dtype=float)
    r = r[np.isfinite(r) & (r > 0)]
    if r.size == 0:
        raise ValueError("No valid radii to bin for volume density.")

    edges = np.logspace(np.log10(r.min()), np.log10(r.max()), nbins + 1)
    counts, _ = np.histogram(r, bins=edges)

    vol = (4.0 / 3.0) * math.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho_hat = counts / vol
    rmid = np.sqrt(edges[1:] * edges[:-1])

    m = counts > 0
    return rmid[m], rho_hat[m]


# -----------------------------
# Velocity projection utilities
# -----------------------------


def project_vel_components_from_xyzv(
    xyzv: np.ndarray, los_axis: str = "z"
) -> Dict[int, np.ndarray]:
    """
    Convert (x,y,z,vx,vy,vz) to projected velocity components corresponding to ScaleFree iproj:
      iproj=1 LOS
      iproj=2 POSR (projected radial)
      iproj=3 POST (projected tangential)

    Default: LOS=z, plane-of-sky=(x,y).
    """
    x, y, z, vx, vy, vz = xyzv.T

    if los_axis == "z":
        X, Y = x, y
        vX, vY = vx, vy
        vlos = vz
    elif los_axis == "y":
        X, Y = x, z
        vX, vY = vx, vz
        vlos = vy
    elif los_axis == "x":
        X, Y = y, z
        vX, vY = vy, vz
        vlos = vx
    else:
        raise ValueError("los_axis must be one of {'x','y','z'}")

    R = np.sqrt(X * X + Y * Y)
    m = R > 0

    vposr = np.full_like(vlos, np.nan, dtype=float)
    vpost = np.full_like(vlos, np.nan, dtype=float)

    vposr[m] = (X[m] * vX[m] + Y[m] * vY[m]) / R[m]
    vpost[m] = (-Y[m] * vX[m] + X[m] * vY[m]) / R[m]

    return {1: vlos, 2: vposr, 3: vpost}


# -----------------------------
# vmoments extraction (robust)
# -----------------------------


def extract_vmoments_products(res) -> Dict[str, Any]:
    """
    Extract key quantities from ScaleFreeRunner.vprofile result blocks.

    Robust to two common structures:
    - vp block as {"columns":[...], "data": (nrow,ncol)}
    - projected_circle_average similarly
    - vp_table as dict iproj -> {"columns":[...], "data":...}
    """
    blocks = getattr(res, "blocks", {}) or {}

    out: Dict[str, Any] = {"vp_rows": {}, "proj_avg": {}, "vp_table": {}}

    # ---- vp summary rows by iproj
    vp = blocks.get("vp")
    if isinstance(vp, dict):
        cols = list(vp.get("columns", []))
        data = vp.get("data")
        if cols and isinstance(data, np.ndarray) and data.ndim == 2 and data.size:
            try:
                i_ip = cols.index("iproj")
            except ValueError:
                i_ip = None

            def _get(col: str, row: np.ndarray) -> float:
                if col in cols:
                    return float(row[cols.index(col)])
                return float("nan")

            if i_ip is not None:
                for row in data:
                    ip = int(row[i_ip])
                    out["vp_rows"][ip] = {
                        "true_V": _get("true_V", row),
                        "true_sig": _get("true_sig", row),
                        "gauss_V": _get("gauss_V", row),
                        "gauss_sig": _get("gauss_sig", row),
                        "h3": _get("h3", row),
                        "h4": _get("h4", row),
                    }

    # ---- projected_circle_average (if present)
    pca = blocks.get("projected_circle_average")
    if isinstance(pca, dict):
        cols = list(pca.get("columns", []))
        data = pca.get("data")
        if cols and isinstance(data, np.ndarray) and data.ndim == 2 and data.size:
            # expected: iproj rho_p v1 v2 v3 v4 ...
            try:
                i_ip = cols.index("iproj")
            except ValueError:
                i_ip = None

            def _get(col: str, row: np.ndarray) -> float:
                if col in cols:
                    return float(row[cols.index(col)])
                return float("nan")

            if i_ip is not None:
                for row in data:
                    ip = int(row[i_ip])
                    out["proj_avg"][ip] = {
                        "rho_p": _get("rho_p", row),
                        "v1": _get("v1", row),
                        "v2": _get("v2", row),
                        "v3": _get("v3", row),
                        "v4": _get("v4", row),
                    }

    # ---- vp_table
    vpt = blocks.get("vp_table")
    if isinstance(vpt, dict):
        for iproj, tbl in vpt.items():
            if not isinstance(tbl, dict):
                continue
            arr = tbl.get("data")
            if (
                isinstance(arr, np.ndarray)
                and arr.ndim == 2
                and arr.size
                and arr.shape[1] >= 2
            ):
                out["vp_table"][int(iproj)] = (
                    arr[:, 0].astype(float),
                    arr[:, 1].astype(float),
                )

    return out


# -----------------------------
# Main
# -----------------------------


def main() -> None:
    p = argparse.ArgumentParser(
        description="Local scalefree mock check (repo import, fast settings)."
    )
    p.add_argument(
        "--n_dens",
        type=int,
        default=8000,
        help="N for density check (<=10000 recommended).",
    )
    p.add_argument(
        "--n_vel",
        type=int,
        default=8000,
        help="N for velocity check (<=10000 recommended).",
    )
    p.add_argument(
        "--nbins",
        type=int,
        default=18,
        help="Angular bins for mock() (smaller = faster).",
    )
    p.add_argument(
        "--outdir",
        type=str,
        default="mock_check_outputs",
        help="Output directory (relative to scripts/).",
    )
    p.add_argument(
        "--seed", type=int, default=42, help="Seed (use non-negative for default_rng)."
    )
    args, _unknown = p.parse_known_args()  # ignore ipykernel flags

    # Force local repo import (script must live in scripts/)
    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner, mock  # noqa: E402

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    outdir = (Path(__file__).resolve().parent / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # -----------------------------
    # Model + view (fast)
    # -----------------------------
    model = dict(
        potential=2,  # logarithmic
        gamma=2.0,
        q=1.0,
        df=1,
        beta=0.189,
        s=0.5,
        t=0.0,
    )
    view = dict(
        inclination=57.1,
        xi=0.0,
        theta=0.0,  # projected-plane angle used by vprofile
    )

    # -----------------------------
    # A) vprofile products (single call)
    # -----------------------------
    runner = ScaleFreeRunner()
    res = runner.vprofile(
        **model,
        **view,
        integration=1,
        ngl_or_eps=0,  # for speed / default
        algorithm=1,  # IMPORTANT: fast
        maxmom=4,
        average=True,
        usevp=True,
        verbose_vp=0,
        output_path=None,
        timeout_s=300,
        debug_prompts=False,
    )
    prod = extract_vmoments_products(res)

    print("\n=== Extracted vprofile quantities ===")
    for ip in (1, 2, 3):
        row = prod["vp_rows"].get(ip, {})
        if not row:
            print(f"iproj={ip}: (no vp row parsed)")
            continue
        print(
            f"iproj={ip}: gauss_V={row['gauss_V']:.6g}, gauss_sig={row['gauss_sig']:.6g}, "
            f"h3={row['h3']:.6g}, h4={row['h4']:.6g} | "
            f"true_V={row['true_V']:.6g}, true_sig={row['true_sig']:.6g}"
        )

    # -----------------------------
    # B) Density check (q=1)
    # -----------------------------
    print(f"\n=== Density check: generating mock positions (N={args.n_dens}) ===")
    xyzv = mock(
        **model,
        nsamples=int(args.n_dens),
        seed=int(args.seed),
        rin=0.5,
        rout=50.0,
        algorithm=1,
        nbins=int(args.nbins),
        debug=True,
    )

    r3d = np.sqrt(xyzv[:, 0] ** 2 + xyzv[:, 1] ** 2 + xyzv[:, 2] ** 2)
    rmid, rho_hat = shell_volume_density(r3d, nbins=25)

    rho_th = volume_density_powerlaw_shape(model["gamma"], rmid)
    k0 = len(rmid) // 2
    scale = rho_hat[k0] / rho_th[k0]
    rho_th_scaled = scale * rho_th

    fig = plt.figure()
    plt.loglog(rmid, rho_hat, marker="o", linestyle="none", label="Mock: shell rho(r)")
    plt.loglog(
        rmid,
        rho_th_scaled,
        linestyle="-",
        label=r"Analytic: $\rho \propto r^{-\gamma}$ (scaled)",
    )
    plt.xlabel("3D radius r")
    plt.ylabel("Volume density (arb.)")
    plt.title(f"Density check (q=1, gamma={model['gamma']})")
    plt.legend()
    fig.savefig(outdir / "density_q1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # -----------------------------
    # C) Flattening check (q<1)
    # -----------------------------
    print("\n=== Flattening check: q<1 ===")
    model_flat = dict(model)
    model_flat["q"] = 0.408

    xyzv_flat = mock(
        **model_flat,
        nsamples=min(8000, int(args.n_dens)),
        seed=int(args.seed) + 1,
        rin=0.5,
        rout=50.0,
        algorithm=1,
        nbins=int(args.nbins),
        debug=False,
    )

    X = xyzv_flat[:, 0]
    Z = xyzv_flat[:, 2]
    q0 = model_flat["q"]

    # Ellipse overlay: R_e = sqrt(x^2 + (z/q)^2)
    Re = np.sqrt(X * X + (Z / q0) ** 2)
    Re0 = np.nanpercentile(Re, 60)

    ang = np.linspace(0, 2 * np.pi, 400)
    ex = Re0 * np.cos(ang)
    ez = q0 * Re0 * np.sin(ang)

    fig = plt.figure()
    idx = np.random.default_rng(int(args.seed)).choice(
        len(X), size=min(6000, len(X)), replace=False
    )
    plt.scatter(X[idx], Z[idx], s=2, alpha=0.25, label="Mock points (x-z)")
    plt.plot(ex, ez, linewidth=2, label=f"Ellipse overlay (axis ratio q={q0})")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"Flattening sanity check (q={q0})")
    plt.legend()
    fig.savefig(outdir / "flattening_overlay.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # -----------------------------
    # D) Velocity checks (single radius)
    # -----------------------------
    print(f"\n=== Velocity check: single radius mock (N={args.n_vel}) ===")
    xyzv_vel = mock(
        **model,
        nsamples=int(args.n_vel),
        seed=int(args.seed) + 2,
        rin=1.0,
        rout=2.0,
        algorithm=1,
        nbins=int(args.nbins),
        debug=True,
    )

    vproj = project_vel_components_from_xyzv(xyzv_vel, los_axis="z")

    for ip in (1, 2, 3):
        v = vproj[ip]
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue

        # trim extreme tails for nicer plots
        lim = np.nanpercentile(np.abs(v), 99.7)
        v = v[np.abs(v) <= lim]

        fig = plt.figure()
        plt.hist(v, bins=80, density=True, alpha=0.35, label="Mock histogram")

        # ScaleFree vp_table overlay (if present)
        if ip in prod["vp_table"]:
            vg, vpg = prod["vp_table"][ip]
            area = np.trapz(vpg, vg)
            if area > 0:
                vpg = vpg / area
            plt.plot(vg, vpg, linewidth=2, label="ScaleFree vp_table (normalized)")

        row = prod["vp_rows"].get(ip, {})
        if row:
            plt.axvline(row["gauss_V"], linewidth=1.2, label="gauss_V (vprofile)")
            plt.text(
                0.05,
                0.85,
                f"gauss_V={row['gauss_V']:.3g}\ngauss_sig={row['gauss_sig']:.3g}\nh3={row['h3']:.3g}, h4={row['h4']:.3g}",
                transform=plt.gca().transAxes,
                va="top",
            )

        label = {1: "LOS", 2: "POSR", 3: "POST"}[ip]
        plt.xlabel(f"Velocity (iproj={ip}, {label})")
        plt.ylabel("PDF")
        plt.title(f"Velocity sanity (algorithm=1, iproj={ip})")
        plt.legend()

        fig.savefig(outdir / f"velocity_iproj{ip}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    print(f"\nDone. Outputs saved in: {outdir}")


if __name__ == "__main__":
    main()
