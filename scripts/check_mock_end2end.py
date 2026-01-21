#!/usr/bin/env python3
"""
Local end-to-end validation of scalefree mock generation (mock-focused).

What this does (high level)
---------------------------
A) Runs ScaleFreeRunner.vprofile once (algorithm=1) and extracts:
   - vp summary rows for iproj=1,2,3 (gauss_V, gauss_sig, h3, h4)

B) Density sanity check (q=1):
   - empirical shell density vs r^{-gamma} shape overlay

C) Flattening sanity check (q<1):
   - scatter in x-z plane with ellipse overlay

D) Velocity sanity check:
   - generate mock at rin..rout ~ [1,2]
   - compute (vlos, vposr, vpost) from (vx,vy,vz)
   - plot histograms
   - overlay analytic curve computed from BALRoGO using:
        vg = linspace(min(v), max(v), 1000)
        fits = [gauss_V, gauss_sig, h3, h4]
        fgrid = dynamics.mom_likelihood_func(fits, vg, zeros, mode="curve")
     If BALRoGO returns scalar inf (non-physical moments), we fall back to
     dynamics.laplace_kernel_pdf / dynamics.uniform_kernel_pdf (same kernel
     family, but bypasses the rejection in mom_likelihood_func).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any, Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Helpers: force repo import
# -----------------------------
def _force_repo_import(repo_root: Path) -> None:
    repo_root = repo_root.resolve()
    sys.path.insert(0, str(repo_root))
    for k in list(sys.modules.keys()):
        if k == "scalefree" or k.startswith("scalefree."):
            del sys.modules[k]


def _import_balrogo_dynamics():
    try:
        from balrogo import dynamics  # type: ignore
        return dynamics
    except Exception:
        import importlib

        return importlib.import_module("dynamics")


# -----------------------------
# Density utilities
# -----------------------------
def volume_density_powerlaw_shape(gamma: float, r: np.ndarray) -> np.ndarray:
    return np.power(r, -gamma)


def shell_volume_density(r: np.ndarray, nbins: int = 25) -> Tuple[np.ndarray, np.ndarray]:
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
def project_vel_components_from_xyzv(xyzv: np.ndarray, los_axis: str = "z") -> Dict[int, np.ndarray]:
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
# vprofile extraction
# -----------------------------
def extract_vmoments_products(res) -> Dict[str, Any]:
    blocks = getattr(res, "blocks", {}) or {}
    out: Dict[str, Any] = {"vp_rows": {}, "vp_table": {}}

    vp = blocks.get("vp")
    if isinstance(vp, dict):
        cols = list(vp.get("columns", []))
        data = vp.get("data")
        if cols and isinstance(data, np.ndarray) and data.ndim == 2 and data.size and "iproj" in cols:
            i_ip = cols.index("iproj")

            def _get(col: str, row: np.ndarray) -> float:
                return float(row[cols.index(col)]) if col in cols else float("nan")

            for row in data:
                ip = int(row[i_ip])
                out["vp_rows"][ip] = {
                    "gauss_V": _get("gauss_V", row),
                    "gauss_sig": _get("gauss_sig", row),
                    "h3": _get("h3", row),
                    "h4": _get("h4", row),
                    "true_V": _get("true_V", row),
                    "true_sig": _get("true_sig", row),
                }

    vpt = blocks.get("vp_table")
    if isinstance(vpt, dict):
        for iproj, tbl in vpt.items():
            if not isinstance(tbl, dict):
                continue
            arr = tbl.get("data")
            if isinstance(arr, np.ndarray) and arr.ndim == 2 and arr.size and arr.shape[1] >= 2:
                out["vp_table"][int(iproj)] = (arr[:, 0].astype(float), arr[:, 1].astype(float))

    return out


# -----------------------------
# BALRoGO curve helpers
# -----------------------------
def balrogo_curve_from_fits(dyn, fits: np.ndarray, vg: np.ndarray) -> Tuple[np.ndarray, str]:
    """
    Compute a pointwise PDF curve on vg using BALRoGO.

    Primary (requested) path:
        fgrid = dyn.mom_likelihood_func(fits, vg, zeros, mode="curve")

    If BALRoGO rejects the moments and returns scalar inf, fall back to:
        dyn.laplace_kernel_pdf(...) or dyn.uniform_kernel_pdf(...)

    Returns
    -------
    (fgrid, tag) where tag indicates which path was used.
    """
    mu, sig, h3, h4 = map(float, fits)
    ex = np.zeros_like(vg, dtype=float)

    ret = dyn.mom_likelihood_func(fits, vg, ex, mode="curve")

    # If mom_likelihood_func returns a curve array, use it
    if isinstance(ret, np.ndarray) and ret.ndim == 1 and ret.shape[0] == vg.shape[0]:
        return ret.astype(float), "mom_likelihood_func(curve)"

    # If it returns scalar inf (common when moments deemed non-physical), fallback
    if np.isscalar(ret) and (not np.isfinite(ret) or float(ret) == float("inf")):
        # Use the same kernel family directly, bypassing the early rejection
        if h4 >= 0:
            fgrid = dyn.laplace_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "laplace_kernel_pdf(fallback)"
        else:
            fgrid = dyn.uniform_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "uniform_kernel_pdf(fallback)"

    # If it's some other scalar, try fallback anyway
    if np.isscalar(ret):
        if h4 >= 0:
            fgrid = dyn.laplace_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "laplace_kernel_pdf(fallback)"
        else:
            fgrid = dyn.uniform_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "uniform_kernel_pdf(fallback)"

    raise RuntimeError(
        "Unexpected return type from dyn.mom_likelihood_func(..., mode='curve'). "
        f"type={type(ret).__name__}"
    )


def normalize_curve(vg: np.ndarray, fgrid: np.ndarray) -> np.ndarray:
    fgrid = np.asarray(fgrid, dtype=float)
    if fgrid.ndim != 1 or fgrid.shape[0] != vg.shape[0]:
        return fgrid

    if hasattr(np, "trapezoid"):
        area = np.trapezoid(fgrid, vg)
    else:
        area = np.trapz(fgrid, vg)

    if np.isfinite(area) and area > 0:
        return fgrid / area
    return fgrid


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    p = argparse.ArgumentParser(description="Local scalefree mock check (fast) with improved velocity overlay.")
    p.add_argument("--n_dens", type=int, default=8000, help="N for density check.")
    p.add_argument("--n_vel", type=int, default=8000, help="N for velocity check.")
    p.add_argument("--nbins", type=int, default=18, help="Angular bins for mock().")
    p.add_argument("--hist_bins", type=int, default=80, help="Histogram bins for velocities.")
    p.add_argument("--outdir", type=str, default="mock_check_outputs", help="Output directory (relative to scripts/).")
    p.add_argument("--seed", type=int, default=42, help="Seed.")
    args, _unknown = p.parse_known_args()

    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner, mock  # noqa: E402

    dyn = _import_balrogo_dynamics()

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    outdir = (Path(__file__).resolve().parent / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    model = dict(
        potential=2,  # logarithmic
        gamma=2.0,
        q=1.0,
        df=1,
        beta=0.189,
        s=0.5,
        t=0.0,
    )
    view = dict(inclination=57.1, xi=0.0, theta=0.0)

    # A) vprofile products
    runner = ScaleFreeRunner()
    res = runner.vprofile(
        **model,
        **view,
        integration=1,
        ngl_or_eps=0,
        algorithm=1,
        maxmom=4,
        average=True,
        usevp=True,
        verbose_vp=0,
        output_path=None,
        timeout_s=300,
        debug_prompts=False,
    )
    prod = extract_vmoments_products(res)

    print("\n=== vprofile summary rows ===")
    for ip in (1, 2, 3):
        row = prod["vp_rows"].get(ip, {})
        if row:
            print(
                f"iproj={ip}: gauss_V={row['gauss_V']:.6g}, gauss_sig={row['gauss_sig']:.6g}, "
                f"h3={row['h3']:.6g}, h4={row['h4']:.6g}"
            )
        else:
            print(f"iproj={ip}: (missing)")

    # B) Density sanity check (q=1)
    print(f"\n=== Density check: N={args.n_dens} ===")
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
    plt.loglog(rmid, rho_th_scaled, linestyle="-", label=r"Analytic: $r^{-\gamma}$ (scaled)")
    plt.xlabel("3D radius r")
    plt.ylabel("Volume density (arb.)")
    plt.title(f"Density check (q=1, gamma={model['gamma']})")
    plt.legend()
    fig.savefig(outdir / "density_q1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # C) Flattening sanity check (q<1)
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
    q0 = float(model_flat["q"])

    Re = np.sqrt(X * X + (Z / q0) ** 2)
    Re0 = float(np.nanpercentile(Re, 60))

    ang = np.linspace(0, 2 * np.pi, 400)
    ex = Re0 * np.cos(ang)
    ez = q0 * Re0 * np.sin(ang)

    fig = plt.figure()
    idx = np.random.default_rng(int(args.seed)).choice(len(X), size=min(6000, len(X)), replace=False)
    plt.scatter(X[idx], Z[idx], s=2, alpha=0.25, label="Mock points (x-z)")
    plt.plot(ex, ez, linewidth=2, label=f"Ellipse overlay (axis ratio q={q0})")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"Flattening sanity check (q={q0})")
    plt.legend()
    fig.savefig(outdir / "flattening_overlay.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # D) Velocity sanity check + BALRoGO overlay
    print(f"\n=== Velocity check: N={args.n_vel} ===")
    model_flat["q"] = 1.0
    model_flat["q"] = 0.5
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
        v = np.asarray(vproj[ip], dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue

        fig = plt.figure()
        plt.hist(v, bins=int(args.hist_bins), density=True, alpha=0.35, label="Mock histogram")

        # Optional: vp_table overlay (if present)
        if ip in prod["vp_table"]:
            vg_tab, vpg_tab = prod["vp_table"][ip]
            vpg_tab = normalize_curve(vg_tab, vpg_tab)
            plt.plot(vg_tab, vpg_tab, linewidth=2, label="ScaleFree vp_table (normalized)")

        # Required: BALRoGO overlay from [gauss_V, gauss_sig, h3, h4]
        row = prod["vp_rows"].get(ip, {})
        if row:
            mu = float(row["gauss_V"])
            sig = float(row["gauss_sig"])
            h3 = float(row["h3"])
            h4 = float(row["h4"])

            if np.all(np.isfinite([mu, sig, h3, h4])) and sig > 0:
                fits = np.array([mu, sig, h3, h4], dtype=float)

                vmin = float(np.min(v))
                vmax = float(np.max(v))
                if vmin == vmax:
                    pad = 3.0 * sig
                    vmin -= pad
                    vmax += pad

                vg = np.linspace(vmin, vmax, 1000)  # EXACT requirement
                fgrid, tag = balrogo_curve_from_fits(dyn, fits, vg)
                fgrid = normalize_curve(vg, fgrid)

                plt.plot(vg, fgrid, linewidth=2, label=f"BALRoGO curve ({tag})")
                plt.axvline(mu, linewidth=1.2, label="gauss_V (vprofile)")
                plt.text(
                    0.05,
                    0.85,
                    f"gauss_V={mu:.3g}\ngauss_sig={sig:.3g}\nh3={h3:.3g}, h4={h4:.3g}",
                    transform=plt.gca().transAxes,
                    va="top",
                )
            else:
                plt.text(
                    0.05,
                    0.85,
                    "vprofile fits non-finite\nor sigma <= 0\n(skipping BALRoGO overlay)",
                    transform=plt.gca().transAxes,
                    va="top",
                )

        label = {1: "LOS", 2: "POSR", 3: "POST"}[ip]
        plt.xlabel(f"Velocity (iproj={ip}, {label})")
        plt.ylabel("PDF")
        plt.title(f"Velocity check (iproj={ip}, hist_bins={args.hist_bins}, vg=1000)")
        plt.legend()
        fig.savefig(outdir / f"velocity_iproj{ip}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    print(f"\nDone. Outputs saved in: {outdir}")


if __name__ == "__main__":
    main()
