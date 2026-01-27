#!/usr/bin/env python3
"""scripts/check_mock_end2end.py

Local end-to-end validation of the *intrinsic* scalefree mock generator.

What is tested
--------------
- Density: recovered shell volume density follows r^{-gamma} (up to scaling).
- Flattening: the (x,z) distribution matches the requested intrinsic axis
ratio q.
- Velocities: the intrinsic velocity components (vr, vtheta, vphi) produced by
  the mock are consistent with the *angle-averaged* intrinsic VP
  information from
  ScaleFree (vprofile with average=True), even though the mock itself uses
  average=False with theta-binning.

Notes
-----
- This script assumes the mock returns *Cartesian* phase-space coordinates.
- No intermediate transformation to observed (LOS/POSr/POSt) coordinates
is used.
- For the velocity goodness-of-fit curves we use the *average=True*
intrinsic VP
  gaussian-fit parameters (gauss_V, gauss_sig) and GH moments (h3,h4).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any, Dict, Tuple

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


def shell_volume_density(
    r: np.ndarray, nbins: int = 25
) -> Tuple[np.ndarray, np.ndarray]:
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
# Spherical <-> Cartesian velocity transforms
# -----------------------------


def sph_vel_from_xyzv(xyzv: np.ndarray) -> Dict[str, np.ndarray]:
    """Return intrinsic spherical velocity components from (x,y,z,vx,vy,vz)."""
    x, y, z, vx, vy, vz = xyzv.T
    r = np.sqrt(x * x + y * y + z * z)
    r = np.where(r > 0, r, np.nan)

    theta = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi = np.arctan2(y, x)

    st = np.sin(theta)
    ct = np.cos(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    # basis vectors
    erx, ery, erz = st * cp, st * sp, ct
    etx, ety, etz = ct * cp, ct * sp, -st
    epx, epy, epz = -sp, cp, 0.0

    vr = vx * erx + vy * ery + vz * erz
    vtheta = vx * etx + vy * ety + vz * etz
    vphi = vx * epx + vy * epy + vz * epz

    return {
        "vr": vr,
        "vtheta": vtheta,
        "vphi": vphi,
        "theta": theta,
        "phi": phi,
        "r": r,
    }


# -----------------------------
# vprofile extraction (intrinsic)
# -----------------------------


def extract_intrinsic_products(res) -> Dict[str, Any]:
    """Extract intrinsic VP fits and VP table (if available)."""
    blocks = getattr(res, "blocks", {}) or {}
    out: Dict[str, Any] = {"vp_rows": {}, "vp_table": {}, "moments": {}}

    # Moments:
    # average=True -> intrinsic_shell_average;
    # average=False -> intrinsic_point
    for key in ("intrinsic_shell_average", "intrinsic_point"):
        blk = blocks.get(key)
        if isinstance(blk, dict):
            cols = list(blk.get("columns", []))
            data = blk.get("data")
            if (
                cols
                and isinstance(
                    data,
                    np.ndarray,
                )
                and data.ndim == 2
                and data.size
            ):
                row = data[0]
                for name in ("rho", "vphi", "vr2", "vth2", "vphi2", "beta"):
                    if name in cols:
                        out["moments"][name] = float(row[cols.index(name)])
            break

    # VP summary: vp_intrinsic
    vp = blocks.get("vp_intrinsic")
    if isinstance(vp, dict):
        cols = list(vp.get("columns", []))
        data = vp.get("data")
        if (
            cols
            and isinstance(data, np.ndarray)
            and data.ndim == 2
            and data.size
            and "icomp" in cols
        ):
            i_ic = cols.index("icomp")

            def _get(col: str, row: np.ndarray) -> float:
                return (
                    float(
                        row[cols.index(col)],
                    )
                    if col in cols
                    else float("nan")
                )

            for row in data:
                ic = int(row[i_ic])
                out["vp_rows"][ic] = {
                    "gauss_V": _get("gauss_V", row),
                    "gauss_sig": _get("gauss_sig", row),
                    "h3": _get("h3", row),
                    "h4": _get("h4", row),
                    "true_V": _get("true_V", row),
                    "true_sig": _get("true_sig", row),
                }

    # vp_table: for intrinsic VPs, parser stores these under vp_table_intrinsic
    vpt = blocks.get("vp_table_intrinsic")
    if isinstance(vpt, dict):
        for icomp, tbl in vpt.items():
            if not isinstance(tbl, dict):
                continue
            arr = tbl.get("data")
            if (
                isinstance(arr, np.ndarray)
                and arr.ndim == 2
                and arr.size
                and arr.shape[1] >= 2
            ):
                out["vp_table"][int(icomp)] = (
                    arr[:, 0].astype(float),
                    arr[:, 1].astype(float),
                )

    return out


# -----------------------------
# BALRoGO curve helpers
# -----------------------------


def balrogo_curve_from_fits(
    dyn, fits: np.ndarray, vg: np.ndarray
) -> Tuple[np.ndarray, str]:
    """Return a PDF curve from BALRoGO given (mean, sigma, h3, h4)."""
    mu, sig, h3, h4 = map(float, fits)
    ex = np.zeros_like(vg, dtype=float)

    ret = dyn.mom_likelihood_func(fits, vg, ex, mode="curve")

    if (
        isinstance(
            ret,
            np.ndarray,
        )
        and ret.ndim == 1
        and ret.shape[0] == vg.shape[0]
    ):
        return ret.astype(float), "mom_likelihood_func(curve)"

    if np.isscalar(ret):
        if h4 >= 0:
            fgrid = dyn.laplace_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return (
                np.asarray(
                    fgrid,
                    dtype=float,
                ),
                "laplace_kernel_pdf(fallback)",
            )
        else:
            fgrid = dyn.uniform_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return (
                np.asarray(
                    fgrid,
                    dtype=float,
                ),
                "uniform_kernel_pdf(fallback)",
            )

    raise RuntimeError(
        "Unexpected return"
        + " type from mom_likelihood_func:"
        + f" {type(ret).__name__}",
    )


def normalize_curve(vg: np.ndarray, fgrid: np.ndarray) -> np.ndarray:
    fgrid = np.asarray(fgrid, dtype=float)
    if fgrid.ndim != 1 or fgrid.shape[0] != vg.shape[0]:
        return fgrid
    area = (
        np.trapezoid(fgrid, vg)
        if hasattr(np, "trapezoid")
        else np.trapz(
            fgrid,
            vg,
        )
    )
    if np.isfinite(area) and area > 0:
        return fgrid / area
    return fgrid


# -----------------------------
# Main
# -----------------------------


def main() -> None:
    p = argparse.ArgumentParser(
        description="Local scalefree"
        + "intrinsic mock check"
        + " (6D Cartesian output)."
    )
    p.add_argument("--n_dens", type=int, default=8000)
    p.add_argument("--n_vel", type=int, default=12000)
    p.add_argument("--nbins", type=int, default=90)
    p.add_argument("--hist_bins", type=int, default=80)
    p.add_argument("--outdir", type=str, default="mock_check_outputs")
    p.add_argument("--seed", type=int, default=42)
    args, _unknown = p.parse_known_args()

    repo_root = Path(__file__).resolve().parents[1]
    _force_repo_import(repo_root)

    import scalefree  # noqa: E402
    from scalefree import ScaleFreeRunner, mock  # noqa: E402

    dyn = _import_balrogo_dynamics()

    print("Imported scalefree from:", Path(scalefree.__file__).resolve())

    outdir = (Path(__file__).resolve().parent / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Model parameters
    # ------------------------------------------------------------------
    model = dict(
        potential=2,
        gamma=2.0,
        q=0.9,
        df=2,
        beta=-0.1,
        s=0.2,
        t=0.0,
    )

    # Prefer the same backend path as scripts/quick_user_run.py
    exe = (repo_root / "fortran_src" / "scalefree_intrvp.e").resolve()

    # ------------------------------------------------------------------
    # A) vprofile reference: intrinsic, average=True
    # ------------------------------------------------------------------
    runner = ScaleFreeRunner(exe_path=exe)
    res = runner.vprofile(
        **model,
        inclination=90.0,
        xi=0.0,
        theta=0.0,
        integration=1,
        ngl_or_eps=0,
        algorithm=3,
        maxmom=30,
        average=True,
        kinematics="intrinsic",
        usevp=True,
        verbose_vp=0,
        output_path=None,
        timeout_s=300,
        debug_prompts=False,
    )
    prod = extract_intrinsic_products(res)

    mom = prod.get("moments", {})
    print("\n=== intrinsic moments (average=True) ===")
    for k in ("rho", "vphi", "vr2", "vth2", "vphi2", "beta"):
        if k in mom:
            print(f"{k:>6} = {mom[k]:.6g}")

    print("\n=== vp_intrinsic summary rows (average=True) ===")
    for ic in (1, 2, 3):
        row = prod["vp_rows"].get(ic, {})
        if row:
            print(
                f"\nicomp={ic}:"
                + f" gauss_V={row['gauss_V']:.6g},"
                + f" gauss_sig={row['gauss_sig']:.6g},"
                + f" h3={row['h3']:.6g},"
                + f" h4={row['h4']:.6g}"
            )
            print(f" true_V={row['true_V']:.6g}," + f" true_sig={row['true_sig']:.6g}")
        else:
            print(f"icomp={ic}: (missing)")

    # Reference mean/sigma from gaussian VP fits
    mu_ref = {
        ic: float(
            prod["vp_rows"].get(ic, {}).get("gauss_V", 0.0),
        )
        for ic in (1, 2, 3)
    }
    sig_ref = {
        ic: float(prod["vp_rows"].get(ic, {}).get("gauss_sig", float("nan")))
        for ic in (1, 2, 3)
    }

    # ------------------------------------------------------------------
    # B) Density sanity check (positions)
    # ------------------------------------------------------------------
    xyzv = mock(
        **model,
        nsamples=int(args.n_dens),
        nbins=int(args.nbins),
        seed=int(args.seed),
        rin=0.5,
        rout=50.0,
        maxmom=8,
        debug=False,
        exe_path=exe,
    )

    r3d = np.sqrt(xyzv[:, 0] ** 2 + xyzv[:, 1] ** 2 + xyzv[:, 2] ** 2)
    rmid, rho_hat = shell_volume_density(r3d, nbins=25)
    rho_th = volume_density_powerlaw_shape(model["gamma"], rmid)
    k0 = len(rmid) // 2
    rho_th_scaled = rho_th * (rho_hat[k0] / rho_th[k0])

    fig = plt.figure()
    plt.loglog(
        rmid,
        rho_hat,
        marker="o",
        linestyle="none",
        label="Mock: shell rho(r)",
    )
    plt.loglog(
        rmid,
        rho_th_scaled,
        linestyle="-",
        label=r"Analytic: $r^{-\gamma}$ (scaled)",
    )
    plt.xlabel("3D radius r")
    plt.ylabel("Volume density (arb.)")
    plt.title(f"Density check (gamma={model['gamma']})")
    plt.legend()
    fig.savefig(outdir / "density_q1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # C) Flattening sanity check
    # ------------------------------------------------------------------
    model_flat = dict(model)
    model_flat["q"] = 0.408
    xyzv_flat = mock(
        **model_flat,
        nsamples=min(8000, int(args.n_dens)),
        nbins=int(args.nbins),
        seed=int(args.seed) + 1,
        rin=0.5,
        rout=50.0,
        maxmom=30,
        debug=False,
        exe_path=exe,
    )

    X = xyzv_flat[:, 0]
    Z = xyzv_flat[:, 2]
    q0 = float(model_flat["q"])

    Re = np.sqrt(X * X + (Z / q0) ** 2)
    Re0 = float(np.nanpercentile(Re, 60))

    ang = np.linspace(0, 2 * np.pi, 400)
    ex_ = Re0 * np.cos(ang)
    ez_ = q0 * Re0 * np.sin(ang)

    fig = plt.figure()
    idx = np.random.default_rng(int(args.seed)).choice(
        len(X), size=min(6000, len(X)), replace=False
    )
    plt.scatter(X[idx], Z[idx], s=2, alpha=0.25, label="Mock points (x-z)")
    plt.plot(ex_, ez_, linewidth=2, label=f"Ellipse overlay (q={q0})")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"Flattening sanity check (q={q0})")
    plt.legend()
    fig.savefig(
        outdir / "flattening_overlay.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # D) Velocity check in intrinsic spherical components
    # ------------------------------------------------------------------
    xyzv_vel = mock(
        **model,
        nsamples=int(args.n_vel),
        nbins=int(args.nbins),
        seed=int(args.seed) + 2,
        rin=1.0,
        rout=2.0,
        maxmom=30,
        debug=False,
        exe_path=exe,
    )

    sv = sph_vel_from_xyzv(xyzv_vel)

    comp_map = {1: "vr", 2: "vtheta", 3: "vphi"}
    comp_label = {1: "Vr", 2: "Vtheta", 3: "Vphi"}

    for ic in (1, 2, 3):
        name = comp_map[ic]
        v = np.asarray(sv[name], dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue

        fig = plt.figure()
        plt.hist(
            v,
            bins=int(args.hist_bins),
            density=True,
            alpha=0.35,
            label="Mock histogram",
        )

        # Optional vp_table overlay
        if ic in prod["vp_table"]:
            vg_tab, vpg_tab = prod["vp_table"][ic]
            vpg_tab = normalize_curve(vg_tab, vpg_tab)
            plt.plot(
                vg_tab,
                vpg_tab,
                linewidth=2,
                label="ScaleFree vp_table",
            )

        # BALRoGO analytic curve from average=True gaussian+GH fits
        row = prod["vp_rows"].get(ic, {})
        if row:
            mu = float(mu_ref[ic])
            sig = float(sig_ref[ic])
            h3 = float(row.get("h3", float("nan")))
            h4 = float(row.get("h4", float("nan")))

            if np.all(np.isfinite([mu, sig, h3, h4])) and sig > 0:
                fits = np.array([mu, sig, h3, h4], dtype=float)

                vmin = float(np.min(v))
                vmax = float(np.max(v))
                if vmin == vmax:
                    pad = 3.0 * sig
                    vmin -= pad
                    vmax += pad

                vg = np.linspace(vmin, vmax, 1000)
                fgrid, tag = balrogo_curve_from_fits(dyn, fits, vg)
                fgrid = normalize_curve(vg, fgrid)

                plt.plot(vg, fgrid, linewidth=2, label="BALRoGO curve")
                plt.axvline(mu, linewidth=1.0, linestyle="--", label="gauss_V")

        plt.xlabel(f"{comp_label[ic]} (icomp={ic})")
        plt.ylabel("PDF")
        plt.xlim(-3, 3)
        plt.title(
            f"Intrinsic velocity check: {comp_label[ic]}"
            + " (theta-binned mock vs avg-fit curve)",
        )
        plt.legend(loc=1)
        fig.savefig(
            outdir / f"velocity_icomp{ic}.png",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(fig)

    print(f"\nDone. Outputs saved in: {outdir}")


if __name__ == "__main__":
    main()
