#!/usr/bin/env python3
"""
Local end-to-end validation of scalefree mock generation (mock-focused),
with apples-to-apples velocity geometry.

Design (per user intent)
------------------------
- The mock is generated star-by-star (average=False): ψ-dependent projected distributions
  are sampled and then rotated back into intrinsic Cartesian (x,y,z,vx,vy,vz) using the
  provided inclination/xi.
- The analytic overlay curve is sourced from vprofile with average=True (sky-averaged
  GH parameters), and compared against the histogram of the star-by-star mock sample.

Critical correctness condition
------------------------------
When building vLOS/vPOSr/vPOSt from (x,y,z,vx,vy,vz) we MUST use the same sky-frame
conventions (inclination, xi) as scalefree.mock / VPOS.md.
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
# Sky-frame rotation utilities
# (match scalefree.mock conventions)
# -----------------------------
def rotate_to_sky_xyz(x, y, z, inc_deg):
    i = np.deg2rad(float(inc_deg))
    ci, si = np.cos(i), np.sin(i)
    xp = y
    yp = -x * ci + z * si
    zp = x * si + z * ci
    return xp, yp, zp


def rotate_to_sky_v(vx, vy, vz, inc_deg):
    # model -> sky, consistent with the above
    i = np.deg2rad(float(inc_deg))
    ci, si = np.cos(i), np.sin(i)
    vxp = vy
    vyp = -vx * ci + vz * si
    vzp = vx * si + vz * ci
    return vxp, vyp, vzp


def rotate_sky_plane_by_xi(
    x: np.ndarray, y: np.ndarray, xi_deg: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Additional rotation in sky plane about z' (LOS) by angle xi.
    If xi=0, identity.
    """
    a = np.deg2rad(float(xi_deg))
    ca = np.cos(a)
    sa = np.sin(a)
    xr = x * ca - y * sa
    yr = x * sa + y * ca
    return xr, yr


def projected_vel_components_from_xyzv(
    xyzv: np.ndarray,
    *,
    inclination_deg: float,
    xi_deg: float,
) -> Dict[int, np.ndarray]:
    """
    Compute projected velocity components in the sky frame:
      iproj=1: vLOS  = vz'
      iproj=2: vPOSr = (x' vx' + y' vy') / R
      iproj=3: vPOSt = (y' vx' - x' vy') / R

    Uses the same inclination convention as scalefree.mock.
    """
    x, y, z, vx, vy, vz = xyzv.T

    x_sky, y_sky, _z_sky = rotate_to_sky_xyz(x, y, z, inclination_deg)
    vx_sky, vy_sky, vz_sky = rotate_to_sky_v(vx, vy, vz, inclination_deg)

    # apply xi rotation in sky plane (positions and velocities)
    x2, y2 = rotate_sky_plane_by_xi(x_sky, y_sky, xi_deg)
    vx2, vy2 = rotate_sky_plane_by_xi(vx_sky, vy_sky, xi_deg)

    R = np.sqrt(x2 * x2 + y2 * y2)
    m = R > 0

    vlos = vz_sky
    vposr = np.full_like(vlos, np.nan, dtype=float)
    vpost = np.full_like(vlos, np.nan, dtype=float)

    vposr[m] = (x2[m] * vx2[m] + y2[m] * vy2[m]) / R[m]
    vpost[m] = (y2[m] * vx2[m] - x2[m] * vy2[m]) / R[m]

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
        if (
            cols
            and isinstance(data, np.ndarray)
            and data.ndim == 2
            and data.size
            and "iproj" in cols
        ):
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
# BALRoGO curve helpers
# -----------------------------
def balrogo_curve_from_fits(
    dyn, fits: np.ndarray, vg: np.ndarray
) -> Tuple[np.ndarray, str]:
    """
    Try mom_likelihood_func(..., mode='curve'); if it returns scalar (often inf),
    fall back to kernel pdf.
    """
    mu, sig, h3, h4 = map(float, fits)
    ex = np.zeros_like(vg, dtype=float)

    ret = dyn.mom_likelihood_func(fits, vg, ex, mode="curve")

    if isinstance(ret, np.ndarray) and ret.ndim == 1 and ret.shape[0] == vg.shape[0]:
        return ret.astype(float), "mom_likelihood_func(curve)"

    if np.isscalar(ret):
        if h4 >= 0:
            fgrid = dyn.laplace_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "laplace_kernel_pdf(fallback)"
        else:
            fgrid = dyn.uniform_kernel_pdf(vg, ex, mu, sig, h3, h4)
            return np.asarray(fgrid, dtype=float), "uniform_kernel_pdf(fallback)"

    raise RuntimeError(
        f"Unexpected return type from mom_likelihood_func: {type(ret).__name__}"
    )


def normalize_curve(vg: np.ndarray, fgrid: np.ndarray) -> np.ndarray:
    fgrid = np.asarray(fgrid, dtype=float)
    if fgrid.ndim != 1 or fgrid.shape[0] != vg.shape[0]:
        return fgrid
    area = np.trapezoid(fgrid, vg) if hasattr(np, "trapezoid") else np.trapz(fgrid, vg)
    if np.isfinite(area) and area > 0:
        return fgrid / area
    return fgrid


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    p = argparse.ArgumentParser(
        description="Local scalefree mock check with correct projected-velocity geometry."
    )
    p.add_argument("--n_dens", type=int, default=8000)
    p.add_argument("--n_vel", type=int, default=8000)
    p.add_argument("--nbins", type=int, default=36)
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

    # Keep model+view consistent across vprofile and mock
    model = dict(
        potential=1,
        gamma=2.0,
        q=1.0,
        df=1,
        beta=0.0,
        s=0.0,
        t=0.0,
    )
    view = dict(inclination=57.1, xi=0.0, theta=0.0)

    # A) vprofile: average=True gives sky-averaged fits (the curve you want)
    runner = ScaleFreeRunner()
    res = runner.vprofile(
        **model,
        **view,
        integration=1,
        ngl_or_eps=0,
        algorithm=3,
        maxmom=20,
        average=True,
        usevp=True,
        verbose_vp=0,
        output_path=None,
        timeout_s=300,
        debug_prompts=False,
    )
    prod = extract_vmoments_products(res)

    print("\n=== vprofile summary rows (average=True; sky-averaged fits) ===")
    for ip in (1, 2, 3):
        row = prod["vp_rows"].get(ip, {})
        if row:
            print(
                f"iproj={ip}: gauss_V={row['gauss_V']:.6g}, gauss_sig={row['gauss_sig']:.6g}, "
                f"h3={row['h3']:.6g}, h4={row['h4']:.6g}"
            )
        else:
            print(f"iproj={ip}: (missing)")

    # B) Density sanity check (positions)
    xyzv = mock(
        **model,
        inclination=view["inclination"],
        xi=view["xi"],
        nsamples=int(args.n_dens),
        seed=int(args.seed),
        rin=0.5,
        rout=50.0,
        algorithm=3,
        maxmom=20,
        average=False,  # star-by-star
        nbins=int(args.nbins),
        debug=True,
    )

    r3d = np.sqrt(xyzv[:, 0] ** 2 + xyzv[:, 1] ** 2 + xyzv[:, 2] ** 2)
    rmid, rho_hat = shell_volume_density(r3d, nbins=25)
    rho_th = volume_density_powerlaw_shape(model["gamma"], rmid)
    k0 = len(rmid) // 2
    rho_th_scaled = rho_th * (rho_hat[k0] / rho_th[k0])

    fig = plt.figure()
    plt.loglog(rmid, rho_hat, marker="o", linestyle="none", label="Mock: shell rho(r)")
    plt.loglog(
        rmid, rho_th_scaled, linestyle="-", label=r"Analytic: $r^{-\gamma}$ (scaled)"
    )
    plt.xlabel("3D radius r")
    plt.ylabel("Volume density (arb.)")
    plt.title(f"Density check (gamma={model['gamma']})")
    plt.legend()
    fig.savefig(outdir / "density_q1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # C) Flattening sanity check
    model_flat = dict(model)
    model_flat["q"] = 0.408
    xyzv_flat = mock(
        **model_flat,
        inclination=view["inclination"],
        xi=view["xi"],
        nsamples=min(8000, int(args.n_dens)),
        seed=int(args.seed) + 1,
        rin=0.5,
        rout=50.0,
        algorithm=3,
        maxmom=20,
        average=False,  # star-by-star
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
    idx = np.random.default_rng(int(args.seed)).choice(
        len(X), size=min(6000, len(X)), replace=False
    )
    plt.scatter(X[idx], Z[idx], s=2, alpha=0.25, label="Mock points (x-z)")
    plt.plot(ex, ez, linewidth=2, label=f"Ellipse overlay (q={q0})")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"Flattening sanity check (q={q0})")
    plt.legend()
    fig.savefig(outdir / "flattening_overlay.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # D) Velocity check:
    # mock: average=False (star-by-star). Projection back uses the same inclination/xi.
    xyzv_vel = mock(
        **model,
        inclination=view["inclination"],
        xi=view["xi"],
        nsamples=int(args.n_vel),
        seed=int(args.seed) + 2,
        rin=1.0,
        rout=2.0,
        algorithm=3,
        maxmom=20,
        average=False,  # IMPORTANT: star-by-star
        nbins=int(args.nbins),
        debug=True,
    )

    vproj = projected_vel_components_from_xyzv(
        xyzv_vel,
        inclination_deg=view["inclination"],
        xi_deg=view["xi"],
    )

    for ip in (1, 2, 3):
        v = np.asarray(vproj[ip], dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue

        fig = plt.figure()
        plt.hist(
            v,
            bins=int(args.hist_bins),
            density=True,
            alpha=0.35,
            label="Mock histogram (star-by-star)",
        )

        # Optional vp_table overlay
        if ip in prod["vp_table"]:
            vg_tab, vpg_tab = prod["vp_table"][ip]
            vpg_tab = normalize_curve(vg_tab, vpg_tab)
            plt.plot(
                vg_tab, vpg_tab, linewidth=2, label="ScaleFree vp_table (normalized)"
            )

        # Analytic overlay from sky-averaged vprofile fits
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

                vg = np.linspace(vmin, vmax, 1000)
                fgrid, tag = balrogo_curve_from_fits(dyn, fits, vg)
                fgrid = normalize_curve(vg, fgrid)

                plt.plot(
                    vg, fgrid, linewidth=2, label=f"BALRoGO curve (avg fits; {tag})"
                )
                # plt.axvline(mu, linewidth=1.2, label="gauss_V (vprofile avg)")
                # plt.text(
                #     0.05,
                #     0.85,
                #     f"gauss_V={mu:.3g}\ngauss_sig={sig:.3g}\nh3={h3:.3g}, h4={h4:.3g}",
                #     transform=plt.gca().transAxes,
                #     va="top",
                # )

        label = {1: "LOS", 2: "POSr", 3: "POSt"}[ip]
        plt.xlabel(f"Velocity (iproj={ip}, {label})")
        plt.ylabel("PDF")
        plt.title(f"Velocity check (iproj={ip}): star-by-star mock vs avg-fit curve")
        plt.legend()
        fig.savefig(outdir / f"velocity_iproj{ip}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    print(f"\nDone. Outputs saved in: {outdir}")


if __name__ == "__main__":
    main()
