#!/usr/bin/env python3
"""scalefree.mock

Intrinsic (6D) mock generator.

This module generates a phase-space sample ``(x,y,z,vx,vy,vz)`` for a ScaleFree
model using the *intrinsic* (local) velocity-profile (VP) formalism.

Requirements imposed by the project request
-------------------------------------------
- The mass model is set by: ``gamma, beta, df, potential, s, t, q``.
- ``algorithm`` is forced to 3 (VP smoothing) while ``maxmom`` remains user-
  configurable.
- The intrinsic polar angle theta is binned into ``nbins`` bins over
  [0, 90] degrees (symmetry about the equatorial plane).
- For each occupied theta-bin we compute intrinsic VP parameters (Gaussian
  mean/sigma and Gauss-Hermite h3/h4) and sample velocities for each intrinsic
  component (Vr, Vtheta, Vphi) using BALRoGO's Sanders–Evans GH PDFs.
- The final output is returned in intrinsic Cartesian coordinates:
  ``(x, y, z, vx, vy, vz)``.

Notes
-----
- This is a Python-only mock generator: it does not rely on any Fortran-side
  mock sampler, but it *does* call the ScaleFree Fortran backend through
  :class:`scalefree.vmoments.ScaleFreeRunner` to obtain VP summaries.
- The position sampler uses the standard axisymmetric ScaleFree density
  ``rho ∝ (R^2 + z^2/q^2)^(-gamma/2)`` with symmetry axis z.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

def _coerce_seed(seed: Optional[int]) -> Optional[int]:
    """NumPy default_rng requires a non-negative seed; map negatives into uint32."""
    if seed is None:
        return None
    s = int(seed)
    return s % (2**32)


def _coerce_potential(potential: Any) -> Callable[[], int]:
    """ScaleFreeRunner accepts many potential encodings; keep backward-compat for int."""
    if callable(potential):
        return potential
    # Allow int-like potentials
    return lambda: int(potential)


# -----------------------------------------------------------------------------
# Density sampling (intrinsic coordinates)
# -----------------------------------------------------------------------------

def _sample_radius_powerlaw(
    rng: np.random.Generator, *, n: int, rin: float, rout: float, gamma: float
) -> np.ndarray:
    """Sample r with p(r) ∝ r^(2-gamma) (CDF ∝ r^(3-gamma))."""
    if rin <= 0 or rout <= rin:
        raise ValueError("Require 0 < rin < rout.")
    a = 3.0 - float(gamma)
    u = rng.random(n)
    if abs(a) > 1e-14:
        return (u * (rout**a - rin**a) + rin**a) ** (1.0 / a)
    # gamma == 3 -> log-uniform
    return rin * np.exp(u * np.log(rout / rin))


def _theta_accept_prob(theta: np.ndarray, *, q: float, gamma: float) -> np.ndarray:
    """Acceptance probability for theta under rho ∝ (R^2 + z^2/q^2)^(-gamma/2)."""
    ct = np.cos(theta)
    st = np.sin(theta)
    denom = st * st + (ct / q) * (ct / q)
    return denom ** (-0.5 * gamma)


def _sample_positions_scalefree(
    *,
    n: int,
    rin: float,
    rout: float,
    gamma: float,
    q: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return xyz (N,3) in the intrinsic model frame."""
    if n < 0:
        raise ValueError("n must be >= 0.")
    if q <= 0:
        raise ValueError("q must be > 0.")

    r = _sample_radius_powerlaw(rng, n=n, rin=rin, rout=rout, gamma=gamma)

    # theta rejection sampling
    theta = np.empty(n, dtype=float)
    filled = 0

    # conservative upper bound; ensures acceptance <= 1 even if q > 1
    prob_max = max(1.0, float(q) ** float(gamma))

    while filled < n:
        m = (n - filled) * 2 + 32
        ct = 2.0 * rng.random(m) - 1.0
        th = np.arccos(ct)

        acc = _theta_accept_prob(th, q=q, gamma=gamma) / prob_max
        keep = rng.random(m) <= acc
        th_keep = th[keep]

        take = min(th_keep.size, n - filled)
        if take > 0:
            theta[filled : filled + take] = th_keep[:take]
            filled += take

    phi = 2.0 * np.pi * rng.random(n)

    st = np.sin(theta)
    ct = np.cos(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    x = r * st * cp
    y = r * st * sp
    z = r * ct

    return np.column_stack([x, y, z])


# -----------------------------------------------------------------------------
# BALRoGO (Sanders–Evans) sampling
# -----------------------------------------------------------------------------

_BALROGO_DYN = None


def _import_balrogo_dynamics():
    """Import BALRoGO dynamics (or local fallback dynamics.py)."""
    global _BALROGO_DYN
    if _BALROGO_DYN is not None:
        return _BALROGO_DYN
    try:
        from balrogo import dynamics  # type: ignore

        _BALROGO_DYN = dynamics
    except Exception:
        import importlib

        _BALROGO_DYN = importlib.import_module("dynamics")
    return _BALROGO_DYN


def _sample_balrogo_gh(
    *,
    mean: float,
    sigma: float,
    h3: float,
    h4: float,
    n: int,
    rng: np.random.Generator,
    nsig: int = 10,
    debug: bool = False,
) -> np.ndarray:
    """Sample N values from BALRoGO's Gauss–Hermite PDF via inverse-CDF sampling."""
    if n <= 0:
        return np.empty((0,), dtype=float)

    if (not np.isfinite(mean)) or (not np.isfinite(sigma)) or (sigma <= 0):
        return rng.normal(loc=float(mean), scale=max(float(sigma), 1e-12), size=n)

    dyn = _import_balrogo_dynamics()

    mom_stats = np.array(
        [
            [float(mean), 0.0],
            [float(sigma), 0.0],
            [float(h3), 0.0],
            [float(h4), 0.0],
        ],
        dtype=float,
    )

    eps = np.zeros(n, dtype=float)

    # BALRoGO uses NumPy global RNG; bridge deterministically from our rng.
    np_seed = int(rng.integers(0, 2**32 - 1, dtype=np.uint32))
    np.random.seed(np_seed)

    samples = dyn.mom_sample_generator(mom_stats, eps=eps, nsig=int(nsig), debug=bool(debug))

    if samples is None:
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    samples = np.asarray(samples, dtype=float)
    if samples.shape != (n,) or not np.all(np.isfinite(samples)):
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    return samples


# -----------------------------------------------------------------------------
# Intrinsic VP extraction (ScaleFreeRunner.vprofile output)
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GHParams:
    mean: float
    sig: float
    h3: float
    h4: float


def _extract_intrinsic_moments(res) -> Tuple[float, float, float, float]:
    """Return (vphi, vr2, vth2, vphi2) from intrinsic moment blocks."""
    blocks = getattr(res, "blocks", {}) or {}

    blk = blocks.get("intrinsic_point")
    if blk is None:
        blk = blocks.get("intrinsic_shell_average")

    if not isinstance(blk, dict):
        raise RuntimeError("vprofile did not return an intrinsic moments block.")

    cols = list(blk.get("columns", []))
    data = np.asarray(blk.get("data"))

    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.ndim != 2 or data.shape[0] < 1:
        raise RuntimeError("Intrinsic moments block has unexpected shape.")

    row = data[0]

    def get(name: str) -> float:
        if name not in cols:
            raise KeyError(f"Missing column '{name}' in intrinsic moments block. cols={cols}")
        return float(row[cols.index(name)])

    vphi = get("vphi")
    vr2 = get("vr2")
    vth2 = get("vth2")
    vphi2 = get("vphi2")
    return vphi, vr2, vth2, vphi2


def _extract_h3_h4(res, icomp: int) -> Tuple[float, float]:
    """Extract (h3, h4) for intrinsic component icomp from vp_intrinsic."""
    blocks = getattr(res, "blocks", {}) or {}
    vp = blocks.get("vp_intrinsic")
    if vp is None:
        raise RuntimeError("vprofile did not return a 'vp_intrinsic' summary block.")

    cols = list(vp.get("columns", []))
    data = np.asarray(vp.get("data"))

    def idx(name: str) -> int:
        try:
            return cols.index(name)
        except ValueError as e:
            raise KeyError(f"Missing column '{name}' in vp_intrinsic summary. cols={cols}") from e

    i_comp = idx("icomp")
    i_h3 = idx("h3")
    i_h4 = idx("h4")

    hit = np.where(data[:, i_comp].astype(int, copy=False) == int(icomp))[0]
    if hit.size != 1:
        raise RuntimeError(
            f"Expected exactly one vp_intrinsic row for icomp={icomp}; found {hit.size}."
        )

    r = data[hit[0]]
    return float(r[i_h3]), float(r[i_h4])


def _scalefree_intrinsic_params_for_theta_bin(
    *,
    runner,
    potential: Any,
    gamma: float,
    q: float,
    df: int,
    beta: float,
    s: float,
    t: float,
    theta_deg: float,
    # numerical
    xi: float,
    integration: int,
    ngl_or_eps: float,
    maxmom: int,
    vp_smooth_eps: float,
    verbose_vp: int,
) -> Dict[int, GHParams]:
    """Compute intrinsic (Vr,Vtheta,Vphi) GH parameters for a single theta-bin."""

    # vprofile still expects inclination/xi/theta inputs.
    # For intrinsic kinematics, inclination is irrelevant; keep a benign value.
    res = runner.vprofile(
        potential=potential,
        gamma=float(gamma),
        q=float(q),
        df=int(df),
        beta=float(beta),
        s=float(s),
        t=float(t),
        inclination=90.0,
        xi=float(xi),
        theta=float(theta_deg),
        integration=int(integration),
        ngl_or_eps=float(ngl_or_eps),
        algorithm=3,
        vp_smooth_eps=float(vp_smooth_eps),
        maxmom=int(maxmom),
        average=False,
        kinematics="intrinsic",
        usevp=True,
        verbose_vp=int(verbose_vp),
        output_path=None,
        debug_prompts=False,
    )

    vphi, vr2, vth2, vphi2 = _extract_intrinsic_moments(res)

    # Means and dispersions to feed BALRoGO.
    mu_r = 0.0
    mu_th = 0.0
    mu_phi = float(vphi)

    sig_r = float(np.sqrt(max(vr2, 0.0)))
    sig_th = float(np.sqrt(max(vth2, 0.0)))
    sig_phi = float(np.sqrt(max(vphi2 - vphi * vphi, 0.0)))

    h3_r, h4_r = _extract_h3_h4(res, 1)
    h3_th, h4_th = _extract_h3_h4(res, 2)
    h3_phi, h4_phi = _extract_h3_h4(res, 3)

    return {
        1: GHParams(mu_r, sig_r, h3_r, h4_r),
        2: GHParams(mu_th, sig_th, h3_th, h4_th),
        3: GHParams(mu_phi, sig_phi, h3_phi, h4_phi),
    }


# -----------------------------------------------------------------------------
# Spherical <-> Cartesian velocity transforms
# -----------------------------------------------------------------------------

def _sph_from_xyz(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.sqrt(x * x + y * y + z * z)
    # theta in [0, pi]
    theta = np.arccos(np.clip(z / np.where(r > 0, r, 1.0), -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def _vxyz_from_sph(
    vr: np.ndarray,
    vtheta: np.ndarray,
    vphi: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    st = np.sin(theta)
    ct = np.cos(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    # unit vectors
    # e_r    = (st*cp, st*sp, ct)
    # e_th   = (ct*cp, ct*sp, -st)
    # e_phi  = (-sp, cp, 0)

    vx = vr * st * cp + vtheta * ct * cp - vphi * sp
    vy = vr * st * sp + vtheta * ct * sp + vphi * cp
    vz = vr * ct - vtheta * st
    return vx, vy, vz


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def mock(
    *,
    # ScaleFree model parameters
    potential: Any = 1,
    gamma: float = 4.0,
    q: float = 1.0,
    df: int = 1,
    beta: float = 0.0,
    s: float = 0.5,
    t: float = 0.0,
    # sampling controls
    nsamples: int = 1000,
    seed: Optional[int] = None,
    rin: float = 1.0,
    rout: float = 1000.0,
    # theta binning
    nbins: int = 36,
    # VP settings (algorithm is forced to 3)
    maxmom: int = 20,
    integration: int = 1,
    ngl_or_eps: float = 0.0,
    vp_smooth_eps: float = 0.0,
    xi: float = 0.0,
    verbose_vp: int = 0,
    # BALRoGO settings
    nsig: int = 10,
    debug: bool = False,
) -> np.ndarray:
    """Generate a 6D intrinsic mock.

    Returns
    -------
    ndarray
        Shape (nsamples, 6) with columns (x, y, z, vx, vy, vz).
    """
    n = int(nsamples)
    if n < 0:
        raise ValueError("nsamples must be >= 0.")
    nb = int(nbins)
    if nb < 1:
        raise ValueError("nbins must be >= 1.")

    rng = np.random.default_rng(_coerce_seed(seed))
    potential = _coerce_potential(potential)

    # 1) sample positions
    xyz = _sample_positions_scalefree(n=n, rin=rin, rout=rout, gamma=gamma, q=q, rng=rng)
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

    # 2) intrinsic spherical angles and symmetric theta in [0, pi/2]
    _r, theta, phi = _sph_from_xyz(x, y, z)
    theta_sym = np.minimum(theta, np.pi - theta)
    theta_deg = np.degrees(theta_sym)

    # 3) bin theta in [0, 90]
    edges = np.linspace(0.0, 90.0, nb + 1)
    bin_id = np.digitize(theta_deg, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, nb - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 4) compute intrinsic VP params per occupied bin
    from scalefree import ScaleFreeRunner  # local import avoids circular issues

    runner = ScaleFreeRunner()

    gh_by_bin: Dict[int, Dict[int, GHParams]] = {}
    occupied = np.unique(bin_id)
    for b in occupied:
        gh_by_bin[int(b)] = _scalefree_intrinsic_params_for_theta_bin(
            runner=runner,
            potential=potential,
            gamma=gamma,
            q=q,
            df=df,
            beta=beta,
            s=s,
            t=t,
            theta_deg=float(centers[int(b)]),
            xi=xi,
            integration=integration,
            ngl_or_eps=ngl_or_eps,
            maxmom=maxmom,
            vp_smooth_eps=vp_smooth_eps,
            verbose_vp=verbose_vp,
        )

    # 5) sample (vr, vtheta, vphi) per bin
    vr = np.empty(n, dtype=float)
    vth = np.empty(n, dtype=float)
    vph = np.empty(n, dtype=float)

    for b in occupied:
        b = int(b)
        sel = bin_id == b
        m = int(np.sum(sel))
        if m == 0:
            continue

        gh = gh_by_bin[b]
        p1, p2, p3 = gh[1], gh[2], gh[3]

        vr[sel] = _sample_balrogo_gh(
            mean=p1.mean,
            sigma=max(p1.sig, 1e-12),
            h3=p1.h3,
            h4=p1.h4,
            n=m,
            rng=rng,
            nsig=nsig,
        )
        vth[sel] = _sample_balrogo_gh(
            mean=p2.mean,
            sigma=max(p2.sig, 1e-12),
            h3=p2.h3,
            h4=p2.h4,
            n=m,
            rng=rng,
            nsig=nsig,
        )
        vph[sel] = _sample_balrogo_gh(
            mean=p3.mean,
            sigma=max(p3.sig, 1e-12),
            h3=p3.h3,
            h4=p3.h4,
            n=m,
            rng=rng,
            nsig=nsig,
        )

    # 6) (vr,vtheta,vphi) -> (vx,vy,vz)
    vx, vy, vz = _vxyz_from_sph(vr, vth, vph, theta=theta, phi=phi)

    out = np.column_stack([x, y, z, vx, vy, vz])

    if debug:
        counts = np.bincount(bin_id, minlength=nb)
        nz = np.count_nonzero(counts)
        print(f"[mock intrinsic] nsamples={n} nbins={nb} occupied_bins={nz} (theta in 0..90 deg)")
        if nz:
            print(f"[mock intrinsic] min_bin_count={counts[counts>0].min()} max_bin_count={counts.max()}")

    return out
