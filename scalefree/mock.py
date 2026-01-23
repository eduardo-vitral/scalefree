#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""scalefree.mock

Intrinsic-θ dependent mock generator (Python-only).

This module generates a synthetic 6D sample in the *model* Cartesian frame:

    (x, y, z, vx, vy, vz)

The generator follows the requested workflow:

1) Sample ``nsamples`` positions from the ScaleFree density distribution.
2) Divide the intrinsic polar angle θ into ``theta_bins`` over [0°, 90°]
   (equatorial symmetry is used: θ -> min(θ, π-θ)).
3) For each occupied θ-bin, run the ScaleFree backend once in *intrinsic*
   mode to obtain:
     - local intrinsic moments (rho, <vphi>, <vr^2>, <vtheta^2>, <vphi^2>)
     - intrinsic VP Gauss–Hermite parameters (h3, h4) for each component
       (icomp = 1: vr, 2: vtheta, 3: vphi)
4) For each star in a θ-bin, sample (vr, vtheta, vphi) using BALRoGO's
   Sanders–Evans-style sampler, with:
     - mean and dispersion anchored to the intrinsic moment block
     - (h3, h4) taken from the intrinsic VP summary
5) Convert (r, θ, φ) and (vr, vtheta, vphi) to Cartesian (x, y, z, vx, vy, vz).

Notes
-----
- The ScaleFree backend is still required (via :class:`scalefree.ScaleFreeRunner`).
- BALRoGO is used for sampling if available; otherwise we fall back to
  importing a local ``dynamics`` module (mirrors prior behavior).
- We *force* VP reconstruction ``algorithm=3`` (per request). ``maxmom`` is
  left free.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

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


def _coerce_potential(potential):
    """vprofile expects a callable; allow int for backward compatibility."""
    if callable(potential):
        return potential
    return lambda: int(potential)


# -----------------------------------------------------------------------------
# Density sampling
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
    """Acceptance ∝ (sin²θ + cos²θ/q²)^(-γ/2) for rho ∝ (R² + z²/q²)^(-γ/2)."""
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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (r, theta, phi, xyz) with xyz shape (n,3) in the model frame."""
    if n < 0:
        raise ValueError("n must be >= 0.")
    if q <= 0:
        raise ValueError("q must be > 0.")

    r = _sample_radius_powerlaw(rng, n=n, rin=rin, rout=rout, gamma=gamma)

    # theta rejection sampling over [0, pi]
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

    return r, theta, phi, np.column_stack([x, y, z])


# -----------------------------------------------------------------------------
# BALRoGO sampling
# -----------------------------------------------------------------------------


_BALROGO_DYN = None


def _import_balrogo_dynamics():
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
    """Sample N values using BALRoGO's Sanders–Evans sampler.

    BALRoGO uses NumPy's *global* RNG. We seed it from our local rng to keep
    ``seed=...`` deterministic.
    """
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

    np_seed = int(rng.integers(0, 2**32 - 1, dtype=np.uint32))
    np.random.seed(np_seed)

    samples = dyn.mom_sample_generator(
        mom_stats,
        eps=eps,
        nsig=int(nsig),
        debug=bool(debug),
    )

    if samples is None:
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    samples = np.asarray(samples, dtype=float)
    if samples.shape != (n,) or not np.all(np.isfinite(samples)):
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    return samples


# -----------------------------------------------------------------------------
# ScaleFree extraction helpers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GHParams:
    mu: float
    sig: float
    h3: float
    h4: float


@dataclass(frozen=True)
class IntrinsicMoments:
    rho: float
    vphi: float
    vr2: float
    vth2: float
    vphi2: float


def _extract_intrinsic_point(res) -> IntrinsicMoments:
    blk = res.blocks.get("intrinsic_point")
    if blk is None:
        raise RuntimeError("vprofile did not return an 'intrinsic_point' block.")

    cols = list(blk.get("columns", []))
    data = np.asarray(blk.get("data"))
    if data.ndim != 2 or data.shape[0] < 1:
        raise RuntimeError("'intrinsic_point' block is empty or malformed.")
    row = data[0]

    def _get(name: str) -> float:
        if name not in cols:
            raise KeyError(f"Missing column '{name}' in intrinsic_point cols={cols}")
        return float(row[cols.index(name)])

    return IntrinsicMoments(
        rho=_get("rho"),
        vphi=_get("vphi"),
        vr2=_get("vr2"),
        vth2=_get("vth2"),
        vphi2=_get("vphi2"),
    )


# Add near top-level in mock.py (module scope)
_WARNED_NO_VP_INTRINSIC = False


def _extract_vp_intrinsic_row(res, icomp: int) -> Dict[str, float]:
    """
    Return the VP summary row for intrinsic component `icomp` (1=vr,2=vtheta,3=vphi).

    Compatibility behavior:
      - Preferred: res.blocks["vp_intrinsic"]
      - Fallback: some older builds might store intrinsic VP summary under "vp"
                  with an "icomp" column.
      - If neither is present: return zeros for h3/h4 (Gaussian fallback).
    """
    global _WARNED_NO_VP_INTRINSIC

    vp = res.blocks.get("vp_intrinsic")

    # Fallback: sometimes intrinsic VP summary may be emitted as kind=vp
    if vp is None:
        cand = res.blocks.get("vp")
        if cand is not None:
            cols = list(cand.get("columns", []))
            if "icomp" in cols:
                vp = cand

    if vp is None:
        if not _WARNED_NO_VP_INTRINSIC:
            _WARNED_NO_VP_INTRINSIC = True
            print(
                "[scalefree.mock] WARNING: backend did not return 'vp_intrinsic'. "
                "Falling back to Gaussian sampling (h3=h4=0). "
                "To enable GH sampling, delete cached backend "
                "'~/Library/Caches/scalefree/scalefree.e' and rerun to rebuild."
            )
        return {"icomp": float(icomp), "h3": 0.0, "h4": 0.0}

    cols = list(vp.get("columns", []))
    data = np.asarray(vp.get("data"))
    if data.ndim != 2 or data.size == 0:
        return {"icomp": float(icomp), "h3": 0.0, "h4": 0.0}

    def idx(name: str) -> int:
        try:
            return cols.index(name)
        except ValueError as e:
            raise KeyError(f"Missing column '{name}' in VP summary cols={cols}") from e

    i_ic = idx("icomp")
    hit = np.where(data[:, i_ic].astype(int, copy=False) == int(icomp))[0]
    if hit.size != 1:
        # Ambiguous; do not crash mock generation
        return {"icomp": float(icomp), "h3": 0.0, "h4": 0.0}

    row = data[hit[0]]
    out = {c: float(row[j]) for j, c in enumerate(cols)}
    # Ensure keys exist for downstream code
    out.setdefault("h3", 0.0)
    out.setdefault("h4", 0.0)
    return out


def _gh_params_for_theta_bin(
    *,
    runner,
    potential: Callable[[], int],
    gamma: float,
    q: float,
    df: int,
    beta: float,
    s: float,
    t: float,
    theta_deg: float,
    integration: int,
    ngl_or_eps: float,
    maxmom: int,
    vp_smooth_eps: float,
    verbose_vp: int,
) -> Dict[int, GHParams]:
    """Run vprofile (intrinsic, point) once for this theta bin center."""
    res = runner.vprofile(
        potential=potential,
        gamma=float(gamma),
        q=float(q),
        df=int(df),
        beta=float(beta),
        s=float(s),
        t=float(t),
        # inclination/xi are irrelevant for intrinsic mode but are always prompted
        # by the Fortran frontend. Use benign defaults.
        inclination=90.0,
        xi=0.0,
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

    mom = _extract_intrinsic_point(res)

    # Component mapping (Fortran): 1=vr, 2=vtheta, 3=vphi
    # Anchor (mu,sig) to intrinsic moments for stability.
    # Use VP-provided (h3,h4) to add non-Gaussianity.
    out: Dict[int, GHParams] = {}
    for icomp in (1, 2, 3):
        vp_row = _extract_vp_intrinsic_row(res, icomp)
        h3 = float(vp_row.get("h3", 0.0))
        h4 = float(vp_row.get("h4", 0.0))

        if icomp == 1:
            mu = 0.0
            sig2 = float(mom.vr2)
        elif icomp == 2:
            mu = 0.0
            sig2 = float(mom.vth2)
        else:
            mu = float(mom.vphi)
            # mom.vphi2 is <vphi^2>, not variance.
            sig2 = float(mom.vphi2) - mu * mu

        sig2 = max(sig2, 1e-14)
        sig = float(np.sqrt(sig2))

        out[icomp] = GHParams(mu=mu, sig=sig, h3=h3, h4=h4)

    return out


# -----------------------------------------------------------------------------
# Coordinate transforms: spherical -> Cartesian
# -----------------------------------------------------------------------------


def _sph_to_cart_v(
    *,
    theta: np.ndarray,
    phi: np.ndarray,
    vr: np.ndarray,
    vtheta: np.ndarray,
    vphi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert (vr, vtheta, vphi) at (theta, phi) to (vx, vy, vz)."""
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    # Unit vectors:
    # e_r     = ( st cp,  st sp,  ct)
    # e_theta = ( ct cp,  ct sp, -st)
    # e_phi   = (   -sp,    cp,   0)
    vx = vr * st * cp + vtheta * ct * cp - vphi * sp
    vy = vr * st * sp + vtheta * ct * sp + vphi * cp
    vz = vr * ct - vtheta * st
    return vx, vy, vz


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------


def mock(
    *,
    # Mass model parameters
    potential: Callable[[], int] = lambda: 1,
    gamma: float = 4.0,
    beta: float = 0.0,
    df: int = 1,
    s: float = 0.5,
    t: float = 0.0,
    q: float = 1.0,
    # Sampling
    nsamples: int = 1000,
    theta_bins: int = 45,
    rin: float = 1.0,
    rout: float = 1000.0,
    seed: Optional[int] = None,
    # vprofile numerical settings
    integration: int = 1,
    ngl_or_eps: float = 0.0,
    maxmom: int = 4,
    vp_smooth_eps: float = 0.0,
    verbose_vp: int = 0,
    # BALRoGO sampling
    nsig: int = 10,
    debug: bool = False,
) -> np.ndarray:
    """Generate a 6D mock in Cartesian coordinates.

    Parameters
    ----------
    theta_bins:
        Number of bins dividing 90° (i.e., [0°, 90°]). Each star is assigned
        to a bin based on its intrinsic θ (using symmetry).
    maxmom:
        Maximum number of intrinsic moments used by the Fortran VP
        reconstruction. ``algorithm`` is forced to 3.
    """
    n = int(nsamples)
    if n < 0:
        raise ValueError("nsamples must be >= 0.")
    nb = int(theta_bins)
    if nb < 1:
        raise ValueError("theta_bins must be >= 1.")
    if q <= 0:
        raise ValueError("q must be > 0.")

    rng = np.random.default_rng(_coerce_seed(seed))
    potential = _coerce_potential(potential)

    # 1) sample positions
    _r, theta, phi, xyz = _sample_positions_scalefree(
        n=n, rin=rin, rout=rout, gamma=gamma, q=q, rng=rng
    )

    # 2) theta symmetry map to [0, 90]
    theta_sym = np.minimum(theta, np.pi - theta)
    theta_sym_deg = np.degrees(theta_sym)

    edges = np.linspace(0.0, 90.0, nb + 1)
    bin_id = np.digitize(theta_sym_deg, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, nb - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 3) compute intrinsic GH params per occupied theta bin
    from scalefree import ScaleFreeRunner  # local import avoids circular issues

    runner = ScaleFreeRunner()
    gh_by_bin: Dict[int, Dict[int, GHParams]] = {}

    occupied = np.unique(bin_id)
    for b in occupied:
        b = int(b)
        gh_by_bin[b] = _gh_params_for_theta_bin(
            runner=runner,
            potential=potential,
            gamma=gamma,
            q=q,
            df=df,
            beta=beta,
            s=s,
            t=t,
            theta_deg=float(centers[b]),
            integration=integration,
            ngl_or_eps=ngl_or_eps,
            maxmom=maxmom,
            vp_smooth_eps=vp_smooth_eps,
            verbose_vp=verbose_vp,
        )

    # 4) sample intrinsic velocities per bin
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
        p_r, p_t, p_p = gh[1], gh[2], gh[3]

        vr[sel] = _sample_balrogo_gh(
            mean=p_r.mu,
            sigma=p_r.sig,
            h3=p_r.h3,
            h4=p_r.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            debug=debug,
        )
        vth[sel] = _sample_balrogo_gh(
            mean=p_t.mu,
            sigma=p_t.sig,
            h3=p_t.h3,
            h4=p_t.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            debug=debug,
        )
        vph[sel] = _sample_balrogo_gh(
            mean=p_p.mu,
            sigma=p_p.sig,
            h3=p_p.h3,
            h4=p_p.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            debug=debug,
        )

    # 5) spherical -> Cartesian velocities
    vx, vy, vz = _sph_to_cart_v(theta=theta, phi=phi, vr=vr, vtheta=vth, vphi=vph)

    out = np.column_stack([xyz, vx, vy, vz])

    if debug:
        counts = np.bincount(bin_id, minlength=nb)
        nz = np.count_nonzero(counts)
        print(
            f"[mock] nsamples={n} theta_bins={nb} occupied_bins={nz} (0..90 deg symmetry)"
        )
        if nz:
            print(
                f"[mock] min_bin_count={counts[counts>0].min()} max_bin_count={counts.max()}"
            )

    return out
