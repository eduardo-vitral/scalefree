#!/usr/bin/env python3
"""scalefree.mock

Intrinsic mock generator based on `scalefree.vmoments`.

This file replaces the previous projected-frame mock generator. It creates
Cartesian 6D samples (x,y,z,vx,vy,vz) by:

1) Sampling positions from the ScaleFree *volume* density with intrinsic
   flattening `q`.
2) Folding the intrinsic meridional angle theta into [0, 90] degrees and
   binning it into `nbins` bins.
3) For each occupied bin, calling `ScaleFreeRunner.vprofile` in intrinsic
   point mode to obtain intrinsic VP diagnostics (`vp_intrinsic`), forcing
   `algorithm=3`.
4) Using the VP Gaussian-fit parameters (`gauss_V`, `gauss_sig`) and GH moments
   (`h3`, `h4`) per intrinsic component (v_r, v_theta, v_phi) to sample
   velocities via BALRoGO.
5) Converting spherical velocities to Cartesian in the model frame.

Notes
-----
- The mass model is specified by: potential, gamma, q, df, beta, s, t.
- `maxmom` is user-controlled; `algorithm` is forced to 3.
- The intrinsic `vprofile` results are at r=1 in dimensionless units, but
  the ScaleFree nature allows scaling; the mock follows the existing package
  convention (as the earlier mock did).

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np


# -----------------------------------------------------------------------------
# Helpers: density sampling (ScaleFree ellipsoidal power-law)
# -----------------------------------------------------------------------------


def _coerce_seed(seed: Optional[int]) -> Optional[int]:
    if seed is None:
        return None
    return int(seed) % (2**32)


def _sample_radius_powerlaw(
    rng: np.random.Generator, *, n: int, rin: float, rout: float, gamma: float
) -> np.ndarray:
    """Sample spherical r with p(r) ∝ r^(2-gamma) over [rin, rout]."""
    if rin <= 0 or rout <= rin:
        raise ValueError("Require 0 < rin < rout.")
    a = 3.0 - float(gamma)
    u = rng.random(n)
    if abs(a) > 1e-14:
        return (u * (rout**a - rin**a) + rin**a) ** (1.0 / a)
    # gamma == 3 -> log-uniform
    return rin * np.exp(u * np.log(rout / rin))


def _theta_accept_prob(
    theta: np.ndarray,
    *,
    q: float,
    gamma: float,
) -> np.ndarray:
    """Acceptance ∝ (sin^2θ + cos^2θ/q^2)^(-gamma/2)."""
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
    """Return xyz (N,3) and spherical (r, theta, phi) in the model frame."""
    if n < 0:
        raise ValueError("n must be >= 0.")
    if q <= 0:
        raise ValueError("q must be > 0.")

    r = _sample_radius_powerlaw(rng, n=n, rin=rin, rout=rout, gamma=gamma)

    # theta rejection sampling over [0, pi]
    theta = np.empty(n, dtype=float)
    filled = 0

    # Conservative bound so acceptance <= 1
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

    return np.column_stack([x, y, z]), r, theta, phi


# -----------------------------------------------------------------------------
# Helpers: BALRoGO GH sampling
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
    nsig: int,
    debug: bool,
) -> np.ndarray:
    if n <= 0:
        return np.empty((0,), dtype=float)

    if (not np.isfinite(mean)) or (not np.isfinite(sigma)) or sigma <= 0:
        return rng.normal(
            loc=float(mean),
            scale=max(float(sigma), 1e-12),
            size=n,
        )

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

    # BALRoGO uses NumPy global RNG; seed deterministically from our RNG.
    np_seed = int(rng.integers(0, 2**32 - 1, dtype=np.uint32))
    np.random.seed(np_seed)

    samples = dyn.mom_sample_generator(
        mom_stats, eps=eps, nsig=int(nsig), debug=bool(debug)
    )
    if samples is None:
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    samples = np.asarray(samples, dtype=float)
    if samples.shape != (n,) or not np.all(np.isfinite(samples)):
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    return samples


# -----------------------------------------------------------------------------
# VP extraction via vmoments
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GHParams:
    mu: float
    sig: float
    h3: float
    h4: float


def _repo_root() -> Path:
    # scalefree/mock.py -> scalefree/ -> repo root
    return Path(__file__).resolve().parents[1]


def _default_intrinsic_exe() -> Path:
    """Mirror scripts/quick_user_run.py default backend path."""
    return _repo_root() / "fortran_src" / "scalefree_intrvp.e"


def _extract_intrinsic_moments(res) -> Tuple[float, float, float, float]:
    """Return (vphi, vr2, vth2, vphi2) from intrinsic_point/shell_average."""
    for k in ("intrinsic_point", "intrinsic_shell_average"):
        if k in res.blocks:
            blk = res.blocks[k]
            cols = [c.lower() for c in blk.get("columns", [])]
            data = blk.get("data")
            if data is None or getattr(data, "size", 0) == 0:
                break
            row = data[0]
            idx = {c: i for i, c in enumerate(cols)}
            vphi = (
                float(
                    row[idx.get("vphi", 1)],
                )
                if "vphi" in idx
                else float(row[1])
            )
            vr2 = (
                float(
                    row[idx.get("vr2", 2)],
                )
                if "vr2" in idx
                else float(row[2])
            )
            vth2 = (
                float(
                    row[idx.get("vth2", 3)],
                )
                if "vth2" in idx
                else float(row[3])
            )
            vphi2 = (
                float(
                    row[idx.get("vphi2", 4)],
                )
                if "vphi2" in idx
                else float(row[4])
            )
            return vphi, vr2, vth2, vphi2
    raise KeyError("No intrinsic moments block found in ScaleFreeResult.")


def _gh_params_for_theta_bin(
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
    debug_prompts: bool,
) -> Dict[int, GHParams]:
    """Compute GH parameters (mu,sig,h3,h4)
    for (v_r, v_th, v_phi) at a theta bin."""

    # print("\n\nChecking quantities:")
    # print(potential, type(potential), 2)
    # print(float(gamma), type(float(gamma)), 2)
    # print(float(q), type(float(q)), 0.608)
    # print(int(df), type(int(df)), 1)
    # print(float(beta), type(float(beta)), 0.189)
    # print(float(s), type(float(s)), 0.5)
    # print(float(t), type(float(t)), 0.0)
    # print(57.1)
    # print(float(xi), float(xi), 0.0)
    # print(float(theta_deg), float(theta_deg), 0.0)
    # print(int(integration), type(int(integration)), 1)
    # print(float(ngl_or_eps), type(float(ngl_or_eps)), 0)
    # print(3)
    # print(int(maxmom), type(int(maxmom)), 8)
    # Intrinsic kinematics: inclination is not used by Fortran for iwhat=0.
    res = runner.vprofile(
        potential=potential,
        gamma=float(gamma),
        q=float(q),
        df=int(df),
        beta=float(beta),
        s=float(s),
        t=float(t),
        inclination=90,
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
        verbose_vp=verbose_vp,
        debug_prompts=debug_prompts,
    )

    # res = runner.vprofile(
    #     potential=2, gamma=2.0, q=0.608, df=1, beta=0.189, s=0.5, t=0.0,
    #     inclination=57.1, xi=0.0, theta=0.0,
    #     integration=1, ngl_or_eps=0, algorithm=3, maxmom=8,
    #     kinematics="intrinsic", average=False,
    #     usevp=True, verbose_vp=0,
    # )

    # Primary source: VP gaussian-fit summary
    vpblk = res.blocks.get("vp_intrinsic")
    by = vpblk.get("by_icomp", {}) if isinstance(vpblk, dict) else {}

    # Fallback (only if VP block is incomplete): intrinsic moments
    vphi_m, vr2_m, vth2_m, vphi2_m = _extract_intrinsic_moments(res)

    out: Dict[int, GHParams] = {}
    for icomp in (1, 2, 3):
        row = by.get(icomp, {}) if isinstance(by, dict) else {}

        mu = float(row.get("gauss_V", np.nan))
        sig = float(row.get("gauss_sig", np.nan))
        h3 = float(row.get("h3", np.nan))
        h4 = float(row.get("h4", np.nan))

        # Robustness: some legacy outputs may omit these keys.
        if not np.isfinite(mu):
            mu = float(row.get("true_V", 0.0))
        if not np.isfinite(sig) or sig <= 0:
            sig = float(row.get("true_sig", np.nan))

        # Final fallback to moment-derived mu/sig if still missing.
        if not np.isfinite(mu) or not np.isfinite(sig) or sig <= 0:
            if icomp == 1:
                mu = 0.0
                sig = float(np.sqrt(max(vr2_m, 0.0)))
            elif icomp == 2:
                mu = 0.0
                sig = float(np.sqrt(max(vth2_m, 0.0)))
            else:
                mu = float(vphi_m)
                sig = float(np.sqrt(max(vphi2_m - vphi_m * vphi_m, 0.0)))

        # h3/h4 fallback (if missing): assume Gaussian (0,0)
        if not np.isfinite(h3):
            h3 = 0.0
        if not np.isfinite(h4):
            h4 = 0.0

        out[icomp] = GHParams(
            mu=float(mu),
            sig=float(sig),
            h3=float(h3),
            h4=float(h4),
        )

    return out


# -----------------------------------------------------------------------------
# Coordinate transforms: spherical -> Cartesian
# -----------------------------------------------------------------------------


def _sph_vel_to_cart(
    *,
    theta: np.ndarray,
    phi: np.ndarray,
    vr: np.ndarray,
    vtheta: np.ndarray,
    vphi: np.ndarray,
) -> np.ndarray:
    """Convert (vr, vtheta, vphi) to (vx, vy, vz)
    for standard physics convention."""
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    # Unit vectors
    # e_r
    erx = st * cp
    ery = st * sp
    erz = ct
    # e_theta
    etx = ct * cp
    ety = ct * sp
    etz = -st
    # e_phi
    epx = -sp
    epy = cp
    epz = 0.0

    vx = vr * erx + vtheta * etx + vphi * epx
    vy = vr * ery + vtheta * ety + vphi * epy
    vz = vr * erz + vtheta * etz + vphi * epz

    return np.column_stack([vx, vy, vz])


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------


def mock(
    *,
    potential: Any,
    gamma: float,
    q: float,
    df: int,
    beta: float,
    s: float,
    t: float,
    nsamples: int,
    nbins: int,
    seed: Optional[int] = None,
    rin: float = 0.5,
    rout: float = 50.0,
    # numerical backend settings
    integration: int = 1,
    ngl_or_eps: float = 0.0,
    vp_smooth_eps: float = 0.0,
    xi: float = 0.0,
    verbose_vp: int = 0,
    maxmom: int = 10,
    # sampling settings
    nsig: int = 10,
    debug: bool = False,
    # advanced: override backend
    exe_path: Optional[Path] = None,
) -> np.ndarray:
    """Generate an intrinsic ScaleFree mock.

    Returns
    -------
    xyzv : ndarray, shape (N, 6)
        Columns: x y z vx vy vz (Cartesian, model frame).
    """

    n = int(nsamples)
    if n < 0:
        raise ValueError("nsamples must be >= 0")

    nb = int(nbins)
    if nb < 1:
        raise ValueError("nbins must be >= 1")

    rng = np.random.default_rng(_coerce_seed(seed))

    # 1) positions
    xyz, r, theta, phi = _sample_positions_scalefree(
        n=n,
        rin=float(rin),
        rout=float(rout),
        gamma=float(gamma),
        q=float(q),
        rng=rng,
    )

    # 2) theta-binning over [0, 90] degrees by symmetry
    theta_fold = np.minimum(theta, np.pi - theta)
    theta_deg = np.degrees(theta_fold)

    edges = np.linspace(0.0, 90.0, nb + 1)
    bin_id = np.digitize(theta_deg, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, nb - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 3) backend runner (prefer the same default as scripts/quick_user_run.py)
    from scalefree import ScaleFreeRunner  # local import avoids circulars

    exe = (
        Path(exe_path).expanduser().resolve()
        if exe_path is not None
        else _default_intrinsic_exe()
    )
    runner = ScaleFreeRunner(exe_path=exe)

    # 4) compute GH params per occupied bin
    gh_by_bin: Dict[int, Dict[int, GHParams]] = {}
    occupied = np.unique(bin_id)

    for b in occupied:
        b = int(b)
        gh_by_bin[b] = _gh_params_for_theta_bin(
            runner=runner,
            potential=potential,
            gamma=float(gamma),
            q=float(q),
            df=int(df),
            beta=float(beta),
            s=float(s),
            t=float(t),
            theta_deg=float(centers[b]),
            xi=float(xi),
            integration=int(integration),
            ngl_or_eps=float(ngl_or_eps),
            maxmom=int(maxmom),
            vp_smooth_eps=float(vp_smooth_eps),
            verbose_vp=int(verbose_vp),
            debug_prompts=bool(debug),
        )

    # 5) sample velocities per bin
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
        p_r = gh[1]
        p_th = gh[2]
        p_ph = gh[3]

        vr[sel] = _sample_balrogo_gh(
            mean=p_r.mu,
            sigma=p_r.sig,
            h3=p_r.h3,
            h4=p_r.h4,
            n=m,
            rng=rng,
            nsig=int(nsig),
            debug=bool(debug),
        )
        vth[sel] = _sample_balrogo_gh(
            mean=p_th.mu,
            sigma=p_th.sig,
            h3=p_th.h3,
            h4=p_th.h4,
            n=m,
            rng=rng,
            nsig=int(nsig),
            debug=bool(debug),
        )
        vph[sel] = _sample_balrogo_gh(
            mean=p_ph.mu,
            sigma=p_ph.sig,
            h3=p_ph.h3,
            h4=p_ph.h4,
            n=m,
            rng=rng,
            nsig=int(nsig),
            debug=bool(debug),
        )

    # 6) spherical -> Cartesian velocities (use original theta, phi)
    vxyz = _sph_vel_to_cart(theta=theta, phi=phi, vr=vr, vtheta=vth, vphi=vph)

    return np.column_stack([xyz, vxyz])


__all__ = ["mock"]
