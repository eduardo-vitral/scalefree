#!/usr/bin/env python3
"""
scalefree.mock

Fast Python-only mock generator that:

1) Samples N points from the ScaleFree density distribution.
2) Bins projected-plane angle into nbins over [0, 180) (symmetry).
3) For each (occupied) bin, runs ScaleFree vprofile once to obtain
   (gauss_V, gauss_sig, h3, h4) for iproj=1,2,3.
4) Assigns each star the parameters of its bin and samples
   (vlos, vposr, vpost) using BALRoGO's Sanders–Evans PDFs.
5) Converts (vlos, vposr, vpost) into (vx, vy, vz) in the model frame.
6) Returns an (N,6) array: (x, y, z, vx, vy, vz).

Notes
-----
- Uses symmetry: only computes bins over [0, 180). Star angles are mapped with psi % 180.
- No Fortran mock code is used. Only the ScaleFree backend (via ScaleFreeRunner.vprofile)
  is used for GH parameters.
- Sampling is fast: loops over bins (<= nbins), not over stars.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple

import numpy as np


def _coerce_seed(seed: Optional[int]) -> Optional[int]:
    """NumPy default_rng requires a non-negative seed; map negatives into uint32 space."""
    if seed is None:
        return None
    s = int(seed)
    return s % (2**32)  # maps negative -> [0, 2**32-1]


def _coerce_potential(potential):
    """vprofile expects a callable; allow int for backward compatibility."""
    if callable(potential):
        return potential
    return lambda: int(potential)


# -----------------------------------------------------------------------------
# Density sampling (matches prior ScaleFree rejection logic)
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
    """
    Acceptance proportional to (sin^2(theta) + cos^2(theta)/q^2)^(-gamma/2),
    corresponding to rho ∝ (R^2 + z^2/q^2)^(-gamma/2) with symmetry axis z.
    """
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
    """Return xyz (N,3) in the model frame."""
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
# Geometry helpers: rotation to/from "sky frame" for a given inclination
# -----------------------------------------------------------------------------

def _rotate_to_sky_xyz(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, inclination_deg: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evans & de Zeeuw (1994) convention as in VPOS.md / scalefree.f:

      x' = y
      y' = -x cos i + z sin i
      z' =  x sin i + z cos i

    z' is LOS; (x', y') is sky plane.
    """
    i = np.deg2rad(float(inclination_deg))
    ci = np.cos(i)
    si = np.sin(i)

    xp = y
    yp = -x * ci + z * si
    zp =  x * si + z * ci
    return xp, yp, zp


def _rotate_from_sky_v(
    vx_sky: np.ndarray, vy_sky: np.ndarray, vz_sky: np.ndarray, inclination_deg: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Inverse mapping back to model frame, consistent with the above:

      vy = vx'
      vx = -vy' cos i + vz' sin i
      vz =  vy' sin i + vz' cos i
    """
    i = np.deg2rad(float(inclination_deg))
    ci = np.cos(i)
    si = np.sin(i)

    vy = vx_sky
    vx = -vy_sky * ci + vz_sky * si
    vz =  vy_sky * si + vz_sky * ci
    return vx, vy, vz



# -----------------------------------------------------------------------------
# BALRoGO (Sanders–Evans) sampling: vectorized inverse-CDF per (mu,sig,h3,h4)
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


def _cdf_from_pdf_grid(x: np.ndarray, pdf: np.ndarray) -> np.ndarray:
    pdf = np.asarray(pdf, dtype=float)
    pdf = np.maximum(pdf, 0.0)
    area = np.trapz(pdf, x)
    if not np.isfinite(area) or area <= 0:
        raise RuntimeError("PDF integral non-positive; parameters likely non-physical.")
    pdf = pdf / area
    dx = np.diff(x)
    c = np.empty_like(x)
    c[0] = 0.0
    c[1:] = np.cumsum(0.5 * (pdf[:-1] + pdf[1:]) * dx)
    c[-1] = 1.0
    return c


def _sample_balrogo_gh(
    *,
    mean: float,
    sigma: float,
    h3: float,
    h4: float,
    n: int,
    rng: np.random.Generator,
    nsig: int = 10,
    grid_n: Optional[int] = None,  # kept for API compatibility; BALRoGO ignores it
    debug: bool = False,
) -> np.ndarray:
    """
    Sample N values from BALRoGO's Sanders–Evans GH sampler.

    Uses dynamics.mom_sample_generator(), which builds the PDF internally and
    inverse-CDF samples. This avoids our own PDF normalization/integration path.

    Notes:
    - dynamics.mom_sample_generator() requires an `eps` array (measurement errors).
      For intrinsic mock sampling we pass eps=0.
    - dynamics.mom_sample_generator() uses NumPy's global RNG. For reproducibility
      with our `rng`, we seed NumPy locally from `rng` before calling.
    """
    if n <= 0:
        return np.empty((0,), dtype=float)

    if (not np.isfinite(mean)) or (not np.isfinite(sigma)) or (sigma <= 0):
        # Fall back safely
        return rng.normal(loc=float(mean), scale=max(float(sigma), 1e-12), size=n)

    dyn = _import_balrogo_dynamics()

    # Build mom_stats shape (4,2): values in col 0, uncertainties in col 1 (unused by sampler)
    mom_stats = np.array(
        [
            [float(mean), 0.0],
            [float(sigma), 0.0],
            [float(h3), 0.0],
            [float(h4), 0.0],
        ],
        dtype=float,
    )

    # BALRoGO requires eps array length n; use 0 (intrinsic, no measurement error)
    eps = np.zeros(n, dtype=float)

    # Deterministic bridge: seed NumPy global RNG from our rng for reproducibility
    np_seed = int(rng.integers(0, 2**32 - 1, dtype=np.uint32))
    np.random.seed(np_seed)

    samples = dyn.mom_sample_generator(
        mom_stats, eps=eps, nsig=int(nsig), debug=bool(debug)
    )

    # BALRoGO returns None for non-physical parameter combinations
    if samples is None:
        # Conservative fallback: pure Gaussian
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    samples = np.asarray(samples, dtype=float)
    if samples.shape != (n,) or not np.all(np.isfinite(samples)):
        return rng.normal(loc=float(mean), scale=float(sigma), size=n)

    return samples


# -----------------------------------------------------------------------------
# ScaleFree GH parameter extraction for a projected-plane angle bin
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GHParams:
    mu: float
    sig: float
    h3: float
    h4: float


def _extract_gh_params(res, iproj: int) -> GHParams:
    vp = res.blocks.get("vp")
    if vp is None:
        raise RuntimeError("vprofile did not return a 'vp' summary block.")

    cols = list(vp.get("columns", []))
    data = np.asarray(vp.get("data"))

    def idx(name: str) -> int:
        try:
            return cols.index(name)
        except ValueError as e:
            raise KeyError(f"Missing column '{name}' in vp summary. cols={cols}") from e

    i_ip = idx("iproj")
    i_mu = idx("gauss_V")
    i_sig = idx("gauss_sig")
    i_h3 = idx("h3")
    i_h4 = idx("h4")

    hit = np.where(data[:, i_ip].astype(int, copy=False) == int(iproj))[0]
    if hit.size != 1:
        raise RuntimeError(
            f"Expected exactly one vp row for iproj={iproj}; found {hit.size}."
        )
    r = data[hit[0]]
    return GHParams(
        mu=float(r[i_mu]), sig=float(r[i_sig]), h3=float(r[i_h3]), h4=float(r[i_h4])
    )


def _scalefree_params_for_bin(
    *,
    runner,
    potential: Callable[[], int],
    gamma: float,
    q: float,
    df: int,
    beta: float,
    s: float,
    t: float,
    inclination: float,
    theta_deg: float,
    xi: float,
    integration: int,
    ngl_or_eps: int,
    algorithm: int,
    maxmom: int,
    average: bool,
    verbose_vp: int,
) -> Dict[int, GHParams]:
    """
    Run vprofile once for this projected-plane angle theta_deg,
    and return GH params for iproj=1,2,3.
    """
    res = runner.vprofile(
        potential=potential,
        gamma=float(gamma),
        q=float(q),
        df=int(df),
        beta=float(beta),
        s=float(s),
        t=float(t),
        inclination=float(inclination),
        xi=float(xi),
        theta=float(theta_deg),
        integration=int(integration),
        ngl_or_eps=int(ngl_or_eps),
        algorithm=int(algorithm),
        maxmom=int(maxmom),
        average=bool(average),
        usevp=True,
        verbose_vp=int(verbose_vp),
        output_path=None,
        debug_prompts=False,
    )
    return {
        1: _extract_gh_params(res, 1),
        2: _extract_gh_params(res, 2),
        3: _extract_gh_params(res, 3),
    }


# -----------------------------------------------------------------------------
# Public API: mock() -> (N,6)
# -----------------------------------------------------------------------------


def mock(
    *,
    # ScaleFree model parameters
    potential: Callable[[], int] = lambda: 1,
    gamma: float = 4.0,
    q: float = 1.0,
    df: int = 1,
    beta: float = 0.0,
    s: float = 0.5,
    t: float = 0.0,
    inclination: float = 90.0,
    xi: float = 0.0,
    # density sampling
    nsamples: int = 1000,
    seed: Optional[int] = None,
    rin: float = 1.0,
    rout: float = 1000.0,
    # vprofile numerical settings
    integration: int = 1,
    ngl_or_eps: int = 0,
    algorithm: int = 3,
    maxmom: int = 4,
    average: bool = False,
    verbose_vp: int = 0,
    # binning / performance
    nbins: int = 180,
    nsig: int = 10,
    grid_n: Optional[int] = None,
    debug: bool = False,
) -> np.ndarray:
    """
    Generate a 6D mock (x,y,z,vx,vy,vz) array.

    Binning symmetry:
      - uses projected-plane angle psi in degrees
      - maps psi -> psi % 180  (so only [0,180) bins are needed)
    """
    n = int(nsamples)
    if n < 0:
        raise ValueError("nsamples must be >= 0.")
    nb = int(nbins)
    if nb < 1:
        raise ValueError("nbins must be >= 1.")

    rng = np.random.default_rng(_coerce_seed(seed))
    potential = _coerce_potential(potential)

    # 1) sample positions in model frame
    xyz = _sample_positions_scalefree(
        n=n, rin=rin, rout=rout, gamma=gamma, q=q, rng=rng
    )
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

    # 2) compute projected-plane angle psi for each star in sky frame
    x_sky, y_sky, z_sky = _rotate_to_sky_xyz(x, y, z, inclination_deg=inclination)
    psi = np.degrees(np.arctan2(y_sky, x_sky))  # [-180, 180]
    psi = (psi + 360.0) % 360.0  # [0, 360)
    psi = psi % 180.0  # [0, 180) symmetry

    # 3) bin in [0,180)
    edges = np.linspace(0.0, 180.0, nb + 1)
    bin_id = np.digitize(psi, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, nb - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 4) compute ScaleFree GH params per occupied bin
    from scalefree import ScaleFreeRunner  # local import avoids circular issues

    runner = ScaleFreeRunner()

    gh_by_bin: Dict[int, Dict[int, GHParams]] = {}
    occupied = np.unique(bin_id)
    for b in occupied:
        gh_by_bin[int(b)] = _scalefree_params_for_bin(
            runner=runner,
            potential=potential,
            gamma=gamma,
            q=q,
            df=df,
            beta=beta,
            s=s,
            t=t,
            inclination=inclination,
            theta_deg=float(centers[int(b)]),
            xi=xi,
            integration=integration,
            ngl_or_eps=ngl_or_eps,
            algorithm=algorithm,
            maxmom=maxmom,
            average=average,
            verbose_vp=verbose_vp,
        )

    # 5) sample (vlos, vposr, vpost) per bin, vectorized per bin
    vlos = np.empty(n, dtype=float)
    vposr = np.empty(n, dtype=float)
    vpost = np.empty(n, dtype=float)

    for b in occupied:
        b = int(b)
        sel = bin_id == b
        m = int(np.sum(sel))
        if m == 0:
            continue

        gh = gh_by_bin[b]
        p1, p2, p3 = gh[1], gh[2], gh[3]

        vlos[sel] = _sample_balrogo_gh(
            mean=p1.mu,
            sigma=p1.sig,
            h3=p1.h3,
            h4=p1.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            grid_n=grid_n,
        )
        vposr[sel] = _sample_balrogo_gh(
            mean=p2.mu,
            sigma=p2.sig,
            h3=p2.h3,
            h4=p2.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            grid_n=grid_n,
        )
        vpost[sel] = _sample_balrogo_gh(
            mean=p3.mu,
            sigma=p3.sig,
            h3=p3.h3,
            h4=p3.h4,
            n=m,
            rng=rng,
            nsig=nsig,
            grid_n=grid_n,
        )

    # 6) convert (vlos, vposr, vpost) -> (vx,vy,vz)
    # In the sky plane, vposr/vpost are polar components relative to psi
    psi_rad = np.radians(psi)
    c = np.cos(psi_rad)
    s_ = np.sin(psi_rad)

    vx_sky = vposr * c - vpost * s_
    vy_sky = vposr * s_ + vpost * c
    vz_sky = vlos

    vx, vy, vz = _rotate_from_sky_v(vx_sky, vy_sky, vz_sky, inclination_deg=inclination)

    out = np.column_stack([x, y, z, vx, vy, vz])

    if debug:
        # lightweight debug: report bin occupancy
        counts = np.bincount(bin_id, minlength=nb)
        nz = np.count_nonzero(counts)
        print(
            f"[mock] nsamples={n} nbins={nb} occupied_bins={nz} (symmetry range 0..180 deg)"
        )
        print(
            f"[mock] min_bin_count={counts[counts>0].min() if nz else 0} max_bin_count={counts.max() if counts.size else 0}"
        )

    return out
