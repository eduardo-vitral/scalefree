from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Optional, Union, Any

import os
import shutil
import subprocess
import re

import numpy as np

try:
    from scipy.optimize import least_squares
except Exception:  # pragma: no cover
    least_squares = None


class HermiteBackendError(RuntimeError):
    """Raised when the ghermite/fitvp backend cannot be found or built."""


# --------------------------
# Backend discovery / build
# --------------------------


def _have_gfortran() -> bool:
    return shutil.which("gfortran") is not None


def _site_packages_root() -> Path:
    # scalefree/hermite.py -> scalefree/ -> (site-packages or repo root)
    return Path(__file__).resolve().parent.parent


def _packaged_fortran_source() -> Path:
    # <root>/fortran_src/fitvp_stdout.f
    return _site_packages_root() / "fortran_src" / "fitvp_stdout.f"


def _user_cache_dir() -> Path:
    if os.name == "nt":
        base = os.environ.get(
            "LOCALAPPDATA",
        ) or str(Path.home() / "AppData" / "Local")
        return Path(base) / "scalefree" / "Cache"

    sysname = ""
    if hasattr(os, "uname"):
        try:
            sysname = os.uname().sysname.lower()
        except Exception:
            sysname = ""
    if sysname == "darwin":
        return Path.home() / "Library" / "Caches" / "scalefree"

    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg) / "scalefree"
    return Path.home() / ".cache" / "scalefree"


def _default_cached_exe() -> Path:
    exe_name = "fitvp.e" if os.name != "nt" else "fitvp.exe"
    return _user_cache_dir() / exe_name


def _compile_backend(*, exe: Path, src: Path) -> None:
    if not _have_gfortran():
        raise HermiteBackendError(
            "gfortran not found. "
            + "Install a Fortran compiler,"
            + " or provide exe_path."
        )
    if not src.exists():
        raise HermiteBackendError(
            "Fortran source not found.\n"
            f"Expected: {src}\n"
            "Ensure fortran_src/fitvp_stdout.f "
            + "is included in your "
            + "package/repo."
        )
    exe.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gfortran",
        "-O2",
        "-std=legacy",
        "-ffixed-line-length-none",
        "-o",
        str(exe),
        str(src),
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise HermiteBackendError(
            "Failed to compile fitvp backend.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{e.stdout}\n\nstderr:\n{e.stderr}"
        ) from e


# --------------------------
# Parsing stdout
# --------------------------

_RE_GAUSS = re.compile(
    r"^#\s*(gauss_moments|gauss_fit)\s+norm\s+(\S+)"
    + r"\s+mean\s+(\S+)\s+dispersion\s+(\S+)\s*$"
)
_RE_GH = re.compile(r"^#\s*gh_moments\s+0-10\s+(.*)$")


def _to_float(tok: Any) -> float:
    t = str(tok).strip().replace("D", "E").replace("d", "e")
    return float(t)


@dataclass(frozen=True)
class HermiteResult:
    gauss_info: Dict[str, float]
    gaussh_info: Dict[str, float]
    h_moments: Dict[str, float]
    stdout: str
    stderr: str
    exe_path: Path


def hermite(
    input_path: Union[str, Path],
    *,
    exe_path: Optional[Union[str, Path]] = None,
    rebuild: bool = False,
    env_var: str = "SCALEFREE_FITVP_EXE",
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
    """
    Compute Gaussian + Gauss-Hermite fits for a
    VP file using the Fortran backend.

    Parameters
    ----------
    input_path:
        Path to the VP data file (2 columns: v, VP(v)),
        optional header allowed.
    exe_path:
        Optional explicit path to a precompiled executable.
    rebuild:
        If True, force recompilation into the cached executable location.
    env_var:
        Name of environment variable that can point to a
        precompiled executable.

    Returns
    -------
    (gauss_info, gaussh_info, h_moments)
        gauss_info: "moment" Gaussian from raw VP moments
        (norm/mean/dispersion)
        gaussh_info: best-fit Gaussian used as GH reference
        (norm/mean/dispersion)
        h_moments: dict h0..h10
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input VP file not found: {input_path}")

    # Resolve executable
    if exe_path is not None:
        exe = Path(exe_path)
    else:
        env_exe = os.environ.get(env_var)
        exe = Path(env_exe) if env_exe else _default_cached_exe()

    src = _packaged_fortran_source()

    if rebuild or (not exe.exists()):
        _compile_backend(exe=exe, src=src)

    # Run backend: pass filepath via stdin
    try:
        p = subprocess.run(
            [str(exe)],
            input=str(input_path) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise HermiteBackendError(
            "fitvp backend failed.\n"
            + f"stdout:\n{e.stdout}\n\n"
            + f"stderr:\n{e.stderr}"
        ) from e

    gauss_info: Dict[str, float] = {}
    gaussh_info: Dict[str, float] = {}
    h_moments: Dict[str, float] = {}

    for line in p.stdout.splitlines():
        m = _RE_GAUSS.match(line.strip())
        if m:
            kind = m.group(1)
            norm = _to_float(m.group(2))
            mean = _to_float(m.group(3))
            disp = _to_float(m.group(4))
            d = {"norm": norm, "mean": mean, "dispersion": disp}
            if kind == "gauss_moments":
                gauss_info = d
            else:
                gaussh_info = d
            continue

        mg = _RE_GH.match(line.strip())
        if mg:
            vals = mg.group(1).split()
            if len(vals) < 11:
                raise HermiteBackendError(
                    "Malformed gh_moments line "
                    + "(expected 11 values). Line:\n"
                    + line
                )
            arr = [_to_float(v) for v in vals[:11]]
            keys = [f"h{i}" for i in range(11)]
            h_moments = {k: float(v) for k, v in zip(keys, arr)}
            continue

    if not gauss_info or not gaussh_info or not h_moments:
        raise HermiteBackendError(
            "Could not parse fitvp output. "
            "Expected # gauss_moments, # gauss_fit and # gh_moments lines.\n\n"
            f"stdout:\n{p.stdout}\n\nstderr:\n{p.stderr}"
        )

    return gauss_info, gaussh_info, h_moments


def hermite_result(
    input_path: Union[str, Path],
    *,
    exe_path: Optional[Union[str, Path]] = None,
    rebuild: bool = False,
    env_var: str = "SCALEFREE_FITVP_EXE",
) -> HermiteResult:
    """
    Same as hermite(), but returns a richer result object
    including stdout/stderr.
    """
    input_path = Path(input_path)

    if exe_path is not None:
        exe = Path(exe_path)
    else:
        env_exe = os.environ.get(env_var)
        exe = Path(env_exe) if env_exe else _default_cached_exe()

    src = _packaged_fortran_source()
    if rebuild or (not exe.exists()):
        _compile_backend(exe=exe, src=src)

    p = subprocess.run(
        [str(exe)],
        input=str(input_path) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )

    gauss_info, gaussh_info, h_moments = hermite(
        input_path, exe_path=exe, rebuild=False, env_var=env_var
    )
    return HermiteResult(
        gauss_info=gauss_info,
        gaussh_info=gaussh_info,
        h_moments=h_moments,
        stdout=p.stdout,
        stderr=p.stderr,
        exe_path=exe,
    )


def gaussh_norm(Ux, ex, mu, sig, hi):
    """
    Vectorized Gauss-Hermite profile used in your fitvp_ghermite-style code.

    Parameters
    ----------
    Ux : array-like
        Coordinates where the VP is evaluated.
    ex : float or array-like
        Measurement uncertainty (same shape as Ux or scalar). This enters as
        a convolution-like broadening term in the analytic expression.
    mu, sig : float
        Mean and dispersion of the underlying Gaussian.
    hi : array-like
        Coefficients array; expects indices at least up to 6 if len(hi) >= 7.
        Uses hi[3], hi[4] and optionally hi[5],
        hi[6] as in your reference code.

    Returns
    -------
    gaussh : ndarray
        Evaluated profile (normalized per the analytic expression).
    """
    Ux = np.asarray(Ux, dtype=float)
    ex = np.asarray(ex, dtype=float)
    if ex.ndim == 0:
        ex = np.full_like(Ux, float(ex))

    mu = float(mu)
    sig = float(sig)
    if sig <= 0:
        # Avoid invalid sqrt/exp; return
        # something that makes optimizer unhappy.
        return np.full_like(Ux, np.inf)

    hi = np.asarray(hi, dtype=float)
    # Ensure indexing exists
    if hi.size < 7:
        hi_pad = np.zeros(7, dtype=float)
        hi_pad[: hi.size] = hi
        hi = hi_pad

    tau = (Ux - mu) / sig

    # NOTE: This is your expression verbatim (with numpy exp).
    if hi.size < 6:
        gaussh = (
            np.sqrt(3)
            * hi[4]
            * (
                3 * ex**8
                + 6 * ex**4 * sig**4 * (-1 + 2 * tau**2)
                + sig**8 * (3 + 4 * tau**2 * (-3 + tau**2))
            )
            + 2
            * np.sqrt(2)
            * (ex**2 + sig**2)
            * (
                3 * ex**6
                + 9 * ex**2 * sig**4
                + 3 * ex**4 * sig**2 * (3 + np.sqrt(3) * hi[3] * tau)
                + sig**6 * (3 + np.sqrt(3) * hi[3] * tau * (-3 + 2 * tau**2))
            )
        ) / (
            3.0
            * np.exp((sig**2 * tau**2) / (2.0 * (ex**2 + sig**2)))
            * (4 + np.sqrt(6) * hi[4])
            * np.sqrt(np.pi)
            * (ex**2 + sig**2) ** 4.5
        )
    else:
        gaussh = (
            15
            * (2 * np.sqrt(3) * hi[4] + np.sqrt(2) * (4 + np.sqrt(5) * hi[6]))
            * ex**12
            + 30
            * ex**10
            * sig**2
            * (
                12 * np.sqrt(2)
                + 2 * np.sqrt(3) * hi[4]
                + 2 * np.sqrt(6) * hi[3] * tau
                + np.sqrt(30) * hi[5] * tau
            )
            + 20
            * ex**6
            * sig**6
            * (
                60 * np.sqrt(2)
                - 3 * np.sqrt(30) * hi[5] * tau
                + 2 * np.sqrt(30) * hi[5] * tau**3
                + 2 * np.sqrt(6) * hi[3] * tau * (3 + tau**2)
                + 6 * np.sqrt(3) * hi[4] * (-1 + 2 * tau**2)
            )
            + 15
            * ex**8
            * sig**4
            * (
                2
                * np.sqrt(2)
                * (
                    30
                    + 6 * np.sqrt(3) * hi[3] * tau
                    + np.sqrt(
                        15,
                    )
                    * hi[5]
                    * tau
                )
                + 3 * np.sqrt(10) * hi[6] * (-1 + 2 * tau**2)
                + 2 * np.sqrt(3) * hi[4] * (-1 + 4 * tau**2)
            )
            + 2
            * ex**2
            * sig**10
            * (
                180 * np.sqrt(2)
                + 15 * np.sqrt(30) * hi[5] * tau
                - 20 * np.sqrt(30) * hi[5] * tau**3
                + 4 * np.sqrt(30) * hi[5] * tau**5
                + 30 * np.sqrt(6) * hi[3] * tau * (-3 + 2 * tau**2)
                + 10 * np.sqrt(3) * hi[4] * (3 - 12 * tau**2 + 4 * tau**4)
            )
            + 5
            * ex**4
            * sig**8
            * (
                2 * np.sqrt(3) * hi[4] * (-3 + 4 * tau**4)
                + 3 * np.sqrt(10) * hi[6] * (3 - 12 * tau**2 + 4 * tau**4)
                + 4
                * np.sqrt(2)
                * (
                    45
                    + 6 * np.sqrt(3) * hi[3] * tau * (-1 + tau**2)
                    + np.sqrt(15) * hi[5] * tau * (-3 + 2 * tau**2)
                )
            )
            + sig**12
            * (
                10 * np.sqrt(3) * hi[4] * (3 - 12 * tau**2 + 4 * tau**4)
                + np.sqrt(2)
                * (
                    60
                    + 30 * np.sqrt(15) * hi[5] * tau
                    - 40 * np.sqrt(15) * hi[5] * tau**3
                    + 8 * np.sqrt(15) * hi[5] * tau**5
                    + 20 * np.sqrt(3) * hi[3] * tau * (-3 + 2 * tau**2)
                    + np.sqrt(5)
                    * hi[6]
                    * (-15 + 90 * tau**2 - 60 * tau**4 + 8 * tau**6)
                )
            )
        ) / (
            30.0
            * np.exp((sig**2 * tau**2) / (2.0 * (ex**2 + sig**2)))
            * (4 + np.sqrt(6) * hi[4] + np.sqrt(5) * hi[6])
            * np.sqrt(np.pi)
            * (ex**2 + sig**2) ** 6.5
        )

    return gaussh


def _trap_weights(x: np.ndarray) -> np.ndarray:
    w = np.empty_like(x, dtype=float)
    w[0] = x[1] - x[0]
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[-1] = x[-1] - x[-2]
    return w


def _raw_moments(v: np.ndarray, vp: np.ndarray) -> Dict[str, float]:
    w = _trap_weights(v)
    sum0 = float(np.sum(vp * w))
    mean = float(np.sum(v * vp * w) / sum0)
    var = float(np.sum((v**2) * vp * w) / sum0 - mean**2)
    disp = float(np.sqrt(max(var, 0.0)))
    return {"norm": sum0, "mean": mean, "dispersion": disp}


def fit_gaussh_vp(
    input_path: Union[str, Path],
    *,
    ex: Union[float, np.ndarray] = 0.0,
    fit_orders: Tuple[int, ...] = (3, 4),
    max_nfev: int = 5000,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
    """
    Fit the 'fitvp_ghermite' analytic Gauss-Hermite model
    (gaussh_norm) to a VP file.

    Parameters
    ----------
    input_path : str|Path
        Two-column file with header allowed: v, VP(v)
    ex : float|array
        Uncertainty term used by gaussh_norm. If scalar, broadcast to v.
    fit_orders : tuple
        Which coefficients to fit among (3,4,5,6). Common: (3,4) or (3,4,5,6).
    max_nfev : int
        Max function evals for least_squares.

    Returns
    -------
    (gauss_info, gaussh_info, h_moments)
        gauss_info: raw-moment Gaussian summary (norm/mean/dispersion)
        gaussh_info: fitted Gaussian summary (norm/mean/dispersion)
        h_moments: dict h0..h10, with fitted h3/h4(/h5/h6) and others set to 0,
                   plus h0 forced to 1.
    """
    if least_squares is None:
        raise RuntimeError(
            "SciPy is required for fit_gaussh_vp ",
            "(scipy.optimize.least_squares). ",
            "Install scipy in your environment.",
        )

    input_path = Path(input_path)
    data = np.loadtxt(input_path, skiprows=1)
    v = np.asarray(data[:, 0], dtype=float)
    vp = np.asarray(data[:, 1], dtype=float)

    gauss_info = _raw_moments(v, vp)

    # Initial guesses
    A0 = gauss_info["norm"] / (
        np.sqrt(2.0 * np.pi) * max(gauss_info["dispersion"], 1e-6)
    )
    mu0 = gauss_info["mean"]
    sig0 = max(gauss_info["dispersion"], 1e-3)

    # Pack parameters: [logA, mu, logSig, h3, h4, (h5), (h6)]
    fit_orders = tuple(fit_orders)
    for o in fit_orders:
        if o not in (3, 4, 5, 6):
            raise ValueError(
                f"fit_orders contains unsupported order {o};",
                " use subset of (3,4,5,6).",
            )

    def pack(logA, mu, logSig, h3=0.0, h4=0.0, h5=0.0, h6=0.0):
        vec = [logA, mu, logSig]
        if 3 in fit_orders:
            vec.append(h3)
        if 4 in fit_orders:
            vec.append(h4)
        if 5 in fit_orders:
            vec.append(h5)
        if 6 in fit_orders:
            vec.append(h6)
        return np.array(vec, dtype=float)

    def unpack(p):
        logA = p[0]
        mu = p[1]
        logSig = p[2]
        idx = 3
        h3 = h4 = h5 = h6 = 0.0
        if 3 in fit_orders:
            h3 = p[idx]
            idx += 1
        if 4 in fit_orders:
            h4 = p[idx]
            idx += 1
        if 5 in fit_orders:
            h5 = p[idx]
            idx += 1
        if 6 in fit_orders:
            h6 = p[idx]
            idx += 1
        return logA, mu, logSig, h3, h4, h5, h6

    p0 = pack(np.log(max(A0, 1e-12)), mu0, np.log(sig0), 0.0, 0.0, 0.0, 0.0)

    # Bounds: must satisfy lb < ub for every parameter
    lo = np.full_like(p0, -np.inf, dtype=float)
    hi = np.full_like(p0, +np.inf, dtype=float)

    # logA bounds
    lo[0] = np.log(1e-20)
    hi[0] = np.log(1e20)

    # mu bounds: allow it to move within a generous window based on data extent
    vmin = float(np.min(v))
    vmax = float(np.max(v))
    vrng = max(vmax - vmin, 1e-6)
    lo[1] = vmin - 0.25 * vrng
    hi[1] = vmax + 0.25 * vrng

    # logSig bounds (sigma positive)
    lo[2] = np.log(1e-6)
    hi[2] = np.log(1e6)

    # coefficient bounds
    for k in range(3, len(p0)):
        lo[k] = -0.5
        hi[k] = 0.5

    wtrap = _trap_weights(v)

    def residuals(p):
        logA, mu, logSig, h3, h4, h5, h6 = unpack(p)
        A = np.exp(logA)
        sig = np.exp(logSig)

        hi_vec = np.zeros(7, dtype=float)
        hi_vec[0] = 1.0
        hi_vec[3] = h3
        hi_vec[4] = h4
        hi_vec[5] = h5
        hi_vec[6] = h6

        model = A * gaussh_norm(v, ex, mu, sig, hi_vec)

        # Weighted residuals: trapz weights act like "dx"; use sqrt(dx) scaling
        # to approximate continuous L2 norm.
        return (vp - model) * np.sqrt(wtrap)

    res = least_squares(
        residuals,
        p0,
        bounds=(lo, hi),
        max_nfev=max_nfev,
        method="trf",
    )

    logA, mu, logSig, h3, h4, h5, h6 = unpack(res.x)
    A = float(np.exp(logA))
    sig = float(np.exp(logSig))

    # Convert amplitude A (peak-like) back to a "norm"
    # comparable to Gaussian integral
    # for reporting consistency. For a pure Gaussian,
    # integral = A * sqrt(2*pi)*sig.
    # For GH this is approximate; still useful as a scale indicator.
    norm_approx = float(A * np.sqrt(2.0 * np.pi) * sig)

    gaussh_info = {
        "norm": norm_approx,
        "mean": float(mu),
        "dispersion": float(sig),
    }

    h_moments = {f"h{i}": 0.0 for i in range(11)}
    h_moments["h0"] = 1.0
    h_moments["h3"] = float(h3)
    h_moments["h4"] = float(h4)
    h_moments["h5"] = float(h5) if 5 in fit_orders else 0.0
    h_moments["h6"] = float(h6) if 6 in fit_orders else 0.0

    return gauss_info, gaussh_info, h_moments


def eval_gaussh_fit(
    v: np.ndarray,
    *,
    ex: Union[float, np.ndarray],
    gaussh_info: Dict[str, float],
    h_moments: Dict[str, float],
) -> np.ndarray:
    """
    Evaluate the fitted GH model (matching fit_gaussh_vp's parameterization)
    on a grid v.

    Uses amplitude A = norm / (sqrt(2*pi)*sigma) with the fitted sigma.
    """
    mu = gaussh_info["mean"]
    sig = gaussh_info["dispersion"]
    A = gaussh_info["norm"] / (np.sqrt(2.0 * np.pi) * sig)

    hi_vec = np.zeros(7, dtype=float)
    hi_vec[0] = 1.0
    for k in (3, 4, 5, 6):
        key = f"h{k}"
        if key in h_moments:
            hi_vec[k] = float(h_moments[key])

    return A * gaussh_norm(v, ex, mu, sig, hi_vec)
