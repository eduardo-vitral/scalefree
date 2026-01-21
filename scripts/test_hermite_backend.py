#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add repo root to sys.path for local testing (no pip install needed)
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

# Fortran-backed fitter (calls the compiled fitvp executable)
from scalefree import hermite
from scalefree import hermite_fit


def write_synthetic_vp(path: Path, *, n: int = 500):
    rng = np.random.default_rng(12345)

    # "Truth" used to generate synthetic VP
    mu = 12.0
    sig = 30.0
    ex = 3.0  # keep nonzero to exercise your definition

    hi_true = np.zeros(7)
    hi_true[0] = 1.0
    hi_true[3] = 0.06  # h3
    hi_true[4] = -0.04  # h4
    # h5,h6 kept at 0

    v = np.linspace(mu - 5 * sig, mu + 5 * sig, n)

    # Scale (overall normalization). Here we treat gaussh_norm as a shape and
    # multiply by A to get VP(v).
    # This matches how you wrote the synthetic generator.
    A = 1.7
    vp_clean = A * hermite.gaussh_norm(v, ex, mu, sig, hi_true)

    # Add small noise; keep positive
    noise = 0.01 * np.max(vp_clean) * rng.normal(size=n)
    vp = np.clip(vp_clean + noise, 1e-14, None)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\tv\tVP(v)\n")
        for vv, ff in zip(v, vp):
            f.write(f"{vv: .10e}\t{ff: .10e}\n")

    return mu, sig, ex, hi_true, v, vp, vp_clean


def trap_weights(x: np.ndarray) -> np.ndarray:
    w = np.empty_like(x, dtype=float)
    w[0] = x[1] - x[0]
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[-1] = x[-1] - x[-2]
    return w


def weighted_metrics(
    y: np.ndarray, yhat: np.ndarray, x: np.ndarray, label: str
) -> None:
    w = trap_weights(x)
    resid = y - yhat
    sse_w = float(np.sum((resid**2) * w))
    mean_w = float(np.sum(y * w) / np.sum(w))
    sst_w = float(np.sum(((y - mean_w) ** 2) * w))
    r2_w = 1.0 - (sse_w / sst_w) if sst_w > 0 else np.nan
    rmse_w = float(np.sqrt(np.sum((resid**2) * w) / np.sum(w)))
    nrmse = rmse_w / (np.max(y) - np.min(y) + 1e-30)

    print(f"\n--- Goodness-of-fit ({label}) ---")
    print(f"Weighted SSE:   {sse_w:.6e}")
    print(f"Weighted RMSE:  {rmse_w:.6e}")
    print(f"NRMSE (range):  {nrmse:.6e}")
    print(f"Weighted R^2:   {r2_w:.6e}")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    vp_path = script_dir / "tmp_vp_gaussh_test.dat"

    mu_t, sig_t, ex, hi_t, v, vp, vp_clean = write_synthetic_vp(vp_path, n=500)

    # ------------------------------------------------------------------
    # FIT: Fortran backend ONLY (Python just calls it and parses stdout)
    # ------------------------------------------------------------------
    gauss_info, gaussh_info, h_moments = hermite_fit(vp_path)

    # ------------------------------------------------------------------
    # ANALYTIC LINES for plotting (Python evaluation ONLY)
    # ------------------------------------------------------------------

    # 1) Gaussian line from Fortran "gaussh_info"
    # (its best-fit Gaussian reference)
    norm_f = gaussh_info["norm"]
    mu_f = gaussh_info["mean"]
    sig_f = gaussh_info["dispersion"]
    y_gauss = (
        norm_f
        / (np.sqrt(2.0 * np.pi) * sig_f)
        * np.exp(-0.5 * ((v - mu_f) / sig_f) ** 2)
    )

    # 2) "Gauss-Hermite" analytic line using YOUR gaussh_norm definition
    #    We build hi from the Fortran-returned h_moments
    # and scale by the Fortran norm.
    hi_fit = np.zeros(7)
    hi_fit[0] = 1.0  # enforce
    hi_fit[3] = float(h_moments.get("h3", 0.0))
    hi_fit[4] = float(h_moments.get("h4", 0.0))
    hi_fit[5] = float(h_moments.get("h5", 0.0))
    hi_fit[6] = float(h_moments.get("h6", 0.0))

    # Key modeling choice:
    # We assume gaussh_norm(...)
    # returns a unit-normalized SHAPE (integral ~ 1),
    # so the fitted "norm" from Fortran provides the overall scale.
    y_gh = norm_f * hermite.gaussh_norm(v, ex, mu_f, sig_f, hi_fit)

    # ------------------------------------------------------------------
    # REPORT
    # ------------------------------------------------------------------
    print("\n--- Truth (synthetic generator) ---")
    print(
        f"mu={mu_t:.6f}",
        f" sig={sig_t:.6f}",
        f" ex={ex:.6f}",
        f" h3={hi_t[3]:.6f}",
        f" h4={hi_t[4]:.6f}",
    )

    print("\n--- Fortran outputs ---")
    print("Gaussian raw moments (gauss_info):")
    print(
        f"  norm={gauss_info['norm']:.6e}",
        f" mean={gauss_info['mean']:.6f}",
        f" disp={gauss_info['dispersion']:.6f}",
    )
    print("Gaussian best-fit reference (gaussh_info):")
    print(
        f"  norm={gaussh_info['norm']:.6e}",
        f" mean={gaussh_info['mean']:.6f}",
        f" disp={gaussh_info['dispersion']:.6f}",
    )

    print(
        "\nGauss-Hermite coefficients used for",
        " gaussh_norm (taken from Fortran stdout):",
    )
    print(
        f"  h3={hi_fit[3]: .6e}",
        f"  h4={hi_fit[4]: .6e}",
        f"  h5={hi_fit[5]: .6e}",
        f"  h6={hi_fit[6]: .6e}",
    )

    # Metrics
    weighted_metrics(
        vp,
        y_gauss,
        v,
        "Gaussian (from Fortran gaussh_info)",
    )
    weighted_metrics(
        vp,
        y_gh,
        v,
        "Gauss-Hermite analytic (gaussh_norm, params from Fortran)",
    )

    # ------------------------------------------------------------------
    # PLOTS
    # ------------------------------------------------------------------
    plt.figure()
    plt.plot(v, vp, label="VP input (noisy)")
    plt.plot(v, vp_clean, label="VP truth (clean)")
    plt.plot(v, y_gauss, label="Gaussian fit (Fortran)")
    plt.plot(v, y_gh, label="Gauss-Hermite line (gaussh_norm; Fortran params)")
    plt.xlabel("v")
    plt.ylabel("VP(v)")
    plt.title(
        "VP: input vs Fortran Gaussian vs analytic GH (using Fortran outputs)",
    )
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Residuals
    plt.figure()
    plt.plot(v, vp - y_gauss, label="Residual (VP - Gaussian)")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("v")
    plt.ylabel("Residual")
    plt.title("Residuals: Gaussian")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.plot(v, vp - y_gh, label="Residual (VP - GH analytic)")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("v")
    plt.ylabel("Residual")
    plt.title("Residuals: Gauss-Hermite analytic (gaussh_norm)")
    plt.legend()
    plt.tight_layout()
    plt.show()
    print("\n\nFinished.")


if __name__ == "__main__":
    main()
