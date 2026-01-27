# Hermite backend playground (`hermite_fit`) — fitting Gauss–Hermite VPs

This note is a hands-on companion to `scripts/test_hermite_backend.py`.

It shows how to:
- generate a **synthetic velocity profile** (VP) from `scalefree.hermite.gaussh_norm`
- run the **Fortran-backed fitter** via `scalefree.hermite_fit(...)`
- reconstruct an analytic Gauss–Hermite VP line in Python from the **Fortran outputs**
- quantify the goodness-of-fit (simple weighted metrics) and visualize residuals

The intent is to give you a copy/paste “playground” for validating and experimenting with the Hermite backend.

---

## 1) What `hermite_fit` does

`scalefree.hermite_fit(path_to_vp_file)` is an alias for `scalefree.hermite.hermite(...)`.

It runs the **Fortran** fitter (`fitvp_stdout`) on a 2-column VP file and parses stdout into:

- `gauss_info` : Gaussian parameters computed from **raw VP moments**  
  Keys: `norm`, `mean`, `dispersion`

- `gaussh_info` : **best-fit Gaussian reference** used for the Gauss–Hermite expansion  
  Keys: `norm`, `mean`, `dispersion`

- `h_moments` : fitted Gauss–Hermite coefficients  
  Keys typically include `h0` … `h10` (you’ll usually care about `h3`–`h6`)

---

## 2) Backend discovery / build behavior (important)

By default, the package will try to find or build a cached executable automatically:

- If you set the environment variable `SCALEFREE_FITVP_EXE`, that path is used.
- Otherwise it uses a cached path under the user cache directory.
- If the executable does not exist (or you pass `rebuild=True`), it will try to compile
  from `fortran_src/fitvp_stdout.f` using `gfortran`.

If you cannot compile on a given machine, point `SCALEFREE_FITVP_EXE` to a precompiled binary.

---

## 3) VP file format

The fitter expects a text file with **two columns**:

1. `v` (velocity grid)
2. `VP(v)` (velocity profile value at that grid)

A header is allowed (the test uses a one-line header). Example:

```
    v    VP(v)
 -1.0000e+02   1.2345e-03
 ...
```

---

## 4) Minimal end-to-end example (inspired by `scripts/test_hermite_backend.py`)

### 4.1 Generate a synthetic VP file

This creates a noisy Gauss–Hermite VP using your Python analytic definition,
then writes it to disk in the format the Fortran backend expects.

```python
from pathlib import Path
import numpy as np

from scalefree import hermite

def write_synthetic_vp(path: Path, *, n: int = 500):
    rng = np.random.default_rng(12345)

    # “Truth” used to generate the synthetic VP
    mu = 12.0
    sig = 30.0
    ex = 3.0          # keep nonzero to exercise your ex-handling
    A = 1.7           # overall amplitude / normalization scale

    # Gauss–Hermite coefficients array (indices match the Fortran convention used here)
    hi_true = np.zeros(7)
    hi_true[0] = 1.0
    hi_true[3] = 0.06   # h3
    hi_true[4] = -0.04  # h4
    # h5,h6 left at 0

    v = np.linspace(mu - 5 * sig, mu + 5 * sig, n)

    # gaussh_norm is treated as a “shape”; we scale by A to get VP(v)
    vp_clean = A * hermite.gaussh_norm(v, ex, mu, sig, hi_true)

    # Add small noise; clip to stay positive
    noise = 0.01 * np.max(vp_clean) * rng.normal(size=n)
    vp = np.clip(vp_clean + noise, 1e-14, None)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\tv\tVP(v)\n")
        for vv, ff in zip(v, vp):
            f.write(f"{vv: .10e}\t{ff: .10e}\n")

    return mu, sig, ex, hi_true, v, vp, vp_clean
```

### 4.2 Run the Fortran fitter and reconstruct analytic curves

```python
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from scalefree import hermite
from scalefree import hermite_fit

vp_path = Path("tmp_vp_gaussh_test.dat")
mu_t, sig_t, ex, hi_t, v, vp, vp_clean = write_synthetic_vp(vp_path, n=500)

# Fortran backend ONLY (Python calls it + parses stdout)
gauss_info, gaussh_info, h_moments = hermite_fit(vp_path)

# 1) Gaussian curve from the Fortran "gaussh_info"
norm_f = gaussh_info["norm"]
mu_f   = gaussh_info["mean"]
sig_f  = gaussh_info["dispersion"]

y_gauss = (
    norm_f
    / (np.sqrt(2.0 * np.pi) * sig_f)
    * np.exp(-0.5 * ((v - mu_f) / sig_f) ** 2)
)

# 2) Gauss–Hermite curve using YOUR gaussh_norm definition + Fortran outputs
hi_fit = np.zeros(7)
hi_fit[0] = 1.0
hi_fit[3] = float(h_moments.get("h3", 0.0))
hi_fit[4] = float(h_moments.get("h4", 0.0))
hi_fit[5] = float(h_moments.get("h5", 0.0))
hi_fit[6] = float(h_moments.get("h6", 0.0))

# Modeling assumption (same as the test script):
# gaussh_norm returns a unit-normalized SHAPE (integral ~ 1),
# so the Fortran “norm” provides the overall VP scale.
y_gh = norm_f * hermite.gaussh_norm(v, ex, mu_f, sig_f, hi_fit)

print("Truth:", mu_t, sig_t, "h3=", hi_t[3], "h4=", hi_t[4])
print("Fortran gauss_info:", gauss_info)
print("Fortran gaussh_info:", gaussh_info)
print("Fortran h_moments (subset):", {k: h_moments[k] for k in ["h3","h4","h5","h6"] if k in h_moments})
```

### 4.3 Plot and residuals

```python
plt.figure()
plt.plot(v, vp, label="VP input (noisy)")
plt.plot(v, vp_clean, label="VP truth (clean)")
plt.plot(v, y_gauss, label="Gaussian fit (Fortran)")
plt.plot(v, y_gh, label="Gauss–Hermite line (gaussh_norm; Fortran params)")
plt.xlabel("v")
plt.ylabel("VP(v)")
plt.title("VP: input vs Fortran Gaussian vs analytic GH (using Fortran outputs)")
plt.legend()
plt.tight_layout()
plt.show()

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
plt.title("Residuals: Gauss–Hermite analytic (gaussh_norm)")
plt.legend()
plt.tight_layout()
plt.show()
```

---

## 6) Common pitfalls / interpretation notes

### 6.1 Normalization (`norm`) conventions
The test assumes:
- `hermite.gaussh_norm(...)` returns a *shape* with integral close to 1 (or at least
  something that becomes a physically meaningful VP once multiplied by an overall scale).
- The Fortran-reported `gaussh_info["norm"]` is the appropriate overall scale factor.

If you change the definition of `gaussh_norm` (or what “norm” means in the backend),
you **must** revisit how you reconstruct `y_gh`.

### 6.2 Role of `ex`
`ex` is carried through the analytic call and is also present in the synthetic generation
specifically to validate that your definition responds sensibly when `ex != 0`.

A good sanity test is to sweep `ex` (e.g. 0, 1, 3, 10) and verify the recovered Gaussian
dispersion behaves as expected (broadening).

### 6.3 Higher-order moments
By default, most workflows focus on `h3` and `h4`. If you start fitting `h5`/`h6`,
check:
- the stability with respect to noise level and grid resolution (`n`)
- whether the fitted coefficients remain within a physically sensible regime

---

## 7) “Playground exercises”

1) **Noise sensitivity**  
Increase the injected noise and track how (h3,h4) drift and how the residuals change.

2) **Grid resolution**  
Repeat with `n = 200, 500, 2000` and check stability of recovered parameters.

3) **Normalization sanity**  
Change the synthetic amplitude `A` and confirm it mostly affects `norm`, not `mean/dispersion`.

4) **Broadening sweep**  
Vary `ex` and confirm the best-fit Gaussian reference updates consistently.

---

## 8) Troubleshooting

- **`HermiteBackendError` / executable not found**  
  Provide an explicit `exe_path=...` or set `SCALEFREE_FITVP_EXE`. If the package attempts
  to compile and fails, confirm `gfortran` is available.

- **Backend compiles but produces nonsense**  
  Inspect the VP file: ensure the second column is positive and the grid is sensible
  (monotonic, reasonable extent).

- **Plot looks “off” but fit numbers look plausible**  
  Double-check the normalization convention used to rebuild `y_gh`.
