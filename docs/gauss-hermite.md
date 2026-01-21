# Gauss–Hermite tools (`scalefree.hermite`)

This module provides two complementary approaches:

1. **Fortran-backed GH fitter** (fast, consistent with legacy workflows)
2. **Pure Python helpers** (`gaussh_norm`, and an optional SciPy fitter)

---

## 1) Fortran-backed fitter: `hermite.hermite(path)`

Use when you have a VP file with two columns:
- velocity `v`
- profile `VP(v)`
(An optional one-line header is tolerated.)

```python
from scalefree import hermite

gauss_info, gaussh_info, h = hermite.hermite("vp.dat")
print(gaussh_info)
print(h["h3"], h["h4"])
```

### Backend resolution

If you do not pass an executable, the function uses:
- `SCALEFREE_FITVP_EXE` if set, otherwise
- a cached executable under the user cache directory,
- otherwise it auto-compiles from packaged Fortran (requires `gfortran`).

---

## 2) Analytic GH profile: `gaussh_norm(...)`

`gaussh_norm(Ux, ex, mu, sig, hi)` evaluates a vectorized Gauss–Hermite profile.

```python
import numpy as np
from scalefree import hermite

v = np.linspace(-10, 10, 200)
ex = 0.0
mu, sig = 0.0, 2.0

hi = np.zeros(7)
hi[3] = 0.05  # h3
hi[4] = 0.02  # h4

vp = hermite.gaussh_norm(v, ex, mu, sig, hi)
```

---

## 3) Optional SciPy fit: `fit_gaussh_vp(...)`

If SciPy is available, you can fit the analytic GH model directly in Python:

```python
from scalefree import hermite

gauss_info, gaussh_info, h = hermite.fit_gaussh_vp(
    "vp.dat",
    ex=0.0,
    fit_orders=(3, 4),  # or (3, 4, 5, 6)
)
```

If SciPy is not installed, this function will raise a clear error.
