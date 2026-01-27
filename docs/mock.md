# mock.py — generating and exploring intrinsic ScaleFree mock datasets

This note is a practical companion to `docs/vmoments.md` and the validation workflow in `scripts/check_mock_end2end.py`.

It shows how to:
- generate intrinsic 6D mock phase-space samples with `scalefree.mock.mock`
- sanity-check the **density slope** and **flattening**
- compare **velocity PDFs** against the intrinsic VP diagnostics returned by `ScaleFreeRunner.vprofile(...)`

The goal is to give you a *playground* you can copy/paste into a notebook or a small script.

---

## 1) What the mock generator produces

`scalefree.mock.mock(...)` returns an `ndarray` with shape `(N, 6)`:

```
[x, y, z, vx, vy, vz]
```

All coordinates are **Cartesian in the model frame** (i.e. **intrinsic**; not LOS/POSR/POST).

Internally, the generator:
1. samples **positions** from the *volume density* ρ(R,z) with intrinsic flattening `q`
2. bins the meridional angle θ and calls `ScaleFreeRunner.vprofile(...)` in intrinsic mode per θ-bin
3. uses the intrinsic VP Gaussian-fit parameters and Gauss–Hermite moments `(gauss_V, gauss_sig, h3, h4)` to sample velocities (via BALRoGO)
4. converts spherical velocities to Cartesian and returns the final 6D sample

---

## 2) Quickstart (minimal working example)

### 2.1 Pick a model

A convenient default (also used in `check_mock_end2end.py`) is:

```python
model = dict(
    potential=2,   # 1=Kepler, 2=Logarithmic
    gamma=2.0,     # density slope
    q=1.0,         # intrinsic axis ratio (1=spherical, <1=oblate)
    df=1,          # DF family
    beta=0.0,      # anisotropy parameter
    s=0.5, t=0.0,  # odd-part parameters
)
```

### 2.2 Generate a dataset

```python
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from scalefree import mock

xyzv = mock(
    **model,
    nsamples=20000,
    nbins=24,
    seed=123,
    rin=0.5,
    rout=50.0,
    maxmom=10,
    debug=False,
    # exe_path=Path(".../fortran_src/scalefree_intrvp.e"),  # optional override
)

print(xyzv.shape)     # (20000, 6)
print(xyzv[:3])
```

### 2.3 Save / reload

```python
np.save("mock_xyzv.npy", xyzv)
# or:
np.savetxt("mock_xyzv.csv", xyzv, delimiter=",", header="x,y,z,vx,vy,vz", comments="")

xyzv2 = np.load("mock_xyzv.npy")
```

---

## 3) Density sanity check (ρ ∝ r^{-γ})

A simple way to verify you are sampling the **volume density** correctly is to compute
a binned shell density estimate and compare it to `r^{-gamma}` up to a scale factor. Notice that the example below is suited only for the case of a spherical model (q = 1), for better visualization.

```python
import numpy as np
import math

def volume_density_powerlaw_shape(gamma, r):
    return np.power(r, -gamma)

def shell_volume_density(r, nbins=25):
    r = np.asarray(r, float)
    r = r[np.isfinite(r) & (r > 0)]
    if r.size == 0:
        raise ValueError("No valid radii.")

    edges = np.logspace(np.log10(r.min()), np.log10(r.max()), nbins + 1)
    counts, _ = np.histogram(r, bins=edges)

    vol = (4.0 / 3.0) * math.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho_hat = counts / vol
    rmid = np.sqrt(edges[1:] * edges[:-1])

    m = counts > 0
    return rmid[m], rho_hat[m]

r3d = np.sqrt(xyzv[:,0]**2 + xyzv[:,1]**2 + xyzv[:,2]**2)
rmid, rho_hat = shell_volume_density(r3d, nbins=25)

rho_th = volume_density_powerlaw_shape(model["gamma"], rmid)
k0 = len(rmid)//2                    # mid bin for scaling
rho_th_scaled = rho_th * (rho_hat[k0] / rho_th[k0])

plt.figure()
plt.loglog(rmid, rho_hat, "o", label="Mock: shell rho(r)")
plt.loglog(rmid, rho_th_scaled, "-", label=r"Analytic: $r^{-\gamma}$ (scaled)")
plt.xlabel("3D radius r")
plt.ylabel("Volume density (arb.)")
plt.legend()
plt.show()
```

If the slope is off:
- increase `nsamples`
- verify `rin/rout` are sensible for your test (too narrow ranges can look noisy)
- check you are not mixing up **volume** and **surface** density expectations

---

## 4) Flattening sanity check (intrinsic q)

A simple geometric check is to plot `(x,z)` and overplot an isodensity ellipse:

```python
import numpy as np
import matplotlib.pyplot as plt

X = xyzv[:,0]
Z = xyzv[:,2]
q0 = float(model["q"])

# "elliptical radius" in the meridional plane
Re = np.sqrt(X*X + (Z/q0)**2)
Re0 = float(np.nanpercentile(Re, 60))   # pick a contour scale

ang = np.linspace(0, 2*np.pi, 400)
ex = Re0 * np.cos(ang)
ez = q0 * Re0 * np.sin(ang)

plt.figure(figsize=(5,5))
plt.scatter(X, Z, s=1, alpha=0.2)
plt.plot(ex, ez, lw=2)
plt.gca().set_aspect("equal", adjustable="box")
plt.xlabel("x")
plt.ylabel("z")
plt.title(f"Flattening sanity check (q={q0})")
plt.show()
```

For `q<1`, the cloud should look visibly squashed in `z` relative to `x`.

---

## 5) Velocity diagnostics vs `vprofile` VP fits

### 5.1 Why this is useful

The mock *samples* velocities using per-θ-bin VP diagnostics from `vprofile` (intrinsic mode).
A good end-to-end check is therefore:
1. run `vprofile(..., kinematics="intrinsic", average=True, usevp=True, algorithm=3)`
2. take the **angle-averaged** intrinsic VP fit parameters for each component
3. compare those reference PDFs to the distribution of `(v_r, v_theta, v_phi)` measured from the mock

### 5.2 Convert Cartesian velocities → intrinsic spherical components

```python
import numpy as np

def sph_vel_from_xyzv(xyzv):
    # Returns intrinsic spherical velocity components from (x,y,z,vx,vy,vz).
    x, y, z, vx, vy, vz = xyzv.T
    r = np.sqrt(x*x + y*y + z*z)
    r = np.where(r > 0, r, np.nan)

    theta = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi = np.arctan2(y, x)

    st = np.sin(theta); ct = np.cos(theta)
    cp = np.cos(phi);   sp = np.sin(phi)

    # basis vectors (er, etheta, ephi)
    erx, ery, erz = st*cp, st*sp, ct
    etx, ety, etz = ct*cp, ct*sp, -st
    epx, epy, epz = -sp, cp, 0.0

    vr     = vx*erx + vy*ery + vz*erz
    vtheta = vx*etx + vy*ety + vz*etz
    vphi   = vx*epx + vy*epy + vz*epz

    return dict(vr=vr, vtheta=vtheta, vphi=vphi, r=r, theta=theta, phi=phi)

sph = sph_vel_from_xyzv(xyzv)
```

### 5.3 Get VP fit parameters from `vprofile(...)`

This mirrors the pattern in `check_mock_end2end.py`. The exact output structure is documented in `docs/vmoments.md`.

```python
from scalefree import ScaleFreeRunner
from pathlib import Path
import numpy as np

exe = Path("fortran_src/scalefree_intrvp.e").resolve()  # adjust to your repo location
runner = ScaleFreeRunner(exe_path=exe)

res = runner.vprofile(
    **model,
    inclination=90.0,
    xi=0.0,
    theta=0.0,
    integration=1,
    ngl_or_eps=0,
    algorithm=3,     # match mock internals
    maxmom=20,
    average=True,    # angle-averaged intrinsic VP
    kinematics="intrinsic",
    usevp=True,
    verbose_vp=0,
    output_path=None,
    timeout_s=300,
    debug_prompts=False,
)

blocks = getattr(res, "blocks", {}) or {}
vp = blocks.get("vp_intrinsic", {})
cols = list(vp.get("columns", []))
data = vp.get("data")

# vp_intrinsic has one row per component with icomp in {1,2,3}
# (v_r, v_theta, v_phi)
vp_rows = {}
if isinstance(data, np.ndarray) and data.ndim == 2 and "icomp" in cols:
    i_ic = cols.index("icomp")
    for row in data:
        ic = int(row[i_ic])
        vp_rows[ic] = {
            "gauss_V":   float(row[cols.index("gauss_V")])   if "gauss_V"   in cols else np.nan,
            "gauss_sig": float(row[cols.index("gauss_sig")]) if "gauss_sig" in cols else np.nan,
            "h3":        float(row[cols.index("h3")])        if "h3"        in cols else np.nan,
            "h4":        float(row[cols.index("h4")])        if "h4"        in cols else np.nan,
        }

print(vp_rows)
```

### 5.4 Plot mock histograms against BALRoGO reference PDFs

If you have `balrogo` installed, you can import `balrogo.dynamics`.
If not, `check_mock_end2end.py` falls back to importing a local `dynamics.py`
module (shipped in this repo) that exposes the needed PDF evaluators.

```python
import numpy as np
import matplotlib.pyplot as plt
from balrogo import dynamics

def balrogo_curve_from_fits(fits, vg):
    # Return a PDF curve from BALRoGO given (mean, sigma, h3, h4).
    ex = np.zeros_like(vg, float)

    # NOTE: BALRoGO API may differ across versions.
    # The repo's dynamics.py is intended as a compatibility shim.
    pdf = dynamics.mom_likelihood_func(fits, vg, ex, mode="curve")
    return pdf

# Grid for plotting
vg = np.linspace(-8, 8, 500)

# Map components: 1=vr, 2=vtheta, 3=vphi
comp = {
    1: ("vr", sph["vr"]),
    2: ("vtheta", sph["vtheta"]),
    3: ("vphi", sph["vphi"]),
}

for ic, (label, v) in comp.items():
    fits = np.array([
        vp_rows[ic]["gauss_V"],
        vp_rows[ic]["gauss_sig"],
        vp_rows[ic]["h3"],
        vp_rows[ic]["h4"],
    ])

    plt.figure()
    plt.hist(v[np.isfinite(v)], bins=60, density=True, alpha=0.4, label="Mock histogram")

    try:
        pdf = balrogo_curve_from_fits(fits, vg)
        plt.plot(vg, pdf, label="Reference PDF from vprofile fits")
    except Exception as e:
        plt.title(f"{label}: could not compute BALRoGO curve ({e})")
    else:
        plt.title(label)

    plt.xlabel("velocity")
    plt.ylabel("PDF")
    plt.legend()
    plt.show()
```

**If your mock histograms do not roughly match:**
- raise `nsamples`
- check you are using `algorithm=3` in `vprofile`
- keep `average=True` for the reference diagnostics (otherwise you are comparing
  against θ-dependent fits rather than an angle-averaged distribution)

---

## 6) Parameter guide (practical)

### 6.1 Sampling controls
- `nsamples`: number of stars (0 is allowed; useful for smoke tests)
- `rin`, `rout`: radial support for position sampling
- `nbins`: number of θ-bins (tradeoff: more bins = more backend calls)

### 6.2 Backend / numerical controls
- `integration`, `ngl_or_eps`: numerical backend settings forwarded into `vprofile`
- `vp_smooth_eps`: smoothing used in VP computations (if enabled in backend)
- `maxmom`: maximum moment order requested from the backend (mock forces `algorithm=3`). The more the merrier, with maxmom > 25 being usually a reliable mock.
- `exe_path`: override the backend executable path

### 6.3 Velocity sampling controls
- `nsig`: velocity grid extent (in σ units) used when sampling from VP models
- `seed`: reproducibility

---

## 7) Suggested “playground exercises”

1) **Reproducibility**  
Generate two mocks with the same `seed` and verify they match bitwise.

2) **Flattening sweep**  
Loop over `q` in `[1.0, 0.8, 0.6, 0.4]` and produce the `(x,z)` scatter plot.

3) **Density slope sweep**  
Loop over `gamma` in `[1.0, 2.0, 3.0]` and compare the shell-density slope.

4) **Velocity anisotropy sweep**  
Vary `beta` and compare the (vr, vtheta, vphi) distributions; validate against `vprofile`. Notice that some `beta` configurations lead to more problematic VP shapes, which end up in less reliable mocks (use the working example above to check it).

---

## 8) Troubleshooting

- **Backend executable not found / compilation issues**  
  Ensure your compiled intrinsic backend exists (commonly `fortran_src/scalefree_intrvp.e`)
  and pass it explicitly via `exe_path=...` or configure `SCALEFREE_EXE`.

- **Very slow runs**  
  Runtime scales with `nbins` because each occupied θ-bin triggers a `vprofile` call.
  Reduce `nbins` or `maxmom` while prototyping.

- **NaNs / non-physical VP fits**  
  If you see `gauss_sig` extremely small or NaNs in h3/h4, focus first on:
  - input physicality (parameter ranges)
  - whether the backend run produced valid VP diagnostics for that θ-bin
  - increasing `maxmom` for stability (if the backend uses higher moments)
