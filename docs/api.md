# API overview

This is a friendly overview of the public API most users need.

> Tip: For advanced control (e.g., VP algorithm knobs), use `ScaleFreeRunner.vprofile()` directly.

---

## `scalefree.vprofile(...)`

**Use when:** you want intrinsic/projected moments and (optionally) velocity profiles and Gauss–Hermite summaries.

Returns a `ScaleFreeResult` with:
- `blocks`: parsed structured blocks as dictionaries and NumPy arrays
- `raw_text`: the structured text that was parsed
- `stdout`, `stderr`: full process outputs (useful for debugging)

Minimal pattern:

```python
from scalefree import vprofile

res = vprofile(
    potential="kepler",
    gamma=2.0, q=0.9, df=1, beta=0.0, s=0.5, t=0.0,
    inclination=60.0, xi=0.0, theta=0.0,
    usevp=True,
    algorithm=3,
)
```

See `docs/vprofile.md` for parameter meanings and output blocks.

---

## `scalefree.ScaleFreeRunner`

**Use when:** you want to reuse a resolved backend executable across many calls, or you want finer control.

```python
from scalefree import ScaleFreeRunner

runner = ScaleFreeRunner(exe_path=None)   # auto-resolve/build
res = runner.vprofile(..., usevp=True)
```

---

## `scalefree.hermite` (module)

### `hermite.hermite(path, ...) -> (gauss_info, gaussh_info, h_moments)`

**Use when:** you have a VP file (two columns: `v`, `VP(v)`) and you want Fortran-backed Gauss–Hermite moments.

Returns three dictionaries:
- `gauss_info`: Gaussian from raw VP moments (norm/mean/dispersion)
- `gaussh_info`: best-fit Gaussian used as GH reference
- `h_moments`: `h0..h10`

### `hermite.gaussh_norm(Ux, ex, mu, sig, hi)`

**Use when:** you want to evaluate the analytic GH profile used by the helper fitter.

### `hermite.fit_gaussh_vp(path, ...)`

**Use when:** you want a pure Python GH fit (SciPy required).

See `docs/gauss-hermite.md`.

---

## `scalefree.mock(...) -> ndarray(N, 6)`

**Use when:** you want a fast synthetic 6D sample `(x,y,z,vx,vy,vz)` from a scale-free model.

Internally:
- samples positions from the scale-free density
- bins stars by projected-plane angle
- runs `vprofile()` per occupied bin to get GH parameters
- samples velocities using Sanders–Evans style GH PDFs (via BALRoGO or fallback import)

See `docs/mock-generator.md`.
