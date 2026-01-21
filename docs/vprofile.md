# `vprofile()` in detail

`vprofile()` is the main entry point for running the ScaleFree Fortran backend and parsing its structured output.

## What it computes

Depending on settings, the backend can output:

- **Projected moments at a point** (block: `projected_point`)
- **Projected moments averaged on a circle** (block: `projected_circle_average`)
- **VP / Gauss–Hermite summary** (block: `vp`)
- **VP tables** for each projection component (block: `vp_table`) if `usevp=True`

The Python wrapper parses these into `res.blocks`.

## Minimal example

```python
from scalefree import vprofile

res = vprofile(
    potential="logarithmic",
    gamma=2.0,
    q=0.608,
    df=1,
    beta=0.189,
    s=0.5,
    t=0.0,
    inclination=57.1,
    xi=0.0,
    theta=0.0,
    average=False,
    usevp=True,
    algorithm=3,
    maxmom=4,
)
```

## Key parameters (practical guide)

### Geometry and model
- `potential`: `"kepler"` / `"logarithmic"` / `1` / `2` / callable returning `1` or `2`
- `gamma`: power-law slope
- `q`: intrinsic axial ratio
- `df`: DF family selector
- `beta`: anisotropy parameter
- `s`, `t`: rotation/odd-part controls (see `docs/theory/rotation.md`)
- `inclination`: viewing inclination in degrees
- `theta`: intrinsic meridional-plane angle (degrees)
- `xi`: projected-plane angle (degrees)

### Output selection
- `average=False`: compute projected quantities at a point
- `average=True`: compute projected quantities averaged on a circle
- `usevp=True`: request velocity-profile output and GH summary
- `maxmom`: number of projected moments used for VP reconstruction

### Numerics and VP reconstruction
- `integration=1`: Gauss–Legendre (default)
- `integration=0`: Romberg
- `ngl_or_eps`:
  - Gauss–Legendre: number of quadrature points (0 lets Fortran choose default)
  - Romberg: fractional accuracy `epsilon`
- `algorithm` (VP shape reconstruction; see `docs/theory/vp-shapes.md`):
  - `1`: direct Vandermonde solve (no regularization)
  - `2`: regularized with user parameter (advanced knob)
  - `3`: automatic regularization / smoothness control (recommended)

## Understanding output blocks

### `projected_point` and `projected_circle_average`

Each is stored as:

```python
blk = res.blocks["projected_point"]  # or "projected_circle_average"
cols = blk["columns"]
data = blk["data"]  # NumPy array
```

Typical columns look like:
- `iproj` plus quantities like `rho_p`, `v1`, `v2`, `v3`, `v4`

### `vp` (VP / GH summary)

This block includes a summary per `iproj` and is the most convenient entry for many workflows:

```python
vp = res.blocks["vp"]
vp_by_iproj = vp["by_iproj"]
print(vp_by_iproj[1])  # dict of columns -> values
```

Common fields include Gaussian and GH summaries such as mean, dispersion, and `h3`, `h4`, ...

### `vp_table` (VP arrays)

If `usevp=True`, you may also receive VP tables for each projection component:

```python
tbl = res.blocks["vp_table"][1]  # iproj=1
v = tbl["data"][:, 0]
vpv = tbl["data"][:, 1]
```

## Reproducibility tip

If you want stable file-based outputs for regression tests, pass `output_path="out.txt"`.
By default, the wrapper parses structured STDOUT instead of leaving files behind.
