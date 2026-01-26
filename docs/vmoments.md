# vmoments.py — computing intrinsic and projected velocity moments (and VP diagnostics)

`scalefree.vmoments` provides a prompt-driven Python wrapper around the ScaleFree Fortran backend. The public entry
point is the `ScaleFreeRunner.vprofile()` method, which can compute:

- **Intrinsic moments** at a point (`iwhat=0`) or as a **mass-weighted shell average** (`iwhat=2`).
- **Projected moments** at a point on the sky (`iwhat=1`) or as a **circle average** (`iwhat=3`).
- Optional **velocity profile (VP) diagnostics** (Gaussian + Gauss–Hermite) and sampled VP tables.

This document explains what is computed, how to request it, and how to interpret the returned arrays and indices.

---

## 1. The main API

### 1.1 `ScaleFreeRunner`

```python
from pathlib import Path
from scalefree.vmoments import ScaleFreeRunner

runner = ScaleFreeRunner()
```

`exe_path` is the path to the Fortran executable. If it does not exist and `gfortran` is available, the package
can compile the backend. You may also use the `SCALEFREE_EXE` environment variable (see the package docs).

### 1.2 `runner.vprofile(...)`

This method runs the interactive Fortran code non-interactively (it answers prompts internally), then parses the output into a structured object:

```python
res = runner.vprofile(
    potential=2,        # 1=Kepler, 2=Logarithmic (see _potential_code)
    gamma=2.0,
    q=0.608,
    df=1,               # distribution function family (Case I / Case II in the Fortran prompts)
    beta=0.189,         # anisotropy parameter (1 - ⟨v_tan²⟩/⟨v_rad²⟩)
    s=0.5, t=0.0,       # odd part parameters
    inclination=57.1,   # degrees
    xi=0.0,             # degrees, projected-plane angle (projected modes only)
    theta=0.0,          # degrees, meridional-plane angle (intrinsic point only)

    # numerical controls
    integration=1,      # 0=Romberg, 1=Gauss–Legendre
    ngl_or_eps=0,       # 0 -> Fortran default
    algorithm=3,        # VP algorithm choice in Fortran (1 default; 2/3 enable regularised VP fits)
    maxmom=8,           # maximum moment order requested from the backend

    # what to compute
    kinematics="projected",   # "intrinsic" | "projected" | "both" | int in {0,1,2,3}
    average=False,            # point vs averaged mode (details below)

    # VP diagnostics
    usevp=True,               # enable VP analysis blocks and VP tables
    verbose_vp=0,             # Fortran verbosity flag for VP intermediate steps

    # output
    output_path=None,         # if set, write a persistent file and parse it
)
```

The return value is a `ScaleFreeResult` with (among other fields):

- `res.blocks`: a dictionary of named blocks (each block contains `columns` and `data`).
- `res.raw_text`: the raw parsed output text (useful for debugging/regression tests).

---

## 2. Kinematics modes and the `average` switch

The ScaleFree backend uses an integer control (`iwhat`) internally. In Python, you select the same behaviour using
`kinematics` and `average`:

| Python request           |    `average` | Fortran `iwhat` | What is computed                                    |
| ------------------------ | -----------: | --------------: | --------------------------------------------------- |
| `kinematics="intrinsic"` |      `False` |               0 | intrinsic moments at a point (set by `theta`)       |
| `kinematics="intrinsic"` |       `True` |               2 | intrinsic **shell average** (mass-weighted)         |
| `kinematics="projected"` |      `False` |               1 | projected moments at a point (set by `xi`)          |
| `kinematics="projected"` |       `True` |               3 | projected **circle average**                        |
| `kinematics="both"`      | `False/True` |  (0→1) or (2→3) | run intrinsic then projected, returning both blocks |

Notes:

- **Angles:** `theta` is only used for intrinsic point calculations (`iwhat=0`).  
  `xi` is only used for projected calculations (`iwhat=1`). In both cases, angles are given in degrees.
- **Scale-free normalisation:** the code works in dimensionless units; results can be rescaled to other radii using the scale-free nature of the models (see the associated theory notes for your workflow).

---

## 3. Output blocks: names, columns, and indices

`res.blocks` is a dictionary. Each block is typically a dict with:
- `columns`: list of column names (strings)
- `data`: a NumPy array of shape `(nrow, ncol)`

### 3.1 Intrinsic moments at a point: `intrinsic_point`

Present when `kinematics="intrinsic"` (or `"both"`) and `average=False`.

**Columns**
- `rho`   : intrinsic density (dimensionless)
- `vphi`  : ⟨v_φ⟩
- `vr2`   : ⟨v_r²⟩
- `vth2`  : ⟨v_θ²⟩
- `vphi2` : ⟨v_φ²⟩

The data contains a single row.

Example access:
```python
blk = res.blocks["intrinsic_point"]
rho, vphi, vr2, vth2, vphi2 = blk["data"][0]
```

### 3.2 Intrinsic shell average: `intrinsic_shell_average`

Present when `kinematics="intrinsic"` (or `"both"`) and `average=True`.

Same as `intrinsic_point`, plus:

- `beta`: anisotropy parameter inferred from the averaged intrinsic moments (mass-weighted).

Example access:
```python
rho, vphi, vr2, vth2, vphi2, beta = res.blocks["intrinsic_shell_average"]["data"][0]
```

### 3.3 Projected moments: `projected_point` and `projected_circle_average`

Present when `kinematics="projected"` (or `"both"`).

Both blocks have the same **row structure**: one row per projected component (`iproj`), with:

- `iproj=1`: LOS (line of sight)
- `iproj=2`: POSR (projected radial component on the sky)
- `iproj=3`: POST (projected tangential component on the sky)

**Columns**
- `iproj` : projection component index (1..3)
- `rho_p` : projected surface density at the chosen sky angle `xi`
- `v1`    : ⟨v⟩ (first raw moment) for the selected component
- `v2`    : ⟨v²⟩ (second raw moment)
- `v3`    : ⟨v³⟩ (third raw moment)
- `v4`    : ⟨v⁴⟩ (fourth raw moment)

Example access (LOS moments at the chosen `xi`):
```python
blk = res.blocks["projected_point"]          # or "projected_circle_average"
rows = blk["data"]                           # shape (3, 6)
los = rows[0]                                # iproj=1 is the first row
iproj, rho_p, v1, v2, v3, v4 = los
```

Derived quantities (using raw moments):
```python
mean = v1
var  = v2 - v1**2
sigma = var**0.5
```

Practical note: depending on parameters and on `maxmom`, some higher-order moments may diverge, in which case the backend may return `Infinity` for `v4` (or higher moments if requested elsewhere). This is expected behaviour for certain scale-free configurations.

---

## 4. VP diagnostics blocks (optional)

VP diagnostics are enabled with `usevp=True`. When enabled, the backend will produce:

- A VP summary block (`vp`) for projected cases.
- A VP summary block (`vp_intrinsic`) for intrinsic components.
- One VP table per component (`vp_table ...` or `vp_table_intrinsic ...`).

### 4.1 Projected VP summary: `vp`

This block contains one row per `iproj` (same mapping as above: 1=LOS, 2=POSR, 3=POST).

**Columns**
- `iproj`
- `true_gam`, `true_V`, `true_sig`  
  Parameters derived directly from the requested raw moments (e.g., mean/dispersion-like quantities).
- `gauss_gam`, `gauss_V`, `gauss_sig`  
  Best-fit Gaussian parameters to the VP.
- `h0 ... h6`  
  Gauss–Hermite coefficients of the VP relative to the best-fit Gaussian.

In most applications, `h3` and `h4` are the most commonly interpreted (asymmetric and symmetric deviations from a
Gaussian, respectively), but the backend provides coefficients through order 6.

Example: extract LOS Gaussian and Gauss–Hermite parameters:
```python
vp = res.blocks["vp"]
cols = vp["columns"]
row_los = vp["by_iproj"][1]   # convenience mapping, if present
# Or: row_los = vp["data"][0]
```

### 4.2 Intrinsic VP summary: `vp_intrinsic`

This is the intrinsic analogue of `vp`. Rows are indexed by `icomp`, where:

- `icomp=1`: intrinsic radial component (v_r)
- `icomp=2`: intrinsic polar component (v_θ)
- `icomp=3`: intrinsic azimuthal component (v_φ)

The columns have the same meaning as in `vp`:
`true_*`, `gauss_*`, and `h0..h6`.

### 4.3 VP tables: `vp_table` and `vp_table_intrinsic`

For each component, the backend also emits a sampled VP table:

- Projected: `res.blocks["vp_table"][iproj]`
- Intrinsic: `res.blocks["vp_table_intrinsic"][icomp]`

Each table has columns:

- `v`: velocity grid point
- `vp`: VP value at that grid point

Example:
```python
tbl = res.blocks["vp_table"][1]        # LOS VP table
v = tbl["data"][:, 0]
vp_vals = tbl["data"][:, 1]
```

---

## 5. Practical patterns

### 5.1 Running “everything” (intrinsic + projected)

```python
res = runner.vprofile(
    potential=2, gamma=2.0, q=0.608, df=1, beta=0.189, s=0.5, t=0.0,
    inclination=57.1, xi=0.0, theta=0.0,
    integration=1, ngl_or_eps=0, algorithm=3, maxmom=8,
    kinematics="both", average=False,
    usevp=True, verbose_vp=0,
)
print(sorted(res.blocks.keys()))
```

### 5.2 Accessing rows reliably via the parsed mappings

For VP blocks, the parser may provide a convenience mapping like `by_iproj` or `by_icomp` so you can access rows by their index without relying on row order:

```python
proj_vp = res.blocks["vp"]["by_iproj"][1]              # LOS row
intr_vp = res.blocks["vp_intrinsic"]["by_icomp"][3]    # v_phi row
```

If the mapping is absent, fall back to row order (`data[0]` is component 1, etc.).

---

## 6. Troubleshooting and tips

- **Unexpected `Infinity` / `nan`:** this typically indicates divergence of a requested moment order for the chosen model parameters. Reduce `maxmom` or interpret only the finite-order results.
- **Reproducibility:** for regression tests and cross-platform comparisons, prefer `output_path=...` so you parse a persistent file written by the backend.
- **Debugging the interactive dialogue:** set `debug_prompts=True` to print the Fortran prompt/response flow.

---

## 7. Minimal “sanity check” script

This mirrors the repository smoke-test behaviour (see `scripts/quick_user_run.py`) but focuses on a single run:

```python
from pathlib import Path
from scalefree.vmoments import ScaleFreeRunner

runner = ScaleFreeRunner(exe_path=Path("fortran_src/scalefree.e"))

res = runner.vprofile(
    potential=2, gamma=2.0, q=0.608, df=1, beta=0.189, s=0.5, t=0.0,
    inclination=57.1, xi=0.0, theta=0.0,
    integration=1, ngl_or_eps=0, algorithm=3, maxmom=8,
    kinematics="projected", average=False,
    usevp=True, verbose_vp=0,
)

print(res.blocks["projected_point"]["columns"])
print(res.blocks["projected_point"]["data"])
print(res.blocks["vp"]["columns"])
print(res.blocks["vp"]["data"])
```
Or similarly:

```python
from pathlib import Path
from scalefree.vmoments import ScaleFreeRunner

runner = ScaleFreeRunner()

res = runner.vprofile(
    potential=2, gamma=2.0, q=0.608, df=1, beta=0.189, s=0.5, t=0.0,
    inclination=57.1, xi=0.0, theta=0.0,
    integration=1, ngl_or_eps=0, algorithm=3, maxmom=8,
    kinematics="intrinsic", average=False,
    usevp=True, verbose_vp=0,
)

print(res.blocks["intrinsic_point"]["columns"])
print(res.blocks["intrinsic_point"]["data"])
print(res.blocks["vp_intrinsic"]["columns"])
print(res.blocks["vp_intrinsic"]["data"])
```