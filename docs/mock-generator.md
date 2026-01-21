# Mock generator (`scalefree.mock`)

`mock()` produces a synthetic 6D sample array:

- output shape: `(N, 6)`
- columns: `(x, y, z, vx, vy, vz)` in the model frame

## What it does (pipeline)

1. Sample `N` positions from the scale-free density.
2. Rotate into the sky frame for the chosen inclination.
3. Bin stars by projected-plane angle in `[0, 180)` (symmetry).
4. For each occupied bin:
   - run `vprofile()` once to get Gauss–Hermite parameters for `iproj=1,2,3`
5. Sample velocities per star using Sanders–Evans style GH PDFs.
6. Convert `(vlos, vposr, vpost)` back into `(vx, vy, vz)`.

## Minimal usage

```python
from scalefree import mock

X = mock(
    potential=lambda: 1,
    gamma=4.0,
    q=0.9,
    beta=0.0,
    s=0.5,
    t=0.0,
    inclination=90.0,
    nsamples=50_000,
    nbins=180,
    usevp=True,
)
```

## Dependencies

- Requires the ScaleFree Fortran backend for `vprofile()`.
- For velocity sampling, the code attempts to import:
  - `from balrogo import dynamics`
  - or falls back to importing a module named `dynamics`

If neither is available, install BALRoGO or ensure `dynamics.py` is importable.

## Performance tips

- `nbins` controls how many times `vprofile()` is called.
  If you need speed, reduce `nbins`.
- For stable results, set `seed=...`.
