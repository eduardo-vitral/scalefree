# Getting started

This package wraps classic Fortran executables with a small Python interface.

## What runs in Fortran vs Python?

- **Fortran backend**
  - `vprofile()` / `ScaleFreeRunner.vprofile()` run the ScaleFree Fortran program.
  - `hermite.hermite()` can run a Fortran Gauss–Hermite fitting helper (`fitvp`-style).

- **Python**
  - Parses structured output from STDOUT (and/or optional output files).
  - Provides a convenience API.
  - `mock()` is Python-only for sampling positions and velocities, but calls `vprofile()` per bin.

## First run: backend compilation

On first use, `vprofile()` typically tries the following, in order:

1. Use `exe_path=...` if you provided it
2. Use `SCALEFREE_EXE` if set
3. Use a cached executable in a user cache directory
4. If none exists, auto-compile from packaged `fortran_src/scalefree.f` (requires `gfortran`)

If `gfortran` is missing, you will get an error explaining how to install it.

### Provide a precompiled executable (recommended on clusters)

If you already compiled the backend yourself:

```bash
export SCALEFREE_EXE=/path/to/scalefree.e
```

Or in Python:

```python
from scalefree import vprofile
res = vprofile(exe_path="/path/to/scalefree.e", ...)
```

## Minimal example

```python
from scalefree import vprofile

res = vprofile(
    potential="logarithmic",
    gamma=2.0, q=0.8, df=1, beta=0.0, s=0.5, t=0.0,
    inclination=60.0, xi=0.0, theta=0.0,
    average=False,
    usevp=True,
    algorithm=3,
)

print(res.blocks.keys())
print(res.blocks["vp"]["by_iproj"][1])  # iproj = 1/2/3
```

## Where do executables get cached?

Both backends cache binaries in a user cache directory (exact location depends on OS).

## Next steps

- Read `docs/vprofile.md` to understand parameters and output blocks.
- Read `docs/gauss-hermite.md` for Gauss–Hermite fitting.
- Read `docs/mock-generator.md` for synthetic data generation.
