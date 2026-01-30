# scalefree

[![pipeline status](https://gitlab.com/eduardo-vitral/scalefree/badges/main/pipeline.svg)](https://gitlab.com/eduardo-vitral/scalefree/-/commits/main)
[![coverage report](https://gitlab.com/eduardo-vitral/scalefree/badges/main/coverage.svg)](https://gitlab.com/eduardo-vitral/scalefree/-/commits/main)
[![pypi](https://img.shields.io/pypi/v/scalefree.svg)](https://pypi.python.org/pypi/scalefree/)
[![python](https://img.shields.io/pypi/pyversions/scalefree.svg)](https://pypi.python.org/pypi/scalefree)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)


<!-- markdownlint-disable-next-line no-inline-html -->
<img alt="logo" align="right" src="https://gitlab.com/eduardo-vitral/scalefree/-/raw/main/images/scalefree_logo.png" width="20%" />


Scale-free dynamical models (Fortran backend) with a small Python interface for computing intrinsic/projected velocity moments and (optionally) reconstructed velocity profiles (VPs) and Gauss–Hermite summaries.

The guiding principles are:
- Keep the user-facing API small and predictable
- Delegate heavy computation to a compiled backend
- Provide structured, parseable outputs for reproducible workflows

This code was developped in the context of the [HSTPROMO](https://www.stsci.edu/~marel/hstpromo.html) collaboration.

---

## What this package provides

### `vprofile(...)`
Runs the ScaleFree Fortran backend and returns intrinsic/projected velocity moments and (optionally) VP / Gauss–Hermite summaries.

Typical uses:
- Line-of-sight (LOS) kinematics
- Proper-motion (plane-of-sky) kinematics (POSr, POSt)
- VP reconstruction from moments (when enabled)

### `hermite(...)`
Utilities to fit Gauss–Hermite coefficients to a velocity-profile file and to evaluate analytic Gauss–Hermite profiles.

### `mock(...)`
A convenience routine for generating synthetic 6D samples `(x, y, z, vx, vy, vz)` from a chosen model configuration (internally uses `vprofile()` to obtain GH parameters across angular bins). Keep in mind that the underlying formalism is scale-free, so please refer to this [link](./docs/theory/main-method.md) for scaling it back to physical values.

---

## Installation

### From source (recommended for development)
```bash
poetry install
```

### With pip
```bash
pip install scalefree --upgrade
```

---

## Fortran requirement

`vprofile()` uses a Fortran executable. If you do not provide a precompiled binary, a local build requires **gfortran**.

Install gfortran:
- Debian/Ubuntu: `sudo apt-get install gfortran`
- Fedora: `sudo dnf install gcc-gfortran`
- macOS (Homebrew): `brew install gcc`
  
---

## Quickstart

### 1) Run `vprofile()`
```python
from scalefree import vprofile

# Minimal projected (point) run + VP/Gauss–Hermite reconstruction
res = vprofile(
    potential="log",   # "log"/"logarithmic" or "kepler" (or 2 / 1)
    gamma=2.0,
    q=0.8,
    df=1,
    beta=0.0,
    s=0.0,
    t=0.0,
    inclination=60.0,  # degrees
    xi=0.0,            # degrees on projected plane (0 = major axis)
    theta=0.0,         # degrees (used only for intrinsic runs; 0 = symmetry axis)
    usevp=True,        # set False if you only want moments
)

# Projected velocity moments (iproj: 1=LOS, 2=POSr, 3=POSt)
proj_los = res.blocks["projected_point"]["by_iproj"][1]
mu = proj_los["v1"]
sigma = (proj_los["v2"] - mu**2) ** 0.5
print(f"LOS: mean={mu:.6g}, sigma={sigma:.6g}")

# Optional: VP / Gauss–Hermite summary (only present when usevp=True)
vp_los = res.blocks["vp"]["by_iproj"][1]
print(f"LOS GH: h3={vp_los['h3']:.6g}, h4={vp_los['h4']:.6g}")
```

### 2) Fit Gauss–Hermite moments from a VP file
```python
from pathlib import Path
import tempfile
import numpy as np

from scalefree import vprofile, hermite

res = vprofile(
    potential="log",
    gamma=2.0,
    q=0.8,
    beta=0.0,
    s=0.0,
    t=0.0,
    inclination=60.0,
    xi=0.0,
    theta=0.0,
    usevp=True,
)

iproj = 1  # 1=LOS, 2=POSr, 3=POSt
vp_data = res.blocks["vp_table"][iproj]["data"]  # (N,2): [v, vp]

# Write to a writable temp location
vp_path = Path(tempfile.gettempdir()) / "my_vp.dat"
np.savetxt(vp_path, vp_data, header="v vp")

gauss_info, gaussh_info, h_moments = hermite.hermite(vp_path)
print(gaussh_info)
print(h_moments["h3"], h_moments["h4"])
```

### 3) Generate a simple mock
```python
from scalefree import mock

X = mock(
    potential=2,
    gamma=2.0,
    q=0.9,
    df=2,
    beta=-0.1,
    s=0.2,
    t=0.0,
    nsamples=10_000,
    nbins=180,
)

print(X.shape)  # (N, 6): (x, y, z, vx, vy, vz)
```

---

## Documentation

To keep this README short, more detailed guides are intended to live under `docs/`:

- `docs/vmoments.md`
- Theory notes:
  - `docs/theory/rotation.md`
  - `docs/theory/vp-shapes.md`
  - `docs/theory/pos-velocity.md`

---

## References and citation guidance

If you use this code in research, please cite the two works below and let us know so you can add your publication to [this list](./docs/used_in.md).

### Core scale-free models (foundational)
- de Bruijne, van der Marel & de Zeeuw (1996), *MNRAS*, 282, 909–925. arXiv: astro-ph/9601044

### Plane of sight moments equations
- Vitral et al. (2024), *ApJ*, 970, 1. DOI: 10.3847/1538-4357/ad571c. arXiv: astro-ph/2407.07769

### BibTeX (copy/paste)
```bibtex
@ARTICLE{1996MNRAS.282..909D,
  author  = {de Bruijne, Jos H. J. and van der Marel, Roeland P. and de Zeeuw, P. Tim},
  title   = {Scale-free dynamical models for galaxies: flattened densities in spherical potentials},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year    = {1996},
  volume  = {282},
  number  = {3},
  pages   = {909--925},
  doi     = {10.1093/mnras/282.3.909},
  eprint  = {astro-ph/9601044},
  archivePrefix = {arXiv}
}

@ARTICLE{2024ApJ...970....1V,
       author = {{Vitral}, Eduardo and {van der Marel}, Roeland P. and {Sohn}, Sangmo Tony and {Libralato}, Mattia and {del Pino}, Andr{\'e}s and {Watkins}, Laura L. and {Bellini}, Andrea and {Walker}, Matthew G. and {Besla}, Gurtina and {Pawlowski}, Marcel S. and {Mamon}, Gary A.},
        title = "{HSTPROMO Internal Proper-motion Kinematics of Dwarf Spheroidal Galaxies. I. Velocity Anisotropy and Dark Matter Cusp Slope of Draco}",
      journal = {\apj},
     keywords = {Dark matter, Dwarf spheroidal galaxies, Astronomy data analysis, Proper motions, Stellar kinematics, Stellar dynamics, Galaxy dynamics, Galaxy structure, 353, 420, 1858, 1295, 1608, 1596, 591, 622, Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2024,
        month = jul,
       volume = {970},
        pages = {1},
          doi = {10.3847/1538-4357/ad571c},
       eprint = {2407.07769},
}
```

---

### Authors:

* Roeland P. van der Marel,

    >1994-1995 :\
        &emsp; development of code\
    address   : Space Telescope Science Institute\
                Research Programs Office (RPO)\
                3700 San Martin Drive\
                Baltimore, MD 21218\
                Tel : (+1) 410 338 4931\
                Fax : (+1) 410 338 4596\
                e-mail   : marel@stsci.edu\
                homepage : https://www.stsci.edu/~marel/
 
* Jos H. J. de Bruijne,

    >1994-1995 :\
        &emsp; testing and application of code\
    address   : Scientific Support Oﬃce\
                Directorate of Science\
                European Space Research and Technology Centre (ESA/ESTEC)\
                Keplerlaan 1\
                2201AZ, Noordwijk\
                The Netherlands
                e-mail   : jos.de.bruijne@cosmos.esa.int \
                homepage : https://www.cosmos.esa.int/web/personal-profiles/jos-de-bruijne

* Eduardo Vitral,

    >2023-present :\
        &emsp; development of the Python interface\
        &emsp; implementation of plane-of-sky routines\
        &emsp; testing and application of code\
    address   : Institute for Astronomy\
                University of Edinburgh\
                Royal Observatory, Blackford Hill\
                Edinburgh EH9 3HJ\
                United Kingdom\
                e-mail   : eduardo.vitral@roe.ac.uk\
                homepage : https://eduardo-vitral.github.io
