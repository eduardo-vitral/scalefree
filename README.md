# scalefree :dizzy:
***

<img alt="logo" align="right" src="./images/hstpromo.jpeg" width="18%" />

### This is a software developped by the [HSTPROMO](https://www.stsci.edu/~marel/hstpromo.html) collaboration

* Scale-free Modeling Software to predict the intrinsic and\
projected kinematical quantities of axisymmetric mass densities\
with anisotropic velocity distributions in spherical potentials.

* *Gauss-Hermite fitting routine.

***

## PART I: GENERAL INFORMATION

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
                homepage : http://sol.stsci.edu/~marel/
 
* Jos H. J. de Bruijne,

    >1994-1995 :\
        &emsp; testing and application of code\
    address   : Sterrewacht Leiden\
                Postbus 9513\
                2300 RA Leiden\
                The Netherlands\
                Tel : (+31) 71 5275878\
                Fax : (+31) 71 5275819\
                e-mail   : debruyne@strw.LeidenUniv.nl\
                homepage : http://www.strw.leidenuniv.nl/~debruyne/

* Eduardo F. Vitral,

    >2023 :\
        &emsp; development of the Python interface\
        &emsp; implementation of plane-of-sky routines\
        &emsp; testing and application of code\
    address   : Space Telescope Science Institute\
                Science Mission Office (SMO)\
                3700 San Martin Drive\
                Baltimore, MD 21218\
                e-mail   : evitral@stsci.edu\
                homepage : http://www.iap.fr/useriap/vitral/index.html

### Summary of the method
 
We consider the case of a scale-free spheroidal mass density with an
anisotropic velocity distribution in a scale-free spherical
potential. The assumption of a spherical potential has the advantage
that all integrals of motion are known explicitly. Two families of
phase-space distribution functions are considered. The **case I**
distribution functions are anisotropic generalizations of the
flattened $f(E,L_z)$ model, which they include as a special case. The
**case II** distribution functions generate flattened
constant-anisotropy models. Free parameters control the radial
power-law slopes of the mass density and potential, the flattening of
the mass distribution, and the velocity dispersion anisotropy. The
models can describe the outer parts of galaxies and the density cusp
structure near a central black hole, but also provide general insight
into the dynamical properties of flattened systems. Because of their
simplicity, they provide a useful complementary approach to the
construction of flattened self-consistent three-integral models for
elliptical galaxies.


### References
 
The method and some applications are discussed in the following paper:
``` 
   `Scale-free dynamical models for galaxies: 
    flattened densities in spherical potentials'
       de Bruijne J., van der Marel R.P., de Zeeuw P.T. 
       MNRAS, 282, 909-925, 1996 
```
which can be retrieved at this [link](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract).
 

### Acknowledgments

If you have found this software useful for your research, we would
appreciate an acknowledgment to *"use of the Scale-free Modeling
Software developed by R.P. van der Marel, J.H.J. de Bruijne and E. Vitral"*, 
along with a citation to the paper from the **References** section above.
The respective `BibTeX` format is:

```bibtex
@ARTICLE{1996MNRAS.282..909D,
       author = {{de Bruijne}, Jos H.~J. and {van der Marel}, Roeland P. and {de Zeeuw}, P. Tim},
        title = "{Scale-free dynamical models for galaxies: flattened densities in spherical potentials}",
      journal = {\mnras},
     keywords = {LINE: PROFILES, GALAXIES: ELLIPTICAL AND LENTICULAR, CD, GALAXIES: INDIVIDUAL: NGC 2434, GALAXIES: INDIVIDUAL: NGC 3706, GALAXIES: KINEMATICS AND DYNAMICS, GALAXIES: STRUCTURE, Astrophysics},
         year = 1996,
        month = oct,
       volume = {282},
       number = {3},
        pages = {909-925},
          doi = {10.1093/mnras/282.3.909},
archivePrefix = {arXiv},
       eprint = {astro-ph/9601044},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


### Bug Reports
 
Please send bug reports and important comments or questions to
marel@stsci.edu
 
 
### Caveats
 
This software comes without guarantee and on an *as is* basis. It is
not being actively maintained or upgraded.
 
All original testing of the code was done on a `UNIX Sparc 2`
workstation with the `SunOS` operating system in the period
1994-1995. In 1997, R.P. van der Marel ensured that the code would 
run on a `Sun Ultra 170` workstation with the `solaris` operating system. 
However, little detailed testing of the code was performed 
on this operating system.
 
It is not guaranteed that the software will work without problems on
other platforms, with different operating systems and different
compilers. However, the code is close to standard fortran, so if any,
only minor revisions will be necessary.
 
R.P. van der Marel will consider all requests for use of this software. 
However, R.P. van der Marel do
not allow people who have received the software from him to distribute
it to third parties. 
Anyone who uses this software is encouraged to
send him an email with their address, so that he can send reports
of bugs and revisions.
 

### Program structure

The following files/folders are present in this package:

- `ghermite/`: Folder containing Gauss-Hermite fitting routines. The files are:
    - `fitvp.e`: Executable.
    - `fitvp.f`: The `Fortran` program to fit Gauss-Hermite moments (up to order 10).
- `scalefree_los/`:
    - `makefile`: makefile that allows compilation of the program `scalefree.f`.
    - `scalefree.e`: Executable.
    - `scalefree.f`: The `Fortran` program that does all the calculations. 
                the file also contains all necessary subroutines, some
                of which are from the [*Numerical Recipes book by Press et
                al.*](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)
- `scalefree_posr/`:
    - `makefile`: makefile that allows compilation of the program `scalefree.f`.
    - `scalefree.e`: Executable.
    - `scalefree.f`: The `Fortran` program that does all the calculations. 
                the file also contains all necessary subroutines, some
                of which are from the [*Numerical Recipes book by Press et
                al.*](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)
- `scalefree_post/`:
    - `makefile`: makefile that allows compilation of the program `scalefree.f`
    and `fitvp.f`.
    - `scalefree.e`: Executable.
    - `scalefree.f`: The `Fortran` program that does all the calculations. 
                the file also contains all necessary subroutines, some
                of which are from the [*Numerical Recipes book by Press et
                al.*](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)
- `examp/`: directory with example input files for the program `scalefree.f`
    and `fitvp.f`.
- `scalefree.py`: The `Python` program that translates the `Fortran` outputs.
- `install.py`: The `Python` program that installs the code.


- `README.md`: This file. Contains instructions and a description of the
                software.


### Installation 

To install the software, run the `install.py` file.\
The code will ask you the path where your scalefree files are installed.\
Once you type it and hit **ENTER**, the following output should be produced:
```console
Provide the path where your scalefree files are installed.
For example: /home/yourname/code-packages/scalefree/ 

/YOUR-PATH/scalefree/

Installation succesfull.
Finished.
```

This will construct the executables *`.e`.


### Examples

The directory `./examp` contains example input files. 

The code `install.py` is itself an example that runs `scalefree`, and
creates executable *`.e` files. Once the software is installed, one
can remove the `exec=True` option for faster outputs.


### Known Bugs

As discussed below, the calculations of the Gauss-Hermite moments of
the projected velocity profiles (VPs) do not work well for certain
parameter ranges. The intrinsic and projected velocity moments should
always be highly accurate.


## PART II: Discussion of the method

Construct galaxy models for flattened systems of test-particles
in spherical potentials, as discussed in [De Bruijne et al (1996)](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract).
 
Let $b$ be a reference length, $\rho_0$ a reference mass density, and
$W$ a reference velocity.
 
Introduce dimensionless units:
 
* $r = r_{true} / b$.
* $R = R_{true} / b$.
* $z = z_{true} / b$.
* $v = v_{true} / (\sqrt{2}  W)$.
* $L = L_{true} / (\sqrt{2}  W  b)$.
* $E = E_{true} / (W^{2})$.
* $\Psi = \Psi_{true} / (W^2)$.
* $\rho = \rho_{true} / \rho_0$.
* $f = f_{true} / [\rho_0 (2 W^2)^{-3/2}]$.
 
We consider two cases for the potential:
 
#### KEPLERIAN POTENTIAL:

Let $M$ be the total mass of the galaxy, such that $\Psi_{true} = G  M / r_{true}$.\
Then choose $W = \sqrt{GM/b}$.\
In these units: $\Psi = 1 / r$.
 
#### LOGARITHMIC POTENTIAL:

Let $V_c$ be the circular velocity, such that $\Psi_{true} = - V_{c}^2  \ln{(r_{true}/b)}$.
Then choose $W = V_c$.\
In these units: $\Psi = - \ln{(r)}$.

***
 
Let the mass density fall off with logarithmic slope $\gamma$, and be 
stratified on spheroids with axial ratio $q$:

$\rho  = r^{-\gamma} (\sin(\theta)^2 + [\cos(\theta)^2/q^2])^{-\gamma/2}$.

where $\theta$ is the polar angle such that:

$R = r \sin(\theta)$ &nbsp; and &nbsp; $z = r \cos(\theta)$.
 
Since the potential is spherical, the quantities $E$, $L^2$, $L_z^2$ are 
integrals of motion. Let $L_{max}(E)$ be the maximum angular momentum that 
can be attained by a star at energy $E = \Psi - v^2$.

* In the **Keplerian** potential: $L^2_{max}(E) = 1 / (4E)$.

* In the **Logarithmic** potential: $L^2_{max}(E) = \exp(-2E -1) / 2$
 
We consider even DFs that are separable functions or quasi-separable 
functions of $E$, $\zeta^2$, $\eta^2$:

* **case I**  : $f_e = g(E)  \zeta^{-2 \beta}  j(e^2  \eta^2)$

* **case II** : $f_e = g(E)  \zeta^{-2 \beta}  h(e^2  \eta^2 / \zeta^2)$

where:
* $\zeta^2 = L^2/L^2_{max}(E)$

* $\eta^2  = L_z^2/L^2_{max}(E)$

* $e^2 = 1 - q^2$

and $\beta$ is a free parameter. With these ansatz's, the functions 
$j$ and $h$ are determined uniquely (see [De Bruijne et al, 1996](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract)).

To allow models with roation, please read this [link](./ROTATION.md).
 
In addition to the calculation of the intrinsic and projected velocity
moments, the program `scalefree.f` also reconstructs the projected
line-of-sight and plane-of-sky
VP shapes (details [here](./VPSHAPES.md)).

### Plane-of-sky additions

The current package is very similar to the original `scalefree` 
software presented in
[De Bruijne et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract),
with the difference that now it is able to solve the velocity
moments of plane of sky coordinates (POSr and POSt). 
The mathematical description of the new models is given [here](./VPOS.md).

### References

>* de Bruijne J., van der Marel R.P., de Zeeuw P.T., 1996, MNRAS, 282, 909
>* Dehnen W., Gerhard O. E., 1994, MNRAS, 268, 1019
>* Gerhard O. E., 1993, MNRAS, 265, 213
>* Press W. H., Teukolsky S. A., Vetterling W. T., Flannery B. P., 
     1992, Numerical Recipes, Second Edition. 
     Cambridge University Press, Cambridge
>* Qian E. E., de Zeeuw P. T., van der Marel R. P., Hunter C., 1995, 
     MNRAS, 274, 602
>* van der Marel R. P., Franx M., 1993, ApJ, 407, 525


## PART III: Description of program in- and output

After the user has entered the required model parameters, and has
answered to the program how it must deal with the various numerical
details, the code is set to run.

There are two main routines in the `scalefree.py` program:

```python
scalefree.hermite(input, exec=False)

    """
    Returns the fits of a Gauss-Hermite adjustment to data.

    Parameters
    ----------
    input: str
        Path where to find the data to be fitted.
        Expected file shape:
        --> First line (header): "	v  VP(v)"
        --> 1st column: x / 2nd column: f(x)
    exec: boolean
        True, if the user wants to generate new .e files.

    Returns
    -------
    fitinfo : dictionaries
        In order:
        gauss_info: real Gaussian fit.
                - norm
                - mean
                - dispersion
        gaussh_info: Gauss-Hermite fit.
            - norm
            - mean
            - dispersion
        h_moments: First 0-10 moments of the
                    Gauss-Hermite fit.
            - hi, with i in [0, 10]
    """
```

```python
scalefree.vprofile(
    potential="1",
    gamma="4",
    q="1",
    df="1",
    beta="0",
    s="0.5",
    t="0",
    inclination="90",
    integration="1",
    ngl="0",
    algorithm="1",
    maxmom="0",
    theta="0",
    xi="0",
    exec=False)

    """
    Returns the fits of a Gauss-Hermite adjustment to data.
    Note: all results are at an
    (intrinsic or projected) radius of 1 in
    dimensionless units. Results can be scaled to other
    radii using the scale-free nature of the models.

    Parameters
    ----------
    potential: str
        Kepler (1) or Logarithmic (2) Potential.
    gamma: str
        Power-law slope gamma of the mass density.
    q: str
        Intrinsic axial ratio q of the mass density.
    df: str
        Case I (1) or Case II (2) DF
    beta: str
        Anisotropy parameter beta of the DF
    s: str
        Odd part parameter s for the DF
    t: str
        Odd part parameter t for the DF
    inclination: str
        Viewing inclination i in degrees (90=edge-on)
    integration: str
        Use Romberg (0) or Gauss-Legendre (1) integration
        for line-of-sight projection
    ngl: str
        Number of quadrature points / fractional accuracy
    algorithm: str
        Algorithm to calculate VPs and Gauss-Hermite moments
            1: Solve VanderMonde matrix directly without
            regularization. Resulting VP will be nonsense
            but the GH moments are generally well
            determined.
            2: Use regularization with a fixed regularization
            parameter.
            3: Use regularization. Increase regularization
            parameter until the VP has no more than 3
            significant local maxima. A local maximum is
            significant if it exceeds the value of its
            neighbors on the grid by eps times the
            absolute VP maximum.
    maxmom: str
        Maximum number of projected moments
        to use (should be an even number)
    theta: str
        Angle theta in the meridional plane
        (in degrees) (0 = symmetry axis)
    xi: str
        Angle on the projected plane
        (in degrees) (0 = major axis)
    exec: boolean
        True, if the user wants to generate new .e files.

    Returns
    -------
    vinfo : objects
        List of fit parameters for each of the three observable
        dimensions: LOS, POSr, POSt (in this order).
        Each dimension is itself a collection of dictionaries.
        In this order:
            intmom: Intrinsic velocity moments:
                - rho
                - <v_ph>
                - <v_r^2>
                - <v_th^2>
                - <v_ph^2>
            projmom: Projected velocity moments:
                - <rho>_p
                - <v>_p
                - <v^2>_p
                - <v^3>_p
                - <v^4>_p
            gauss_info: real Gaussian fit.
                - norm
                - mean
                - dispersion
            gaussh_info: Gauss-Hermite fit.
                - norm
                - mean
                - dispersion
            h_moments: First 0-10 moments of the
                       Gauss-Hermite fit.
                - hi, with i in [0, 6]
            vinfo: Velocity distribution function
                - x
                - f(x)
    """
```

***