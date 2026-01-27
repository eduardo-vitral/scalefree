
## Discussion of the method

Construct galaxy models for flattened systems of test-particles
in spherical potentials, as discussed in [De Bruijne et al (1996)](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract).
 
Let $b$ be a reference length, $\rho_0$ a reference mass density, and
$W$ a reference velocity.
 
Introduce dimensionless units:
 
* $r = r_{\mathrm{true}} / b$.
* $R = R_{\mathrm{true}} / b$.
* $z = z_{\mathrm{true}} / b$.
* $v = v_{\mathrm{true}} / (\sqrt{2}  W)$.
* $L = L_{\mathrm{true}} / (\sqrt{2}  W  b)$.
* $E = E_{\mathrm{true}} / (W^{2})$.
* $\Psi = \Psi_{\mathrm{true}} / (W^2)$.
* $\rho = \rho_{\mathrm{true}} / \rho_0$.
* $f = f_{\mathrm{true}} / [\rho_0 (2 W^2)^{-3/2}]$.
 
We consider two cases for the potential:
 
#### KEPLERIAN POTENTIAL:

Let $M$ be the total mass of the galaxy, such that $\Psi_{\mathrm{true}} = G  M / r_{\mathrm{true}}$.\
Then choose $W = \sqrt{GM/b}$.\
In these units: $\Psi = 1 / r$.
 
#### LOGARITHMIC POTENTIAL:

Let $V_c$ be the circular velocity, such that $\Psi_{\mathrm{true}} = - V_{c}^2  \ln{(r_{\mathrm{true}}/b)}$.
Then choose $W = V_c$.\
In these units: $\Psi = - \ln{(r)}$.

***
 
Let the mass density fall off with logarithmic slope $\gamma$, and be 
stratified on spheroids with axial ratio $q$:

$\rho  = r^{-\gamma} (\sin(\theta)^2 + [\cos(\theta)^2/q^2])^{-\gamma/2}$.

where $\theta$ is the polar angle such that:

$R = r \sin(\theta)$ &nbsp; and &nbsp; $z = r \cos(\theta)$.
 
Since the potential is spherical, the quantities $E$, $L^2$, $L_z^2$ are 
integrals of motion. Let $L_{\mathrm{max}}(E)$ be the maximum angular momentum that 
can be attained by a star at energy $E = \Psi - v^2$.

* In the **Keplerian** potential: $L^2_{\mathrm{max}}(E) = 1 / (4E)$.

* In the **Logarithmic** potential: $L^2_{\mathrm{max}}(E) = \exp(-2E -1) / 2$
 
We consider even DFs that are separable functions or quasi-separable 
functions of $E$, $\zeta^2$, $\eta^2$:

* **case I**  : $f_e = g(E)  \zeta^{-2 \beta}  j(e^2  \eta^2)$

* **case II** : $f_e = g(E)  \zeta^{-2 \beta}  h(e^2  \eta^2 / \zeta^2)$

where:
* $\zeta^2 = L^2/L^2_{\mathrm{max}}(E)$

* $\eta^2  = L_z^2/L^2_{\mathrm{max}}(E)$

* $e^2 = 1 - q^2$

and $\beta$ is a free parameter. With these ansatz's, the functions 
$j$ and $h$ are determined uniquely (see [De Bruijne et al, 1996](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract)). For the **case I** DF, notice that while the $\beta$ parameter above does not necessarily match the definition in ([Vitral et al. 2024](https://iopscience.iop.org/article/10.3847/1538-4357/ad571c), Eq. [5]), the current software handles this internally such that the provided $\beta$ argument *ALWAYS* follows the latter definition. In practice, the code fits the necessary case I-$\beta$ that reproduces the desired $\beta_{\mathrm B}$.

To allow models with rotation, please read this [link](./rotation.md).
 
In addition to the calculation of the intrinsic and projected velocity
moments, the program `scalefree.f` also reconstructs the projected
line-of-sight and plane-of-sky
VP shapes (details [here](./vp-shapes.md)).

### Plane-of-sky additions

The current package is very similar to the original `scalefree` 
software presented in
[De Bruijne et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract),
with the difference that now it is able to solve the velocity
moments of plane of sky coordinates (POSr and POSt). 
The mathematical description of the new models is given [here](./pos-velocities$\beta$.md).

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
>* Vitral E., et al., 2024, ApJ, 970, 1