# VELOCITY PROFILE SHAPES
***

In addition to the calculation of the intrinsic and projected velocity
moments, the program `scalefree.f` also reconstructs the projected
line-of-sight and plane-of-sky
VP shapes. This however is a rather more difficult
problem than the calculation of the projected velocity moments
themselves. The approach that we have used is to characterize each
projected VP as a set of $N$ values on an array. These values are then
related to the (known) first N projected velocity moments through the
solution of a so-called **VanderMonde matrix**. Once a VP is known on an
array of velocity values, characteristic quantities such as the
Gauss-Hermite moments 
([van der Marel & Franx, 1993](https://ui.adsabs.harvard.edu/abs/1993ApJ...407..525V/abstract)) 
can be evaluated in
straightforward fashion.
 
The difficulty lies in finding the solution of the VanderMonde matrix
equation, which is an ill-conditioned problem. We have implemented two
different approaches:

1. The first approach is to solve the VanderMonde matrix equation
directly, without any regularization or other tampering with the
equation itself. This yields solutions that are oscillatory on the
scale of the velocity array spacing, due to the ill-conditioned nature
of the problem (noise amplification). Although the VP itself is
ridiculous, it turns out that quantities like the Gauss-Hermite
moments (which are integrals over the VP) are often well and correctly
determined. The question remains how many ($N$) projected velocity
moments should be used to get the most accurate Gauss-Hermite moment
determinations. The program chooses $N$ by starting with a low value of
$N$, and then solving the VanderMonde matrix equation repeatedly until
the Gauss-Hermite moments have converged to within a preset
tolerance. The reason that this approach works is because a
Gauss-Hermite series bears a strong resemblance to a Fourier series.
Gauss-Hermite moments $h_i$ of higher order measure power in the VP on
higher frequencies 
([Gerhard, 1993](https://ui.adsabs.harvard.edu/abs/1993MNRAS.265..213G/abstract)). 
For high $N$, numerical errors cause
the VanderMonde solution to become oscillatory on the scale of the
velocity array spacing. This, however, does not influence the
observables (${\gamma}$, ${V}$, ${\sigma}$, ${h_3}$, ${h_4}$,...), 
because these only measure
power on relatively low frequencies.

2. The second approach is to add a regularization term to the
VanderMonde matrix equation, and then solve it in a least-squares
sense (cf. Section 18.5 of [*Numerical Recipes*](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)). This works well, but
leaves open the question how the regularization parameter must be
chosen. The program gives two options: 
* A. Let the user give the
regularization parameter by hand; 
* B. First start with the case of
essentially no regularization. This yields an oscillatory solution
with many local peaks. Then increase the amount of regularization in
steps until the solution has become such that there are no more than 3
significant local maxima (a local maximum is considered significant if
it exceeds the value of its neighbors on the grid by a preset
tolerance eps times the absolute VP maximum).

The accuracy of algorithm (1) was tested as described in Appendix A of
[De Bruijne et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996MNRAS.282..909D/abstract). 
To recapitulate, the results are as follows. A
program was written that calculates the VP for the spherical ($q=1$)
case of the models by direct evaluation of the three-dimensional
integral along the line of sight that defines the VP. The
(${\gamma}$, ${V}$, ${\sigma}$, ${h_3}$, ${h_4}$,...)
 calculated with the program were found to
be accurate, except for the case of a logarithmic potential and large
anisotropy (${\beta < -4.0}$ or ${\beta > 0.8}$). This has two reasons: 
* (i) the logarithmic potential has no escape velocity so that it is more
difficult to represent the VP on a finite velocity array;
* (ii) the VP becomes discontinuous in the limit ${\beta \rightarrow -\infty}$ (only
circular orbits), and singular in the limit ${\beta \rightarrow 1}$ (only radial
orbits) (e.g., 
[van der Marel & Franx, 1993](https://ui.adsabs.harvard.edu/abs/1993ApJ...407..525V/abstract)).

A second test was
provided by the ${f(E,L_z)}$ case (${\beta=0}$), for which the
three-dimensional VP integral could be calculated as in [Qian et
al. (1995)](https://ui.adsabs.harvard.edu/abs/1995MNRAS.274..602Q/abstract). Again, the 
(${\gamma}$, ${V}$, ${\sigma}$, ${h_3}$, ${h_4}$,...)
calculated with the
VanderMonde algorithm were generally found to be accurate, with the
exception of rather flattened models ${q < 0.6}$ in a logarithmic
potential. This is understood from the fact that the VPs of these
models become contrived and double-peaked for large flattening ([Dehnen
& Gerhard, 1994](https://ui.adsabs.harvard.edu/abs/1994MNRAS.268.1019D/abstract)). 
In both test cases inaccurate results were
accompanied by non-zero values of ${h_1}$ and ${h_2}$, which should be zero by
definition. For the general case one may therefore take the values of
${h_1}$ and ${h_2}$ as an indicator of the numerical accuracy of the
VanderMonde algorithm. From this it was found that the algorithm (1)
fails only for the case of a logarithmic potential, large flattening
and large anisotropy.

If knowledge of the actual shape of the VP is required (rather than
just a measurement of Gauss-Hermite moments), then algorithm (2) is
the only possible solution. The algorithm (2) was implemented into the
software after the de Bruijne et al. paper was written, and no
detailed study of its accuracy was made. It is likely that the
algorithm can be optimized further, in particular the choice of the
optimal regularization parameter. Users may want work on this
themselves. It appears useful for many application to compare the
results of approach (1) to those of approach (2), as a test of the
accuracy of the results. Note that algorithm (2) is likely to be
inaccurate for the same set of parameter ranges for which algorithm
(1) is inaccurate (see above).