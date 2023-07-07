# ROTATING MODELS
***

To allow the construction of rotating models we adopt an odd part of
the form:
 
$f_o = f_e  (2s-1)  (\eta^2)^t$
 
The parameters $s$ and $t$ determine the amount of streaming in 
the model. The maximum streaming model has $t=0$ and $s = 0$ or $1$.
The total DF is ${f = f_e + f_o}$.
 
All the intrinsic velocity moments for the models can be expressed as
a power series in ${[e^2 \sin^2(\theta)]}$. The program `scalefree.f`
evaluates the velocity moments. These are then used to calculate the
projected line-of-sight and plane-of-sky velocity moments on the sky 
as one-dimensional
integrals. This can be done with essentially any given accuracy,
either by Gauss-Legendre quadrature or Romberg integration (to be
specified by the user). These parts of the program should work
flawlessly for all allowed input parameter values.