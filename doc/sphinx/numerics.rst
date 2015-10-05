Numerical technique
###################

MagIC is a pseudo-spectral MHD code. This numerical technique was originally
developed by P. Gilman and G. Glatzmaier for the spherical geometry.  In this
approach the unknowns are expanded into complete sets of functions in radial
and angular directions: Chebyshev polynomials in the radial directions and
spherical harmonic functions in the azimuthal and latitudinal directions.  This
allows to express all partial derivatives analytically.  Employing
orthogonality relations of spherical harmonic functions and using collocation
in radius then lead to algebraic equations that are integrated in time with a
mixed implicit/explicit time stepping scheme.  The nonlinear terms and the
Coriolis force are evaluated in the physical (or grid) space rather than in
spectral space.  Although this approach requires costly numerical
transformations between the two representations (from spatial to spectral using
Legendre and Fourier transforms), the resulting decoupling of all spherical
harmonic modes leads to a net gain in computational speed.  Before explaining
these methods in more detail, we introduce the poloidal/toroidal decomposition.


Poloidal/toroidal decomposition
===============================

Any vector :math:`\vec{v}` that fulfills  :math:`\vec{\nabla}\cdot\vec{v}=0`
can be decomposed in a poloidal and toroidal part :math:`W` and :math:`Z`,
respectively

.. math::
  \vec{v} = \vec{\nabla}\times\left(\vec{\nabla}\times W\,\vec{e_r}\right) +
  \vec{\nabla}\times Z\,\vec{e_r}.

Three unknown vector components are thus replaced by two scalar fields,
the poloidal potential :math:`W` and the toroidal potential :math:`Z`.
This decomposition is unique, aside from an arbitrary radial function  :math:`f(r)`
that can be added to :math:`W` or :math:`Z` without affecting :math:`\vec{v}`.

In the anelastic approximation, such a decomposition can be used for the
mass flux :math:`\tilde{\rho}\vec{u}` and the magnetic field :math:`\vec{B}`. This yields

.. math::
  \begin{aligned}
  \tilde{\rho} \vec{u} & = \vec{\nabla}\times(\vec{\nabla}\times W\,\vec{e_r}) +
  \vec{\nabla}\times Z\,\vec{e_r} \\
  \vec{B} & = \vec{\nabla}\times(\vec{\nabla}\times g\,\vec{e_r}) +
  \vec{\nabla}\times h\,\vec{e_r}.
  \end{aligned}
  :label: eqToroPolo

The two scalar potentials of a divergence free vector field can be extracted
from its radial component and the radial component of its curl:


.. math::
 \begin{aligned}
  \vec{e_r}\cdot \vec{v} &=  - \Delta_H W, \\
  \vec{e_r}\cdot\left(\vec{\nabla}\times\vec{B}\right) & = - \Delta_H Z,
 \end{aligned}

where the operator :math:`\Delta_H` denotes the horizontal part of the Laplacian:

.. math::
  \Delta_H= \frac{1}{r^{2}\sin{\theta}}
  \frac{\partial}{\partial\theta}\left(\sin{\theta}\frac{\partial}{\partial\theta}\right)
  + \frac{1}{r^{2}\sin^2{\theta}} \frac{\partial^{2}}{\partial^{2}\phi}.



Spherical harmonic representation
=================================

Spherical harmonic functions :math:`Y_{\ell m}` are a natural choice for the
horizontal expansion in colatitude :math:`\theta` and longitude :math:`\phi`:

.. math::
  Y_{\ell m}(\theta,\phi) = P_{\ell m}(\cos{\theta})\,e^{i m \phi},

where :math:`\ell` and :math:`m` denote spherical harmonic degree and order, respectively,
:math:`P_{\ell m}` is an associated Legendre function.  Different normalization are in
use. Here we adopt a complete normalization so that the orthogonality relation
reads 

.. math::
   \int_{0}^{2\pi} d\,\phi \int_{0}^{\pi}
  \sin{\theta}\, d\theta\; Y_{\ell m}(\theta,\phi)\,Y_{\ell^\prime
  m^\prime}(\theta,\phi) \; =  \; \delta_{\ell \ell^\prime}\delta_{m m^\prime}.

This means that

.. math::
  Y_{\ell m}(\theta,\phi) = \sqrt{\dfrac{1}{2\pi}}\dfrac{(2\ell+1)(\ell-|m|)!}{2(\ell+|m|)!}
  P_{\ell m}(\cos{\theta})\,e^{i m \phi}\,(-1)^m,

For example, the spherical harmonic representation of the
magnetic poloidal potential :math:`g(r,\theta,\phi)`, truncated at degree and order
:math:`\ell_{max}`, then reads

.. math::
  g(r,\theta,\phi) = \sum_{\ell=0}^{\ell_{max}}\sum_{m=-\ell}^{\ell} g_{\ell m}(r)\,Y_{\ell m}(\theta,\phi),
  :label: eqSpatSpec

with

.. math::
  g_{\ell m}(r) = \frac{1}{\pi}\,\int_{0}^{\pi} d \theta \sin{\theta}\; g_m(r,\theta)\;
  P_{\ell m}(\cos{\theta}),
  :label: eqLegTF1

.. math:: 
  g_{m}(r,\theta) = \frac{1}{2\pi}\,\int_{0}^{2\pi} d \phi\; g(r,\theta,\phi)\; e^{- i m \phi} .
  :label: eqLegTF2

The potential :math:`g(r,\theta,\phi)` is a real function so that
:math:`g_{\ell m}^\star(r)=g_{\ell,-m}(r)`, where the asterisk denotes the complex conjugate.
Thus, only coefficients with :math:`m \ge 0` have to be considered. The same kind of
expansion is made for the toroidal magnetic potential, the mass flux potentials,
pressure and entropy (or temperature).

The equations :eq:`eqLegTF1` and :eq:`eqLegTF2` define a two-step transform
from the longitude/latitude representation to the spherical harmonic
representation :math:`(r,\theta,\phi)\longrightarrow(r,\ell,m)`.  The equation
:eq:`eqSpatSpec` formulates the inverse procedure
:math:`(r,\ell,m)\longrightarrow(r,\theta,\phi)`. Fast-Fourier transforms are
employed in the longitudinal direction, requiring (at least) :math:`N_\phi = 2 \ell_{max}+1`
evenly spaced grid points :math:`\phi_i`.  
MagIC relies on the Gauss-Legendre quadrature for evaluating the integral
:eq:`eqLegTF1`


.. math::
   g_{\ell m}(r) = \frac{1}{N_{\theta}}
  \sum_{j=1}^{N_{\theta}}\,w_j\,g_m(r,\theta_j)\; P_{\ell m}(\cos{\theta_j}),

where :math:`\theta_j` are the :math:`N_{\theta}` Gaussian quadrature points
defining the latitudinal grid, and :math:`w_j` are the respective weights.  Pre-stored
values of the associated Legendre functions at grid points :math:`\theta_j` in
combination with a FFT in :math:`\phi` provide the inverse transform :eq:`eqSpatSpec`.
Generally, :math:`N_\phi=  2 N_\theta` is chosen, which provides
isotropic resolution in the equatorial region.  Choosing
:math:`\ell_{max}= [ \min(2 N_\theta,N_\phi)-1]/3` prevents aliasing errors.

Radial representation
=====================

In MagIC, the radial dependencies are expanded into complete sets of functions: the 
Chebyshev polynomials :math:`{\cal C}(x)`.  The polynomial of degree :math:`n` is defined by


.. math::
  {\cal C}_n(x)=\cos\left[n\,\arccos(x)\right]\quad -1\leq x \leq 1\,.

When truncating at degree :math:`N`, the radial expansion of the poloidal
magnetic potential reads

.. math::
  g_{\ell m}(r) = \sum_{n=0}^{N} g_{\ell mn}\;{\cal C}_n(r) ,
  :label: eqGridCheb

with

.. math::
   g_{\ell mn} = \frac{2-\delta_{n0}}{\pi}\int_{-1}^{1} 
   \frac{d x\, g_{\ell m}(r(x))\;{\cal C}_n(x)}{\sqrt{1-x^2}} .
  :label: eqSpecCheb

The Chebyshev definition space :math:`(-1\leq x\leq 1)` is then linearly mapped
onto a radius range :math:`(r_i\leq r \leq r_o)` by

.. math::
   x(r)=  2 \frac{r-r_i}{r_o-r_i} - 1 .

In addition, nonlinear mapping can be defined to modify the radial dependence of the
grid-point density.

When choosing the :math:`N_r` extrema of :math:`{\cal C}_{N_r-1}`  as radial grid points,

.. math::
  x_k=\cos{\left(\pi \frac{(k-1)}{N_r-1}\right)}\;\;\;,\;\;\; k=1,2,\ldots,N_r ,

the values of the Chebyshev polynomials at these points are simply given by
the cosine functions:

.. math::
  {\cal C}_{nk} = {\cal C}_n(x_k)=\cos{\left(\pi \frac{ n (k-1)}{N_r-1}\right)} .

This particular choice has two advantages.
For one, the grid points become denser toward the inner and outer
radius and better resolve potential thermal and viscous boundary layers.
In addition, FFTs can be employed to switch between
grid representation :eq:`eqGridCheb` and Chebyshev representations :eq:`eqSpecCheb`,
rendering this procedure a fast-Chebyshev transform.
Choosing :math:`N_r>N` provides radial dealiasing.
