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
  :label: eqDeltaH

where the operator :math:`\Delta_H` denotes the horizontal part of the Laplacian:

.. math::
  \Delta_H= \frac{1}{r^{2}\sin{\theta}}
  \frac{\partial}{\partial\theta}\left(\sin{\theta}\frac{\partial}{\partial\theta}\right)
  + \frac{1}{r^{2}\sin^2{\theta}} \frac{\partial^{2}}{\partial^{2}\phi}.
  :label: eqLaplaceH



Spherical harmonic representation
=================================

Spherical harmonic functions :math:`Y_\ell^m` are a natural choice for the
horizontal expansion in colatitude :math:`\theta` and longitude :math:`\phi`:

.. math::
  Y_\ell^m(\theta,\phi) = P_{\ell}^m(\cos{\theta})\,e^{i m \phi},

where :math:`\ell` and :math:`m` denote spherical harmonic degree and order, respectively,
:math:`P_\ell^m` is an associated Legendre function.  Different normalization are in
use. Here we adopt a complete normalization so that the orthogonality relation
reads 

.. math::
   \int_{0}^{2\pi} d\,\phi \int_{0}^{\pi}
   \sin{\theta}\, d\theta\; Y_\ell^m(\theta,\phi)\,Y_{\ell^\prime}^{m^\prime}
   (\theta,\phi) \; =  \; \delta_{\ell \ell^\prime}\delta^{m m^\prime}.
   :label: eqOrthoYlm

This means that

.. math::
  Y_{\ell}^{m}(\theta,\phi) = \sqrt{\dfrac{1}{2\pi}}\dfrac{(2\ell+1)(\ell-|m|)!}{2(\ell+|m|)!}
  P_\ell^m(\cos{\theta})\,e^{i m \phi}\,(-1)^m,

For example, the spherical harmonic representation of the
magnetic poloidal potential :math:`g(r,\theta,\phi)`, truncated at degree and order
:math:`\ell_{max}`, then reads

.. math::
  g(r,\theta,\phi) = \sum_{\ell=0}^{\ell_{max}}\sum_{m=-\ell}^{\ell} g_{\ell m}(r)\,Y_{\ell}^{m}(\theta,\phi),
  :label: eqSpatSpec

with

.. math::
  g_{\ell m}(r) = \frac{1}{\pi}\,\int_{0}^{\pi} d \theta \sin{\theta}\; g_m(r,\theta)\;
  P_\ell^m(\cos{\theta}),
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
  \sum_{j=1}^{N_{\theta}}\,w_j\,g_m(r,\theta_j)\; P_\ell^m(\cos{\theta_j}),

where :math:`\theta_j` are the :math:`N_{\theta}` Gaussian quadrature points
defining the latitudinal grid, and :math:`w_j` are the respective weights.  Pre-stored
values of the associated Legendre functions at grid points :math:`\theta_j` in
combination with a FFT in :math:`\phi` provide the inverse transform :eq:`eqSpatSpec`.
Generally, :math:`N_\phi=  2 N_\theta` is chosen, which provides
isotropic resolution in the equatorial region.  Choosing
:math:`\ell_{max}= [ \min(2 N_\theta,N_\phi)-1]/3` prevents aliasing errors.

Special relations
-----------------

The action of a horizontal Laplacian :eq:`eqLaplaceH` on spherical harmonics can be
analytically expressed by

.. math::
   \Delta_H Y_{\ell}^{m} = -\dfrac{\ell(\ell+1)}{r^2}\,Y_{\ell}^{m}\,.
   :label: eqHorizLaplYlm

They are several useful recurrence relations for the Legendre polynomials that will
be further employed to compute Coriolis forces and the :math:`\theta` and :math:`\phi`
derivatives of advection and Lorentz forces:

.. math::
   \cos\theta\,\partial P_\ell^m = c_{\ell+1}^m\,P_{\ell+1}^m
   +c_{\ell}^m\,P_{\ell-1}^m\,,
   :label: eqCosPlm

and

.. math::
   \sin\theta\,\dfrac{\partial P_\ell^m}{\partial \theta} = (\ell-1)\,c_{\ell+1}^m\,P_{\ell+1}^m
   -(\ell+2)\,c_{\ell}^m\,P_{\ell-1}^m \,,
   :label: eqSindPlm

where :math:`c_\ell^m` is defined by

.. math::
   c_\ell^m = \sqrt{\dfrac{(\ell+m)(\ell-m)}{(2\ell-1)(2\ell+1)}}\,.
   :label: eqClmOp



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
   :label: eqChebMap

In addition, nonlinear mapping can be defined to modify the radial dependence of the
grid-point density.

When choosing the :math:`N_r` extrema of :math:`{\cal C}_{N_r-1}`  as radial grid points,

.. math::
   x_k=\cos{\left(\pi \frac{(k-1)}{N_r-1}\right)}\;\;\;,\;\;\; k=1,2,\ldots,N_r ,
   :label: eqChebGrid

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

Spectral equations
==================

We have now introduced the necessary tools for deriving the
spectral equations.
Taking the **radial components** of the Navier-Stokes equation :eq:`eqNSNd`
and the induction equation :eq:`eqIndNd` provides the equations
for the poloidal potentials :math:`W(r,\theta,\phi)` and :math:`g(r,\theta,\phi)`.
The **radial component of the curl** of these equations provides
the equations for the toroidal counterparts
:math:`Z(r,\theta,\phi)` and :math:`h(r,\theta,\phi)`.
The pressure remains an additional unknown. Hence one more equation 
involving :math:`W_{\ell mn}` and :math:`p_{\ell mn}`
is required. It is obtained by taking the
**horizontal divergence** of the Navier-Stokes equation :eq:`eqNSNd`.

Expanding all potentials in spherical harmonics and Chebyshev polynomials,
multiplying with :math:`{Y_{\ell}^{m}}^\star`, and integrating over spherical surfaces
(while making use of
the orthogonality relation :eq:`eqOrthoYlm` results in equations for the
coefficients :math:`W_{\ell mn}`, :math:`Z_{\ell mn}`, :math:`g_{\ell mn}`, 
:math:`h_{\ell mn}`, :math:`p_{\ell mn}` and :math:`T_{\ell mn}`,
respectively.


Equation for the poloidal potential :math:`W`
---------------------------------------------

The temporal evolution of :math:`W` is obtained by taking :math:`\vec{e_r}\cdot` of each
term entering the Navier-Stokes equation :eq:`eqNSNd`. For the
time-derivative, one gets using :eq:`eqDeltaH`:

.. math::
   \tilde{\rho}\vec{e_r}\cdot\left(\dfrac{\partial \vec{u}}{\partial t}\right) =
   \dfrac{\partial}{\partial t}(\vec{e_r}\cdot\tilde{\rho}\vec{u})=-\Delta_H\dfrac{\partial
   W}{\partial t}.

For the viscosity term, one gets

.. math::
   \begin{aligned}
   \vec{e_r}\cdot\vec{\nabla}\cdot \mathsf{S} = & -\nu\,\Delta_H\left[\dfrac{\partial^2 W}
   {\partial r^2}
   +\left\lbrace 2\dfrac{d\ln\nu}{dr}-\dfrac{1}{3}\dfrac{d\ln\tilde{\rho}}{dr}\right\rbrace
   \dfrac{\partial W}{\partial r} \right. \\
   & -\left. \left\lbrace -\Delta_H + \dfrac{4}{3}\left(\dfrac{d^2\ln\tilde{\rho}}{dr^2}
   +\dfrac{d\ln\nu}{dr} \dfrac{d\ln\tilde{\rho}}{dr}  +
   \dfrac{1}{r}\left[3\dfrac{d\ln\nu}{dr}+
   \dfrac{d\ln\tilde{\rho}}{dr}\right] \right) \right\rbrace W\right],
   \end{aligned}

.. note:: In case of a constant kinematic viscosity, the :math:`d\ln\nu/dr`
          terms vanish. If in addition,the background density is constant, the
          :math:`d\ln\tilde{\rho}/dr` terms also vanish. In that Boussinesq limit, this
          viscosity term would then be simplified as

          .. math::
            \vec{e_r}\cdot\Delta \vec{u} = -\Delta_H\left[\dfrac{\partial^2 W}{\partial r^2}
            +\Delta_H\,W\right]

Using Eq. :eq:`eqHorizLaplYlm` then allows to finally write the time-evolution equation
for the poloidal potential :math:`W_{\ell m n}`:

.. math::
   \boxed{
   \begin{aligned}
   E\,\dfrac{\ell(\ell+1)}{r^2}\left[\left\lbrace\dfrac{\partial}{\partial t} + 
   \nu\,\dfrac{\ell(\ell+1)}{r^2} + \dfrac{4}{3}\,\nu\,\left(\dfrac{d^2\ln\tilde{\rho}}{dr^2}
   +\dfrac{d\ln\nu}{dr} \dfrac{d\ln\tilde{\rho}}{dr}  +
   \dfrac{1}{r}\left[3\dfrac{d\ln\nu}{dr}+
   \dfrac{d\ln\tilde{\rho}}{dr}\right] \right)\right\rbrace\right. & \,{\cal C}_n  & \\
   -\nu\,\left\lbrace 2\dfrac{d\ln\nu}{dr}-\dfrac{1}{3}\dfrac{d\ln\tilde{\rho}}{dr}\right\rbrace
   &\,{\cal C}'_n & \\
   -\nu & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& W_{\ell m n} \\
   + \left[{\cal C}'_n -\dfrac{d\ln\tilde{\rho}}{dr}{\cal C}_n\right] & & P_{\ell m n} \\
   - \left[\dfrac{Ra\,E}{Pr}\,\tilde{\rho}\,g(r)\right] & \,{\cal C}_n & s_{\ell m n} \\
   = {\cal N}^W = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \vec{F} & &
   \end{aligned}}
   :label: eqSpecW

Here, :math:`d\Omega` is the spherical surface element. We use the summation convention
for the Chebyshev index :math:`n`. The radial derivatives of Chebyshev
polynomials are denoted by primes.
   

Equation for the toroidal potential :math:`Z`
---------------------------------------------

The temporal evolution of :math:`Z` is obtained by taking the radial component of the
curl of the Navier-Stokes equation (i.e.  :math:`\vec{e_r}\cdot\vec{\nabla}\times`). For 
the time derivative, one gets using :eq:`eqDeltaH`:

.. math::
   \vec{e_r}\cdot\vec\nabla\times\left(\dfrac{\partial\tilde{\rho}\vec{u}}{\partial t}\right)=
   \dfrac{\partial}{\partial t}(\vec{e_r}\cdot\vec{\nabla}\times\tilde{\rho}
   \vec{u})=-\dfrac{\partial}{\partial t}(\Delta_H Z) = -\Delta_H\dfrac{\partial Z}{\partial t}

The pressure gradient, one has

.. math::
   \vec{\nabla}\times \left[\tilde{\rho}\vec{\nabla}\left(\dfrac{p'}{\tilde{\rho}}\right)\right] = 
   \vec{\nabla} \tilde{\rho} \times \vec{\nabla}\left(\dfrac{p'}{\tilde{\rho}}\right) + 
   \underbrace{\tilde{\rho} \vec{\nabla} \times \left[\vec{\nabla}\left( \dfrac{p'}{\tilde{\rho}}
   \right)\right]}_{=0}.

This expression has no component along :math:`\vec{e_r}`, as a consequence, there is
no pressure gradient contribution here. The
gravity term also vanishes as :math:`\vec{\nabla}\times(\tilde{\rho}g(r)\vec{e_r})` has no
radial component. 

.. math::
   \begin{aligned}
   \vec{e_r}\cdot\vec{\nabla}\times\left[\vec{\nabla}\cdot\mathsf{S}\right] = &
   -\nu\,\Delta_H\left[\dfrac{\partial^2 Z}{\partial r^2}
   +\left(\dfrac{d\ln\nu}{dr}-\dfrac{d\ln\tilde{\rho}}{dr}\right)\,\dfrac{\partial Z}{\partial r}  \right.\\
   & \left. - \left(\dfrac{d\ln\nu}{dr}\dfrac{d\ln\tilde{\rho}}{dr}+
     \dfrac{2}{r}\dfrac{d\ln\nu}{dr}+
     \dfrac{d^2\ln\tilde{\rho}}{dr^2}+\dfrac{2}{r}
   \dfrac{d\ln\tilde{\rho}}{dr}-\Delta_H\right) Z \right].
   \end{aligned}

.. note:: Once again, this viscous term can be greatly simplified in the Boussinesq limit:

          .. math::
            \vec{e_r}\cdot\vec{\nabla}\times\left(\Delta \vec{u}\right) = 
            -\Delta_H\left[\dfrac{\partial^2 Z}{\partial r^2}
            +\Delta_H\,Z\right]

Using Eq. :eq:`eqHorizLaplYlm` then allows to finally write the time-evolution equation
for the poloidal potential :math:`Z_{\ell m n}`:

.. math::
   \boxed{
   \begin{aligned}
   E\,\dfrac{\ell(\ell+1)}{r^2}\left[\left\lbrace\dfrac{\partial}{\partial t} + 
   \nu\,\dfrac{\ell(\ell+1)}{r^2} + \nu\,\left(\dfrac{d\ln\nu}{dr}\dfrac{d\ln\tilde{\rho}}{dr}+
   \dfrac{2}{r}\dfrac{d\ln\nu}{dr}+ \dfrac{d^2\ln\tilde{\rho}}{dr^2}+\dfrac{2}{r}
   \dfrac{d\ln\tilde{\rho}}{dr}\right)\right\rbrace\right. & \,{\cal C}_n  & \\
   -\nu\,\left(\dfrac{d\ln\nu}{dr}-\dfrac{d\ln\tilde{\rho}}{dr}\right) &\,{\cal C}'_n & \\
   -\nu & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& Z_{\ell m n} \\
   = {\cal N}^Z = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \left(\vec{\nabla}\times\vec{F}\right) & &
   \end{aligned}}
   :label: eqSpecZ


Equation for pressure :math:`P`
-------------------------------

The evolution of equation for pressure is obtained by taking the horizontal
divergence (i.e. :math:`\vec{\nabla}_H\cdot`)
of the Navier-Stokes equation. This operator is defined such
that

.. math::
   \vec{\nabla}_H\cdot\vec{a} = r\sin \dfrac{\partial (\sin\theta\,a_\theta)}{\partial \theta}
   +r\sin \dfrac{\partial a_\phi}{\partial \phi}.

This relates to the total divergence via:

.. math::
   \vec{\nabla}\cdot\vec{a}= \dfrac{1}{r^2}\dfrac{\partial(r^2 a_r)}{\partial r}+ 
   \vec{\nabla}_H\cdot\vec{a}.

The time-derivative term is thus expressed by

.. math::
   \begin{aligned} 
   \vec{\nabla}_H\cdot\left(\tilde{\rho}\dfrac{\partial \vec{u}}{\partial t}\right) 
   &= \dfrac{\partial}{\partial t}\left[\vec{\nabla}_H\cdot(\tilde{\rho}\vec{u}
   )\right] \\
   & =  \dfrac{\partial}{\partial t}\left[\vec{\nabla}\cdot(\tilde{\rho}\vec{u})
   -\dfrac{1}{r^2}\dfrac{\partial(r^2\tilde{\rho} u_r)}{\partial r}\right] \\
   & = -\dfrac{\partial}{\partial t}\left[\dfrac{\partial (\tilde{\rho} u_r)}{\partial r}
   +\dfrac{2\tilde{\rho} u_r}{r}\right] \\
   & = \dfrac{\partial}{\partial t}\left[\dfrac{\partial (\Delta_H W)}{\partial r}
   +\dfrac{2}{r}\Delta_H W\right] \\
   & = \Delta_H\dfrac{\partial}{\partial t}\left(\dfrac{\partial W}{\partial r}\right)
   \end{aligned}

We note that the gravity term vanishes since :math:`\vec{\nabla}_H\cdot(\tilde{\rho}
g(r)\vec{e_r}) = 0`. Concerning the pressure gradient, one has

.. math::
   -\vec{\nabla}_H\cdot\left[\tilde{\rho} \vec{\nabla}\left(\dfrac{p'}{\tilde{\rho}}
   \right)\right] = -\left\lbrace\vec{\nabla}\cdot\left[\tilde{\rho} \vec{\nabla}
   \left(\dfrac{p'}{\tilde{\rho}}\right)\right]-
   \dfrac{1}{r^2}\dfrac{\partial}{\partial r}\left[ r^2 \tilde{\rho} 
   \dfrac{\partial}{\partial r}\left(\dfrac{p'}{\tilde{\rho}}\right)\right] \right\rbrace = 
   -\Delta_H \, p'.

The viscosity term then reads

.. math::
  \begin{aligned}
  \vec{\nabla}_H\cdot \left( \vec{\nabla}\cdot\mathsf{S} \right) = & \nu\,\Delta_H\left[
  \dfrac{\partial^3 W}{\partial r^3} + \left(\dfrac{d\ln\nu}{dr}-
  \dfrac{d\ln\tilde{\rho}}{dr}\right) \dfrac{\partial^2 W}{\partial r^2} \right. \\
  & - \left[\dfrac{d^2\ln\tilde{\rho}}{dr^2} + \dfrac{d\ln\nu}{dr}\dfrac{d\ln\tilde{\rho}}{dr}+
  \dfrac{2}{r}\left(\dfrac{d\ln\nu}{dr}+\dfrac{d\ln\tilde{\rho}}{dr}\right)
  -\Delta_H \right]\dfrac{\partial W}{\partial r} \\
  & \left. -\left( \dfrac{2}{3}\dfrac{d\ln\tilde{\rho}}{dr}+\dfrac{2}{r}+\dfrac{d\ln\nu}{dr}
  \right)\Delta_H\,W \right].
  \end{aligned}

.. note:: Once again, this viscous term can be greatly simplified in the Boussinesq limit:

          .. math::
            \vec{\nabla}_H\cdot\left(\Delta \vec{u}\right) = 
            -\Delta_H\left[\dfrac{\partial^3 W}{\partial r^3}
            +\Delta_H\,\dfrac{\partial W}{\partial r}-\dfrac{2}{r}\Delta_H\,W\right]

Using Eq. :eq:`eqHorizLaplYlm` then allows to finally write the equation for the pressure
:math:`P_{\ell m n}`:

.. math::
   \boxed{
   \begin{aligned}
   E\,\dfrac{\ell(\ell+1)}{r^2}\left[
   -\nu\,\left( \dfrac{2}{3}\dfrac{d\ln\tilde{\rho}}{dr}+\dfrac{2}{r}+\dfrac{d\ln\nu}{dr}
   \right)\dfrac{\ell(\ell+1)}{r^2} \right.
   & \,{\cal C}_n  & \\
   \left\lbrace\dfrac{\partial}{\partial t} + 
   \nu\,\dfrac{\ell(\ell+1)}{r^2} + \nu\,\left[\dfrac{d^2\ln\tilde{\rho}}{dr^2}+
    \dfrac{d\ln\nu}{dr}\dfrac{d\ln\tilde{\rho}}{dr}+
   \dfrac{2}{r}\left(\dfrac{d\ln\nu}{dr}+\dfrac{d\ln\tilde{\rho}}{dr}\right)\right]\right\rbrace
   & \,{\cal C}'_n  & \\
   -\nu\,\left(  \dfrac{d\ln\nu}{dr}-\dfrac{d\ln\tilde{\rho}}{dr}
   \right) &\,{\cal C}''_n & \\
   -\nu & \,{\cal C}'''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& W_{\ell m n} \\
   + \left[\dfrac{\ell(\ell+1)}{r^2}\right] & \,{\cal C}_n & P_{\ell m n} \\
   = {\cal N}^P = -\int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{\nabla}_H\cdot\vec{F} & &
   \end{aligned}}
   :label: eqSpecP


.. note:: We note that the terms on the left hand side of :eq:`eqSpecW`, :eq:`eqSpecZ` and
          :eq:`eqSpecP` resulting from the viscous term, the pressure gradient,
          the buoyancy term, and the explicit time derivative completely decouple 
          in spherical harmonic degree and order.
          
          The terms that do not decouple, namely Coriolis force, Lorentz force and 
          advection of momentum, are collected on the right-hand side
          of :eq:`eqSpecW`, :eq:`eqSpecZ` and :eq:`eqSpecP` into the forcing term
          :math:`\vec{F}`:

          .. math::
             \vec{F}=-2\,\tilde{\rho}\,\vec{e_z}\times\vec{u} - E\,\tilde{\rho}\,
             \vec{u}\cdot\vec{\nabla}\,\vec{u} 
             +\frac{1}{Pm}\left(\vec{\nabla}\times\vec{B}\right)\times\vec{B}
             :label: eqForcing

Resolving :math:`\vec{F}` into potential functions is not required. Its
numerical evaluation is discussed :ref:`below <secNonlinearEqs>`.


Equation for entropy :math:`s`
------------------------------

The equation for the entropy (or temperature in the Boussinesq limit) is given by

.. math::
   \boxed{
   \begin{aligned}
   \dfrac{1}{Pr}\left[\left(Pr\dfrac{\partial}{\partial t} + 
   \kappa\,\dfrac{\ell(\ell+1)}{r^2} 
   \right)\right. & \,{\cal C}_n  & \\
   -\kappa\,\left(\dfrac{d\ln\kappa}{dr}+\dfrac{d\ln\tilde{\rho}}{dr}+
   +\dfrac{dln\tilde{T}}{dr}+\dfrac{2}{r}\right) 
   &\,{\cal C}'_n & \\
   -\kappa & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& s_{\ell m n} \\
   = {\cal N}^S = -\int d\Omega\,{Y_{\ell}^{m}}^\star\,\left[\vec{u}\cdot\vec{\nabla}s+
   \dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}\left(\Phi_\nu+
   \dfrac{\lambda}{Pm^2\,E}\,j^2\right) \right] & &
   \end{aligned}}
   :label: eqSpecS

In this expression, :math:`j=\vec{\nabla}\times\vec{B}` is the current. Once again,
the numerical evaluation of the right-hand-side (i.e. the non-linear terms) is
discussed :ref:`below <secNonlinearEqs>`.


Equation for the poloidal magnetic potential :math:`g`
------------------------------------------------------

The equation for the poloidal magnetic field coefficient reads


.. math::
   \boxed{
   \begin{aligned}
   \dfrac{\ell(\ell+1)}{r^2}\left[\left(\dfrac{\partial}{\partial t} + 
   \dfrac{1}{Pm}\lambda\,\dfrac{\ell(\ell+1)}{r^2} 
   \right)\right. & \,{\cal C}_n  & \\
   -\dfrac{1}{Pm}\,\lambda & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& g_{\ell m n} \\
   = {\cal N}^g = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \vec{D} & &
   \end{aligned}}
   :label: eqSpecG



Equation for the toroidal magnetic potential :math:`h`
------------------------------------------------------

The equation for the toroidal magnetic field coefficient reads

.. math::
   \boxed{
   \begin{aligned}
   \dfrac{\ell(\ell+1)}{r^2}\left[\left(\dfrac{\partial}{\partial t} + 
   \dfrac{1}{Pm}\lambda\,\dfrac{\ell(\ell+1)}{r^2} 
   \right)\right. & \,{\cal C}_n  & \\
   -\dfrac{1}{Pm}\,\dfrac{d\lambda}{dr} &\,{\cal C}'_n & \\
   -\dfrac{1}{Pm}\,\lambda & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& h_{\ell m n} \\
   = {\cal N}^h = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \left(\vec{\nabla}\times \vec{D}\right) & &
   \end{aligned}}
   :label: eqSpecH


We have now derived a full set of equations
:eq:`eqSpecW`, :eq:`eqSpecZ`, :eq:`eqSpecP`, :eq:`eqSpecS`, :eq:`eqSpecG` and
:eq:`eqSpecH`,
each describing the evolution of a single spherical harmonic mode of the
six unknown fields (assuming that the terms on the right hand side
are given). Each equation couples :math:`N+1` Chebyshev coefficients
for a given spherical harmonic mode :math:`(\ell,m)`.
Typically, a collocation method is employed to solve for the Chebyshev coefficients.
This means that the equations are required to be exactly satisfied at :math:`N-1`
grid points defined by the equations :eq:`eqChebMap` and :eq:`eqChebGrid`.
Excluded are the points :math:`r=r_i` and :math:`r=r_o`, where the 
:ref:`boundary conditions <secBoundaryConditions>` provide
additional constraints on the set of Chebyshev coefficients.


Time-stepping schemes
=====================

Implicit time stepping schemes theoretically offer increased stability and
allow for larger time steps.
However, fully implicit approaches have the disadvantage that
the nonlinear-terms couple all spherical harmonic modes.
The potential gain in computational speed is therefore lost at
higher resolution, where one very large matrix has to be dealt with
rather than a set of much smaller ones.
Similar considerations hold for the Coriolis force, one of
the dominating forces in the system and therefore a prime candidate for
implicit treatment. However, the Coriolis term couples modes :math:`(\ell,m,n)` with
:math:`(\ell+1,m,n)` and :math:`(\ell-1,m,n)` and also couples poloidal and
toroidal flow potentials. An implicit treatment of the Coriolis term therefore
also results in a much larger (albeit sparse) inversion matrix.

W consequently adopt in **MagIC** a mixed implicit/explicit algorithm.
Nonlinear and Coriolis terms, collected on the right hand side of equations
:eq:`eqSpecW`, :eq:`eqSpecZ`, :eq:`eqSpecP`, :eq:`eqSpecS`, :eq:`eqSpecG`
and :eq:`eqSpecH` are treated explicitly with a second order
`Adams-Bashforth <https://en.wikipedia.org/wiki/Linear_multistep_method>`_ . 
Terms collected on the left hand side are
time-stepped with an implicit modified `Crank-Nicolson
<https://en.wikipedia.org/wiki/Crank–Nicolson_method>`_ algorithm.
While the equations are coupled radially, they decouple for all spherical
harmonic modes. 

.. note::  The poloidal flow potential :eq:`eqSpecW` and the pressure :eq:`eqSpecP`
           are nevertheless coupled for a given spherical harmonic mode.

As an example, we derive the time stepping equation for the poloidal
magnetic potential of degree :math:`\ell` and order :math:`m`,
denoting the explicit nonlinear term at radial grid point :math:`r_k` with

.. math::
  D_{k\ell m}(t)= \int d\Omega\; {Y_{\ell}^{m}}^\star\; \vec{e_r} \cdot \vec{D}(t,r_k,\theta,\phi)\;\; .

After discretization of the partial time derivative,
:math:`\partial g_{\ell mn}/\partial t = [g_{\ell mn}(t+\delta t) - g_{\ell mn}(t)]/\delta t`
where :math:`\delta t` is the time step, we can formulate the left hand side
of :eq:`eqSpecG` as a matrix multiplication. The matrices :math:`\mathsf{A}`  and 
:math:`\mathsf{G}` are defined by

.. math::
    {A}_{kn} = \dfrac{\ell (\ell+1)}{r_k^2}\,\dfrac{1}{\delta t} {\cal C}_{nk}\;\;

and

.. math::
    {G}_{kn}=\dfrac{\ell(\ell+1)}{r_k^2}\,\dfrac{1}{Pm}\left( \dfrac{\ell(\ell+1}{r_k^2} 
    {\cal C}_{nk}-{\cal C}''_{nk} \right)\;\;,

where :math:`{\cal C}_{nk}={\cal C}_n(r_k)`. The matrices depend on :math:`\ell` 
but not on :math:`m`.  Advancing time from :math:`t` to :math:`t+\delta t` is 
then a matter of solving

.. math::
      \left( {A}_{kn} + \alpha {G}_{kn}\right)\;g_{\ell mn}(t+\delta t) =
      \left( {A}_{kn} - (1 - \alpha) {G}_{kn} \right)\;g_{\ell mn}(t) +
      \frac{3}{2} D_{k\ell m}(t) - \frac{1}{2} D_{k\ell m}(t-\delta t)\;\;.

The classical Crank-Nicholson scheme is recovered for :math:`\alpha=0.5`, but
it seems that a slightly larger weight of :math:`\alpha=0.6` helps to stabilize
the time integration.  Since the stability requirements limiting :math:`\delta
t` will usually change during a computational run, the time step should be
adjusted accordingly.  The matrix :math:`\mathsf{G}` remains unchanged, but
:math:`\mathsf{A}` has to be updated whenever :math:`\delta t` is changed.
This, in turn, requires a new triangulation of matrix :math:`A_{kn}+\alpha G_{kn}`,
which is then stored for subsequent time steps until the next adjustment of
:math:`\delta t` is in order. 

**Courant's condition** offers a guideline
concerning the value of :math:`\delta t`, demanding that :math:`\delta t` should be smaller
than the advection time between two grid points.  Strong Lorentz forces require
an additional stability criterion that is obtained by replacing the flow speed
by Alfvén's velocity in a modified Courant criterion.
The explicit treatment of the Coriolis force requires that the time step is
limited to a fraction of the rotation period, which may be the relevant
criterion at low Ekman number when flow and magnetic field remain weak.
Non-homogeneous grids and other numerical effects generally require an
additional safety factor in the choice of :math:`\delta t`.


.. _secNonlinearEqs:

Coriolis force and non-linear terms
===================================

.. math::
   {\cal A}_r=
   -\tilde{\rho}\left(
   u_r\dfrac{\partial u_r}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_r}{\partial \theta}+
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_r}{\partial \phi}
  -\dfrac{u_\theta^2+u_\phi^2}{r}\right)+
   \dfrac{1}{Pm}\left(j_\theta\,B_\phi-j_\phi\,B_\theta\right)\, ,
   :label: eqAdvTheta

.. math::
   {\cal A}_\theta=
   -\tilde{\rho}\left(
   u_r\dfrac{\partial u_\theta}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_\theta}{\partial \theta} +
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_\theta}{\partial \phi}+
   \dfrac{u_r u_\theta}{r}-\dfrac{\cos\theta}{r\sin\theta}u_\phi^2\right)+
   \dfrac{1}{Pm}\left(j_\phi\,B_r-j_r\,B_\phi\right)\, ,
   :label: eqAdvTheta

.. math::
   {\cal A}_\phi=
   -\tilde{\rho}\left(
   u_r\dfrac{\partial u_\phi}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_\phi}{\partial \theta} +
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_\phi}{\partial \phi}+
   \dfrac{u_r u_\phi}{r} +\dfrac{\cos\theta}{r\sin\theta}u_\theta u_\phi\right)+
   \dfrac{1}{Pm}\left(j_r\,B_\theta-j_\theta\,B_r\right)\, ,
   :label: eqAdvPhi

.. _secBoundaryConditions:

Boundary conditions and inner core
==================================

Mechanical boundary conditions
------------------------------

Since the system of equations is formulated on a radial grid, boundary
conditions can simply be satisfied by replacing the collocation equation
at grid points :math:`r_i` and :math:`r_o` with appropriate expressions.
The condition of zero radial flow on the boundaries implies

.. math::
  {\cal C}_n(r) W_{\ell mn} = 0 \;\;\mbox{at}\;\; r=r_i,r_o\;\;.
  :label: eqBcRigid1

Note that the summation convention with respect to
radial modes :math:`n` is used again.
**The no-slip** condition further requires that the
horizontal flow components also have to vanish, provided
the two boundaries are at rest. This condition is fulfilled when

.. math::
   {\cal C}'_n(r) W_{\ell mn} = 0\;\;\mbox{at}\;\; r=r_i,r_o
  :label: eqBcRigid2

and

.. math::
   {\cal C}_n(r) Z_{\ell mn} = 0\;\;\mbox{at}\;\; r=r_i,r_o
  :label: eqBcRigid3

for all spherical harmonic modes :math:`(\ell,m)`.
The conditions :eq:`eqBcRigid1`-:eq:`eqBcRigid3`
replace the poloidal flow potential equations :eq:`eqSpecW`
and the pressure equation :eq:`eqSpecP`, respectively, at
the collocation points :math:`r_i` and :math:`r_o`.

If the inner-core and/or the mantle are allowed to react to torques,
a condition based on the conservation of
angular momentum replaces condition :eq:`eqBcRigid3` for the mode
:math:`(\ell =1,m=0)`:

.. math::
  \mathsf{I} \dfrac{\partial\vec{\omega}}{\partial t}= \vec{\Gamma}\;\;.

The tensor :math:`\mathsf{I}` denotes the moment of inertia of inner core or mantle,
respectively, :math:`\vec{\omega}` is the mantle or inner-core rotation rate relative
to that of the reference frame, and :math:`\vec{\Gamma}` is the respective torque.

**Free-slip boundary conditions** require that the viscous stress vanishes, which
in turn implies that the non-diagonal components :math:`\mathsf{Sr}_{r\phi}` and
:math:`\mathsf{S}_{r\theta}` of the rate-of-strain tensor vanish.
Translated to the spectral representation this requires

.. math::
  \left[{\cal C}''_n(r) -\left(\frac{2}{r}+\dfrac{d\ln\tilde{\rho}}{dr}\right)\,{\cal C}'_n(r)
  \right] W_{\ell mn} = 0 \;\;\mbox{and}\;\;
  \left[{\cal C}'_n(r) -\left(\frac{2}{r}+\dfrac{d\ln\tilde{\rho}}{dr}\right)\,{\cal C}_n(r)
  \right] z_{\ell mn} = 0\;.


Magnetic boundary conditions and inner core
-------------------------------------------

Magnetic boundary conditions at the interface with an insulating mantle
or insulating inner core are similarly implemented.
The toroidal magnetic field cannot enter any
insulator and therefore has to vanish at the boundary

.. math::
  {\cal C}_n(r) h_{\ell mn} = 0\;\;\mbox{at}\;\; r=r_i\;\;\mbox{and/or}\;\; r=r_o\;\;.

Matching conditions for the poloidal magnetic field with a source-free
external potential field require that the following equations are satisfied at
the boundary grid points:

.. math::
  {\cal C}_n^\prime(r) g_{\ell mn} - {\cal C}_n(r)\frac{\ell+1}{r} g_{\ell mn} = 0 \;\;\;\mbox{at}\;\;\; r=r_i,

.. math::
  {\cal C}_n^\prime(r) g_{\ell mn} + {\cal C}_n(r)\frac{\ell}{r} g_{\ell mn} = 0 \;\;\;\mbox{at}\;\;\; r=r_o.

If the inner core is modeled as an electrical conductor, a simplified dynamo
equation has to be solved in which the fluid flow is replaced by the
solid-body rotation of the inner core. The latter is described by a single toroidal
flow mode :math:`(\ell=1,m=0)`. The resulting nonlinear terms can be expressed by a simple
spherical harmonic expansion, where the superscript :math:`I` denotes values in the
inner core and :math:`\omega_I` its differential rotation rate:

.. math::
   \int d\Omega\; {Y_{\ell}^{m}}^\star\;\vec{e_r}\cdot\left[\vec{\nabla}\times
   \left(\vec{u^I}\times\vec{B^I}\right)\right] =
   - i\,\omega_I\,m\,\dfrac{\ell(\ell+1)}{r^2}\;g_{\ell m}^I(r)\; ,
   :label: glmnI

.. math::
   \int d\Omega\; {Y_{\ell}^{m}}^\star\;\vec{e_r}\cdot\left[\vec{\nabla}\times
   \vec{\nabla}\times\left(\vec{u^I}\times\vec{B^I}\right)\right] =
   - i\,\omega_I\,m\,\dfrac{\ell(\ell+1)}{r^2} \;h_{\ell m}^I(r)\;
   :label: hlmnI

The expensive back and forth transformations between spherical
harmonic and grid representations are therefore not required for advancing the
inner-core magnetic field in time.

In the inner core the magnetic potentials are again conveniently 
expanded into Chebyshev polynomials. The Chebyshev variable :math:`x`
spans the whole diameter of the inner core, so that grid points are dense 
near the inner-core boundary but sparse in the center. The mapping is given by:

.. math::
   x(r)=  \dfrac{r}{r_i}\;\;,\;\;-r_i\le r\le r_i\;\;.
   :label: mapRIC

Each point in the inner core is thus represented twice, by
grid points :math:`({r,\theta,\phi})` and :math:`({-r,\pi-\theta,\phi+\pi})`.
Since both representations must be identical, this imposes a symmetry
constraint that can be fulfilled when the radial expansion
comprises only polynomials of even order:

.. math::
  g_{\ell m}^I(r) = \left(\dfrac{r}{r_i}\right)^{\ell+1}
                    \,\sum_{i=0}^{M-1} g_{\ell m\,2i}^I\;{\cal C}_{2i}(r)\;\;.
  :label: radI

An equivalent expression holds for the toroidal potential in the inner core.
FFTs can again by employed efficiently
for the radial transformation, using the :math:`M` extrema of
:math:`{\cal C}_{2 M-1}(r)` with :math:`x>0` as grid points.

The sets of spectral magnetic field equations for the inner and the outer core
are coupled via continuity equations for the magnetic field and the
horizontal electric field.
Continuity of the magnetic field is assured by (**i**) continuity of the
toroidal potential, (**ii**) continuity of the poloidal potential, and
(**iii**) continuity of the radial derivative of
the latter. Continuity of the horizontal electric field demands (**iv**) that
the radial derivative of the toroidal potential is continuous, provided
that the horizontal flow and the electrical conductivity are continuous
at the interface.
These four conditions replace the spectral equations
:eq:`eqSpecG`, :eq:`eqSpecH` on the outer-core side and
equations :eq:`glmnI`, :eq:`hlmnI` on the inner-core side.
Employing free-slip conditions or allowing for electrical conductivity differences
between inner and outer core leads to more complicated and even non-linear matching
conditions.

