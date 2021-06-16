.. _secNumerics:

Numerical technique
###################

**MagIC** is a pseudo-spectral MHD code. This numerical technique was
originally developed by P. Gilman and G. Glatzmaier for the spherical geometry.
In this approach the unknowns are expanded into complete sets of functions in
radial and angular directions: **Chebyshev polynomials or Finite differences in
the radial direction** and **spherical harmonic functions in the azimuthal and
latitudinal directions**. This allows to express all partial derivatives
analytically.  Employing orthogonality relations of spherical harmonic
functions and using collocation in radius then lead to algebraic equations that
are integrated in time with a **mixed implicit/explicit time stepping scheme**.
The nonlinear terms and the Coriolis force are evaluated in the physical (or
grid) space rather than in spectral space. Although this approach requires
costly numerical transformations between the two representations (from spatial
to spectral using Legendre and Fourier transforms), the resulting decoupling of
all spherical harmonic modes leads to a net gain in computational speed.
Before explaining these methods in more detail, we introduce the
poloidal/toroidal decomposition.


Poloidal/toroidal decomposition
===============================

Any vector :math:`\vec{v}` that fulfills  :math:`\vec{\nabla}\cdot\vec{v}=0`, i.e.
a so-called *solenoidal field*,
can be decomposed in a poloidal and a toroidal part :math:`W` and :math:`Z`,
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
  \boxed{
  \begin{aligned}
  \tilde{\rho} \vec{u} & = \vec{\nabla}\times(\vec{\nabla}\times W\,\vec{e_r}) +
  \vec{\nabla}\times Z\,\vec{e_r}\,, \\
  \vec{B} & = \vec{\nabla}\times(\vec{\nabla}\times g\,\vec{e_r}) +
  \vec{\nabla}\times h\,\vec{e_r}\,.
  \end{aligned}}
  :label: eqToroPolo

The two scalar potentials of a divergence free vector field can be extracted
from its radial component and the radial component of its curl using the fact that
the toroidal field has not radial component: 

.. math::
  \begin{aligned}
  \vec{e_r}\cdot \tilde{\rho}\vec{u} &=  - \Delta_H\,W, \\
  \vec{e_r}\cdot\left(\vec{\nabla}\times\vec{u}\right) & = - \Delta_H\,Z,
  \end{aligned}
  :label: eqDeltaH

where the operator :math:`\Delta_H` denotes the horizontal part of the Laplacian:

.. math::
  \Delta_H= \frac{1}{r^{2}\sin{\theta}}
  \frac{\partial}{\partial\theta}\left(\sin{\theta}\frac{\partial}{\partial\theta}\right)
  + \frac{1}{r^{2}\sin^2{\theta}} \frac{\partial^{2}}{\partial^{2}\phi}.
  :label: eqLaplaceH


The equation :eq:`eqToroPolo` can be expanded in spherical coordinates.
The three components of :math:`\tilde{\rho}\vec{u}`
are given by

.. math::
   \tilde{\rho}\vec{u} = -(\Delta_H W)\,\vec{e_r} + \left( \dfrac{1}{r}
   \dfrac{\partial^2 W}{\partial r \partial \theta} + 
   \dfrac{1}{r\sin\theta}\dfrac{\partial Z}{\partial \phi}\right)\,\vec{e_\theta} 
   +\left(\dfrac{1}{r\sin\theta}\dfrac{\partial^2 W}{\partial r\partial \phi}-
   \dfrac{1}{r}\dfrac{\partial Z}{\partial\theta} \right)\,\vec{e_\phi},
   :label: eqToroPolo1

while the curl of :math:`\tilde{\rho}\vec{u}` is expressed by

.. math::
   \begin{aligned}
   \vec{\nabla}\times\tilde{\rho}\vec{u} = &-(\Delta_H Z)\,\vec{e_r}+
   \left[-\dfrac{1}{r\sin\theta}\dfrac{\partial}{\partial\phi}\left(
   \dfrac{\partial^2}{\partial r^2}+\Delta_H \right) W + \dfrac{1}{r}
   \dfrac{\partial^2 Z}{\partial r\partial\theta}\right]\vec{e_\theta} \\
   &+\left[\dfrac{1}{r}\dfrac{\partial}{\partial\theta}\left(
   \dfrac{\partial^2}{\partial r^2}+\Delta_H\right) W + \dfrac{1}{r \sin\theta} 
   \dfrac{\partial^2 Z}{\partial r\partial\phi}\right]\vec{e_\phi},
   \end{aligned}
 :label: eqToroPolo2
 
Using the horizontal part of the divergence operator

.. math::
   \vec{\nabla}_H = \dfrac{1}{r\sin\theta} \dfrac{\partial}{\partial \theta}\sin\theta\;\vec{e}_\theta 
   + \dfrac{1}{r\sin\theta} \dfrac{\partial}{\partial \phi}\;\vec{e}_\phi

above expressions can be simplified to 

.. math::
   \tilde{\rho}\vec{u} = -\Delta_H\;\vec{e_r}\; W + \vec{\nabla}_H \dfrac{\partial}{\partial r}\;W 
                         + \vec{\nabla}_H\times\vec{e}_r\;Z
                         
and

.. math::
   \nabla\times\tilde{\rho}\vec{u} = -\Delta_H\;\vec{e}_r\;Z + \vec{\nabla}_H \dfrac{\partial}{\partial r}\;Z 
                         - \vec{\nabla}_H\times\Delta_H\vec{e}_r\;W\;\;.

Below we will use the fact that the horizontal components of the poloidal field depend 
on the radial derivative of the poloidal potential. 

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
  Y_{\ell}^{m}(\theta,\phi) = \left(\dfrac{(2\ell+1)}{4\pi}\dfrac{(\ell-|m|)!}{(\ell+|m|)!}\right)^{1/2}
  P_\ell^m(\cos{\theta})\,e^{i m \phi},

As an example, the spherical harmonic representation of the
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
pressure, entropy (or temperature) and chemical composition.

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

.. seealso:: In MagIC, the Legendre functions are defined in the subroutine
             :f:subr:`plm_theta <plms_theta/plm_theta()>`. The Legendre transforms
             from spectral to grid space are computed in the module
             :f:mod:`legendre_spec_to_grid`, while the backward transform (from grid
             space to spectral space) are computed in the module
             :f:mod:`legendre_grid_to_spec`. The fast Fourier transforms are computed
             in the module :f:mod:`fft`.

Special recurrence relations
----------------------------

The action of a horizontal Laplacian :eq:`eqLaplaceH` on spherical harmonics can be
analytically expressed by

.. math::
   \boxed{
   \Delta_H Y_{\ell}^{m} = -\dfrac{\ell(\ell+1)}{r^2}\,Y_{\ell}^{m}\,.
   }
   :label: eqHorizLaplYlm


They are several useful recurrence relations for the Legendre polynomials that will
be further employed to compute Coriolis forces and the :math:`\theta` and :math:`\phi`
derivatives of advection and Lorentz forces.
Four of these operators are used in **MagIC**. The first one is defined by

.. math::
   \vartheta_1 = \dfrac{1}{\sin\theta}\dfrac{\partial}{\partial\theta}\sin^2\theta
   =\sin\theta\dfrac{\partial}{\partial\theta}+2\cos\theta\,.

The action of this operator on a Legendre polynomials is given by

.. math::
   \vartheta_1 = (\ell+2)\,c_{\ell+1}^m\,P_{\ell+1}^m(\cos\theta)
   -(\ell-1)\,c_\ell^m\,P_{\ell-1}^m(\cos\theta),

where :math:`c_\ell^m` is defined by

.. math::
   c_\ell^m = \sqrt{\dfrac{(\ell+m)(\ell-m)}{(2\ell-1)(2\ell+1)}}\,.
   :label: eqClmOp

*How is that implemented in the code?* Let's assume we want the spherical harmonic contribution
of degree :math:`\ell` and order `m` for the expression

.. math::
   \dfrac{1}{\sin\theta}\dfrac{\partial}{\partial\theta}(\sin\theta\,f(\theta))\,.

In order to employ the operator :math:`\vartheta_1` for the derivative, we thus define a
new function

.. math::
   F(\theta)=f(\theta)/\sin\theta\,,

so that

.. math::
   \dfrac{1}{\sin\theta}\dfrac{\partial}{\partial\theta}[\sin\theta\,f(\theta)]
   =\vartheta_1 F(\theta)\,.

Expanding :math:`F(\theta)` in Legendre polynomials and using the respective
orthogonality relation we can then map out the required contribution in the following way:

.. math::
  \boxed{
  \int_0^\pi d\theta\,\sin\theta\,P_\ell^m\vartheta_1\sum_{\ell'}F_{\ell'}^m P_{\ell'}^m
  =(\ell+1)\,c_{\ell}^m\,F_{\ell-1}^m-\ell\,c_{\ell+1}^m\,F_{\ell+1}^m\,.}
  :label: eqOpTheta1

Here, we have assumed that the Legendre functions are completely normalized such that

.. math::
   \int_0^\pi d\theta\,\sin\theta\,P_\ell^m P_{\ell'}^m = \delta_{\ell \ell'}\,.

.. seealso:: This operator is defined in the module :f:mod:`horizontal_data` by the variables
   :f:var:`dTheta1S <dtheta1s>` for the first part of the right-hand 
   side of :eq:`eqOpTheta1` and :f:var:`dTheta1A <dtheta1a>` for the 
   second part.

The second operator used to formulate colatitude derivatives is

.. math::
   \vartheta_2 = \sin\theta\dfrac{\partial}{\partial\theta}\,.

The action of this operator on the Legendre polynomials reads

.. math::
   \vartheta_2 P_\ell^m(\cos\theta)=\ell\,c_{\ell+1}^m\,P_{\ell+1}^m(\cos\theta)
   -(\ell+1)\,c_\ell^m\,P_{\ell-1}^m(\cos\theta)\,,

so that

.. math::
   \boxed{
   \int_0^\pi d\theta\,\sin\theta \,P_\ell^m\vartheta_2\sum_{\ell'}f_{\ell'}^m P_{\ell'}^m
   =(\ell-1)\,c_{\ell}^m\,f_{\ell-1}^m-(\ell+2)\,c_{\ell+1}^m\,f_{\ell+1}^m\,.}
  :label: eqOpTheta2

.. seealso:: This operator is defined in the module :f:mod:`horizontal_data` by the variables
   :f:var:`dTheta2S <dtheta2s>` for the first part of the right-hand 
   side of :eq:`eqOpTheta2` and :f:var:`dTheta2A <dtheta2a>` for the 
   second part.


The third combined operator is defined by:

.. math::
   \vartheta_3 = \sin\theta\dfrac{\partial}{\partial\theta}+\cos\theta\,L_H\,,

where :math:`-L_H/r^2=\Delta_H`.

Acting with :math:`\vartheta_3` on a Legendre function gives:

.. math::
   \vartheta_3 P_\ell^m(\cos\theta)=\ell(\ell+1)\,c_{\ell+1}^m\,P_{\ell+1}^m(\cos\theta)
   +(\ell-1)(\ell+1)\,c_\ell^m\,P_{\ell-1}^m(\cos\theta)\,,

which results into:

.. math::
  \boxed{
  \int_0^\pi d\theta\,\sin\theta\,P_\ell^m\vartheta_3\sum_{\ell'}f_{\ell'}^m P_{\ell'}^m
  =(\ell-1)(\ell+1)\,c_{\ell}^m\,f_{\ell-1}^m+\ell(\ell+2)\,c_{\ell+1}^m\,f_{\ell+1}^m\,.}
  :label: eqOpTheta3

.. seealso:: This operator is defined in the module :f:mod:`horizontal_data` by the variables
   :f:var:`dTheta3S <dtheta3s>` for the first part of the right-hand 
   side of :eq:`eqOpTheta3` and :f:var:`dTheta3A <dtheta3a>` for the 
   second part.


The fourth (and last) combined operator is defined by:

.. math::
   \vartheta_4 = \dfrac{1}{\sin\theta}\dfrac{\partial}{\partial\theta}\sin^2\theta\,L_H
   =\vartheta1\,L_H\,,

Acting with :math:`\vartheta_3` on a Legendre function gives:

.. math::
   \vartheta_4 P_\ell^m(\cos\theta)=\ell(\ell+1)(\ell+2)\,c_{\ell+1}^m\,P_{\ell+1}^m(\cos\theta)
   -\ell(\ell-1)(\ell+1)\,c_\ell^m\,P_{\ell-1}^m(\cos\theta)\,,

which results into:

.. math::
  \boxed{
  \int_0^\pi d\theta\,\sin\theta\,P_\ell^m\vartheta_4\sum_{\ell'}f_{\ell'}^m P_{\ell'}^m
  =\ell(\ell-1)(\ell+1)\,c_{\ell}^m\,f_{\ell-1}^m-\ell(\ell+1)(\ell+2)\,c_{\ell+1}^m\,f_{\ell+1}^m\,.}
  :label: eqOpTheta4

.. seealso:: This operator is defined in the module :f:mod:`horizontal_data` by the variables
   :f:var:`dTheta4S <dtheta4s>` for the first part of the right-hand 
   side of :eq:`eqOpTheta4` and :f:var:`dTheta4A <dtheta4a>` for the 
   second part.



Radial representation
=====================

In MagIC, the radial dependencies are either expanded into complete sets of functions, the 
Chebyshev polynomials :math:`{\cal C}(x)`; or discretized using finite differences. For the former approach, the Chebyshev polynomial of degree :math:`n` is defined by


.. math::
  {\cal C}_n(x)\approx\cos\left[n\,\arccos(x)\right]\quad -1\leq x \leq 1\,.

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
   x_k=\cos{\left[\pi \frac{(k-1)}{N_r-1}\right]}\;\;\;,\;\;\; k=1,2,\ldots,N_r ,
   :label: eqChebGrid

the values of the Chebyshev polynomials at these points are simply given by
the cosine functions:

.. math::
  {\cal C}_{nk} = {\cal C}_n(x_k)=\cos{\left[\pi \frac{ n (k-1)}{N_r-1}\right]} .

This particular choice has two advantages.
For one, the grid points become denser toward the inner and outer
radius and better resolve potential thermal and viscous boundary layers.
In addition, type I Discrete Cosine Transforms (DCTs) can be employed to switch between
grid representation :eq:`eqGridCheb` and Chebyshev representations :eq:`eqSpecCheb`,
rendering this procedure a fast-Chebyshev transform.

.. seealso:: The Chebyshev (Gauss-Lobatto) grid is defined in the module
             :f:mod:`chebyshev_polynoms_mod`. The cosine transforms are computed in the
             modules :f:mod:`cosine_transform_even` and :f:mod:`fft_fac_mod`.

Spectral equations
==================

We have now introduced the necessary tools for deriving the
spectral equations.
Taking the **radial components** of the Navier-Stokes equation
and the induction equation provides the equations
for the poloidal potentials :math:`W(r,\theta,\phi)` and :math:`g(r,\theta,\phi)`.
The **radial component of the curl** of these equations provides
the equations for the toroidal counterparts
:math:`Z(r,\theta,\phi)` and :math:`h(r,\theta,\phi)`.
The pressure remains an additional unknown. Hence one more equation 
involving :math:`W_{\ell mn}` and :math:`p_{\ell mn}`
is required. It is obtained by taking the
**horizontal divergence** of the Navier-Stokes equation.

Expanding all potentials in spherical harmonics and Chebyshev polynomials,
multiplying with :math:`{Y_{\ell}^{m}}^\star`, and integrating over spherical surfaces
(while making use of
the orthogonality relation :eq:`eqOrthoYlm` results in equations for the
coefficients :math:`W_{\ell mn}`, :math:`Z_{\ell mn}`, :math:`g_{\ell mn}`, 
:math:`h_{\ell mn}`, :math:`P_{\ell mn}` and :math:`s_{\ell mn}`,
respectively.


Equation for the poloidal potential :math:`W`
---------------------------------------------

The temporal evolution of :math:`W` is obtained by taking :math:`\vec{e_r}\cdot` of each
term entering the Navier-Stokes equation. For the
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
            +\Delta_H\,W\right]\,.

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
   - \left[\dfrac{Ra_\xi\,E}{Sc}\,\tilde{\rho}\,g(r)\right] & \,{\cal C}_n & \xi_{\ell m n} \\
   = {\cal N}^W_{\ell m} = \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^W =\int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \vec{F}\,. & &
   \end{aligned}}
   :label: eqSpecW

Here, :math:`d\Omega` is the spherical surface element. We use the summation convention
for the Chebyshev index :math:`n`. The radial derivatives of Chebyshev
polynomials are denoted by primes.

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecW` are coded in
             the subroutines :f:subr:`get_wpMat <updatewp_mod/get_wpmat()>`
   

Equation for the toroidal potential :math:`Z`
---------------------------------------------

The temporal evolution of :math:`Z` is obtained by taking the radial component of the
curl of the Navier-Stokes equation (i.e.  :math:`\vec{e_r}\cdot\vec{\nabla}\times`). For 
the time derivative, one gets using :eq:`eqDeltaH`:

.. math::
   \vec{e_r}\cdot\vec\nabla\times\left(\dfrac{\partial\tilde{\rho}\vec{u}}{\partial t}\right)=
   \dfrac{\partial}{\partial t}(\vec{e_r}\cdot\vec{\nabla}\times\tilde{\rho}
   \vec{u})=-\dfrac{\partial}{\partial t}(\Delta_H Z) =
   -\Delta_H\dfrac{\partial Z}{\partial t}\,.

The pressure gradient, one has

.. math::
   \vec{\nabla}\times \left[\tilde{\rho}\vec{\nabla}\left(\dfrac{p'}{\tilde{\rho}}\right)\right] = 
   \vec{\nabla} \tilde{\rho} \times \vec{\nabla}\left(\dfrac{p'}{\tilde{\rho}}\right) + 
   \underbrace{\tilde{\rho} \vec{\nabla} \times \left[\vec{\nabla}\left( \dfrac{p'}{\tilde{\rho}}
   \right)\right]}_{=0}\,.

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
            +\Delta_H\,Z\right]\,.

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
   = {\cal N}^Z_{\ell m} = \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^Z = 
   \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \left(\vec{\nabla}
   \times\vec{F}\right)\,. & &
   \end{aligned}}
   :label: eqSpecZ

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecZ` are coded in
             the subroutines :f:subr:`get_zMat <updatez_mod/get_zmat()>`


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
   )\right], \\
   & =  \dfrac{\partial}{\partial t}\left[\vec{\nabla}\cdot(\tilde{\rho}\vec{u})
   -\dfrac{1}{r^2}\dfrac{\partial(r^2\tilde{\rho} u_r)}{\partial r}\right], \\
   & = -\dfrac{\partial}{\partial t}\left[\dfrac{\partial (\tilde{\rho} u_r)}{\partial r}
   +\dfrac{2\tilde{\rho} u_r}{r}\right], \\
   & = \dfrac{\partial}{\partial t}\left[\dfrac{\partial (\Delta_H W)}{\partial r}
   +\dfrac{2}{r}\Delta_H W\right], \\
   & = \Delta_H\dfrac{\partial}{\partial t}\left(\dfrac{\partial W}{\partial r}\right).
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
            +\Delta_H\,\dfrac{\partial W}{\partial r}-\dfrac{2}{r}\Delta_H\,W\right]\,.

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
   = {\cal N}^P_{\ell m} = -\int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^P=-\int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{\nabla}_H\cdot\vec{F}\,. & &
   \end{aligned}}
   :label: eqSpecP

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecP` are coded in
             the subroutines :f:subr:`get_wpMat <updatez_mod/get_wpmat()>`


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
             +\frac{1}{Pm}\left(\vec{\nabla}\times\vec{B}\right)\times\vec{B}\,.
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
   = {\cal N}^S_{\ell m} = \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^S = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\left[-\vec{u}\cdot\vec{\nabla}s+
   \dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}\left(\Phi_\nu+
   \dfrac{\lambda}{Pm^2\,E}\,j^2\right) \right]\,. & &
   \end{aligned}}
   :label: eqSpecS

In this expression, :math:`j=\vec{\nabla}\times\vec{B}` is the current. Once again,
the numerical evaluation of the right-hand-side (i.e. the non-linear terms) is
discussed :ref:`below <secNonLinearS>`.

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecS` are coded in
             the subroutines :f:subr:`get_sMat <updatez_mod/get_smat()>`

Equation for chemical composition :math:`\xi`
---------------------------------------------

The equation for the chemical composition is given by

.. math::
   \boxed{
   \begin{aligned}
   \dfrac{1}{Sc}\left[\left(Sc\dfrac{\partial}{\partial t} + 
   \kappa_\xi\,\dfrac{\ell(\ell+1)}{r^2} 
   \right)\right. & \,{\cal C}_n  & \\
   -\kappa_\xi\,\left(\dfrac{d\ln\kappa_\xi}{dr}+\dfrac{d\ln\tilde{\rho}}{dr}+
   +\dfrac{2}{r}\right) 
   &\,{\cal C}'_n & \\
   -\kappa_\xi & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& \xi_{\ell m n} \\
   = {\cal N}^\xi_{\ell m} = \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^\xi = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\left[-\vec{u}\cdot\vec{\nabla}\xi
   \right]\,. & &
   \end{aligned}}
   :label: eqSpecXi

Once again, the numerical evaluation of the right-hand-side (i.e. the
non-linear term) is discussed :ref:`below <secNonLinearXi>`.

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecXi` are coded in
             the subroutines :f:subr:`get_xiMat <updatexi_mod/get_ximat()>`



Equation for the poloidal magnetic potential :math:`g`
------------------------------------------------------

The equation for the poloidal magnetic potential is the radial 
component of the dynamo equation since 

.. math::
  \vec{e_r}\cdot\left(\dfrac{\partial \vec{B}}{\partial t}\right) =
   \dfrac{\partial}{\partial t}(\vec{e_r}\cdot\vec{B})=-\Delta_H\dfrac{\partial
   g}{\partial t}.

The spectral form then reads 

.. math::
   \boxed{
   \begin{aligned}
   \dfrac{\ell(\ell+1)}{r^2}\left[\left(\dfrac{\partial}{\partial t} + 
   \dfrac{1}{Pm}\lambda\,\dfrac{\ell(\ell+1)}{r^2} 
   \right)\right. & \,{\cal C}_n  & \\
   -\dfrac{1}{Pm}\,\lambda & \,{\cal C}''_n \left. \phantom{\dfrac{d\nu}{dr}}\right]& g_{\ell m n} \\
   = {\cal N}^g_{\ell m} = \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^g=\int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \vec{D}\,. & &
   \end{aligned}}
   :label: eqSpecG

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecG` are coded in
             the subroutines :f:subr:`get_bMat <updateb_mod/get_bmat()>`



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
   = {\cal N}^h_{\ell m}= \int d\Omega\,{Y_{\ell}^{m}}^\star\,{\cal N}^h = \int d\Omega\,{Y_{\ell}^{m}}^\star\,\vec{e_r}\cdot \left(\vec{\nabla}\times \vec{D}\right)\,. & &
   \end{aligned}}
   :label: eqSpecH

.. seealso:: The exact computation of the linear terms of :eq:`eqSpecH` are coded in
             the subroutines :f:subr:`get_bMat <updateb_mod/get_bmat()>`

.. note:: We note that the terms on the left hand side of :eq:`eqSpecG` and :eq:`eqSpecH`
          resulting from the magnetic diffusion term
          and the explicit time derivative completely decouple 
          in spherical harmonic degree and order.
          
          The dynamo term does not decouple:

          .. math::
             \vec{D}=\vec{\nabla}\times\left(\vec{u}\times\vec{B}\right)\,.
             :label: eqDynamoTerm


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

We consequently adopt in **MagIC** a mixed implicit/explicit algorithm.
The general differential equation in time can be written in the form

.. math:: \dfrac{\partial }{\partial t} y = \mathcal{I}(y,t) + \mathcal{E}(y,t),\quad y(t_o)=y_o\;.
  :label: eqTstep

where :math:`\mathcal{I}` denotes the terms treated in an implicit time step 
and :math:`\mathcal{E}` the terms treated explicitly, i.e. the nonlinear and Coriolis contributions.  In MagIC, two families of time-stepping schemes are supported: IMEX multistep and IMEX multistage methods.

First of all, the IMEX multistep methods correpond to time schemes where the solution results from the combination of several previous steps (such as for instance the Crank-Nicolson/Adams-Bashforth scheme). In this case, a general :math:`k`-step IMEX multistep scheme reads

.. math:: \left(I-b_o^{\mathcal I} \delta t\,\mathcal{I}\right)y_{n+1}=\sum_{j=1}^k a_j y_{n+1-j}+\delta t\sum_{j=1}^k \left(b_j^\mathcal{E} \mathcal{E}_{n+1-j}+b_{j}^\mathcal{I}\mathcal{I}_{n+1-j}\right)\,,

where :math:`I` is the identity matrix. The vectors :math:`\vec{a}`, :math:`\vec{b}^\mathcal{E}` and :math:`\vec{b}^\mathcal{I}` correspond to the weighting factors of an IMEX multistep scheme.  For instance, the commonly-used second-order scheme assembled from the combination of a Crank-Nicolson for the implicit terms and a second-order Adams-Bashforth for the explicit terms (hereafter CNAB2) corresponds to the following vectors: :math:`\vec{a}=(1,0)`, :math:`\vec{b}^{\mathcal{I}}=(1/2,1/2)` and :math:`\vec{b}^{\mathcal{E}}=(3/2,-1/2)` for a constant timestep size :math:`\delta t`.  

In addition to CNAB2, MagIC supports several semi-implicit backward differentiation schemes of second, third and fourth order that are known to have good stability properties (heareafter SBDF2, SBDF3 and SBDF4), a modified CNAB2 from `Ascher et al. (1995) <https://doi.org/10.1137/0732037>`_ (termed MODCNAB) and the CNLF scheme (combination of Crank-Nicolson and Leap-Frog for the explicit terms).

MagIC also supports several IMEX Runge-Kutta multistage methods, frequently
called DIRK, an acronym that stands for *Diagonally Implicit Runge Kutta*. For
such schemes, the equation :eq:`eqTstep` is time-advanced from :math:`t_n` to :math:`t_{n+1}` by solving :math:`\nu` sub-stages

.. math:: (I-a_{ii}^\mathcal{I}\delta t \mathcal{I})y_{i}=y_n+\delta t \sum_{j=1}^{i-1}\left(a_{i,j}^{\mathcal{E}}\mathcal{E}_j+a_{i,j}^\mathcal{I}\mathcal{I}_j\right), \quad 1\leq i\leq \nu,

where :math:`y_i` is the intermediate solution at the stage :math:`i`. The matrices :math:`a_{i,j}^\mathcal{E}` and :math:`a_{i,j}^\mathcal{I}` constitute the so-called Butcher tables that correspond to a DIRK scheme. MagIC supports several second and third order schemes: ARS222 and ARS443 from `Ascher et al. (1997) <https://doi.org/10.1016/S0168-9274(97)00056-1>`_, LZ232 from `Liu and Zou (2006) <https://doi.org/10.1016/j.cam.2005.02.020>`_, PC2 from `Jameson et al. (1981) <https://doi.org/10.2514/6.1981-1259>`_ and  BPR353 from `Boscarino et al. (2013) <https://doi.org/10.1137/110842855>`_.


In the code the equation :eq:`eqTstep` is formulated for each unknown spectral coefficient  
(expect pressure) of spherical harmonic degree :math:`\ell` and order :math:`m` 
and for each radial grid point :math:`r_k`. 
Because non-linear terms and the Coriolis force are treated explicitly, 
the equations decouple for all spherical harmonic modes.
The different radial grid points, however, couple via the 
Chebychev modes and form a linear algebraic system of equations that can 
be solved with standard methods for the different spectral contributions. 

For example the respective system of equations for the poloidal magnetic potential :math:`g` time advanced with a CNAB2 reads

.. math::
      \left( \mathcal{A}_{kn} - \dfrac{1}{2}\,\delta t\,\mathcal{I}_{kn}\right)\;g_{\ell mn}(t+\delta t) =
      \left( \mathcal{A}_{kn} + \dfrac{1}{2}\,\delta t\,\mathcal{I}_{kn} \right)\;g_{\ell mn}(t) +
      \frac{3}{2}\,\delta t\,\mathcal{E}_{k\ell m}(t) - \frac{1}{2}\,\delta t\,\mathcal{E}_{k\ell m}(t-\delta t) 
      :label: imex

with 

.. math::
    \mathcal{A}_{kn} = \dfrac{\ell (\ell+1)}{r_k^2}\,{\cal C}_{nk}\;,

.. math::
    \mathcal{I}_{kn}=\dfrac{\ell(\ell+1)}{r_k^2}\,\dfrac{1}{Pm}\left({\cal C}''_{nk} - \dfrac{\ell(\ell+1)}{r_k^2}\; 
    {\cal C}_{nk}\right)\;,

and :math:`{\cal C}_{nk}={\cal C}_n(r_k)`.
:math:`\mathcal{A}_{kn}` is a matrix that converts the poloidal field modes :math:`g_{\ell mn}` 
to the radial magnetic field :math:`B_r(r_k,\ell,m)` for a given spherical harmonic contribution.

Here :math:`k` and :math:`n` number the radial grid points and the Chebychev coefficients, respectively. 
Note that the Einstein sum convention is used for Chebychev modes :math:`n`.

:math:`\mathcal{I}_{kn}` is the matrix describing the implicit contribution which is purely diffusive here. 
Neither  :math:`\mathcal{A}_{kn}` nor :math:`\mathcal{I}_{kn}` depend on time but the former 
needs to be updated when the time step :math:`\delta t` is changed. 
The only explicit contribution is the nonlinear dynamo term 

.. math:: \mathcal{E}_{k\ell m}(t)= {\cal N}_{k\ell m}^g = \int d\Omega\; {Y_{\ell}^{m}}^\star\; 
          \vec{e_r} \cdot \vec{D}(t,r_k,\theta,\phi)\;\; .  

:math:`\mathcal{E}_{k\ell m}` is a one dimensional vector for all spherical harmonic combinations 
:math:`\ell m`.
  

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

Coriolis force and nonlinear terms
==================================

.. _secNonLinearW:

Nonlinear terms entering the equation for :math:`W`
---------------------------------------------------

The nonlinear term :math:`{\cal N}^W` that enters the equation for the poloidal potential
:eq:`eqSpecW` contains the radial component of the advection, the Lorentz force 
and Coriolis force. In spherical coordinate, the two first contributions read:

.. math::
   \tilde{\rho}\left(\vec{u}\cdot\vec{\nabla}\vec{u}\right)=
   \left\lbrace
   \begin{aligned}
   {\cal A}_r \\
   {\cal A}_\theta \\
   {\cal A}_\phi
   \end{aligned}
   \right\rbrace
   =
   \left\lbrace
   \begin{aligned}
   -\tilde{\rho}\,E\,\left(
   u_r\dfrac{\partial u_r}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_r}{\partial \theta}+
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_r}{\partial \phi}
   -\dfrac{u_\theta^2+u_\phi^2}{r}\right)+
   \dfrac{1}{Pm}\left(j_\theta\,B_\phi-j_\phi\,B_\theta\right)\, , \\
   -\tilde{\rho}\,E\,\left(
   u_r\dfrac{\partial u_\theta}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_\theta}{\partial \theta} +
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_\theta}{\partial \phi}+
   \dfrac{u_r u_\theta}{r}-\dfrac{\cos\theta}{r\sin\theta}u_\phi^2\right)+
   \dfrac{1}{Pm}\left(j_\phi\,B_r-j_r\,B_\phi\right)\, ,\\
   -\tilde{\rho}\,E\,\left(
   u_r\dfrac{\partial u_\phi}{\partial r}+
   \dfrac{u_\theta}{r}\dfrac{\partial u_\phi}{\partial \theta} +
   \dfrac{u_\phi}{r\sin\theta}\dfrac{\partial u_\phi}{\partial \phi}+
   \dfrac{u_r u_\phi}{r} +\dfrac{\cos\theta}{r\sin\theta}u_\theta u_\phi\right)+
   \dfrac{1}{Pm}\left(j_r\,B_\theta-j_\theta\,B_r\right)\, ,
   \end{aligned}
   \right\rbrace
   :label: eqAdvection

The Coriolis force can be expressed as a function of the potentials :math:`W` and
:math:`Z` using :eq:`eqToroPolo1`

.. math::
   2\tilde{\rho} \vec{e_r}\cdot(\vec{u}\times\vec{e_z})=2\sin\theta\,\tilde{\rho}
   u_\phi=\dfrac{2}{r}\left(\dfrac{\partial^2 W}{\partial r\partial \phi}-\sin\theta
   \dfrac{\partial Z}{\partial \theta}\right)\,.

The nonlinear terms that enter the equation for the poloidal potential :eq:`eqSpecW` thus 
reads:

.. math::
   {\cal N}^W = \dfrac{2}{r}\left(\dfrac{\partial^2 W}{\partial r\partial \phi}-\sin\theta
   \dfrac{\partial Z}{\partial \theta}\right)+{\cal A}_r\,.

The :math:`\theta`-derivative entering the radial component of the Coriolis force is thus the
operator :math:`\vartheta_2` defined in :eq:`eqOpTheta1`. Using the recurrence
relation, one thus finally gets in spherical harmonic space:

.. math::
   \boxed{
   {\cal N}^W_{\ell m}  = \dfrac{2}{r}\left[i m \dfrac{\partial W_\ell^m}{\partial r}-(\ell-1)c_\ell^m
   Z_{\ell-1}^m+(\ell+2)c_{\ell+1}^m Z_{\ell+1}^m\right]
   +{{\cal A}_r}_\ell^m\, .
   }
   :label: eqNLW

To get this expression, we need to first compute :math:`{\cal A}_r` in the physical space. This
term is computed in the subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in
the module :f:mod:`grid_space_arrays_mod`. :math:`{\cal A}_r` is then transformed to the
spectral space by using a Legendre and a Fourier transform to produce :math:`{{\cal A}_r}_\ell^m`.

.. seealso:: The final calculations of :eq:`eqNLW` are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`.

.. _secNonLinearZ:

Nonlinear terms entering the equation for :math:`Z`
---------------------------------------------------

The nonlinear term :math:`{\cal N}^Z` that enters the equation for the toroidal potential
:eq:`eqSpecZ` contains the radial component of the curl of the advection and Coriolis force.
The Coriolis force can be rewritten as a function of :math:`W` and :math:`Z`:

.. math::
    \begin{aligned}
    \vec{e_r}\cdot\vec{\nabla}\times\left[(2\tilde{\rho}\vec{u})\times
    \vec{e_z}\right] & =2\vec{e_r}\cdot\left[(\vec{e_z}\cdot\vec{\nabla})(\tilde{\rho}
    \vec{u})\right], \\
    & = 2\left[\cos\theta\dfrac{\partial (\tilde{\rho} u_r)}{\partial r}
    -\dfrac{\sin\theta}{r}\dfrac{\partial (\tilde{\rho}
    u_r)}{\partial \theta}+\dfrac{\tilde{\rho} u_\theta\sin\theta}{r}\right], \\
    & = 2\left[-\cos\theta\dfrac{\partial}{\partial r}(\Delta_H W)+
    \dfrac{\sin\theta}{r}\dfrac{\partial}{\partial \theta}(\Delta_H
    W)+\dfrac{\sin\theta}{r^2}\dfrac{\partial^2 W}{\partial r\partial \theta}+
    \dfrac{1}{r^2}\dfrac{\partial Z}{\partial \phi}\right].
    \end{aligned}

Using the :math:`\vartheta` operators defined in :eq:`eqOpTheta1`-:eq:`eqOpTheta4` then
allows to rewrite the Coriolis force in the following way:

.. math::
   \vec{e_r}\cdot\vec{\nabla}\times\left[(2\tilde{\rho}\vec{u})\times
   \vec{e_z}\right]=\dfrac{2}{r^2}\left(\vartheta_3\,\dfrac{\partial W}{\partial r}
   -\dfrac{1}{r}\,\vartheta_4\,W+ \dfrac{\partial Z}{\partial \phi} \right)\,.
   :label: eqCorZNL

The contributions of nonlinear advection and Lorentz forces that enter the equation
for the toroidal potential are written this way:

.. math::
   \dfrac{1}{r\sin\theta}\left[
   \dfrac{\partial (\sin\theta{\cal A}_\phi)}{\partial \theta} -
   \dfrac{\partial {\cal A}_\theta}{
   \partial\phi}\right]\,.

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, the actual
strategy is to follow the following steps:

1. Compute the quantities :math:`{\cal A}_\phi/r\sin\theta`
   and :math:`{\cal A}_\theta/r\sin\theta` in the physical space. In the code, this step
   is computed in the subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in 
   the module :f:mod:`grid_space_arrays_mod`. 

2. Transform :math:`{\cal A}_\phi/r\sin\theta` and :math:`{\cal A}_\theta/r\sin\theta` to
   the spectral space (thanks to a Legendre and a Fourier transform). In MagIC, this step
   is computed in the modules :f:mod:`legendre_grid_to_spec` and :f:mod:`fft`. After
   this step :math:`{{\cal A}t}_{\ell}^m` and :math:`{{\cal A}p}_{\ell}^m` are defined.

3. Calculate the colatitude and theta derivatives using the recurrence relations:

   .. math::
      \vartheta_1\,{{\cal A}p}_{\ell}^m-\dfrac{\partial {{\cal A}t}_{\ell}^m}{\partial \phi}\,.
      :label: eqAdvZNL

Using :eq:`eqCorZNL` and :eq:`eqAdvZNL`, one thus finally gets

.. math::
   \boxed{
   \begin{aligned}
   {\cal N}^Z_{\ell m}  = & \dfrac{2}{r^2}\left[(\ell-1)(\ell+1)\,c_\ell^m\,
   \dfrac{\partial W_{\ell-1}^m}{\partial r}+\ell(\ell+2)\,c_{\ell+1}^m\,
   \dfrac{\partial W_{\ell+1}^m}{\partial r} \right. \\
   & \left. -\dfrac{\ell(\ell-1)(\ell+1)}{r}\,c_\ell^m\,W_{\ell-1}^m+
   \dfrac{\ell(\ell+1)(\ell+2)}{r}\,c_{\ell+1}^m\,W_{\ell+1}^m+
   im\,Z_\ell^m\right] \\
   & + (\ell+1)\,c_\ell^m\,{{\cal A}p}_{\ell-1}^m-
   \ell\,c_{\ell+1}^m\,{{\cal A}p}_{\ell+1}^m
   -im\,{{\cal A}t}_{\ell}^m\,.
   \end{aligned}
   }
   :label: eqNLZ

.. seealso:: The final calculations of :eq:`eqNLZ` are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`.

.. _secNonLinearP:

Nonlinear terms entering the equation for :math:`P`
---------------------------------------------------

The nonlinear term :math:`{\cal N}^P` that enters the equation for the pressure
:eq:`eqSpecP` contains the horizontal divergence of the advection and Coriolis force.
The Coriolis force can be rewritten as a function of :math:`W` and :math:`Z`:

.. math::
    \begin{aligned}
    \vec{\nabla}_H\cdot\left[(2\tilde{\rho}\vec{u})\times
    \vec{e_z}\right] & =2\vec{e_z}\cdot\left[\vec{\nabla}\times(\tilde{\rho}
    \vec{u})\right] -\left(\dfrac{\partial}{\partial r}+\dfrac{2}{r}\right)
    \left[\vec{e_r}\cdot(2\tilde{\rho}\vec{u}\times\vec{e_z})\right],\\
    & = -2\cos\theta\,\Delta_H Z-2\sin\theta\left[-\dfrac{1}{r\sin\theta}
    \dfrac{\partial}{\partial\phi}\left(
    \dfrac{\partial^2}{\partial r^2}+\Delta_H \right) W +
    \dfrac{1}{r}\dfrac{\partial^2 Z}{\partial r\partial\theta}\right]
    \\
    & \phantom{=\cos\theta} -\left(\dfrac{\partial}{\partial r}+\dfrac{2}{r}\right)
    \left[2\sin\theta\tilde{\rho}u_\phi\right], \\
    & = 2\left[\dfrac{1}{r}\left(\Delta_H+\dfrac{\partial^2}{\partial r^2}\right)
    \dfrac{\partial W}{\partial \phi}-\cos\theta\Delta_H Z -\dfrac{\sin\theta}{r}
    \dfrac{\partial^2 Z}{\partial r \partial \theta}\right] \\
    & \phantom{=\cos\theta} -\left(\dfrac{\partial}{\partial r}+\dfrac{2}{r}\right)
    \left[\dfrac{2}{r}\left(\dfrac{\partial^2 W}{\partial r\partial\phi}-\sin\theta
    \dfrac{\partial Z}{\partial \theta}\right)\right], \\
    & = 2\left(\dfrac{\Delta_H}{r}\dfrac{\partial W}{\partial \phi}-\dfrac{1}{r^2}
    \dfrac{\partial^2 W}{\partial\phi\partial r} -\cos\theta\Delta_H\,Z
    +\dfrac{\sin\theta}{r^2}\dfrac{\partial Z}{\partial \theta}\right).
    \end{aligned}

Using the :math:`\vartheta` operators defined in :eq:`eqOpTheta3`-:eq:`eqOpTheta4` then
allows to rewrite the Coriolis force in the following way:

.. math::
   \vec{\nabla}_H\cdot\left[(2\tilde{\rho}\vec{u})\times
   \vec{e_z}\right]=\dfrac{2}{r^2}\left(-\dfrac{L_H}{r}\,\dfrac{\partial W}{\partial \phi}
   -\dfrac{\partial^2 W}{\partial\phi\partial r}+\vartheta_3\, Z
   \right)\,.
   :label: eqCorPNL

The contributions of nonlinear advection and Lorentz forces that enter the equation
for pressure are written this way:

.. math::
   \dfrac{1}{r\sin\theta}\left[
   \dfrac{\partial (\sin\theta{\cal A}_\theta)}{\partial \theta} +
   \dfrac{\partial {\cal A}_\phi}{
   \partial\phi}\right]\,.

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, we then follow
the same three steps as for the advection term entering the equation for :math:`Z`.

.. math::
   \vartheta_1\,{{\cal A}t}_{\ell}^m+\dfrac{\partial {{\cal A}p}_{\ell}^m}{\partial \phi}\,.
   :label: eqAdvPNL

Using :eq:`eqCorPNL` and :eq:`eqAdvPNL`, one thus finally gets

.. math::
   \boxed{
   \begin{aligned}
   {\cal N}^P_{\ell m}  = & \dfrac{2}{r^2}\left[-im\,\dfrac{\ell(\ell+1)}{r}\,W_\ell^m
   -im\,\dfrac{\partial W_\ell^m}{\partial r}+(\ell-1)(\ell+1)\,c_\ell^m\,
   Z_{\ell-1}^m+\ell(\ell+2)\,c_{\ell+1}^m\,
   Z_{\ell+1}^m \right] \\
   & + (\ell+1)\,c_\ell^m\,{{\cal A}t}_{\ell-1}^m-
   \ell\,c_{\ell+1}^m\,{{\cal A}t}_{\ell+1}^m
   +im\,{{\cal A}p}_{\ell}^m\,.
   \end{aligned}
   }
   :label: eqNLP

.. seealso:: The final calculations of :eq:`eqNLP` are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`.

.. _secNonLinearS:

Nonlinear terms entering the equation for :math:`s`
---------------------------------------------------

The nonlinear terms that enter the equation for entropy/temperature
:eq:`eqSpecS` are twofold: (i) the advection term, (ii) the viscous and Ohmic
heating terms (that vanish in the Boussinesq limit of the Navier Stokes equations).

Viscous and Ohmic heating are directly calculated in the physical space by the
subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in
the module :f:mod:`grid_space_arrays_mod`. Let's introduce :math:`{\cal H}`, the sum
of the viscous and Ohmic heating terms.

.. math::
   {\cal H} = \dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}\left(\Phi_\nu+
   \dfrac{\lambda}{Pm^2\,E}\,j^2\right)\,.

Expanding this term leads to:

.. math::
   \begin{aligned}
   {\cal H}=& \dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}\left[
   \tilde{\rho}\nu\left\lbrace 2\left(\dfrac{\partial u_r}{ \partial r}\right)^2
   +2\left(\dfrac{1}{r}\dfrac{\partial u_\theta}{\partial\theta}+\dfrac{u_r}{r}
   \right)^2+2\left( \dfrac{1}{r\sin\theta}\dfrac{\partial u_\phi}{\partial\phi}
   + \dfrac{u_r}{r}+\dfrac{\cos\theta}{r\sin\theta}u_\theta \right)^2\right.\right. \\
   & \phantom{\dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}}
   +\left(r\dfrac{\partial}{\partial r}\left(\dfrac{u_\theta}{r}
   \right)+\dfrac{1}{r}\dfrac{\partial u_r}{\partial\theta}\right)^2+
   \left(r\dfrac{\partial}{\partial r}\left(\dfrac{u_\phi}{r}\right)+
   \dfrac{1}{r\sin\theta}\dfrac{\partial u_r}{\partial\phi}  \right)^2 \\
   & \phantom{\dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}}\left.
   + \left(\dfrac{\sin\theta}{r}\dfrac{\partial}{\partial\theta}\left(
   \dfrac{u_\phi}{\sin\theta}\right)+\dfrac{1}{r\sin\theta}
   \dfrac{\partial u_\theta}{\partial\phi}\right)^2 
   -\dfrac{2}{3}\,\left(\dfrac{d\ln\tilde{\rho}}{dr}\,u_r\right)^2 \right\rbrace \\
   & \phantom{\dfrac{Pr\,Di}{Ra}\dfrac{1}{\tilde{\rho}\tilde{T}}}\left.
   +  \dfrac{\lambda}{Pm^2\,E}\,\left\lbrace 
   j_r^2+j_\theta^2+j_\phi^2\right\rbrace\right]\,.
   \end{aligned}
   :label: eqHeatingEntropy

This term is then transformed to the spectral space with a Legendre and a Fourier
transform to produce :math:`{\cal H}_\ell^m`.

The treatment of the advection term :math:`-\vec{u}\cdot\vec{\nabla}s` is a bit different.
It is in a first step rearranged as follows

.. math::
   -\vec{u}\cdot\vec{\nabla}s = -\dfrac{1}{\tilde{\rho}}\left[
   \vec{\nabla}\cdot\left(\tilde{\rho}s\vec{u} \right)-
   s\underbrace{\vec{\nabla}\cdot\left(\tilde{\rho}\vec{u} \right)}_{=0}\right]\,.

The quantities that are calculated in the physical space are thus simply the product of
entropy/temperature :math:`s` by the velocity components. This defines three variables
defined in the grid space that are computed in the subroutine :f:subr:`get_nl 
<grid_space_arrays_mod/get_nl()>`:

.. math::
   \mathcal{US}_r = \tilde{\rho}s u_r,\quad  \mathcal{US}_\theta = \tilde{\rho}s u_\theta,
   \quad \mathcal{US}_\phi = \tilde{\rho}s u_\phi,

To get the actual advection term, one must then apply the divergence operator to get:

.. math::
   -\vec{u}\cdot\vec{\nabla}s = -\dfrac{1}{\tilde{\rho}}\left[
   \dfrac{1}{r^2}\dfrac{\partial}{\partial r}\left(r^2\,\mathcal{US}_r\right)+
   \dfrac{1}{r\sin\theta}\dfrac{\partial}{\partial\theta}\left(\sin\theta\,\mathcal{US}_\theta
   \right)+\dfrac{1}{r\sin\theta}\dfrac{\partial\,\mathcal{US}_\phi}{\partial\phi}\right]\,.

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, the actual
strategy is then to follow the following steps:

1. Compute the quantities :math:`r^2\,\mathcal{US}_r`, :math:`\mathcal{US}_\phi/r\sin\theta`
   and :math:`\mathcal{US}_\theta/r\sin\theta` in the physical space. In the code, this step
   is computed in the subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in 
   the module :f:mod:`grid_space_arrays_mod`. 

2. Transform :math:`r^2\,\mathcal{US}_r`, :math:`\mathcal{US}_\phi/r\sin\theta` 
   and :math:`\mathcal{US}_\theta/r\sin\theta` to
   the spectral space (thanks to a Legendre and a Fourier transform). In MagIC, this step
   is computed in the modules :f:mod:`legendre_grid_to_spec` and :f:mod:`fft`. After
   this step :math:`{\mathcal{US}r}_{\ell}^m`, :math:`{\mathcal{US}t}_{\ell}^m` 
   and :math:`{\mathcal{US}p}_{\ell}^m` are defined.

3. Calculate the colatitude and theta derivatives using the recurrence relations:

   .. math::
      -\dfrac{1}{\tilde{\rho}}\left[
      \dfrac{1}{r^2}\dfrac{\partial\, {\mathcal{US}r}_\ell^m}{\partial r}+
      \vartheta_1\,{\mathcal{US}t}_\ell^m+
      \dfrac{\partial\,{\mathcal{US}p}_\ell^m}{\partial \phi}\right]\,.
      :label: eqAdvSNL

Using :eq:`eqHeatingEntropy` and :eq:`eqAdvSNL`, one thus finally gets

.. math::
   \boxed{
   {\cal N}^S_{\ell m}  = -\dfrac{1}{\tilde{\rho}}\left[
   \dfrac{1}{r^2}\dfrac{\partial\, {\mathcal{US}r}_\ell^m}{\partial r}
   + (\ell+1)\,c_\ell^m\,{\mathcal{US}t}_{\ell-1}^m-
   \ell\,c_{\ell+1}^m\,{\mathcal{US}t}_{\ell+1}^m
   +im\,{\mathcal{US}p}_\ell^m\right]+{\cal H}_\ell^m\,.
   }
   :label: eqNLS

.. seealso:: The :math:`\theta` and :math:`\phi` derivatives that enter :eq:`eqNLS` 
             are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`. The radial derivative
             is computed afterwards at the very beginning of
             :f:subr:`updateS <updates_mod/updates()>`.

.. _secNonLinearXi:

Nonlinear terms entering the equation for :math:`\xi`
-----------------------------------------------------

The nonlinear term that enters the equation for chemical composition
:eq:`eqSpecXi` is the advection term. This term is treated the same way
as the advection term that enters the entropy equation.
It is in a first step rearranged as follows

.. math::
   -\vec{u}\cdot\vec{\nabla}\xi = -\dfrac{1}{\tilde{\rho}}\left[
   \vec{\nabla}\cdot\left(\tilde{\rho}\xi\vec{u} \right)-
   \xi\underbrace{\vec{\nabla}\cdot\left(\tilde{\rho}\vec{u} \right)}_{=0}\right]\,.

The quantities that are calculated in the physical space are thus simply the
product of composition :math:`\xi` by the velocity components. This
defines three variables defined in the grid space that are computed in the
subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>`:

.. math::
   \mathcal{UX}_r = \tilde{\rho}\xi u_r,\quad  \mathcal{US}_\theta = \tilde{\rho}\xi u_\theta,
   \quad \mathcal{UX}_\phi = \tilde{\rho}\xi u_\phi,

To get the actual advection term, one must then apply the divergence operator
to get:

.. math::
   -\vec{u}\cdot\vec{\nabla}\xi = -\dfrac{1}{\tilde{\rho}}\left[
   \dfrac{1}{r^2}\dfrac{\partial}{\partial r}\left(r^2\,\mathcal{UX}_r\right)+
   \dfrac{1}{r\sin\theta}\dfrac{\partial}{\partial\theta}\left(\sin\theta\,\mathcal{UX}_\theta
   \right)+\dfrac{1}{r\sin\theta}\dfrac{\partial\,\mathcal{UX}_\phi}{\partial\phi}\right]\,.

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, the actual
strategy is then to follow the following steps:

1. Compute the quantities :math:`r^2\,\mathcal{UX}_r`, :math:`\mathcal{UX}_\phi/r\sin\theta`
   and :math:`\mathcal{UX}_\theta/r\sin\theta` in the physical space. In the code, this step
   is computed in the subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in 
   the module :f:mod:`grid_space_arrays_mod`. 

2. Transform :math:`r^2\,\mathcal{UX}_r`, :math:`\mathcal{UX}_\phi/r\sin\theta` 
   and :math:`\mathcal{UX}_\theta/r\sin\theta` to
   the spectral space (thanks to a Legendre and a Fourier transform). In MagIC, this step
   is computed in the modules :f:mod:`legendre_grid_to_spec` and :f:mod:`fft`. After
   this step :math:`{\mathcal{UX}r}_{\ell}^m`, :math:`{\mathcal{UX}t}_{\ell}^m` 
   and :math:`{\mathcal{UX}p}_{\ell}^m` are defined.

3. Calculate the colatitude and theta derivatives using the recurrence relations:

   .. math::
      -\dfrac{1}{\tilde{\rho}}\left[
      \dfrac{1}{r^2}\dfrac{\partial\, {\mathcal{UX}r}_\ell^m}{\partial r}+
      \vartheta_1\,{\mathcal{UX}t}_\ell^m+
      \dfrac{\partial\,{\mathcal{UX}p}_\ell^m}{\partial \phi}\right]\,.

One thus finally gets

.. math::
   \boxed{
   {\cal N}^\xi_{\ell m}  = -\dfrac{1}{\tilde{\rho}}\left[
   \dfrac{1}{r^2}\dfrac{\partial\, {\mathcal{UX}r}_\ell^m}{\partial r}
   + (\ell+1)\,c_\ell^m\,{\mathcal{UX}t}_{\ell-1}^m-
   \ell\,c_{\ell+1}^m\,{\mathcal{UX}t}_{\ell+1}^m
   +im\,{\mathcal{UX}p}_\ell^m\right]\,.
   }
   :label: eqNLXi

.. seealso:: The :math:`\theta` and :math:`\phi` derivatives that enter :eq:`eqNLXi` 
             are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`. The radial derivative
             is computed afterwards at the very beginning of
             :f:subr:`updateXi <updatexi_mod/updatexi()>`.

.. _secNonLinearG:

Nonlinear terms entering the equation for :math:`g`
---------------------------------------------------

The nonlinear term that enters the equation for the poloidal potential of the magnetic
field :eq:`eqSpecG` is the radial component of the induction term :eq:`eqDynamoTerm`.
In the following we introduce the electromotive force 
:math:`{\cal F} = \vec{u}\times\vec{B}` with its three components 

.. math::
   {\cal F}_r=u_\theta B_\phi-u_\phi B_\theta,\quad
   {\cal F}_\theta=u_\phi B_r-u_r B_\phi,\quad
   {\cal F}_\phi=u_r B_\theta-u_\theta B_r\,.

The radial component of the induction term then reads:

.. math::
  {\cal N}^g = \vec{e_r}\cdot\left[\vec{\nabla}\times\left(\vec{u}\times\vec{B}\right)\right]
   =\dfrac{1}{r\sin\theta}\left[\dfrac{\partial\,(\sin\theta {\cal F}_\phi)}{\partial\theta}
   -\dfrac{\partial {\cal F}_\theta}{\partial \phi}\right]\,.

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, we then
follow the usual following steps:

1. Compute the quantities :math:`r^2\,\mathcal{F}_r`, :math:`\mathcal{F}_\phi/r\sin\theta`
   and :math:`\mathcal{F}_\theta/r\sin\theta` in the physical space. In the code, this step
   is computed in the subroutine :f:subr:`get_nl <grid_space_arrays_mod/get_nl()>` in 
   the module :f:mod:`grid_space_arrays_mod`. 

2. Transform :math:`r^2\,\mathcal{F}_r`, :math:`\mathcal{F}_\phi/r\sin\theta` 
   and :math:`\mathcal{F}_\theta/r\sin\theta` to
   the spectral space (thanks to a Legendre and a Fourier transform). In MagIC, this step
   is computed in the modules :f:mod:`legendre_grid_to_spec` and :f:mod:`fft`. After
   this step :math:`{\mathcal{F}_r}_{\ell}^m`, :math:`{\mathcal{F}_\theta}_{\ell}^m` 
   and :math:`{\mathcal{F}_\phi}_{\ell}^m` are defined.

3. Calculate the colatitude and theta derivatives using the recurrence relations:

   .. math::
      \vartheta_1\,{\mathcal{F}_\phi}_\ell^m-
      \dfrac{\partial\,{\mathcal{F}_\theta}_\ell^m}{\partial \phi}\,.

We thus finally get

.. math::
   \boxed{
   {\cal N}^g_{\ell m}  = 
   (\ell+1)\,c_\ell^m\,{\mathcal{F}_\phi}_{\ell-1}^m-\ell\,c_{\ell+1}^m\,
   {\mathcal{F}_\phi}_{\ell+1}^m -im\,{\mathcal{F}_\theta}_{\ell}^m\,.
   }
   :label: eqNLG

.. seealso:: The final calculations of :eq:`eqNLG` are done in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`.

.. _secNonLinearH:

Nonlinear terms entering the equation for :math:`h`
---------------------------------------------------

The nonlinear term that enters the equation for the toroidal potential of the magnetic
field :eq:`eqSpecH` is the radial component of the curl of the 
induction term :eq:`eqDynamoTerm`:

.. math::
   \begin{aligned}
   {\cal N}^h = \vec{e_r}\cdot\left[\vec{\nabla}\times\vec{\nabla}\times\left(\vec{u}\times\vec{B}\right)
   \right]
   & =\vec{e_r}\cdot\left[\vec{\nabla}\left(\vec{\nabla}\cdot\vec{\mathcal{F}}\right)
   -\Delta\vec{\mathcal{F}}\right], \\
   & = \dfrac{\partial}{\partial r}\left[\dfrac{1}{r^2}
   \dfrac{\partial(r^2 {\mathcal{F}}_r)}{\partial r} + \dfrac{1}{r\sin\theta}
   \dfrac{\partial(\sin\theta\,{\mathcal{F}}_\theta)}{\partial\theta}+\dfrac{1}{r\sin\theta}
   \dfrac{\partial{\mathcal{F}}_\phi}{\partial\phi} \right] \\
   & \phantom{=\ }-
   \Delta {\mathcal{F}}_r+\dfrac{2}{r^2}\left[{\mathcal{F}}_r +\dfrac{1}{\sin\theta}
   \dfrac{\partial(\sin\theta\,{\mathcal{F}}_\theta)}{\partial\theta}+
   \dfrac{1}{\sin\theta}\dfrac{\partial {\mathcal{F}}_\phi}{\partial \phi}\right], \\
   & = \dfrac{1}{r^2}\dfrac{\partial}{\partial r}\left[\dfrac{r}{\sin\theta}\left(
   \dfrac{\partial(\sin\theta\,{\mathcal{F}}_\theta)}{\partial\theta}+
   \dfrac{\partial{\mathcal{F}}_\phi}{\partial\phi} \right)\right]-\Delta_H\,{\mathcal{F}}_r\,.
   \end{aligned}

To make use of the recurrence relations :eq:`eqOpTheta1`-:eq:`eqOpTheta4`, we then follow
the same steps than for the nonlinear terms that enter the equation for poloidal potential
of the magnetic field :math:`g`:

.. math::
   \dfrac{1}{r^2}\dfrac{\partial }{\partial r}\left[r^2\left(\vartheta_1\,
   {\mathcal{F}t}_\ell^m+\dfrac{\partial\,{\mathcal{F}p}_\ell^m}{\partial \phi}\right)\right]
   +L_H\, {\mathcal{F}r}_\ell^m\,.

We thus finally get

.. math::
   \boxed{
   {\cal N}^h_{\ell m}  =\ell(\ell+1)\,{\mathcal{F}r}_{\ell}^m+
   \dfrac{1}{r^2}\dfrac{\partial}{\partial r}\left[r^2\left\lbrace
   (\ell+1)\,c_\ell^m\,{\mathcal{F}t}_{\ell-1}^m-\ell\,c_{\ell+1}^m\,
   {\mathcal{F}t}_{\ell+1}^m +im\,{\mathcal{F}p}_{\ell}^m\right\rbrace
   \right]\,.
   }
   :label: eqNLH

.. seealso:: The :math:`\theta` and :math:`\phi` derivatives that enter :eq:`eqNLH` 
             are computed in the subroutine 
             :f:subr:`get_td <nonlinear_lm_mod/get_td()>`. The remaining radial derivative
             is computed afterwards at the very beginning of
             :f:subr:`updateB <updateb_mod/updateb()>`.


.. _secBoundaryConditions:

Boundary conditions and inner core
==================================

Mechanical boundary conditions
------------------------------

Since the system of equations is formulated on a radial grid, boundary
conditions can simply be satisfied by replacing the collocation equation
at grid points :math:`r_i` and :math:`r_o` with appropriate expressions.
The condition of zero radial flow on the boundaries implies that the poloidal 
potential has to vanish, i.e. :math:`W(r_o)=0` and :math:`W(r_i)=0`. 
In Chebychev representation this implies 

.. math::
  {\cal C}_n(r) W_{\ell mn} = 0 \;\;\mbox{at}\;\; r=r_i,r_o\;\;.
  :label: eqBcRigid1

Note that the summation convention with respect to
radial modes :math:`n` is used again.
**The no-slip** condition further requires that the
horizontal flow components also have to vanish, provided
the two boundaries are at rest. This condition is fulfilled for

.. math:: 
   \dfrac{\partial W}{\partial r}=0\;\;\mbox{and}\;\; Z=0,
   
at the respective boundary. In spectral space these conditions read 

.. math::
   {\cal C}'_n(r) W_{\ell mn} = 0\;\;\mbox{at}\;\; r=r_i,r_o\,,
  :label: eqBcRigid2

and

.. math::
   {\cal C}_n(r) Z_{\ell mn} = 0\;\;\mbox{at}\;\; r=r_i,r_o\,,
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
   \mathsf{I} \dfrac{\partial\omega}{\partial t}= \Gamma_L+\Gamma_\nu\,.

The tensor :math:`\mathsf{I}` denotes the moment of inertia of inner core or mantle,
respectively, :math:`\omega` is the mantle or inner-core rotation rate relative
to that of the reference frame, and :math:`\Gamma_{L,\nu}` are the respective torques
associated with Lorentz or viscous forces. The torques are expressed by

.. math::
   \Gamma_L = \dfrac{1}{E\,Pm}\oint B_r B_\phi\,r\sin\theta\,\mathrm{d}S\,,

and

.. math::
   \Gamma_\nu = \oint \tilde{\rho} \tilde{\nu} r\dfrac{\partial}{\partial r}\left(\dfrac{u_\phi}{r}\right) r\sin\theta\,\mathrm{d}S\,,

where :math:`\mathrm{d}S = r^2\sin\theta \mathrm{d}\theta\mathrm{d}\phi` and :math:`r\in[r_i,r_o]` in the above expressions. Using the following equality

.. math::
   \oint \tilde{\rho} r\sin\theta u_\phi\,\mathrm{d} S=4\sqrt{\dfrac{\pi}{3}}Z_{10}r^2,
   
the viscous torques can be expressed by

.. math::
   \Gamma_\nu = \pm4\sqrt{\dfrac{\pi}{3}}\tilde{\nu}r^2\left[\dfrac{\partial Z_{10}}{\partial r}-\left(\dfrac{2}{r}+\beta\right)Z_{10}\right]\,,

where the sign in front depends whether :math:`r=r_o` or :math:`r=r_i`.


**Free-slip boundary conditions** require that the viscous stress vanishes, which
in turn implies that the non-diagonal components :math:`\mathsf{Sr}_{r\phi}` and
:math:`\mathsf{S}_{r\theta}` of the rate-of-strain tensor vanish. 
Translated to the spectral representation this requires

.. math::
  \left[{\cal C}''_n(r) -\left(\frac{2}{r}+\dfrac{d\ln\tilde{\rho}}{dr}\right)\,{\cal C}'_n(r)
  \right] W_{\ell mn} = 0 \;\;\mbox{and}\;\;
  \left[{\cal C}'_n(r) -\left(\frac{2}{r}+\dfrac{d\ln\tilde{\rho}}{dr}\right)\,{\cal C}_n(r)
  \right] Z_{\ell mn} = 0\;.
  
We show the derivation for the somewhat simpler Boussinesq approximation which yields the condition 

.. math::
   \dfrac{\partial}{\partial r} \dfrac{\vec{u}_H}{r} = 0
   
where the index H denotes the horizonal flow components. 
In terms of poloidal and toroidal components this implies 

.. math::
   \dfrac{\partial}{\partial r} \dfrac{1}{r} \left( \vec{\nabla}_H \dfrac{\partial W}{\partial r}\right) =
   \vec{\nabla}_H \dfrac{1}{r} \left( \dfrac{\partial^2}{\partial r^2} - \dfrac{2}{r} \dfrac{\partial}{\partial r} \right) W = 0
   
and

.. math::
   \dfrac{\partial}{\partial r} \dfrac{1}{r} \nabla\times \vec{e}_r Z = 
   \nabla\times \vec{e}_r \dfrac{1}{r} \left( \dfrac{\partial}{\partial r} - \dfrac{2}{r} \right) Z = 0
   
which can be fulfilled with 

.. math:: 
   \left( \dfrac{\partial^2}{\partial r^2} - \dfrac{2}{r} \dfrac{\partial}{\partial r} \right) W = 0
   
and 
   
.. math::
   \left( \dfrac{\partial}{\partial r} - \dfrac{2}{r} \right) Z = 0\;.
   
In spectral representation this then reads

.. math:: 
   \left({\cal C}''_n - \dfrac{2}{r}{\cal C}'_n
   \right) W_{\ell mn} = 0 \;\;\mbox{and}\;\;
   \left({\cal C}'_n - \frac{2}{r}{\cal C}_n
   \right) Z_{\ell mn} = 0\;.
   
   
Thermal boundary conditions
---------------------------

For Entropy or temperature in the Boussinesq approximation either fixed value of fixed flux conditions are used. 
The former implies 

.. math:: s=\mbox{const.}\;\;\mbox{or}\;\;T=\mbox{const.}

at :math:`r_i` and/or :math:`r_o`, while the latter means

.. math:: \dfrac{\partial}{\partial r} s=\mbox{const.}\;\;\mbox{or}\;\;\dfrac{\partial}{\partial r}  T=\mbox{const.}

In spectral representation for example the respective entropy condition read

.. math:: {\cal C}_n s_{\ell mn}=\mbox{const.}\;\;\mbox{or}\;\;{\cal C'}_n s_{\ell mn}=\mbox{const.}

Appropriate constant values need to be chosen and are instrumental in driving the dynamo when 
flux conditions are imposed. 


Boundary conditions for chemical composition
--------------------------------------------

For the chemical composition, either the value or the flux is imposed at the
boundaries. The former implies:

.. math:: \xi=\mbox{const.}

at :math:`r_i` and/or :math:`r_o`, while the latter means

.. math:: \dfrac{\partial}{\partial r} \xi=\mbox{const.}

In spectral representation, this then reads

.. math:: {\cal C}_n \xi_{\ell mn}=\mbox{const.}\;\;\mbox{or}\;\;{\cal C'}_n \xi_{\ell mn}=\mbox{const.}


Magnetic boundary conditions and inner core
-------------------------------------------

Three different magnetic boundary conditions are implemented in MagIC. 
The most simple one is the conceptual condition at the boundary to an infinite conductor. 
Surface current in this conductor will prevent the internally produced magnetic 
field from penetrating so that the field has to vanish at the boundary. 
The condition are thus the same as for a rigid flow (with boundaries at rest). 
We only provide the spectral representation here: 

.. math::
  {\cal C}_n(r) W_{\ell mn} = 0 \;\;\mbox{at}\;\; r=r_i,r_o\;\;.
  :label: eqBcMagRigid

Note that the summation convention with respect to
radial modes :math:`n` is used again.
**The no-slip** condition further requires that the
horizontal flow components also have to vanish, provided
the two boundaries are at rest. This condition is fulfilled for
   
.. math::
   {\cal C}_n(r) g_{\ell mn} = 0\;\;,\;\;{\cal C}'_n(r) g_{\ell mn} = 0
   \;\;\mbox{and}\;\;{\cal C}_n(r) h_{\ell mn} = 0.
  :label: eqBcMag0

More complex are the conditions to an electrical insulator. 
Here we actually use matching condition to a potential field condition
that are formulated like boundary conditions.  
Since the electrical currents have to vanish in the insulator we have :math:`\nabla\times\vec{B}`, 
which means that the magnetic field is a potential field :math:`\vec{B}^I=-\vec{\nabla} V` 
with :math:`\Delta V=0`. This Laplace equation implies a coupling between radial and 
horizontal derivatives which is best solved in spectral space. Two potential contributions 
have to be considered depending whether the field is produced above the interface radius 
:math:`r_{BC}` or below. We distinguish these contributions with upper indices I for internal 
or below and E for external or above. The total potential then has the form:

.. math:: V_{\ell m}(r) = r_{BC} V_{\ell m}^I \left(\dfrac{r_{BC}}{r}\right)^{\ell+1} + 
          r_{BC} V_{\ell m}^E \left(\dfrac{r}{r_{BC}}\right)^{\ell}.

with the two spectral potential representations :math:`V_{\ell m}^I` and  :math:`V_{\ell m}^E`. 
This provides well defined radial derivative for both field contributions. 
For boundary :math:`r_o` we have to use the first contribution and match the respective field as well 
as its radial derivative to the dynamo solution. The toroidal field cannot penetrate the 
insulator and thus simply vanishes which yields :math:`h=0` or 

.. math:: {\cal C}_n h_{\ell mn} = 0

in spectral space. The poloidal field then has to match the potential field which implies 

.. math:: \nabla_H \dfrac{\partial}{\partial r} g = -\nabla_H V^I

for the horizontal components and 

.. math:: \dfrac{\nabla_H^2}{r^2} g = \dfrac{\partial}{\partial r} V^I

for the radial. In spectral space these condition can be reduce to

.. math:: {\cal C'}_n(r) g_{\ell mn} = V^I_{lm}\;\;\mbox{and}\;\;
          \dfrac{\ell(\ell+1)}{r^2} {\cal C}_n g_{\ell mn} = - \dfrac{\ell+1}{r} V^I_{lm}.
          
Combining both allows to eliminate the potential and finally leads to the spectral condition used in MagIC:

.. math:: \left( {\cal C'}_n(r_o) + \frac{\ell}{r_o} {\cal C}_n(r_o) \right) g_{\ell mn} = 0

Analogous consideration lead to the respective condition at the interface to an insulating inner core:

.. math:: \left( {\cal C'}_n(r_i) - \frac{\ell+1}{r_i} {\cal C}_n(r_i) \right) g_{\ell mn} = 0.

If the inner core is modelled as an electrical conductor, a simplified dynamo
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
   - i\,\omega_I\,m\,\dfrac{\ell(\ell+1)}{r^2} \;h_{\ell m}^I(r)\;.
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

