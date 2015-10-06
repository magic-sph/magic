Formulation of the (magneto)-hydrodynamics problem
##################################################

The general equations describing thermal convection and dynamo action of a
rotating compressible fluid are the starting point from which the Boussinesq or
the anelastic approximations are developed.  In MagIC, we consider a spherical
shell rotating about the vertical :math:`z` axis with a constant angular
velocity :math:`\Omega`. The conservation of mass is expressed by the
continuity equation:

.. math::
   \dfrac{\partial \rho}{\partial t} + \vec{\nabla} \cdot \rho \vec{u} = 0,
   :label: eqContinuity

The conservation of momentum by

.. math::
   \rho\left(\dfrac{\partial \vec{u}}{\partial t}+ \vec{u}\cdot\vec{\nabla}\
   \vec{u} \right) =-\vec{\nabla} p +
   \dfrac{1}{\mu_0}(\vec{\nabla}\times\vec{B})\times\vec{B} +\rho
   \vec{g}-2\rho\vec{\Omega}\times\vec{u}+ \vec{\nabla}\cdot\mathsf{S},
   :label: eqNS

where :math:`\mathsf{S}` corresponds to the rate-of-strain tensor given by:

.. math::
   \begin{aligned}
   S_{ij} & = 2\nu\rho\left[e_{ij}-\dfrac{1}{3}\delta_{ij}\,\vec{\nabla}\cdot\vec{u} \right], \\
   e_{ij} & =\dfrac{1}{2}\left(\dfrac{\partial u_i}{\partial x_j}+\dfrac{\partial
   u_j}{\partial x_i}\right).
   \end{aligned}

Concerning the energy equation, several forms are possible (using internal energy, temperature
or entropy). Here we use entropy :math:`s` as the main variable, which leads to:

.. math::
  \rho T\left(\dfrac{\partial s}{\partial t}+\vec{u}\cdot \vec{\nabla} s \right) = \vec{\nabla}\cdot (K\vec{\nabla} T) + \Phi_\nu +\lambda\left(\vec{\nabla}\times\vec{B}\right)^2,
  :label: eqEntropy

where :math:`\Phi_\nu` corresponds to the viscous heating expressed by

.. math::
   \Phi_\nu = 2\rho\left[e_{ij}e_{ji}-\dfrac{1}{3}\left(\vec{\nabla}\cdot\vec{u}\right)^2\right].

The induction equation is obtained from the Maxwell equations (ignoring displacement current)
and Ohm's law (neglecting Hall effect):

.. math::
  \dfrac{\partial \vec{B}}{\partial t} = \vec{\nabla} \times \left( \vec{u}\times\vec{B}-\lambda\,\vec{\nabla}\times\vec{B}\right).
  :label: eqInduction

In those equations, the symbols :math:`\vec{u}`, :math:`\vec{B}`, :math:`p` and
:math:`s` correspond to velocity, magnetic field, pressure and entropy.
:math:`\vec{g}` corresponds to gravity and :math:`\rho` to density. :math:`\lambda`
is the magnetic diffusivity, :math:`\mu_0` the magnetic permeability, :math:`\nu` the kinematic
viscosity and :math:`K` the thermal conductivity.

.. figure:: figs/shell.png
   :width: 500px
   :align: center
   :alt: caption

   Sketch of the spherical shell model and its system of coordinate.

The anelastic approximation
===========================

MagIC relies on the anelastic approximation of the Navier Stokes equations. There are
different flavours of this approximation but they all assume that:

   1. the departures from the thermodynamic state variables due to convection are small
      with respect to the reference state,

   2. the short-period acoustic waves are filtered-out.

This second assumption implies that larger numerical timesteps will be allowed, since
the typical timescale of the fast acoustic waves is typically much slower than the convective
turnover time. The strict elimination of the acoustic waves formally means

.. math::

   \frac{\partial \rho'}{\partial t}=0

in the continuity equation, where :math:`\rho'` corresponds here to the density perturbations
of the reference state.

The reference state is a background state against which perturbations are
described. In MagIC, this background state is assumed to only depends on one
spatial variable (radius), but in general it may as well be explicitly
time-dependent. Each thermodynamic variable :math:`f` is then expressed as a sum
of a spherically-symmetric time-independent quantity :math:`\tilde{f}` and a fluctuating
quantity :math:`f'`:

.. math::

  f(r,\theta,\phi,t) = \tilde{f}(r)+f'(r,\theta,\phi,t)


This separation of variables is then introduced in the set of equations
:eq:`eqContinuity`-:eq:`eqInduction` to perform a formal **scale analysis**, taking
into account that :math:`f'/\tilde{f} = \epsilon << 1`. The variables are then expanded
in power series of :math:`\epsilon` and only the highest order terms are retained.


An adiabatic reference state
============================

In a vigorously convecting astrophysical of geophysical system (like the convecting interior
of a planet or a star), the super-adiabaticity  of the fluid is extremely small, since the
transport of heat by convective motions is highly efficient. Therefore, the reference background state can be assumed to be perfectly adiabatic and obey to the following equations

.. math::
   \dfrac{d \tilde{T}}{dr} = -\dfrac{\alpha g \tilde{T}}{c_p},
   :label: eqAdiabatRef
   

where :math:`c_p` is the heat capacity and :math:`\alpha` expressed by

.. math::
   c_p = T\left(\dfrac{\partial s}{\partial T}\right)_p, \quad\text{and}\quad
  \alpha = -\dfrac{1}{\rho}\left(\dfrac{\partial\rho}{\partial T}\right)_p.


At this stage, it becomes convenient to start introducing non-dimensionalised quantities.
The background quantities (density temperature and transport properties) are non-dimensionalised using their values at the spherical shell outer boundary :math:`r_o`. The shell thickness
:math:`d=r_o-r_i` is used as the reference lenghtscale. The dimensionless form of Eq. :eq:`eqAdiabatRef` then reads:

.. math::
  \dfrac{d \tilde{T}}{d r} = -Di\,\alpha(r) g(r) \tilde{T}(r),
  :label: eqAdiaRefNd

where :math:`\tilde{T}`, :math:`\alpha` and :math:`g` have been non-dimensionalised using their values at the outer boundary. :math:`Di` is the dissipation number expressed by

.. math::
   Di = \dfrac{\alpha_o g_o d}{c_p}.
   :label: eqDissipNb

:math:`Di` is a measure of the thermal effects due to compressibility, namely
viscous and ohmic heating. :math:`Di` is also the ratio between two length
scales: the thickness of the spherical shell :math:`d` divided by the
temperature scale heights :math:`H_T=-(d\ln T/dr)^{-1}`.

When :math:`Di \ll 1`, the region where convection develops has a nearly constant reference
temperature. Since there is no basic temperature stratification, viscous heating (which is
the consequence of the thermal stratification due to compressibility) becomes negligible.

.. note:: The Boussinesq limit can thus be recovered by using :math:`Di \rightarrow 0`.

Provided an equation of state is given, it is then possible to integrate Eq. :eq:`eqAdiaRefNd`,to obtain the adiabatic background state.

Analytical solution in the limit of an ideal gas
------------------------------------------------

In the limit of an ideal gas which follows :math:`\tilde{p}=\tilde{\rho}\tilde{T}` and has
:math:`\alpha=1/\tilde{T}`, one directly gets:

.. math::
   \begin{aligned}
   \dfrac{d \tilde{T}}{dr}  & = -Di\,g(r), \\
   \tilde{\rho} & = \tilde{T}^{1/(\gamma-1)},
   \end{aligned}

where :math:`\gamma=c_p/c_v`. If we now in addition make the assumption of a
centrally-condensed mass in the center of the spherical shell of radius
:math:`r\in[r_i,r_o]`, i.e. :math:`g\propto, 1/r^2`, this leads to

.. math::
   \begin{aligned}
    \tilde{T}(r) & =Di\frac{r_o^2}{r}+(1-Di\,r_o), \\
    \tilde{\rho}(r) & = \tilde{T}^m, \\
    Di & = \dfrac{r_i}{r_o}\left(\exp\dfrac{N_\rho}{m}-1\right),
   \end{aligned}

where :math:`N_\rho=\ln(\rho_i/\rho_o)` is the number of density scale heights of the reference
state and :math:`m=1/(\gamma-1)` is the polytropic index.
   

.. warning:: * The relationship between :math:`N_\rho` and the dissipation number
               :math:`Di` directly depends on the gravity profile. The formula above
               is only valid when :math:`g\propto 1/r^2`.
             * In this formulation, when you change the polytropic index :math:`m`, you
               also change the nature of the fluid you're modelling since you accordingly
               modify :math:`\gamma=c_p/c_v`.


MHD equations
=============

One of the assumptions of the anelastic approximation is that the fluctuations due to convection are much smaller than the reference state:

.. math::
   \epsilon \simeq \dfrac{\rho'}{\tilde{\rho}}\simeq \dfrac{T'}{\tilde{T}}\simeq \dfrac{p'}{\tilde{p}}\simeq s' \ll 1.

In the following, we will treat the equations
:eq:`eqContinuity`-:eq:`eqInduction` in nondimensional form. There is no unique
way to scale the equations and as a consequence different sets of
non-dimensional numbers are employed. For convection-driven dynamos, there is four
independent control parameters.

We use here the viscous diffusion time :math:`d^2/\nu_o` (where :math:`\nu_o` is the kinematic
viscosity at the outer boundary as a time unit and :math:`\nu_o/d` as the reference velocity.
Magnetic field is expressed in units of :math:`\sqrt{\rho_o\mu_0\lambda_i\Omega}`, where
:math:`\rho_o` is the density at the outer boundary and :math:`\lambda_i` is the magnetic
diffusivity at the **inner** boundary.

.. note:: All the transport properties except the magnetic diffusivity are normalised to their
          values at the outer boundary. The motivation to rather base the reference magnetic
          diffusivity to the **inner** boundary is twofold: (i) it allows an easier control
          of the possible continuous conductivity value in the inner core; (ii) it is a more
          natural choice when modelling gas giants planets which exhibit strong electrical 
          conductivity decays in the outer layer.

This leads to the following sets of dimensionless equations:

.. math::
   E\left(\dfrac{\partial \vec{u}}{\partial t}+\vec{u}\cdot\vec{\nabla}\vec{u}\right)
   +2\vec{e_z}\times\vec{u}= -\vec{\nabla}\left({\dfrac{p'}{\tilde{\rho}}}\right)+\dfrac{Ra\,E}{Pr}g(r)
   \,s'\,\vec{e_r} + \dfrac{1}{Pm\,\tilde{\rho}}\left(\vec{\nabla}\times \vec{B} 
   \right)\times \vec{B}+ \dfrac{E}{\tilde{\rho}} \vec{\nabla}\cdot \mathsf{S},
   :label: eqNSNd

.. math::
   \vec{\nabla}\cdot\tilde{\rho}\vec{u}=0,
   :label: eqContNd

.. math::
   \vec{\nabla}\cdot\vec{B}=0,
   :label: eqMagNd

.. math::
   \dfrac{\partial \vec{B}}{\partial t} = \vec{\nabla} \times \left( \vec{u}\times\vec{B}\right)-\dfrac{1}{Pm}\vec{\nabla}\times\left(\lambda(r)\,\vec{\nabla}\times\vec{B}\right).
   :label: eqIndNd

Entropy equation and turbulent diffusion
----------------------------------------

The entropy equation usually requires an additional assumption in most of the existing 
anelastic approximations. Indeed, if one simply expands Eq. :eq:`eqEntropy` with the classical
temperature diffusion an operator of the form:

.. math::
   \epsilon\,\vec{\nabla}\cdot \left( K \vec{\nabla} T'\right)+\vec{\nabla}\cdot \left( K \vec{\nabla} \tilde{T}\right),

will remain the right-hand side of the equation. At first glance, there seems
to be a :math:`1/\epsilon` factor between the first term and the second one,
which would suggest to keep only the second term in this expansion. However,
for astrophysical objects which exhibit strong convective driving (and hence
large Rayleigh numbers), the diffusion of the adiabatic background is actually
very small and may be comparable or even smaller in magnitude as the :math:`\epsilon`
terms representing the usual convective perturbations. For the Earth core for instance,
if one assumes that the typical temperature fluctuations are of the order of 1 mK and
the temperature contrast between the inner and outer core is of the order of 1000 K, then
:math:`\epsilon \sim 10^{-6}`. The ratio of the two terms can thus be estimated as

.. math:: \epsilon \dfrac{T'/\delta^2}{T/d^2},
   :label: eqEpsRatio

where :math:`d` is the thickness of the inner core and :math:`\delta` is the typical thermal
boundary layer thickness. This ratio is exactly one when :math:`\delta =1\text{ m}`, a
plausible value for the Earth inner core. 

In numerical simulations however, the over-estimated diffusivities restrict the computational
capabilities to much lower Rayleigh numbers. As a consequence, the actual boundary layers
in a global DNS will be much thicker and the ratio :eq:`eqEpsRatio` will be much smaller than
unity. The second terms will thus effectively acts  as a radial-dependent heat source or sink
that will drive or hinder convection. This is one of the physical motivation to rather introduce a **turbulent diffusivity** that will be approximated by

.. math:: \kappa \tilde{\rho}\tilde{T} \vec{\nabla} s,

where :math:`\kappa` is the turbulent diffusivity. Entropy diffusion is assumed to dominate
over temperature diffusion in turbulent flows.

The choice of the entropy scale to non-dimensionalise Eq. :eq:`eqEntropy` also depends on
the nature of the boundary conditions: it can be simply the entropy contrast over the layer
:math:`\Delta s` when the entropy is held constant at both boundaries, or :math:`d\,(ds /dr)`
when flux-based boundary conditions are employed. We will restrict to the first option in
the following, but keep in mind that depending on your setup, the entropy reference scale
(and thus the Rayleigh number definition) might change.


.. math::
  \tilde{\rho}\tilde{T}\left(\dfrac{\partial s'}{\partial t} + 
  \vec{u}\cdot\vec{\nabla} s'\right) =
  \dfrac{1}{Pr}\vec{\nabla}\cdot\left(\kappa(r)\tilde{\rho}\tilde{T}\vec{\nabla} s'\right) +
  \dfrac{Pr\,Di}{Ra}\Phi_\nu +
  \dfrac{Pr\,Di}{Pm^2\,E\,Ra}\lambda(r)\left(\vec{\nabla}
  \times\vec{B}\right)^2,
  :label: eqEntropyNd

The Boussinesq limits of the equation :math:`Di \rightarrow 0`
--------------------------------------------------------------

When the dissipation number :math:`Di\rightarrow 0` then :math:`\tilde{T}=\text{cst.}`.
If in addition to that if :math:`\gamma \neq 1`, the density background :math:`\tilde{\rho}`
is also constant. 

.. note:: The peculiar configuration of :math:`\gamma=1` corresponds to the so-called
          zero-Grüneisen limit of the Navier-Stokes equation (or isothermal) and is 
          a special case in which :math:`Di=0` but a density background (controlled 
          by :math:`N_\rho`) is still allowed. 

A brief look at Eq. :eq:`eqEntropyNd` then shows than viscous and Ohmic heating will disappear
from the entropy equation. Furthermore, temperature and entropy fluctuations become equivalent
quantities. If in addition to that we also neglect the possible radial-dependence of the
transport properties (electrical conductivity, viscosity and thermal diffusivity),
the set of equations :eq:`eqNSNd`-:eq:`eqEntropyNd` thus simplifies to the classical
Boussinesq set of equations:

.. math::
   E\left(\dfrac{\partial \vec{u}}{\partial t}+\vec{u}\cdot\vec{\nabla}\vec{u}\right)
   +2\vec{e_z}\times\vec{u}= -\vec{\nabla}p'+\dfrac{Ra\,E}{Pr}g(r)
   \,T'\,\vec{e_r} + \dfrac{1}{Pm}\left(\vec{\nabla}\times \vec{B} 
   \right)\times \vec{B}+ E\,\Delta \vec{u},

.. math::
   \vec{\nabla}\cdot\vec{u}=0,

.. math::
   \vec{\nabla}\cdot\vec{B}=0,

.. math::
   \dfrac{\partial \vec{B}}{\partial t} = \vec{\nabla} \times \left( \vec{u}\times\vec{B}\right)+\dfrac{1}{Pm}\Delta\vec{B}.

.. math::
  \dfrac{\partial T'}{\partial t} + 
  \vec{u}\cdot\vec{\nabla} T' =
  \dfrac{1}{Pr}\Delta T'.




Dimensionless control parameters
--------------------------------

The equations :eq:`eqNSNd`-:eq:`eqEntropyNd` are governed by four dimensionless numbers: the
Ekman number

.. math::
   E = \frac{\nu}{\Omega d^2},
   :label: eqEkman

the Rayleigh number

.. math::
   Ra = \frac{\alpha_o g_o T_o d^3 \Delta s}{c_p \kappa_o \nu_o},
   :label: eqRayleigh

the Prandtl number

.. math::
   Pr = \frac{\nu_o}{\kappa_o},
   :label: eqPrandtl

and the magnetic Prandtl number

.. math::
   Pm = \frac{\nu_o}{\lambda_i}.
   :label: eqmaPrandtl

In addition to these four numbers, the reference state is controlled by the geometry of
the spherical shell given by its radius ratio

.. math::
   \eta = \frac{r_i}{r_o},
   :label: eqRadratio

and the background density and temperature profiles, either controlled by :math:`Di` or
by :math:`N_\rho` and :math:`m`.

Variants of the non-dimensional equations and control parameters result from
different choices for the fundamental scales. For the length scale often
:math:`r_o` is chosen instead of :math:`d`. Other natural scales for time are the
magnetic or the thermal diffusion time, or the rotation period.
There are also different options for scaling the magnetic field strength.
The prefactor of two, which is retained in the
Coriolis term in :eq:`eqNSNd`, is often incorporated into the definition of the
Ekman number.


Usual diagnostic quantities
---------------------------

Characteristic properties of the solution are usually expressed in terms
of non-dimensional diagnostic parameters.
In the context of the geodynamo for instance, the two
most important ones are the magnetic Reynolds number :math:`Rm` and
the Elsasser number :math:`\Lambda`. Usually the rms-values of the velocity
:math:`u_{rms}` and of the magnetic field :math:`B_{rms}` inside the spherical shell
are taken as characteristic values. The magnetic Reynolds number

.. math::
   Rm =  \frac{u_{rms}d}{\lambda_i}

can be considered as a measure for the flow velocity and describes
the ratio of advection of the magnetic field to magnetic diffusion.
Other characteristic non-dimensional numbers related to the flow velocity are
the (hydrodynamic) Reynolds number

.. math::
   Re = \frac{u_{rms} d}{\nu_o},

which measures the ratio of inertial forces to viscous forces,
and the Rossby number

.. math::
   Ro = \frac{u_{rms}}{\Omega d} ,

a measure for the ratio of inertial to Coriolis forces.

.. math::
   \Lambda = \frac{B_{rms}^2}{\mu_0\lambda_i\rho_o\Omega}

measures the ratio of Lorentz to Coriolis forces and is
equivalent to the square of the non-dimensional magnetic field strength
in the scaling chosen here.



Boundary conditions and treatment of inner core
===============================================

Mechanical conditions
---------------------

In its simplest form, when modelling the geodynamo, the fluid shell is treated
as a container with rigid, impenetrable, and co-rotating walls. This implies
that within the rotating frame of reference all velocity components vanish at
:math:`r_o` and :math:`r_i`.  In case of modelling the free surface of a gas
giant planets or a star, it is preferable to rather replace the condition of
zero horizontal velocity by one of vanishing viscous shear stresses (the
so-called free-slip condition).

Furthermore, even in case of modelling the liquid iron core of a terrestrial
planet, there is no a priori reason why the inner core should necessarily
co-rotate with the mantle. Some models for instance allow for differential
rotation of the inner core and mantle with respect to the reference frame.  The
change of rotation rate is determined from the net torque. Viscous,
electromagnetic, and torques due to gravitational coupling between density
heterogeneities in the mantle and in the inner core contribute.

Magnetic boundary conditions and inner core conductivity
--------------------------------------------------------

When assuming that the fluid shell is surrounded by electrically insulating  regions
(inner core and external part),
the magnetic field inside the fluid shell matches continuously
to a potential field in both the exterior and the interior regions. Alternative
magnetic boundary conditions (like cancellation of the horizontal component of the field
) are also possible.

Depending on the physical problem you want to model, treating the inner core as an 
insulator is not realistic either, and it might instead be more appropriate to
assume that it has the same electrical conductivity as
the fluid shell. In this case, an equation equivalent to :eq:`eqIndNd` must
be solved for the inner core, where the velocity field simply
describes the solid body rotation of the inner core with respect to the reference frame.
At the inner core boundary a continuity condition for the magnetic field and the
horizontal component of the electrical field apply.

Thermal boundary conditions and distribution of buoyancy sources
----------------------------------------------------------------

In many dynamo models, convection is simply driven by an imposed fixed
super-adiabatic entropy contrast between the inner and outer boundaries.  This
approximation is however not necessarily the best choice, since for instance,
in the present Earth,  convection is thought to be driven by a combination of
thermal and compositional buoyancy.  Sources of heat are the release of latent
heat of inner core solidification and the secular cooling of the outer and
inner core, which can effectively be treated like a heat source.  The heat loss
from the core is controlled by the convecting mantle, which effectively imposes
a condition of fixed heat flux at the core-mantle boundary on the dynamo. The
heat flux is in that case spatially and temporally variable. 
