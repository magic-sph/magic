.. _secPhysNml:

Physical parameters namelist
============================

This namelist contains all the appropriate relevant control physical parameters.

Dimensionless control parameters
--------------------------------

* **ra** (default :f:var:`ra=0.0 <ra>`) is a real. This the thermal Rayleigh number expressed by

  .. math::
     Ra = \frac{\alpha g_o \Delta T d^3}{\kappa\nu}

* **raxi** (default :f:var:`raxi=0.0 <raxi>`) is a real. This the compositional Rayleigh number expressed by

  .. math::
     Ra_\xi = \frac{\alpha g_o \Delta \xi d^3}{\kappa_\xi\nu}

* **ek** (default :f:var:`ek=1e-3 <ek>`) is a real. This is the Ekman number expressed by

  .. math::
     E = \frac{\nu}{\Omega d^2}

* **pr** (default :f:var:`pr=1.0 <pr>`) is a real. This is the Prandtl number expressed by

  .. math::
     Pr = \frac{\nu}{\kappa}

* **sc** (default :f:var:`sc=10.0 <sc>`) is a real. This is the Schmidt number expressed by

  .. math::
     Sc = \frac{\nu}{\kappa_\xi}

* **prmag** (default :f:var:`prmag=5.0 <prmag>`) is a real. This is the magnetic Prandtl number expressed by

  .. math::
     Pm = \frac{\nu}{\lambda}

* **po** (default :f:var:`po=0.0 <po>`) is a real. This is the Poincaré number expressed by

  .. math::
     Po = \frac{\Omega_p}{\Omega}

* **prec_angle** (default :f:var:`prec_angle=23.5 <prec_angle>`) is a real. This is the angle between the precession and the rotation axes expressed in degrees.

* **radratio** (default :f:var:`radratio=0.35 <radratio>`) is a real. This is the ratio of the inner core radius :math:`r_i` to the outer core radius :math:`r_o`:

  .. math::
     \eta = \frac{r_i}{r_o}

.. _varstrat:

* **strat** (default :f:var:`strat=0.0 <strat>`) is a real. This is the number of density scale heights of the reference state:

  .. math::
     N_\rho = \ln \frac{\tilde{\rho}(r_i)}{\tilde{\rho}(r_o)}

* **DissNb** (default :f:var:`DissNb=0.0 <dissnb>`) is a real. This is the dissipation number:

  .. math::
     Di = \frac{\alpha_o g_o d}{c_p}

  .. warning:: This can only be provided  as a **replacement** input of **strat**.
               I.E., when one wants to define a reference state, one has to specify
               **either** ``strat`` **or** ``DissNb`` in the input namelist.

* **polind** (default :f:var:`polind=1.5 <polind>`) is a real. This is the polytropic index, which relates the background temperature to the background density:

  .. math::
     \tilde{\rho} = \tilde{T}^m

  ..

  .. warning:: Be careful: in its current version the code only handles **adiabatic** backgrounds, therefore changing ``polind`` physically means that the nature of the fluid (in particular its Grüneisen parameter) will change. For an ideal gas, it actually always follows :math:`m+1=\frac{\gamma -1}{\gamma}`

  ..

* **l_isothermal** (default :f:var:`l_isothermal=.false.`) is a logical. When set to ``.true.``, makes the temperature background isothermal (i.e. :math:`\tilde{T}=cst.`). In that case, the dissipation number :math:`Di` vanishes and there is no viscous and Ohmic heating left. The only difference with the Boussinesq set of equations are thus restricted to the density background :math:`\tilde{\rho}` and its radial derivatives that enters the viscous stress. This approximation is also called the **zero Grüneisen parameter** and was extensively explored by Denise Tortorella during her `PhD <http://www.mps.mpg.de/3183008/Dissertation_2005_Tortorella__Denise_Aida1.pdf>`_.


.. _varepsc0:

Heat sources and sinks
----------------------

* **epsc0** (default :f:var:`epsc0=0.0 <epsc0>`) is a real. This is the volumetric heat source :math:`\epsilon_0` that enters the thermal equilibrium relation:

  .. math::
     -\nabla\cdot\left(\tilde{\rho}\tilde{T}\nabla s\right) + \epsilon_0\,f(r)=0
     :label: heatEq

  ..

  The radial function :math:`f(r)` can be modified with the variable ``nVarEps`` that enters the same input namelist.

* **epscxi0** (default :f:var:`epscxi0=0.0 <epscxi0>`) is a real. This is the volumetric source :math:`\epsilon_\xi` that enters the compositional equilibrium relation:

  .. math::
     -\nabla\cdot\left(\tilde{\rho}\nabla \xi\right) + \epsilon_\xi=0
     :label: xiEq

  ..

* **nVarEps** (default :f:var:`nVarEps=0 <nvareps>`) is an integer. This is used to modify the radial-dependence ofthe volumetric heat source, i.e. :math:`f(r)` that enters equation :eq:`heatEq`.

  +---------------+-------------------------------------------------------------+
  | ``nVarEps=0`` | Constant, i.e. :math:`f(r)=\hbox{cst.}`.                    |
  +---------------+-------------------------------------------------------------+
  | ``nVarEps=1`` | Proportional to density, i.e. :math:`f(r)=\tilde{\rho}(r)`. |
  +---------------+-------------------------------------------------------------+
  | ``nVarEps=2`` | Proportional to density times temperature,                  |
  |               | i.e. :math:`f(r)=\tilde{\rho}(r)\tilde{T}`.                 |
  +---------------+-------------------------------------------------------------+

.. _varinterior_model:

Realistic interior models
-------------------------

* **interior_model** (default :f:var:`interior_model="None" <interior_model>`) is a character string. This defines a polynomial fit of the density profile of the interior structure of several astrophysical objects. Possible options are ``"earth"``, ``"jupiter"``, ``"saturn"`` and ``"sun"`` (the naming is **not** case sensitive).

  .. warning:: When ``interior_model`` is defined the variables ``strat``, ``polind``, ``g0``, ``g1`` and ``g2`` are not interpreted.

  ..

  The subroutine :f:subr:`radial <radial_functions/radial>` gives the exact details of the implementation.

* **r_cut_model** (default :f:var:`r_cut_model=0.98 <r_cut_model>`) is a real. This defines the cut-off radius of the reference model, i.e. the fluid domain is restricted to radii with :math:`r\leq r_{cut}`.

The following input parameters will thus define a polynomial fit to the expected interior structure of Jupiter until 99% of Jupiter's radius (assumed here at the 1 bar level)

   .. code-block:: fortran

           interior_model="JUP",
	   r_cut_model   =0.99e0,


Gravity
-------

The radial dependence of the gravity profile can be adjusted following

.. math::
   g(r)=g_0+g_1\,\frac{r}{r_o}+g_2\,\left(\frac{r_o}{r}\right)^2
   :label: eqGravity

The three following parameters are used to set this profile

* **g0** (default :f:var:`g0=0 <g0>`) is the pre-factor of the constant part of the gravity profile, i.e. :math:`g_0` in equation :eq:`eqGravity`.

* **g1** (default :f:var:`g1=1 <g1>`) is the pre-factor of the linear part of the gravity profile, i.e. :math:`g_1` in equation :eq:`eqGravity`.

* **g2** (default :f:var:`g2=0 <g2>`) is the pre-factor of the :math:`1/r^2` part of the gravity profile, i.e. :math:`g_2` in equation :eq:`eqGravity`.


Centrifugal acceleration
------------------------

The centrifugal acceleration can be computed for a polytropic background

* **dilution_fac** (default :f:var:`dilution_fac=0.0 <dilution_fac>`) is the ratio of the centrifugal acceleration at the equator to the surface gravitational acceleration.

  .. math::
     m=\frac{\Omega^2 d}{g_o}
     :label: dilution_fac

  ..

Phase field
------------

* **stef** (default :f:var:`stef=1.0 <stef>` is a real. This is the Stefan number (used when the phase field is plugged in). It is expressed by the ratio of the latent heat per unit mass associated with the solid-liquid transition and the specific heat:

  .. math::
     St = \frac{\mathcal{L}}{c_p\Delta T}

* **tmelt** (default :f:var:`tmelt=0.0 <tmelt>`) is a real. This is the dimensionless melting temperature.

* **epsPhase** (default :f:var:`epsPhase=0.01 <epsPhase>`) is a real. This is the dimensionless interface thickness between the solid and the liquid phase (sometimes known as the Cahn number).

* **phaseDiffFac** (default :f:var:`phaseDiffFac=1.0 <phaseDiffFac>`) is a real. This is a coefficient that goes in front of the diffusion term in the phase field equation.

* **penaltyFac** (default :f:var:`penaltyFac=1.0 <penaltyFac>`) is a real. This is coefficient used for the penalisation of the velocity field in the solid phase. The smaller the coefficient, the stronger the penalisation. Since this is a nonlinear term, it is handled explicitly and the time step size should be decreased with the square of :f:var:`penaltyFac`.


Transport properties
--------------------

* **difExp** (default :f:var:`difExp=-0.5 <difexp>`) is a real. This is the exponent that is used when :f:var:`nVarVisc=2 <nvarvisc>`, :f:var:`nVarDiff=2 <nvardiff>` or :f:var:`nVarCond=4 <nvarcond>`.


.. _varnVarCond:

Electrical conductivity
+++++++++++++++++++++++

There are several electrical conductivity profiles implemented in the code that
can be chosen with the :f:var:`nVarCond <nvarcond>` input variable. The following one
corresponds to a constant electrical conductivity in the deep interior
(:math:`r<r_m`) and an exponential decay in the outer layer.

.. math::
  \sigma(r)=1+ (\sigma_m-1)\left(\frac{r-r_i}{r_m-r_i}\right)^a \quad \hbox{for}\quad r<r_m, \\
  \sigma(r)=\sigma_m \exp \left[a \left(\frac{r-r_m}{r_m-r_i}\right)\frac{\sigma_m-1}{\sigma_m}\right]
  \quad\hbox{for}\quad r\geq r_m.
  :label: eqElecCond

* **nVarCond** (default :f:var:`nVarCond=0 <nvarcond>`) is an integer. This is used to modify the radial-dependence of the electrical conductivity.

  +----------------+-----------------------------------------------------------------------+
  | ``nVarCond=0`` | Constant electrical conductivity, i.e. :math:`\sigma=\hbox{cst.}`     |
  +----------------+-----------------------------------------------------------------------+
  | ``nVarCond=1`` | :math:`\sigma\propto\tanh[a(r-r_m)]`                                  |
  +----------------+-----------------------------------------------------------------------+
  | ``nVarCond=2`` | See equation :eq:`eqElecCond`.                                        |
  +----------------+-----------------------------------------------------------------------+
  | ``nVarCond=3`` | Magnetic diffusivity proportional to :math:`1/\tilde{\rho}`, i.e.     |
  |                |   .. math::                                                           |
  |		   |      \lambda=\frac{\tilde{\rho}_i}{\tilde{\rho}}                      |
  +----------------+-----------------------------------------------------------------------+
  | ``nVarCond=4`` | Radial profile of the form:                                           |
  |                |   .. math::                                                           |
  |                |      \lambda=\left(\frac{\tilde{\rho}(r)}                             |
  |                |       {\tilde{\rho}_i}\right)^{\alpha}                                |
  +----------------+-----------------------------------------------------------------------+

* **con_RadRatio**  (default :f:var:`con_RadRatio=0.75 <con_radratio>`) is a real. This defines the transition radius :math:`r_m` that enters equation :eq:`eqElecCond`.

* **con_DecRate** (default :f:var:`con_DecRate=9 <con_decrate>`) is an integer. This defines the decay rate :math:`a` that enters equation :eq:`eqElecCond`.

* **con_LambdaMatch** (default :f:var:`con_LambdaMatch=0.6 <con_lambdamatch>`) is a real. This is the value of the conductivity at the transition point :math:`\sigma_m` that enters equation :eq:`eqElecCond`.

* **con_LambdaOut** (default :f:var:`con_LambdaOut=0.1 <con_lambdaout>`) is a real. This is the value of the conduvity at the outer boundary. This parameter is only used when ``nVarCond=1``.

* **con_FuncWidth** (default :f:var:`con_FuncWidth=0.25 <con_funcwidth>`) is a real. This parameter is only used when ``nVarCond=1``.

* **r_LCR**  (default :f:var:`r_LCR=2.0 <r_lcr>`) is a real. ``r_LCR`` possibly defines a low-conductivity region for :math:`r\geq r_{LCR}`, in which the electrical conductivity vanishes, i.e. :math:`\lambda=0`.

.. _varnVarDiff:

Thermal diffusivity
+++++++++++++++++++

* **nVarDiff** (default :f:var:`nVarDiff=0 <nvardiff>`) is an integer. This is used to change the radial-dependence of the thermal diffusivity:

  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=0`` | Constant thermal diffusivity :math:`\kappa`                                |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=1`` | Constant thermal conductivity, i.e.                                        |
  |                |    .. math:: \kappa =\frac{\tilde{\rho}_i}{\tilde{\rho}(r)}                |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=2`` | Radial profile of the form:                                                |
  |                |    .. math:: \kappa=\left(\frac{\tilde{\rho}(r)}                           |
  |                |              {\tilde{\rho}_i}\right)^{\alpha}                              |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=3`` | polynomial-fit to an interior model of Jupiter                             |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=4`` | polynomial-fit to an interior model of the Earth liquid core               |
  +----------------+----------------------------------------------------------------------------+

.. _varnVarVisc:

Viscosity
+++++++++

* **nVarVisc** (default :f:var:`nVarVisc=0 <nvarvisc>`) is an integer. This is used to change the radial-dependence of the viscosity:

  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=0`` | Constant kinematic viscosity :math:`\nu`                                |
  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=1`` | Constant dynamic viscosity, i.e.                                        |
  |                |    .. math:: \nu =\frac{\tilde{\rho}_o}{\tilde{\rho}(r)}                |
  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=2`` | Radial profile of the form:                                             |
  |                |    .. math:: \nu=\left(\frac{\tilde{\rho}(r)}                           |
  |                |              {\tilde{\rho}_i}\right)^{\alpha}                           |
  +----------------+-------------------------------------------------------------------------+

  where :math:`\alpha` is an exponent set by the namelist input variable ``difExp``.


Anelastic liquid equations
--------------------------

.. warning:: This part is still work in progress. The input parameters here are likely to
             be changed in the future.

* **epsS** (default :f:var:`epsS=0.0 <epss>`) is a real. It controls the deviation to the adiabat. It can be related to the small parameter :math:`\epsilon`:

  .. math:: \epsilon \simeq \frac{\Delta T}{T} \simeq \frac{\Delta s}{c_p}

* **cmbHflux** (default :f:var:`cmbHflux=0.0 <cmbhflux>`) is a real. This is the CMB heat flux that enters the calculation of the reference state of the liquid core of the Earth, when the anelastic liquid approximation is employed.

* **slopeStrat** (default :f:var:`slopeStrat=20.0 <slopestrat>`) is a real. This parameter controls the transition between the convective layer and the stably-stratified layer below the CMB.


Boundary conditions
-------------------

.. _secThermalBcs:

Thermal boundary conditions
+++++++++++++++++++++++++++

* **ktops** (default :f:var:`ktops=1 <ktops>`) is an  integer to specify the outer boundary entropy (or temperature) boundary condition:

  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=1`` | Fixed entropy at outer boundary: :math:`s(r_o)=s_{top}`                             |
  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=2`` | Fixed entropy flux at outer boundary: :math:`\partial s(r_o)/\partial r = q_t`      |
  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=3`` | Fixed temperature at outer boundary: :math:`T(r_o)=T_{top}`                         |
  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=4`` | Fixed temperature flux at outer boundary: :math:`\partial T(r_o)/\partial r = q_t`  |
  +-------------+-------------------------------------------------------------------------------------+

* **kbots** (default :f:var:`ktops=1 <kbots>`) is an  integer to specify the inner boundary entropy (or temperature) boundary condition.

* **s_top** (default :f:var:`s_top= 0 0 0.0 0.0 <s_top>`) is a real array of lateraly varying outer heat boundary conditions. Each four consecutive numbers are interpreted as follows:

  1. Spherical harmonic degree :math:`\ell`

  2. Spherical harmonic order :math:`m`

  3. Real amplitude (:math:`\cos` contribution)

  4. Imaginary amplitude (:math:`\sin` contribution)

  For example, if the boundary condition should be a combination of an :math:`(\ell=1,m=0)` sherical harmonic with the amplitude 1 and an :math:`(\ell=2,m=1)` spherical harmonic with the amplitude (0.5,0.5) the respective namelist entry could read:


  .. code-block:: fortran

     s_top = 1, 0, 1.0, 0.0, 2, 1, 0.5, 0.5, ! The comas could be left away.

* **s_bot** (default :f:var:`s_bot=0 0 0.0 0.0 <s_bot>`) is a real array. This is the same as ``s_top`` but for the bottom boundary.

* **impS** (default :f:var:`impS=0 <imps>`) is an integer. This is a  flag to indicate if there is a localized entropy disturbance, imposed at the CMB. The number of these input boundary conditions is stored in ``n_impS`` (the maximum allowed is 20), and it's given by the number of ``sCMB`` defined in the same namelist. The default value of ``impS`` is zero (no entropy disturbance). If it is set in the namelist for an integer greater than zero, then ``sCMB`` has to be also defined in the namelist, as shown below.

* **sCMB** (default :f:var:`sCMB=0.0 0.0 0.0 0.0 <scmb>`) is a real array of CMB heat boundary conditions (similar to the case of ``s_bot`` and ``s_top``). Each four consecutive numbers are interpreted as follows:

  1. Highest amplitude value of the entropy boundary condition, stored in array :f:var:`peakS(20) <peaks>`. When ``impS<0``, ``peakS`` is a relative amplitude in comparison to the :math:`(\ell=0,m=0)` contribution (for example, the case ``s_top= 0 0 -1 0``).

  2. :math:`\theta` coordinate (input has to be given in degrees), stored in array :f:var:`thetaS(20) <thetas>`.

  3. :math:`\phi` coordinate (input has to be given in degrees), stored in array :f:var:`phiS(20) <phis>`.

  4. Angular width (input has to be given in degrees), stored in array :f:var:`widthS(20) <widths>`.

.. _secCompBcs:

Boundary conditions for chemical composition
++++++++++++++++++++++++++++++++++++++++++++

* **ktopxi** (default :f:var:`ktopxi=1 <ktopxi>`) is an  integer to specify the outer boundary chemical composition boundary condition:

  +--------------+-------------------------------------------------------------------------------------+
  | ``ktopxi=1`` | Fixed composition at outer boundary: :math:`\xi(r_o)=\xi_{top}`                     |
  +--------------+-------------------------------------------------------------------------------------+
  | ``ktopxi=2`` | Fixed composition flux at outer boundary: :math:`\partial \xi(r_o)/\partial r = q_t`|
  +--------------+-------------------------------------------------------------------------------------+

* **kbotxi** (default :f:var:`ktopxi=1 <kbotxi>`) is an  integer to specify the inner boundary chemical composition boundary condition.

* **xi_top** (default :f:var:`xi_top= 0 0 0.0 0.0 <xi_top>`) is a real array of lateraly varying outer chemical composition boundary conditions. Each four consecutive numbers are interpreted as follows:

  1. Spherical harmonic degree :math:`\ell`

  2. Spherical harmonic order :math:`m`

  3. Real amplitude (:math:`\cos` contribution)

  4. Imaginary amplitude (:math:`\sin` contribution)

  For example, if the boundary condition should be a combination of an :math:`(\ell=1,m=0)` sherical harmonic with the amplitude 1 and an :math:`(\ell=2,m=1)` spherical harmonic with the amplitude (0.5,0.5) the respective namelist entry could read:


  .. code-block:: fortran

     xi_top = 1, 0, 1.0, 0.0, 2, 1, 0.5, 0.5, ! The comas could be left away.

* **xi_bot** (default :f:var:`xi_bot=0 0 0.0 0.0 <xi_bot>`) is a real array. This is the same as ``xi_top`` but for the bottom boundary.

* **impXi** (default :f:var:`impXi=0 <impxi>`) is an integer. This is a  flag to indicate if there is a localized chemical composition disturbance, imposed at the CMB. The number of these input boundary conditions is stored in ``n_impXi`` (the maximum allowed is 20), and it's given by the number of ``xiCMB`` defined in the same namelist. The default value of ``impXi`` is zero (no chemical composiiton disturbance). If it is set in the namelist for an integer greater than zero, then ``xiCMB`` has to be also defined in the namelist, as shown below.

* **xiCMB** (default :f:var:`xiCMB=0.0 0.0 0.0 0.0 <xicmb>`) is a real array of CMB chemical composition boundary conditions (similar to the case of ``xi_bot`` and ``xi_top``). Each four consecutive numbers are interpreted as follows:

  1. Highest amplitude value of the chemical composition boundary condition, stored in the array :f:var:`peakXi(20) <peakxi>`. When ``impXi<0``, ``peakXi`` is a relative amplitude in comparison to the :math:`(\ell=0,m=0)` contribution (for example, the case ``xi_top= 0 0 -1 0``).

  2. :math:`\theta` coordinate (input has to be given in degrees), stored in array :f:var:`thetaXi(20) <thetaxi>`.

  3. :math:`\phi` coordinate (input has to be given in degrees), stored in array :f:var:`phiXi(20) <phixi>`.

  4. Angular width (input has to be given in degrees), stored in array :f:var:`widthXi(20) <widthxi>`.


Boundary conditions for phase field
+++++++++++++++++++++++++++++++++++

.. _secPhaseFieldBcs:

* **ktopphi** (default :f:var:`ktopphi=1 <ktopphi>`) is an  integer to specify the boundary condition of the phase field at the outer boundary.

  +---------------+-------------------------------------------------------------------------------------+
  | ``ktopphi=1`` | Fixed phase field at outer boundary: :math:`\phi(r_o)=\phi_{top}`                   |
  +---------------+-------------------------------------------------------------------------------------+
  | ``ktopphi=2`` | Fixed phase field gradient : :math:`\partial \phi(r_o)/\partial r = 0`              |
  +---------------+-------------------------------------------------------------------------------------+

* **kbotphi** (default :f:var:`kbotphi=1 <kbotphi>`) is an  integer to specify the boundary condition of the phase field at the inner boundary.

.. _secMechanicalBcs:

Mechanical boundary conditions
++++++++++++++++++++++++++++++

* **ktopv** (default :f:var:`ktopv=2 <ktopv>`) is an integer, which corresponds to the mechanical boundary condition for :math:`r=r_o`.

  +-------------+--------------------------------------------------------------------+
  | ``ktopv=1`` | Stress-free outer boundary for :math:`r=r_o`:                      |
  |             |   .. math::                                                        |
  |             |      W_{\ell m}(r=r_o)=0, \quad                                    |
  |             |      \frac{\partial}{\partial r}\left(\frac{1}{r^2\tilde{\rho}}    |
  |             |      \frac{\partial W_{\ell m}}{\partial r}\right)=0 \\            |
  |             |      \frac{\partial}{\partial r}\left(\frac{1}{r^2\tilde{\rho}}    |
  |             |       Z_{\ell m}\right)=0                                          |
  +-------------+--------------------------------------------------------------------+
  | ``ktopv=2`` | Rigid outer boundary for :math:`r=r_o`:                            |
  |             |    .. math::                                                       |
  |             |       W_{\ell m}=0,\quad                                           |
  |             |       \frac{\partial W_{\ell m}}{\partial r}=0, \\                 |
  |             |       Z_{\ell m}=0                                                 |
  +-------------+--------------------------------------------------------------------+


* **kbotv** (default :f:var:`kbotv=2 <kbotv>`) is an integer, which corresponds to the mechanical boundary condition for :math:`r=r_i`.

.. _secMagneticBcs:

Magnetic boundary conditions
++++++++++++++++++++++++++++


* **ktopb** (default :f:var:`ktopb=1 <ktopb>`) is an integer, which corresponds to the magnetic boundary condition for :math:`r=r_o`.

  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=1`` | Insulating outer boundary:                                                      |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial g_{\ell m}}{\partial r}+\frac{\ell}{r}\,g_{\ell m}=0,\quad |
  |             |       \frac{\partial h_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=2`` | Perfect condutor:                                                               |
  |             |    .. math::                                                                    |
  |             |       g_{\ell m} = \frac{\partial^2 g_{\ell m}}{\partial r^2}=0,\quad           |
  |             |       \frac{\partial h_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=3`` | Finitely conducting mantle                                                      |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=4`` | Pseudo-vacuum outer boundary:                                                   |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial g_{\ell m}}{\partial r}=0,\quad  h_{\ell m}=0              |
  +-------------+---------------------------------------------------------------------------------+

* **kbotb** (default :f:var:`kbotb=1 <kbotb>`) is an integer, which corresponds to the magnetic boundary condition for :math:`r=r_i`.

  +-------------+---------------------------------------------------------------------------------+
  | ``kbotb=1`` | Insulating inner boundary:                                                      |
  |             |    .. math::                                                                    |
  |             |     \frac{\partial g_{\ell m}}{\partial r}-\frac{\ell+1}{r}\,g_{\ell m}=0,\quad |
  |             |       \frac{\partial h_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``kbotb=2`` | Perfectly-conducting inner core:                                                |
  |             |    .. math::                                                                    |
  |             |       g_{\ell m} = \frac{\partial^2 g_{\ell m}}{\partial r^2}=0,\quad           |
  |             |       \frac{\partial h_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``kbotb=3`` | Finitely conducting inner core                                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``kbotb=4`` | Pseudo-vacuum outer boundary:                                                   |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial g_{\ell m}}{\partial r}=0,\quad  h_{\ell m}=0              |
  +-------------+---------------------------------------------------------------------------------+

.. _secPressureBcs:

Boundary condition for spherically-symmetric pressure
+++++++++++++++++++++++++++++++++++++++++++++++++++++

* **ktopp** (default :f:var:`ktopp=1 <ktopp>`) is an integer, which corresponds to the boundary condition for the spherically-symmetric pressure at :math:`r=r_o`.

  +-------------+--------------------------------------------------------------------+
  | ``ktopp=1`` | The integral of the spherically-symmetric density perturbation     |
  |             | vanishes.                                                          |
  +-------------+--------------------------------------------------------------------+
  | ``ktopp=2`` | The spherically-symmetric pressure fluctuation vanishes at the     |
  |             | outer boundary.                                                    |
  +-------------+--------------------------------------------------------------------+
