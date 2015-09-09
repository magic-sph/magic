Physical parameters namelist
============================

This namelist contains all the appropriate relevant control physical parameters.

Dimensionless control parameters
--------------------------------

* **ra** (default ``ra=1.1e5``) is a real. This the Rayleigh number expressed by
  
  .. math::
     Ra = \frac{\alpha g_o \Delta T d^3}{\kappa\nu}

* **ek** (default ``ek=1e-3``) is a real. This is the Ekman number expressed by

  .. math::
     E = \frac{\nu}{\Omega d^2}

* **pr** (default ``pr=1.0``) is a real. This is the Prandtl number expressed by

  .. math::
     Pr = \frac{\nu}{\kappa}

* **prmag** (default ``prmag=5.0``) is a real. This is the magnetic Prandtl number expressed by

  .. math::
     Pm = \frac{\nu}{\lambda}

* **radratio** (default ``radratio=0.35``) is a real. This is the ratio of the inner core radius :math:`r_i` to the outer core radius :math:`r_o`:

  .. math::
     \eta = \frac{r_i}{r_o}

* **strat** (default ``strat=0.0``) is a real. This is the number of density scale heights of the reference state:

  .. math::
     N_\rho = \ln \frac{\tilde{\rho}(r_i)}{\tilde{\rho}(r_o)}

* **polind** (default ``polind=1.5``) is a real. This is the polytropic index, which relates the background temperature to the background density:

  .. math::
     \tilde{\rho} = \tilde{T}^m

Gravity
-------

The radial dependence of the gravity profile can be adjusted following

.. math::
   g(r)=g_0+g1\,\frac{r}{r_o}+g2\,\left(\frac{r_o}{r}\right)^2

The three following parameters are used to set this profile

* **g0** (default ``g0=0``) is the pre-factor of the constant part of the gravity profile.

* **g1** (default ``g1=1``) is the pre-factor of the linear part of the gravity profile.

* **g2** (default ``g2=0``) is the pre-factor of the :math:`1/r^2` part of the gravity profile.
     
Transport properties
--------------------

* **difExp** (default ``difExp=-0.5``) is a real. This is the exponent that is used when ``nVarVisc=2``, ``nVarDiff=2`` or ``nVarCond=4``.


Electrical conductivity
+++++++++++++++++++++++

* **nVarCond** (default ``nVarCond=0``) is an integer. This is used to modify the radial-dependence of the electrical conductivity.

Thermal diffusivity
+++++++++++++++++++

* **nVarDiff** (default ``nVarCond=0``) is an integer. This is used to change the radial-dependence of the thermal diffusivity:

  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=0`` | Constant thermal diffusivity :math:`\kappa`                                |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=1`` | Constant thermal conductivity, i.e.                                        |
  |                | :math:`\kappa =\frac{\tilde{\rho}_i}{\tilde{\rho}(r)}`                     |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=2`` | Radial profile of the form                                                 |
  |                | :math:`\kappa=\left(\frac{\tilde{\rho}(r)}{\tilde{\rho}_i}\right)^{\alpha}`|
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=3`` | polynomial-fit to an interior model of Jupiter                             |
  +----------------+----------------------------------------------------------------------------+
  | ``nVarDiff=4`` | polynomial-fit to an interior model of the Earth liquid core               |
  +----------------+----------------------------------------------------------------------------+

Viscosity
+++++++++

* **nVarVisc** (default ``nVarVisc=0``) is an integer. This is used to change the radial-dependence of the viscosity:

  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=0`` | Constant kinematic viscosity :math:`\nu`                                |
  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=1`` | Constant dynamic viscosity, i.e.                                        |
  |                | :math:`\nu =\frac{\tilde{\rho}_o}{\tilde{\rho}(r)}`                     |
  +----------------+-------------------------------------------------------------------------+
  | ``nVarVisc=2`` | Radial profile of the form                                              |
  |                | :math:`\nu=\left(\frac{\tilde{\rho}(r)}{\tilde{\rho}_i}\right)^{\alpha}`|
  +----------------+-------------------------------------------------------------------------+

  where :math:`alpha` is an exponent set by the namelist input variable ``difExp``.


Boundary conditions
-------------------

Thermal boundary conditions
+++++++++++++++++++++++++++

* **ktops** (default ``ktops=1``) is an  integer to specify the outer boundary entropy (or temperature) boundary condition:

  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=1`` | Fixed entropy at outer boundary: :math:`s(r_o)=s_{top}`                             |
  +-------------+-------------------------------------------------------------------------------------+
  | ``ktops=2`` | Fixed entropy flux at outer boundary: :math:`\partial s(r_o)/\partial r = s_{top}`  |
  +-------------+-------------------------------------------------------------------------------------+

* **kbots** (default ``ktops=1``) is an  integer to specify the inner boundary entropy (or temperature) boundary condition.

* **s_top** (default ``s_top= 0 0 0.0 0.0``) is a real array of lateraly varying outer heat boundary conditions. Each four consecutive numbers are interpreted as follows:

  1. Spherical harmonic degree :math:`\ell`

  2. Spherical harmonic order :math:`m`

  3. Real amplitude (:math:`\cos` contribution)

  4. Imaginary amplitude (:math:`\sin` contribution)

  For example, if the boundary condition should be a combination of an :math:`(\ell=1,m=0)` sherical harmonic with the amplitude 1 and an :math:`(\ell=2,m=1)` spherical harmonic with the amplitude (0.5,0.5) the respective namelist entry could read: 
  
  
  .. code:: fortran
   
     s_top = 1, 0, 1.0, 0.0, 2, 1, 0.5, 0.5, !The comas could be left away.

* **s_bot** (default ``s_bot=0 0 0.0 0.0``) is a real array. This is the same as ``s_top`` but for the bottom boundary.

* **impS** (default ``impS=0``) is an integer. This is a  flag to indicate if there is a localized entropy disturbance, imposed at the CMB. The number of these input boundary conditions is stored in ``n_impS`` (the maximum allowed is 20), and it's given by the number of ``sCMB`` defined in the same namelist. The default value of ``impS`` is zero (no entropy disturbance). If it is set in the namelist for an integer greater than zero, then ``sCMB`` has to be also defined in the namelist, as shown below.

* **sCMB** (default ``sCMB=0.0 0.0 0.0 0.0``) is a real array of CMB heat boundary conditions (similar to the case of ``s_bot`` and ``s_top``). Each four consecutive numbers are interpreted as follows:

  1. Highest amplitude value of the entropy boundary condition, stored in array ``peakS(20)``. When ``impS<0``, ``peakS`` is a relative amplitude in comparison to the :math:`(\ell=0,m=0)` contribution (for example, the case ``s_top= 0 0 -1 0``).

  2. :math:`\theta` coordinate (input has to be given in degrees), stored in array ``thetaS(20)``.

  3. :math:`\phi` coordinate (input has to be given in degrees), stored in array ``phiS(20)``.

  4. Angular width (input has to be given in degrees), stored in array ``widthS(20)``.


Mechanical boundary conditions
++++++++++++++++++++++++++++++

* **ktopv** (default ``ktopv=2``) is an integer, which corresponds to the mechanical boundary condition for :math:`r=r_o`.

  +-------------+--------------------------------------------------------------------+
  | ``ktopv=1`` | Stress-free outer boundary for :math:`r=r_o`:                      |
  |             |   .. math::                                                        |
  |             |      w_{\ell m}(r=r_o)=0, \quad                                    |
  |             |      \frac{\partial}{\partial r}\left(\frac{1}{r^2\tilde{\rho}}    |
  |             |      \frac{\partial w_{\ell m}}{\partial r}\right)=0 \\            |
  |             |      \frac{\partial}{\partial r}\left(\frac{1}{r^2\tilde{\rho}}    |
  |             |       z_{\ell m}\right)=0                                          |
  +-------------+--------------------------------------------------------------------+
  | ``ktopv=2`` | Rigid outer boundary for :math:`r=r_o`:                            |
  |             |    .. math::                                                       |
  |             |       w_{\ell m}=0,\quad                                           |
  |             |       \frac{\partial w_{\ell m}}{\partial r}=0, \\                 |
  |             |       z_{\ell m}=0                                                 |
  +-------------+--------------------------------------------------------------------+


* **kbotv** (default ``kbotv=2``) is an integer, which corresponds to the mechanical boundary condition for :math:`r=r_i`.

Magnetic boundary conditions
++++++++++++++++++++++++++++


* **ktopb** (default ``ktopb=1``) is an integer, which corresponds to the magnetic boundary condition for :math:`r=r_o`.

  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=1`` | Insulating outer boundary:                                                      |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial b_{\ell m}}{\partial r}+\frac{\ell}{r}\,b_{\ell m}=0,\quad |
  |             |       \frac{\partial j_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=3`` | Finitely conducting mantle                                                      |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=4`` | Pseudo-vacuum outer boundary:                                                   |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial b_{\ell m}}{\partial r}=0,\quad  j_{\ell m}=0              |
  +-------------+---------------------------------------------------------------------------------+

* **kbotb** (default ``kbotb=1``) is an integer, which corresponds to the magnetic boundary condition for :math:`r=r_i`.

  +-------------+---------------------------------------------------------------------------------+
  | ``kbotb=1`` | Insulating inner boundary:                                                      |
  |             |    .. math::                                                                    |
  |             |     \frac{\partial b_{\ell m}}{\partial r}-\frac{\ell+1}{r}\,b_{\ell m}=0,\quad |
  |             |       \frac{\partial j_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=2`` | Perfectly-conducting innner core:                                               |
  |             |    .. math::                                                                    |
  |             |       b_{\ell m} = \frac{\partial b_{\ell m}}{\partial r}=0,\quad               |
  |             |       \frac{\partial j_{\ell m}}{\partial r}=0                                  |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=3`` | Finitely conducting innner core                                                 |
  +-------------+---------------------------------------------------------------------------------+
  | ``ktopb=4`` | Pseudo-vacuum outer boundary:                                                   |
  |             |    .. math::                                                                    |
  |             |       \frac{\partial b_{\ell m}}{\partial r}=0,\quad  j_{\ell m}=0              |
  +-------------+---------------------------------------------------------------------------------+
