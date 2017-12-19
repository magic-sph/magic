.. _secControlNml:

Control namelist
================

This namelist defines the numerical parameters of the problem plus the
variables that control and organize the run.

.. _varmode:

* **mode** (default :f:var:`mode=0 <mode>`) is an integer which controls the type of calculation performed.

  +---------+--------------------------------------------------------+
  | mode=0  | Self-consistent dynamo                                 |
  +---------+--------------------------------------------------------+
  | mode=1  | Convection                                             |
  +---------+--------------------------------------------------------+
  | mode=2  | Kinematic dynamo                                       |
  +---------+--------------------------------------------------------+
  | mode=3  | Magnetic decay modes                                   |
  +---------+--------------------------------------------------------+
  | mode=4  | Magneto convection                                     |
  +---------+--------------------------------------------------------+
  | mode=5  | Linear onset of convection                             |
  +---------+--------------------------------------------------------+
  | mode=6  | Self-consistent dynamo, but with no Lorentz force      |
  +---------+--------------------------------------------------------+
  | mode=7  | Super-rotating inner core or mantle, no convection and |
  |         | no magnetic field                                      |
  +---------+--------------------------------------------------------+
  | mode=8  | Super-rotating inner core or mantle, no convection     |
  +---------+--------------------------------------------------------+
  | mode=9  | Super-rotating inner core or mantle, no convection     |
  |         | and no Lorentz force                                   |
  +---------+--------------------------------------------------------+
  | mode=10 | Super-rotating inner core or mantle, no convection,    |
  |         | no magnetic field, no Lorentz force and no advection   |
  +---------+--------------------------------------------------------+

.. _varTAG:

* **tag** (default :f:var:`tag="default" <tag>`) is a character string, used as an extension for all output files.

* **n_time_steps** (default :f:var:`n_time_steps=100 <n_time_steps>`) is an integer, the number of time steps to be performed.

* **tEND** (default :f:var:`tEND=0.0 <tend>`) is a real, which can be used to force the code to stop when :math:``t=tEND``. This is only used when ``t/=tEND``.

* **alpha** (default :f:var:`alpha=0.5 <alpha>`) is a real. This is the weight used for current time step in implicit time step.

Default scales
--------------

* **n_tScale** (default :f:var:`n_tScale=0 <n_tscale>`) is an integer, which determines the time scaling

  +-------------+----------------------------+---------------------+
  | n_tScale=0  | Use viscous time scale.    | :math:`d^2/\nu`     |
  +-------------+----------------------------+---------------------+
  | n_tScale=1  | Use magnetic time scale.   | :math:`d^2/\eta`    |
  +-------------+----------------------------+---------------------+
  | n_tScale=2  | Use thermal time scale.    | :math:`d^2/\kappa`  |
  +-------------+----------------------------+---------------------+
  | n_tScale=3  | Use rotational time scale. | :math:`\Omega^{-1}` |
  +-------------+----------------------------+---------------------+

* **n_lScale** (default :f:var:`n_lScale=0 <n_lscale>`) is an integer which determines the reference length scale.

  +-------------+------------------------------------------+
  | n_lScale=0  | Use outer core.                          |
  +-------------+------------------------------------------+
  | n_lScale=1  | Use total core.                          |
  +-------------+------------------------------------------+


* **enscale** (default :f:var:`enscale=1.0 <enscale>`) is a real. This is the scaling for energies.

Update control
--------------

* **l_update_v** (default :f:var:`l_update_v=.true. <l_update_v>`) is a logical that specifies whether the velocity field should be time-stepped or not.

* **l_update_b** (default :f:var:`l_update_b=.true. <l_update_b>`) is a logical that specifies whether the magnetic field should be time-stepped or not.

* **l_update_s** (default :f:var:`l_update_s=.true. <l_update_s>`) is a logical that specifies whether the entropy/temperature should be time-stepped or not.

Time step control
-----------------

A modified courant criteria including a modified Alfven-velocity is used to
account for the magnetic field. The relative and absolute importance of flow
and Alfven-velocity can be controled by **courfac** and **alffac** respectively.
The parameter **l_cour_alf_damp** allows to choose wheter the actual Alven speed
is used to estimate the Courant condition or if damping (see Christensen et al.,
GJI, 1999) is included.

* **dtstart** (default :f:var:`dtstart=0.0 <dtstart>`) is a real, which is used as the initial time step if the starting solution is initialized (see below) and :math:`\hbox{dtstart}>0`.

* **dtMax** (default :f:var:`dtMax=1e-4 <dtmax>`) is a  real. This is the maximum allowed time step :math:`\delta t`. If :math:`\delta t > \hbox{dtmax}`, the time step is decreased to at least dtmax (See routine `dt_courant`). Run is stopped if :math:`\delta t < \hbox{dtmin}` and :math:`\hbox{dtmin}=10^{-6}\,\hbox{dtmax}`.

* **courfac** (default :f:var:`courfac=2.5 <courfac>`) is a real used to scale velocity in courant criteria.

* **alffac** (default :f:var:`alffac=1.0 <alffac>`) is a  real, used to scale Alfven-velocity in courant criteria.

* **l_cour_alf_damp** (default :f:var:`l_cour_alf_damp=.true. <l_cour_alf_damp>`) is a logical. This is used to decide whether the damping of the Alven waves is taken into account when estimating the Courant condition (see Christensen et al., GJI, 1999). At low Ekman numbers, this criterion might actually lead to spurious oscillations/instabilities of the code.

* **n_cour_step** (default :f:var:`n_cour_step=10 <n_cour_step>`) is an integer. This is the number of time steps before consecutive checking of courant criteria. Note: the courant criteria is checked always after the time step has been changed if ``n_cour_step>0``.


Run time
--------

The total desired runtime (in human units and not in CPU units) can be specified with the three variables **runHours**, **runMinutes** and **runSeconds**.

* **runHours** (default :f:var:`runHours=0 <runhours>`) is an integer that controls the number of run hours. 

* **runMinutes** (default :f:var:`runMinutes=0 <runminutes>`) is an integer that controls the .

* **runSeconds** (default :f:var:`runSeconds=0 <runseconds>`) is an integer that controls the number of run hours.


Here is an example for a run of 23h30:

.. code-block:: fortran

   runHours   = 23,
   runMinutes = 30,


Hyperdiffusivity
----------------

Hyperdiffusion can be applied by multiplying the diffusion operators by a factor of the form

.. math::
   d(\ell)=1+D\left[\frac{\ell+1-\ell_{hd}}{\ell_{max}+1-\ell_{hd}} \right]^{\beta}

for the spherical harmonic degrees :math:`\ell \geq \ell_{hd}`.

* **difnu** (default :f:var:`difnu=0.0 <difnu>`) is a real. This is the amplitude :math:`D` of the viscous hyperdiffusion.

* **difkappa** (default :f:var:`difkappa=0.0 <difkappa>`) is a real. This is the amplitude :math:`D` of the thermal hyperdiffusion.

* **difeta** (default :f:var:`difeta=0.0 <difeta>`) is a real. This is the amplitude :math:`D` of the magnetic hyperdiffusion.

* **ldif** (default :f:var:`ldif=1 <ldif>`) is an integer. This is the degree :math:`\ell_{hd}` where hyperdiffusion starts to act.

* **ldifexp** (default :f:var:`ldifexp=-1 <ldifexp>`) is an integer. This is the exponent :math:`\beta` of hyperdiffusion.


Angular momentum correction
---------------------------

In case of the use of stress-free boundary conditions at both boundaries, it is safer to ensure
that the angular momentum is correctly conserved. This can be enforced through the following
input variables:

* **l_correct_AMe** (default :f:var:`l_correct_AMe=.false. <l_correct_ame>`) is a logical. This is used to correct the equatorial angular momentum.

* **l_correct_AMz** (default :f:var:`l_correct_AMz=.false. <l_correct_amz>`) is a logical. This is used to correct the axial angular momentum.


.. _varl_newmap:

Radial scheme and mapping of the Gauss-Lobatto grid
---------------------------------------------------

In MagIC, one can either use finite differences or Chebyshev polynomials for the radial integration scheme. This choice is controlled by the following input parameter:

* **radial_scheme** (default :f:var:`radial_scheme='CHEB' <radial_scheme>`) is a character string.

  +-----------------------+--------------------------------+
  | radial_scheme='CHEB'  | Use Chebyshev polynomials      |
  +-----------------------+--------------------------------+
  | radial_scheme='FD'    | Use finite differences         |
  +-----------------------+--------------------------------+

When Chebyshev polynomials are used, it is also possible to use a non-linear
mapping function to concentrate/diperse grid points around a point inside the
domain. 


* **l_newmap** (default :f:var:`l_newmap=.false. <l_newmap>`) is a logical. A radial mapping can be applied to the Chebyshev grid when ``l_newmap`` is set to ``.true.``. The radial profile of the mapping function is then stored during the initialisation of the code in the file :ref:`rNM.TAG <secMappingFile>`.

* **map_function** (default :f:var:`map_function='arcsin' <map_function>`) is a character string. This allows to select which mapping function is used:

  +-----------------------+-----------------------------------------------------------------------------------------------------------+
  | map_function='TAN'    | Use a tangent mapping  (see `Bayliss and Turkel 1992 <https://doi.org/10.1016/0021-9991(92)90012-N>`_)    |
  +-----------------------+-----------------------------------------------------------------------------------------------------------+
  | map_function='ARCSIN' | Use finite differences (see `Kosloff and Tal-Ezer 1993 <https://doi.org/10.1006/jcph.1993.1044>`_)        |
  +-----------------------+-----------------------------------------------------------------------------------------------------------+

If the tangent mapping is used, the function that re-distributes the collocation 
points is expressed by

.. math::
   r=\frac{1}{2}\left(\alpha_2+\frac{\textrm{tan}\left[\lambda(r_{cheb}-x_0)\right]}{\alpha_1}\right) + \frac{r_i+r_o}{2} \textrm{ ,}

where the Gauss-Lobatto collocation points are

.. math::
   r_{cheb}&=\textrm{cos}\left( \frac{\pi(k-1)}{N_r} \right) \textrm{ , }\;\; k=1,2,...,n_r \textrm{ , }\; n_r=n\_r\_max

and :math:`r\!\in\![r_i,r_o]`, :math:`r_{cheb}\!\in\![-1.0,1.0]`. The parameters to calculate :math:`r` are

.. math::
   \lambda&=\frac{\textrm{tan}^{-1}\left(\alpha_1(1-\alpha_2)\right)}{1-x_0} \\
   x_0&=\frac{K-1}{K+1} \\
   K&=\frac{\textrm{tan}^{-1}\left(\alpha_1(1+\alpha_2)\right)}{\textrm{tan}^{-1}\left(\alpha_1(1-\alpha_2)\right)} \textrm{ .}

The coefficient :math:`\alpha_1` determines the degree of concentration/dispersion of the grid points around :math:`r_{cheb}\!=\!\alpha_2`. If :math:`\alpha_1` is too high, the :math:`r` function becomes nearly discontinuous. To avoid numerical problems, :math:`\alpha_1` should remain close to unity.

If the arcsin mapping is used, the function that re-distributes the collocation points
is given by

.. math::
   r=\frac{1}{2}\left[ \frac{\textrm{arcin}\left(\alpha_1 r_{cheb}\right)}{\textrm{arcsin} \alpha_1} \right]+\frac{r_i+r_o}{2} \textrm{ ,}

In the Kosloff and Tal-Ezer mapping, :math:`\alpha_1` transforms the Gauss-Lobatto
grid into a more regularly-spaced grid. When :math:`\alpha_1 \rightarrow 0` one 
recovers the Gauss-Lobatto grid, while :math:`\alpha_1 \rightarrow 1` yields a
regular grid. 

.. warning:: The Kosloff-Tal-Ezer mapping becomes singular when :math:`\alpha_1=1`.
             Acceptable values are :math:`0<\alpha_1<1`. Note that the error increases
	     as :math:`\epsilon=\left(\frac{1-\sqrt{1-\alpha_1^2}}{\alpha_1}\right)^{N_r}`.

..


* **alph1** (default :f:var:`alph1=0.8 <alph1>`) is a real. This is a control parameter of the mapping function.

* **alph2** (default :f:var:`alph2=0.0 <alph2>`) is a real. This is a control parameter of the mapping function.


Miscellaneous
-------------

* **l_non_rot** (default :f:var:`l_non_rot=.false. <l_non_rot>`) is a logical. Use it when you want to do non-rotating numerical simulations.

* **anelastic_flavour** (default :f:var:`anelastic_flavour="None" <anelastic_flavour>`) is a character string. This allows to change the thermal diffusion operator used within the anelastic approximation. Possible values are:

   +---------------------------+------------------------------------+
   | anelastic_flavour='LBR'   | Entropy diffusion                  |
   +---------------------------+------------------------------------+
   | anelastic_flavour='ENT'   | Entropy diffusion                  |
   +---------------------------+------------------------------------+
   | anelastic_flavour='ALA'   | Anelastic liquid approximation     |
   +---------------------------+------------------------------------+
   | anelastic_flavour='TDIFF' | Temperature diffusion              |
   +---------------------------+------------------------------------+
   | anelastic_flavour='TEMP'  | Temperature diffusion              |
   +---------------------------+------------------------------------+

* **thermo_variable** (default :f:var:`thermo_variable="S" <thermo_variable>`) is a character string. This allows to change the default thermodynamic variable (and hence change the entropy/temperature equation used). This switch only matters when one wants to run an anelastic model. Possible values are:

   +-----------------------+-----------------------------------------+
   | thermo_variable='S'   | Use entropy as a primitive variable     |
   +-----------------------+-----------------------------------------+
   | thermo_variable='T'   | Use temperature as a primitive variable |
   +-----------------------+-----------------------------------------+

* **polo_flow_eq** (default :f:var:`polo_flow_eq="WP" <polo_flow_eq>`) is a character string. This allows to change how the equation for the poloidal flow potential is constructed. One can either use the radial component of the Navier-Stokes equation and hence keep a coupled system that involve the poloidal potential :math:`W` and the pressure :math:`p`, or take the radial component of the double-curl of the Navier-Stokes equation to suppress pressure.

   +---------------------+-----------------------------------------+
   | polo_flow_eq='WP'   | Use the pressure formulation            |
   +---------------------+-----------------------------------------+
   | polo_flow_eq='DC'   | Use the double-curl formulation         |
   +---------------------+-----------------------------------------+

