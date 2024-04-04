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

* **l_update_xi** (default :f:var:`l_update_xi=.true. <l_update_xi>`) is a logical that specifies whether the chemical composition should be time-stepped or not.

* **l_update_phi** (default :f:var:`l_update_phi=.true. <l_update_phi>`) is a logical that specifies whether the phase field should be time-stepped or not.


Time step control
-----------------

A modified Courant criterion including a modified Alfven-velocity is used to
account for the magnetic field. The relative and absolute importance of flow
and Alfven-velocity can be controled by **courfac** and **alffac** respectively.
The parameter **l_cour_alf_damp** allows to choose whether the actual Alven speed
is used to estimate the Courant condition or if damping is included. Practically,
the timestep size is controlled as follows

.. math::
   \delta t < \min_{V}\left( c_I\,E,\, \dfrac{\delta r}{|u_r|},\, \dfrac{\delta h}{u_h} \right)

where :math:`u_h=(u_\theta^2+u_\phi^2)^{1/2}`, :math:`\delta h = \dfrac{r}{\sqrt{\ell(\ell+1)}}`, and :math:`\delta r` is the radial grid interval. The first term in the left hand side accounts for the explicit treatment of the Coriolis term.

.. math::
   {|u_r|}=c_F{|u_{F,r}|}+c_A\dfrac{u_{A,r}^2}{\left[u_{A,r}^2+\left(\frac{1+Pm^{-1}}{2\delta r}\right)^2\right]^{1/2}}\,,

where :math:`u_{F,r}` is the radial component of the fluid velocity and :math:`u_{A,r}=Br/\sqrt{E\,Pm}` is the radial Alven velocity. The denominator of the rightmost term accounts for the damping of the Alven waves.

* **dtMax** (default :f:var:`dtMax=1e-4 <dtmax>`) is a  real. This is the maximum allowed time step :math:`\delta t`. If :math:`\delta t > \hbox{dtmax}`, the time step is decreased to at least dtMax (See routine `dt_courant`). Run is stopped if :math:`\delta t < \hbox{dtmin}` and :math:`\hbox{dtmin}=10^{-6}\,\hbox{dtmax}`.

* **courfac** (default :f:var:`courfac=2.5 <courfac>`) is a real used to scale velocity in Courant criteria. This parameter corresponds to :math:`c_F` in the above equation.

* **alffac** (default :f:var:`alffac=1.0 <alffac>`) is a  real, used to scale Alfven-velocity in Courant criteria. This parameter corresponds to :math:`c_A` in the above equation.

* **intfac** (default :f:var:`intfac=0.15 <intfac>`) is a  real, used to scale Coriolis factor in Courant criteria. This parameter corresponds to :math:`c_I` in the above equation.

* **l_cour_alf_damp** (default :f:var:`l_cour_alf_damp=.true. <l_cour_alf_damp>`) is a logical. This is used to decide whether the damping of the Alven waves is taken into account when estimating the Courant condition (see Christensen et al., GJI, 1999). At low Ekman numbers, this criterion might actually lead to spurious oscillations/instabilities of the code. When turn to False, :math:`{|u_r|}=c_F{|u_{F,r}|}+c_A{|u_{A,r}|}`.

* **time_scheme** (default :f:var:`time_scheme='CNAB2' <time_scheme>`) is a character string. This is used to choose the time step integrator used in the code among the following implicit-explicit time schemes:

  +-----------------------+-------------------------------------------------------+
  | time_scheme='CNAB2'   | Crank-Nicolson and 2nd order Adams-Bashforth scheme   |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='CNLF'    | Crank-Nicolson and Leap-Frog scheme                   |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='MODCNAB' | Modified CN/AB2                                       |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='SBDF2'   | Semi-implicit backward difference scheme of 2nd order |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='SBDF3'   | Semi-implicit backward difference scheme of 3rd order |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='SBDF4'   | Semi-implicit backward difference scheme of 4th order |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='ARS222'  | Semi-implicit S-DIRK of 2nd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='ARS232'  | Semi-implicit S-DIRK of 2nd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='CK232'   | Semi-implicit S-DIRK of 2nd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='LZ232'   | Semi-implicit S-DIRK of 2nd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='PC2'     | Semi-implicit S-DIRK of 2nd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='CB3'     | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='ARS343'  | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='MARS343' | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='ARS443'  | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='BPR353'  | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='BHR553'  | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='DBM453'  | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='LZ453'   | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='KC343'   | Semi-implicit S-DIRK of 3rd order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='KC564'   | Semi-implicit S-DIRK of 4th order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='KC674'   | Semi-implicit S-DIRK of 4th order                     |
  +-----------------------+-------------------------------------------------------+
  | time_scheme='KC785'   | Semi-implicit S-DIRK of 5th order                     |
  +-----------------------+-------------------------------------------------------+


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

* **difchem** (default :f:var:`difchem=0.0 <difchem>`) is a real. This is the amplitude :math:`D` of the hyperdiffusion applied to chemical composition.

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
  | map_function='ARCSIN' | Use an arcsin mapping  (see `Kosloff and Tal-Ezer 1993 <https://doi.org/10.1006/jcph.1993.1044>`_)        |
  +-----------------------+-----------------------------------------------------------------------------------------------------------+
  | map_function='TT'     | Use the mapping by `Tee and Trefethen 2006 <https://doi.org/10.1137/050641296>`_                          |
  +-----------------------+-----------------------------------------------------------------------------------------------------------+
  | map_function='JAFARI' | Use the mapping by `Jafari-Varzaneh and Hosseini 2014 <https://doi.org/10.1007/s11075-014-9883-3>`_       |
  +-----------------------+-----------------------------------------------------------------------------------------------------------+

If the tangent mapping is used, the function that re-distributes the collocation 
points is expressed by

.. math::
   r=\frac{1}{2}\left(\alpha_2+\frac{\textrm{tan}\left[\lambda(x_{cheb}-x_0)\right]}{\alpha_1}\right) + \frac{r_i+r_o}{2} \textrm{ ,}

where the Gauss-Lobatto collocation points are

.. math::
   x_{cheb}=\textrm{cos}\left( \frac{\pi(k-1)}{N_r} \right) \textrm{ , }\;\; k=1,2,...,n_r \textrm{ , }\; n_r=n\_r\_max

and :math:`r\!\in\![r_i,r_o]`, :math:`x_{cheb}\!\in\![-1.0,1.0]`. The parameters to calculate :math:`r` are

.. math::
   \lambda&=\frac{\textrm{tan}^{-1}\left(\alpha_1(1-\alpha_2)\right)}{1-x_0} \\
   x_0&=\frac{K-1}{K+1} \\
   K&=\frac{\textrm{tan}^{-1}\left(\alpha_1(1+\alpha_2)\right)}{\textrm{tan}^{-1}\left(\alpha_1(1-\alpha_2)\right)} \textrm{ .}

The coefficient :math:`\alpha_1` determines the degree of concentration/dispersion of the grid points around :math:`x_{cheb}\!=\!\alpha_2`. If :math:`\alpha_1` is too high, the :math:`r` function becomes nearly discontinuous. To avoid numerical problems, :math:`\alpha_1` should remain close to unity.

If the arcsin mapping is used, the function that re-distributes the collocation points
is given by

.. math::
   r=\frac{1}{2}\left[ \frac{\textrm{arcin}\left(\alpha_1 x_{cheb}\right)}{\textrm{arcsin} \alpha_1} \right]+\frac{r_i+r_o}{2} \textrm{ ,}

In the Kosloff and Tal-Ezer mapping, :math:`\alpha_1` transforms the Gauss-Lobatto
grid into a more regularly-spaced grid. When :math:`\alpha_1 \rightarrow 0` one 
recovers the Gauss-Lobatto grid, while :math:`\alpha_1 \rightarrow 1` yields a
regular grid. 

.. warning:: The Kosloff-Tal-Ezer mapping becomes singular when :math:`\alpha_1=1`.
             Acceptable values are :math:`0<\alpha_1<1`. Note that the error increases
	     as :math:`\epsilon=\left(\frac{1-\sqrt{1-\alpha_1^2}}{\alpha_1}\right)^{N_r}`.

..

If the Tee and Trefethen sinh mapping is employed, the grid points are redistributed in the following manner

.. math::
   r=\frac{1}{2}\left(\alpha_2+\frac{\textrm{sinh}\left[A(x_{cheb}-1)+B\right]}{\alpha_1}\right) + \frac{r_i+r_o}{2} \textrm{ ,}

where

.. math::
   A=\frac{1}{2}\left[\textrm{sinh}(\alpha_1(1-\alpha_2))+\textrm{sinh}(\alpha_1(1+\alpha_2)) \right], \quad B = \textrm{sinh}(\alpha_1(1-\alpha_2))

With this mapping, :math:`\alpha_1` is directly related to the stiffness of the transition.


If the Jafari-Varzaneh and Hosseini mapping is employed, similarly to the tangent mapping, :math:`\alpha_1` determines the degree of concentration of the grid points around :math:`x_{cheb}\!=\!\alpha_2`. This is expected to do a better job than the tangent mapping, both in terms of matrix conditioning and in terms of reducing the Gibbs phenomenon around a steep change (Allen-Cahn type of equations involved in the phase field model comes to mind).


* **alph1** (default :f:var:`alph1=0.8 <alph1>`) is a real. This is a control parameter of the mapping function.

* **alph2** (default :f:var:`alph2=0.0 <alph2>`) is a real. This is a control parameter of the mapping function. The default value of :math:`0` corresponds to the center of the grid.


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

* **polo_flow_eq** (default :f:var:`polo_flow_eq="WP" <polo_flow_eq>`) is a character string. This allows to change how the equation for the poloidal flow potential is constructed. One can either use the radial component of the Navier-Stokes equation and hence keep a coupled system that involve the poloidal potential :math:`W` and the pressure :math:`p`, or take the radial component of the double-curl of the Navier-Stokes equation to suppress pressure.

   +---------------------+-----------------------------------------+
   | polo_flow_eq='WP'   | Use the pressure formulation            |
   +---------------------+-----------------------------------------+
   | polo_flow_eq='DC'   | Use the double-curl formulation         |
   +---------------------+-----------------------------------------+

* **mpi_transp** (default :f:var:`mpi_transp="auto" <mpi_tansp>`) is a character string. It allows to change the way the global MPI transposes are handled by the code. By default, the code tries to determine by itself the fastest method. One can nevertheless force the code to use local communicators (such as Isend/Irecv/waitall), make use of the native alltoallv MPI variant or choose the alltoallw variant instead.

   +--------------------+--------------------------------------------------+
   | mpi_transp='auto'  | Automatic determination of the fastest transpose |
   +--------------------+--------------------------------------------------+
   | mpi_transp='p2p'   | Use Isend/Irecv/Waitall communicators            |
   +--------------------+--------------------------------------------------+
   | mpi_transp='a2av'  | Use alltoallv communicators                      |
   +--------------------+--------------------------------------------------+
   | mpi_transp='a2aw'  | Use alltoallw communicators                      |
   +--------------------+--------------------------------------------------+

* **mpi_packing** (default :f:var:`mpi_packing="packed" <mpi_packing>`) is a character string. It allows to change the size of the global MPI transposes. One can choose between some packing of the fields into buffers (default) or a sequence of single field transposes. There is a possible automatic detection but testing unfortunately reveals frequent false detection.

   +------------------------+--------------------------------------------------+
   | mpi_packing='auto'     | Automatic determination of the fastest transpose |
   +------------------------+--------------------------------------------------+
   | mpi_packing='packed'   | Pack some fields into buffers                    |
   +------------------------+--------------------------------------------------+
   | mpi_packing='single'   | Transpose each field individually                |
   +------------------------+--------------------------------------------------+

* **l_adv_curl** (default :f:var:`l_adv_curl=.false. <l_adv_curl>`) is a logical. When set to True, the advection term is treated as :math:`\vec{u}\times\vec{\omega}` instead of :math:`\vec{u}\vec{\nabla}\vec{u}`. The practical consequence of that is to reduce the number of spectral/spatial Spherical Harmonic Transforms and hence to speed-up the code. Because of the treatment of the viscous heating term in the anelastic approximation, this is only an option when considering Boussinesq models.
