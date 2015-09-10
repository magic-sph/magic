Control namelist
================

This namelist defines the numerical parameters of the problem plus the
variables that control and organize the run.

* **mode** (default ``mode=0``) is an integer which controls the type of calculation performed.

  +---------+--------------------------------------------------------+
  | mode=0  | Selfconsistent dynamo                                  |
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
  | mode=6  | Selfconsistent dynamo, but with no Lorentz force       |
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

* **tag** (default ``tag=default``) is a character string, used as an extension for all output files.

* **n_time_steps** (default ``n_time_steps=100``) is an integer, the number of time steps to be performed.

* **tEND** (default ``tEND=0.0``) is a real, which can be used to force the code to stop when :math:``t=tEND``. This is only used when ``t/=tEND``.

* **alpha** (default ``alpha=0.5``) is a real. This is the weight used for current time step in implicit time step.

Default scales
--------------

* **n_tScale** (default ``n_tScale=0``) is an integer, which determines the time scaling

  +-------------+---------------------------+
  | n_tScale=0  | Use viscous time scale.   |
  +-------------+---------------------------+
  | n_tScale=1  | Use magnetic time scale.  |
  +-------------+---------------------------+
  | n_tScale=2  | Use thermal time scale.   |
  +-------------+---------------------------+

* **n_lScale** (default ``n_lScale=0``) is an integer which determines the reference length scale.

  +-------------+------------------------------------------+
  | n_lScale=0  | Use outer core.                          |
  +-------------+------------------------------------------+
  | n_lScale=1  | Use total core.                          |
  +-------------+------------------------------------------+


* **enscale** (default ``enscale=1.0``) is a real. This is the scaling for energies.

Update control
--------------

* **l_update_v** (default ``l_update_v=.true.``) is a logical that specifies whether the velocity field should be time-stepped or not.

* **l_update_b** (default ``l_update_b=.true.``) is a logical that specifies whether the magnetic field should be time-stepped or not.

* **l_update_s** (default ``l_update_s=.true.``) is a logical that specifies whether the entropy/temperature should be time-stepped or not.

Time step control
-----------------

A modified courant criteria including a modified Alfven-velocity is used to
account for the magnetic field. The relative and absolute importance of flow
and Alfven-velocity can be controled by **courfac** and **alffac** respectively.

* **dtstart** (default ``dtstart=0.0``) is a real, which is used as the initial time step if the starting solution is initialized (see below) and :math:`\hbox{dtstart}>0`.

* **dtMax** (default ``dtMax=1e-4``) is a  real. This is the maximum allowed time step :math:`\delta t`. If :math:`\delta t > \hbox{dtmax}`, the time step is decreased to at least dtmax (See routine `dt_courant`). Run is stopped if :math:`\delta t < \hbox{dtmin}` and :math:`\hbox{dtmin}=10^{-6}\,\hbox{dtmax}`.

* **courfac** (default ``courfac=2.5``) is a real used to scale velocity in courant criteria.

* **alffac** (default ``alffac=1.0``) is a  real, used to scale Alfven-velocity in courant criteria.

* **n_cour_step** (default ``n_cour_step=10``) is an integer. This is the number of time steps before consecutive checking of courant criteria. Note: the courant criteria is checked always after the time step has been changed if ``n_cour_step}>0``.


Run time
--------

The total desired runtime (in human units and not in CPU units) can be specified with the three variables **runHours**, **runMinutes** and **runSeconds**.

* **runHours** (default ``runHours=0``) is an integer that controls the number of run hours. 

* **runMinutes** (default ``runMinutes=0``) is an integer that controls the .

* **runSeconds** (default ``runSeconds=0``) is an integer that controls the number of run hours.


Here is an example for a run of 23h30:

.. code:: fortran

   runHours   = 23,
   runMinutes = 30,


Hyperdiffusivity
----------------

Hyperdiffusion can be applied by multiplying the diffusion operators by a factor of the form

.. math::
   d(\ell)=1+D\left[\frac{\ell+1-\ell_{hd}}{\ell_{max}+1-\ell_{hd}} \right]^{\beta}

for the spherical harmonic degrees :math:`\ell \geq \ell_{hd}`.

* **difnu** (default ``difnu=0.0``) is a real. This is the amplitude :math:`D` of the viscous hyperdiffusion.

* **difkappa** (default ``difkappa=0.0``) is a real. This is the amplitude :math:`D` of the thermal hyperdiffusion.

* **difeta** (default ``difeta=0.0``) is a real. This is the amplitude :math:`D` of the magnetic hyperdiffusion.

* **ldif** (default ``ldif=1``) is an integer. This is the degree :math:`\ell_{hd}` where hyperdiffusion starts to act.

* **ldifexp** (default ``ldifexp=-1``) is an integer. This is the exponent :math:`\beta` of hyperdiffusion.


Angular momentum correction
---------------------------

In case of the use of stress-free boundary conditions at both boundaries, it is safer to ensure
that the angular momentum is correctly conserved. This can be enforced through the following
input variables:

* **l_correct_AMe** (default ``l_correct_AMe=.false.``) is a logical. This is used to correct the equatorial angular momentum.

* **l_correct_AMz** (default ``l_correct_AMz=.false.``) is a logical. This is used to correct the axial angular momentum.


Mapping of the Gauss-Lobatto grid
---------------------------------

* **l_newmap** (default ``l_newmap=.false.``) is a logical. A radial mapping can be applied to the Chebyshev grid.

* **alph1** (default ``alph1=2.0``) is a real. This is a control parameter of the mapping function.

* **alph2** (default ``alph2=0.0``) is a real. This is a control parameter of the mapping function.


Miscellaneous
-------------

* **l_non_rot** (default ``l_non_rot=.false.``) is a logical. Use it when you want to do non-rotating numerical simulations.

