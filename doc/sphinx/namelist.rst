Input parameters
################

True runtime input parameters are read from STDIN as namelists, a Fortran
feature. A namelist is identified by its unique name `&name`. The
name-statement is followed by the parameters that are part of the namelist in
the format `parameter=value,`. The namelist is closed by a backslash. The
subroutine `defaultNamelists` (in the module `Namelist.f90`) defines a default
value for each parameter. Only the parameters whose value should differ from
its default have to be stated in the namelist.

An example for the short namelist defining inner core parameters is::

   &inner_core
   sigma_ratio = 1.0,
   nRotIc = 1

Comas can be used to seperate namelist entries since they are not interpreted by the code.

Magic uses six namelists:

1. `&grid` for resolution
2. `&control` for control parameters and numerical parameters.
3. `&phys_param` for the physical parameters.
4. `&start_field` to define the starting fields.
5. `&output_control` for defining the output.
6. `&mantle` for setting mantle parameters.
7. `&inner_core` for setting inner core parameters.
8. `&B_external` for setting up an external field contribution

Grid namelist
=============

Control namelist
================

This namelist defines the numerical parameters of the problem plus the
variables that control and organize the run. An exception are the grid-size and
truncation parameters, that have to be set before compiling the code.

* **mode** is an integer which controls the type of calculation performed.

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

* **tag** is a character string, used as an extension for all output files.

* **n_time_steps** is an integer, the number of time steps to be performed.

* **n_tScale** is an integer, which determines the time scaling

  +-------------+---------------------------+
  | n_tScale=0  | Use viscous time scale.   |
  +-------------+---------------------------+
  | n_tScale=1  | Use magnetic time scale.  |
  +-------------+---------------------------+
  | n_tScale=2  | Use thermal time scale.   |
  +-------------+---------------------------+

  The default is n_tScale=0

* **n_lScale** is an integer which determines the reference length scale.

  +-------------+------------------------------------------+
  | n_lScale=0  | Use outer core.                          |
  +-------------+------------------------------------------+
  | n_lScale=1  | Use total core.                          |
  +-------------+------------------------------------------+

* **alpha** is a real. This is the weight used for current time step in implicit time step.

* **enscale** is a real. This is the scaling for energies.

* **l_update_v** is a logical that specifies whether the velocity field should be time-stepped or not.

* **l_update_b** is a logical that specifies whether the magnetic field should be time-stepped or not.

* **l_update_s** is a logical that specifies whether the entropy/temperature should be time-stepped or not.

Time step control
-----------------

A modified courant criteria including a modified Alfven-velocity is used to
account for the magnetic field. The relative and absolute importance of flow
and Alfven-velocity can be controled by **courfac** and **alffac** respectively.

* **dtstart** is a real, which is used as the initial time step if the starting solution is initialized (see below) and :math:`\hbox{dtstart}>0`.

* **dtmax** is a  real. This is the maximum allowed time step :math:`\delta t`. If :math:`\delta t > \hbox{dtmax}`, the time step is decreased to at least dtmax (See routine `dt_courant`). Run is stopped if :math:`\delta t < \hbox{dtmin}` and :math:`\hbox{dtmin}=10^{-6}\,\hbox{dtmax}`.

* **courfac** is a real used to scale velocity in courant criteria.

* **alffac** is a  real, used to scale Alfven-velocity in courant criteria.

* **n_cour_step** is an integer. This is the number of time steps before consecutive checking of courant criteria. Note: the courant criteria is checked always after the time step has been changed if :math:`\hbox{n\_cour\_step}>0`.

Hyperdiffusivity
----------------

Hyperdiffusion is applied by multiplying a factor of the form

.. math::
   d(\ell)=1+D\left[\frac{\ell+1-\ell_{hd}}{\ell_{max}+1-\ell_{hd}} \right]^{\beta}
to the diffusion terms for degrees :math:`\ell \geq \ell_{hd}`.

* **difnu** is a real. This is the amplitude :math:`D` of the viscous hyperdiffusion.

* **difkappa** is a real. This is the amplitude :math:`D` of the thermal hyperdiffusion.

* **difeta** is a real. This is the amplitude :math:`D` of the magnetic hyperdiffusion.

* **ldif** is an integer. This is the degree :math:`\ell_{hd}` where hyperdiffusion starts to act.

* **ldifexp** is an integer. This is the exponent :math:`\beta` of hyperdiffusion.
