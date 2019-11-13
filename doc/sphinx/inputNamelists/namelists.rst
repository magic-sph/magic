.. _secNamelists:

Input parameters
################

True runtime input parameters are read from STDIN as namelists, a Fortran
feature. A namelist is identified by its unique name `&name`. The
name-statement is followed by the parameters that are part of the namelist in
the format `parameter=value,`. The namelist is closed by a backslash. The
subroutine `defaultNamelists` (in the module `Namelist.f90`) defines a default
value for each parameter. Only the parameters whose value should differ from
its default have to be stated in the namelist.

An example for the short namelist defining inner core parameters is

.. code-block:: fortran

   &inner_core
     sigma_ratio = 1.0,
     nRotIc      = 1

Comas can be used to seperate namelist entries since they are not interpreted by the code.

Magic uses the following **eight namelists** :

.. topic:: Namelists

     1. :ref:`&grid <secGridNml>` for resolution

     2. :ref:`&control <secControlNml>` for control parameters and numerical parameters.

     3. :ref:`&phys_param <secPhysNml>` for the physical parameters.

     4. :ref:`&B_external <secBextnml>` for setting up an external field contribution

     5. :ref:`&start_field <secStartNml>` to define the starting fields.

     6. :ref:`&output_control <secOutputNml>` for defining the output.

     7. :ref:`&mantle <secMantle>` for setting mantle parameters.

     8. :ref:`&inner_core <secInnerCore>` for setting inner core parameters.

The number of possible input parameters has grown to more than 100/150. **Don't be confused
by all the possible options though, since all parameters are internally set to a useful 
default value!** 

Practically, in a production run, the number of parameters you may want
to adjust is thus much smaller. As an example, the following namelist shows you how
to initiate and quickly run one of the anelastic benchmarks by (`Jones et al., 2011 
<http://dx.doi.org/10.1016/j.icarus.2011.08.014>`_):

.. code-block:: fortran

     &grid
      n_r_max     =97,           ! 97 radial grid points
      n_cheb_max  =95,
      n_phi_tot   =288,          ! 288 points in the azimuthal direction
      n_r_ic_max  =17,
     n_cheb_ic_max=15,
      minc        =1,            ! No azimuthal symmetry
     /
     &control
      mode        =1,            ! This is a non-magnetic case
      tag         ="test",       ! Trailing name of the outputs produced by the code
      n_time_steps=50000,        ! Number of time steps
      courfac     =2.5D0,        ! Courant factor (flow)
      alffac      =1.0D0,        ! Courant factor (magnetic field)
      dtmax       =1.0D-4,       ! Maximum allowed time-step
      alpha       =0.6D0,
      runHours    =23,           ! Run time (hours)
      runMinutes  =30,           ! Run time (minutes)
      time_scheme ='CNAB2',      ! Name of the time stepper
     /
     &phys_param
      ra          =1.48638035D5, ! Rayleigh number
      ek          =1.0D-3,       ! Ekman number
      pr          =1.0D0,        ! Prandtl number
      strat       =5.D0,         ! Density contrast
      polind      =2.0D0,        ! Polytropic index
      radratio    =0.35D0,       ! Aspect ratio of the spherical shell
      g0          =0.D0,         ! Gravity profile
      g1          =0.D0,
      g2          =1.D0,
      ktops       =1,            ! Entropy boundary condition
      kbots       =1,
      ktopv       =1,            ! Mechanical boundary condition
      kbotv       =1,
     /
     &start_field
      l_start_file=.false.,
      start_file  ="checkpoint_end.CJ3",
      init_s1     =1919,         ! Initial entropy perturbation pattern
      amp_s1      =0.01,         ! Amplitude of the initial perturbation
     /
     &output_control
      n_log_step  =50,           ! Store time series every 50 time steps
      n_graphs    =1,            ! 1 G_#.TAG file produced at the end of the run
      n_specs     =5,            ! 5 spectra produced during the run
      n_rsts      =1,            ! 1 checkpoint_end.TAG file produced at the end of the run
      runid       ="C.Jones bench", 
     /
     &mantle
      nRotMa      =0             ! Non-rotating mantle
     /
     &inner_core
      sigma_ratio =0.d0,         ! Non-conducting inner core
      nRotIC      =0,            ! Non-rotating inner core
     /

This example might then be easily adapted to your desired configuration.


.. toctree::
   :hidden:
   :maxdepth: 1

   gridNamelist.rst
   controlNamelist.rst
   physNamelist.rst
   BextNamelist.rst
   startNamelist.rst
   outNamelist.rst
   mantle_icNamelist.rst



