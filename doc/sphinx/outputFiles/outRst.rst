.. _secRestartFile:

Restart files ``checkpoint_*.TAG``
===================================

.. note:: These frequency of writing these files are determined by the standard
 inputs mentioned in the section on :ref:`restart files <secRstIn>` in the
 :ref:`output control namelist <secOutputNml>`. If nothing is specified then, by
 default one restart file is written at the end of the run.
 
.. note:: A restart file is read **only** when :ref:`l_start = .true. <varl_start>`

These are unformatted fortran files containing a snapshot of information about
spectral coefficients and physical and grid parameters. As the name suggests,
these files are used to 'restart' a run from a specific time. One such file is
read by the code at the beginning and are used as initial conditions for the
run. These are very useful for continuing a simulation for a long time on
computing clusters where the time for a single run is limited.

The file to be read at the beginning is specified by the input parameter
:ref:`start_file <varstart_file>` which takes in a string providing path to the
file.

These files are written by the subroutine :f:subr:`store <storeCheckPoints/store>`.

The following notations will be used for the coefficients of potentials
(note that scalar fields like temperature and pressure do not have a poloidal/toroidal
decomposition):

    +----------------+-----------+------------+
    | Field          | Poloidal  | Toroidal   |
    +================+===========+============+
    | Magnetic       | :f:var:`b`| :f:var:`aj`|
    +----------------+-----------+------------+
    | Velocity       | :f:var:`w`| :f:var:`z` |
    +----------------+-----------+------------+
    | Temperature    |      :f:var:`s`        |
    +----------------+-----------+------------+
    | Pressure       |      :f:var:`p`        |
    +----------------+------------------------+

Time derivatives are denoted with a self-explanatory notation. e.g, :f:var:`dbdt` is the first derivative of :f:var:`b`.

The word ``Last`` appended to a variable name denotes that the value is of the
time-step previous to the one during which the file is being written. They are
needed for the time-stepping schemes.

``_ic`` with a variable name says that it belongs to the Inner Core.
  
  .. code-block:: fortran
    
    !--------
    ! Line 1
    !--------

    time*tScale, dt*tScale, ra, pr, prmag, ek, radratio, inform, n_r_max,
    n_theta_max, n_phi_tot, minc, nalias, n_r_ic_max, sigma_ratio

    if (l_heat):                    !Run involving heat transport
                                    !(Convection)
    !---------
    ! Line 2
    !---------

       w,z,p,s                      
    
    !---------
    ! Line 3
    !---------

  
       dsdtLast,dwdtLast,dzdtLast,dpdtLast
    
    else:                           
    
    !---------
    ! Line 2
    !---------

  
       w,z,p
   
    !---------
    ! Line 3
    !---------

      dwdtLast,dzdtLast,dpdtLast

    if (l_mag):                                    !If magnetic run
    
    !---------
    ! Line 4
    !---------


      b, aj, dbdtLast, djdtLast

     if(l_mag .and. l_cond_ic):                    !If magnetic run
                                                   !and conducting inner core
    !---------
    ! Line 5
    !---------

      b_ic, aj_ic, dbdt_icLast, djdt_icLast
    
    !--------------------------------------------------
    ! Line 4 or 5 or 6 depending on l_mag and l_cond_ic
    !--------------------------------------------------

    lorentz_torque_icLast, lorentz_torque_maLast, !Information about torques,
    omega_ic1, omegaOsz_ic1, tOmega_ic1,          !prescribed rotation and
    omega_ic2, omegaOsz_ic2, tOmega_ic2,          !oscillation rates,
    omega_ma1, omegaOsz_ma1, tOmega_ma1,          !and the time step-size
    omega_ma2, omegaOsz_ma2, tOmega_ma2,
    dtNew

The checkpoint files can be read using the python class :py:class:`MagicCheckpoint <magic.MagicCheckpoint>`.

    >>> chk = MagicCheckpoint(filename='checkpoint_end.test')
    >>> # print size of poloidal and l_max
    >>> print(chk.wpol.shape, chk.l_max)
    >>> # convert from cheb to FD using 96 grid points
    >>> chk.cheb2fd(96)
    >>> write new file
    >>> chk.write('checkpoint_fd.test')
