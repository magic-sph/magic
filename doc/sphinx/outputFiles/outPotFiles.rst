
Potential files ``[V|B|T]_lmr_#.TAG``
=====================================

.. _secVpotFile:

``V_lmr_#.TAG`` and ``V_lmr_ave.TAG``
-------------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeVpot=.true. <varl_storeVpot>`.


These files contain a snapshot of the poloidal and toroidal flow potentials
:f:var:`w` and :f:var:`z` in the radial space for all spherical harmonic
degrees and orders. They basically contain two arrays of dimension
(:f:var:`lm_max`, :f:var:`n_r_max`).  The detailed calculations are done in
the subroutine :f:subr:`write_Pot <out_coeff/write_pot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_r_max, n_r_ic_max, minc, lm_max           ! Header (truncation)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (physics)

       !------------
       ! Line 3
       !-------------

       time                                               ! Header (time)

       !------------
       ! Line 4
       !-------------

       r(1), ..., r(n_r_max),rho0(1), ..., rho0(n_r_max)  ! Radius, Density

       !------------
       ! Line 5
       !-------------
 
                                                          ! Poloidal potential
       w(lm=1,n_r=1), w(lm=2, n_r=1), ..., w(lm=lm_max, n_r=1),
       ...
       w(lm=1,n_r=n_r_max), ..., w(lm=lm_max,n_r=n_r_max)

       !------------
       ! Line 6
       !-------------
 
                                                          ! Toroidal potential
       z(lm=1,n_r=1), z(lm=2, n_r=1), ..., z(lm=lm_max, n_r=1),
       ...
       z(lm=1,n_r=n_r_max), ..., z(lm=lm_max,n_r=n_r_max)


.. _secBpotFile:

``B_lmr_#.TAG``, ``B_lmr_ave.TAG``
----------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeBpot=.true. <varl_storeBpot>`.

These files contain a snapshot of the poloidal and toroidal magnetic potentials
:f:var:`b` and :f:var:`aj` in the radial space for all spherical harmonic
degrees and orders.  The detailed calculations are done in
the subroutine :f:subr:`write_Pot <out_coeff/write_pot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_r_max, n_r_ic_max, minc, lm_max           ! Header (truncation)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (physics)

       !------------
       ! Line 3
       !-------------

       time                                               ! Header (time)

       !------------
       ! Line 4
       !-------------

       r(1), ..., r(n_r_max),rho0(1), ..., rho0(n_r_max)  ! Radius, Density

       !------------
       ! Line 5
       !-------------
 
                                                          ! Poloidal potential
       b(lm=1,n_r=1), b(lm=2, n_r=1), ..., b(lm=lm_max, n_r=1),
       ...
       b(lm=1,n_r=n_r_max), ..., b(lm=lm_max,n_r=n_r_max)

       !------------
       ! Line 6
       !-------------
 
                                                          ! Toroidal potential
       aj(lm=1,n_r=1), aj(lm=2, n_r=1), ..., aj(lm=lm_max, n_r=1),
       ...
       aj(lm=1,n_r=n_r_max), ..., aj(lm=lm_max,n_r=n_r_max)

       !**************************************************************************!
       ! The two following lines are optional and are only written when there is  !
       ! an electrically-conducting inner-core                                    !
       !**************************************************************************!

       !------------
       ! Line 7
       !-------------
 
       time,                                              ! Time and poloidal potential
       b_ic(lm=1,n_r=1), b_ic(lm=2, n_r=1), ..., b_ic(lm=lm_max, n_r=1),
       ...
       b_ic(lm=1,n_r=n_r_max), ..., b_ic(lm=lm_max,n_r=n_r_max)

       !------------
       ! Line 8
       !-------------
 
       time,                                              ! Time and toroidal potential
       aj_ic(lm=1,n_r=1), aj_ic(lm=2, n_r=1), ..., aj_ic(lm=lm_max, n_r=1),
       ...
       aj_ic(lm=1,n_r=n_r_max), ..., aj_ic(lm=lm_max,n_r=n_r_max)



.. _secTpotFile:

``T_lmr_#.TAG``, ``T_lmr_ave.TAG``
----------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeTpot=.true. <varl_storeTpot>`.

These files contain a snapshot of the temperature/entropy
:f:var:`s` in the spectral and radial spaces for all spherical harmonic
degrees and orders. They basically contain one array of dimension
(:f:var:`lm_max`, :f:var:`n_r_max`).  The detailed calculations are done in
the subroutine :f:subr:`write_Pot <out_coeff/write_pot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_r_max, n_r_ic_max, minc, lm_max           ! Header (truncation)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (physics)

       !------------
       ! Line 3
       !-------------

       time                                               ! Header (time)

       !------------
       ! Line 4
       !-------------

       r(1), ..., r(n_r_max),rho0(1), ..., rho0(n_r_max)  ! Radius, Density

       !------------
       ! Line 5
       !-------------
 
                                                          ! Temperature/entropy
       s(lm=1,n_r=1), s(lm=2, n_r=1), ..., s(lm=lm_max, n_r=1),
       ...
       s(lm=1,n_r=n_r_max), ..., s(lm=lm_max,n_r=n_r_max)
