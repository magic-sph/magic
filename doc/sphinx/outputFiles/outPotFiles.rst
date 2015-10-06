
Potential files ``[V|B|T]pot_#.TAG``
====================================

.. _secVpotFile:

``Vpot_#.TAG`` and ``Vpot_ave.TAG``
-----------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeVpot=.true. <varl_storeVpot>`.


These files contain a snapshot of the poloidal and toroidal flow potentials
:f:var:`w` and :f:var:`z` in the Chebyshev space for all spherical harmonic
degrees and orders. They basically contain two arrays of dimension
(:f:var:`lm_max`, :f:var:`n_cheb_max`).  The detailed calculations are done in
the subroutine :f:subr:`storePot <store_pot_mod/storepot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_cheb_max, n_cheb_ic_max, minc, lm_max     ! Header (truncation informations)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (informations about phyics)

       !------------
       ! Line 3
       !-------------
 
       time,                                              ! Time and poloidal potential
       w(lm=1,n_cheb=1), w(lm=2, n_cheb=1), ..., w(lm=lm_max, n_cheb=1),
       ...
       w(lm=1,n_cheb=n_cheb_max, ..., w(lm=lm_max,n_cheb=n_cheb_max)

       !------------
       ! Line 4
       !-------------
 
       time,                                              ! Time and toroidal potential
       z(lm=1,n_cheb=1), z(lm=2, n_cheb=1), ..., z(lm=lm_max, n_cheb=1),
       ...
       z(lm=1,n_cheb=n_cheb_max, ..., z(lm=lm_max,n_cheb=n_cheb_max)


.. _secBpotFile:

``Bpot_#.TAG``, ``Bpot_ave.TAG``
--------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeBpot=.true. <varl_storeBpot>`.

These files contain a snapshot of the poloidal and toroidal magnetic potentials
:f:var:`b` and :f:var:`aj` in the Chebyshev space for all spherical harmonic
degrees and orders.  The detailed calculations are done in
the subroutine :f:subr:`storePot <store_pot_mod/storepot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_cheb_max, n_cheb_ic_max, minc, lm_max     ! Header (truncation informations)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (informations about phyics)

       !------------
       ! Line 3
       !-------------
 
       time,                                              ! Time and poloidal potential
       b(lm=1,n_cheb=1), b(lm=2, n_cheb=1), ..., b(lm=lm_max, n_cheb=1),
       ...
       b(lm=1,n_cheb=n_cheb_max, ..., b(lm=lm_max,n_cheb=n_cheb_max)

       !------------
       ! Line 4
       !-------------
 
       time,                                              ! Time and toroidal potential
       aj(lm=1,n_cheb=1), aj(lm=2, n_cheb=1), ..., aj(lm=lm_max, n_cheb=1),
       ...
       aj(lm=1,n_cheb=n_cheb_max, ..., aj(lm=lm_max,n_cheb=n_cheb_max)

       !**************************************************************************!
       ! The two following lines are optional and are only written when there is  !
       ! an electrically-conducting inner-core                                    !
       !**************************************************************************!

       !------------
       ! Line 5
       !-------------
 
       time,                                              ! Time and poloidal potential
       b_ic(lm=1,n_cheb=1), b_ic(lm=2, n_cheb=1), ..., b_ic(lm=lm_max, n_cheb=1),
       ...
       b_ic(lm=1,n_cheb=n_cheb_max, ..., b_ic(lm=lm_max,n_cheb=n_cheb_max)

       !------------
       ! Line 6
       !-------------
 
       time,                                              ! Time and toroidal potential
       aj_ic(lm=1,n_cheb=1), aj_ic(lm=2, n_cheb=1), ..., aj_ic(lm=lm_max, n_cheb=1),
       ...
       aj_ic(lm=1,n_cheb=n_cheb_max, ..., aj_ic(lm=lm_max,n_cheb=n_cheb_max)



.. _secTpotFile:

``Tpot_#.TAG``, ``Tpot_ave.TAG``
--------------------------------

.. note:: These output files are **only** written when either when :ref:`l_storePot=.true. <varl_storePot>` or when :ref:`l_storeTpot=.true. <varl_storeTpot>`.

These files contain a snapshot of the temperature/entropy
:f:var:`s` in the spectral and Chebyshev spaces for all spherical harmonic
degrees and orders. They basically contain one array of dimension
(:f:var:`lm_max`, :f:var:`n_cheb_max`).  The detailed calculations are done in
the subroutine :f:subr:`storePot <store_pot_mod/storepot()>`. The outputs are
stored as a fortran unformatted file which follows the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       l_max, n_cheb_max, n_cheb_ic_max, minc, lm_max     ! Header (truncation informations)

       !-------------
       ! Line 2
       !-------------

       ra, ek, pr, prmag, sigma_ratio, omega_ma, omega_ic ! Header (informations about phyics)

       !------------
       ! Line 3
       !-------------
 
       time,                                              ! Time and temperature/entropy
       s(lm=1,n_cheb=1), s(lm=2, n_cheb=1), ..., s(lm=lm_max, n_cheb=1),
       ...
       s(lm=1,n_cheb=n_cheb_max, ..., s(lm=lm_max,n_cheb=n_cheb_max)
