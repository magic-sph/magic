.. _secPotFiles:

Potential files ``[V|B|T|Xi]_lmr_#.TAG``
=========================================

Those files contain a snapshot of either poloidal/toroidal potentials
``V_lmr_#.TAG`` and ``B_lmr_#.TAG`` or a scalar like temperature/entropy
or chemical composition (``T_lmr_#.TAG`` or ``Xi_lmr_#.TAG``) in the
radial space for all spherical harmonic degrees and orders.
The detailed calculations are done in
the subroutine :f:subr:`write_Pot <out_coeff/write_pot()>`. The outputs are
stored as a fortran unformatted file with a stream access. It has the
following structure

   .. code-block:: fortran

       !-------------
       ! Header
       !-------------
       version
       time, ra, pr, raxi, sc, prmag, ek, radratio, sigma_ration ! Parameters
       n_r_max, n_r_ic_max, l_max, minc, lm_max                  ! Truncation
       omega_ic, omega_ma                                        ! Rotation rates
       r(1), r(2), ..., r(n_r_max)                               ! Radius
       rho0(1), rho0(2), ..., rho0(n_r_max)                      ! Background density

       !-------------
       ! Poloidal potential or scalar
       !-------------
       w(lm=1,n_r=1), w(lm=2, n_r=1), ..., w(lm=lm_max, n_r=1),
       ...
       w(lm=1,n_r=n_r_max), ..., w(lm=lm_max,n_r=n_r_max)

       !-------------
       ! If stored: toroidal potential
       !-------------
       z(lm=1,n_r=1), z(lm=2, n_r=1), ..., z(lm=lm_max, n_r=1),
       ...
       z(lm=1,n_r=n_r_max), ..., z(lm=lm_max,n_r=n_r_max)

       !**************************************************************************!
       ! This last part is optional and are is written when there is              !
       ! an electrically-conducting inner-core                                    !
       !**************************************************************************!
       b_ic(lm=1,n_r=1), b_ic(lm=2, n_r=1), ..., b_ic(lm=lm_max, n_r=1),
       ...
       b_ic(lm=1,n_r=n_r_max), ..., b_ic(lm=lm_max,n_r=n_r_max)
       aj_ic(lm=1,n_r=1), aj_ic(lm=2, n_r=1), ..., aj_ic(lm=lm_max, n_r=1),
       ...
       aj_ic(lm=1,n_r=n_r_max), ..., aj_ic(lm=lm_max,n_r=n_r_max)


The potential files can be read and transformed to the physical space using the python class :py:class:`MagicPotential <magic.MagicPotential>`.

    >>> p = MagicPotential(field='V')
    >>> print(p.pol[p.idex[3,2], 32]) # print w(l=3,m=2,n_r=32)

Once transformed to the physical space using a Fourier and a Legendre transform, they
can be displayed:

    >>> p.equat(field='vr', cm='jet', levels=50) # Equatorial cut of vr
    >>> p.avg(field='vp') # Azimuthal average of vphi
    >>> p.surf(field='vt', r=0.8) # Radial cut of vtheta at r=0.8r_o
