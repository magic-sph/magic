!$Id$
!***********************************************************************
    SUBROUTINE initB(b,aj,b_ic,aj_ic,                     &
                     lorentz_torque_ic,lorentz_torque_ma, &
                     lmStart,lmStop)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to initialize the magnetic field   |
!  |  according to the control parameters imagcon and init_b1/2.       |
!  |  In addition CMB and ICB peak values are calculated for           |
!  |  magneto convection.                                              |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
     
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE init_fields
    USE Grenoble
    USE blocking
    USE horizontal_data
    USE logic
    USE const
    USE usefull, only: random

    IMPLICIT NONE

    INTEGER :: lmStart,lmStop

!-- output:
    REAL(kind=8) :: lorentz_torque_ic
    REAL(kind=8) :: lorentz_torque_ma

    COMPLEX(kind=8) :: b(lm_max,n_r_max)
    COMPLEX(kind=8) :: aj(lm_max,n_r_max)
    COMPLEX(kind=8) :: b_ic(lm_max,n_r_ic_max)
    COMPLEX(kind=8) :: aj_ic(lm_max,n_r_ic_max)

!-- local:
    INTEGER :: lm,lm0,lmStart00
    INTEGER :: n_r
    REAL(kind=8) :: b_pol,b_tor
    COMPLEX(kind=8) :: aj0(n_r_max+1)
    COMPLEX(kind=8) :: aj0_ic(n_r_ic_max)
    REAL(kind=8) :: arg,aj_ic1,aj_ic2

    REAL(kind=8) :: b1(n_r_max)
    REAL(kind=8) :: b1_ic(n_r_ic_max)
    REAL(kind=8) :: bR,bI,rr
    REAL(kind=8) :: aVarCon,bVarCon
    INTEGER :: bExp

    INTEGER :: l1m0,l2m0,l3m0,l1m1

!-- end of declaration
!-------------------------------------------------------------------
     
    lmStart00=MAX(lmStart,1)

!        IF ( lGrenoble ) CALL getB0

    l1m0=lm2(1,0)
    l2m0=lm2(2,0)
    l3m0=lm2(3,0)
    l1m1=lm2(1,1)

    lm0=l2m0 ! Default quadrupole field

    IF ( imagcon == -1 ) THEN

    !----- impose l=1,m=0 poloidal field at ICB:
        lm0=l1m0
        bpeakbot=-dsqrt(pi/3.d0)*r_icb**2*amp_b1
        bpeaktop=0.d0

    ELSE IF ( imagcon == -2 ) THEN

    !----- impose l=1,m=0 poloidal field at CMB:
        lm0=l1m0
        bpeakbot=0.d0
        bpeaktop=-dsqrt(pi/3.d0)*r_cmb**2*amp_b1

    ELSE IF ( imagcon == 1 ) THEN

    !----- impose l=2,m=0 toroidal field at ICB:
        lm0=l2m0
        bpeakbot=4.d0/3.d0*dsqrt(pi/5.d0)*r_icb*amp_b1
        bpeaktop=0.d0

    ELSE IF ( imagcon == 10 ) THEN

    !----- impose l=2,m=0 toroidal field at ICB and CMB:
        lm0=l2m0
        bpeakbot=4.d0/3.d0*dsqrt(pi/5.d0)*r_icb*amp_b1
        bpeaktop=4.d0/3.d0*dsqrt(pi/5.d0)*r_cmb*amp_b1

    ELSE IF ( imagcon == 11 ) THEN

    !----- same as imagcon.eq.10 but opposite sign at CMB:
        lm0=l2m0
        bpeakbot=4.d0/3.d0*dsqrt(pi/5.d0)*r_icb*amp_b1
        bpeaktop=-4.d0/3.d0*dsqrt(pi/5.d0)*r_cmb*amp_b1

    ELSE IF ( imagcon == 12 ) THEN

    !----- impose l=1,m=0 toroidal field at ICB and CMB:
        lm0=l1m0
        bpeakbot=2.d0*dsqrt(pi/3.d0)*r_icb*amp_b1
        bpeaktop=2.d0*dsqrt(pi/3.d0)*r_cmb*amp_b1

    ELSE IF ( imagcon == 0 ) THEN

        lm0=l2m0
        bpeakbot=0.d0
        bpeaktop=0.d0

    ELSE IF ( imagcon == -10 ) THEN

    !----- Test of variable conductivity case with analytical solution:
    !      Assume the magnetic diffusivity is lambda=r**5, that the aspect ratio
    !      is 0.5, and that there is no flow.
    !      The analytical stationary solution for the (l=3,m=0) toroidal field
    !      with bounday condition aj(r=r_ICB)=1, aj(r=r_CMB)=0 is then
    !      given by jVarCon(r)!
    !      A disturbed solution is used to initialize aj,
    !      the disturbance should decay with time.
    !      The solution is stored in file testVarCond.TAG at the end of the run,
    !      where the first column denotes radius, the second is aj(l=3,m=0,r) and
    !      the third is jVarCon(r). Second and third should be identical when
    !      the stationary solution has been reached.
        lm0=l3m0  ! This is l=3,m=0
        bpeakbot=1.D0
        bpeaktop=0.D0
        aVarCon =-1.D0/255.D0
        bVarCon =256.D0/255.D0
        IF ( lmStart <= lm0 .AND. lmStop >= lm0 ) THEN ! select processor
            DO n_r=1,n_r_max             ! Diffusive toroidal field
                jVarCon(n_r)=aVarCon*r(n_r)**2 + bVarCon/(r(n_r)**6)
                aj(lm0,n_r) =jVarCon(n_r) + 0.1D0*DSIN((r(n_r)-r_ICB)*pi)
            END DO
        END IF

    END IF

    IF ( init_b1 == 1 .OR. imagcon > 0 ) THEN

    !----- Conductive solution for toroidal field,
    !      diffusion equation solved in j_cond, amplitude defined
    !      by bpeaktop and bpeakbot respectively.
    !      bpeakbot is only used for insulating inner core !

         
        IF ( lmStart <= lm0 .AND. lmStop >= lm0 ) THEN ! select processor
            CALL j_cond(lm0,aj0,aj0_ic)
            DO n_r=1,n_r_max             ! Diffusive toroidal field
                aj(lm0,n_r)=aj0(n_r)
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    aj_ic(lm0,n_r)=aj0_ic(n_r)
                END DO
            END IF
        END IF
         
    ELSE IF ( init_b1 == 2 ) THEN  ! l=1,m=0 analytical toroidal field
    ! with a maximum of amp_b1 at mid-radius
    ! between r_icb and r_cmb for an insulating
    ! inner core and at r_cmb/2 for a conducting
    ! inner core

        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_tor=-2.d0*amp_b1*dsqrt(pi/3.d0)  ! minus sign makes phi comp. > 0
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_max
                    aj(l1m0,n_r)=aj(l1m0,n_r) + &
                                 b_tor*r(n_r)*dsin(pi*r(n_r)/r_cmb)
                END DO
                DO n_r=1,n_r_ic_max
                    aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r) + &
                                    b_tor*r_ic(n_r)*dsin(pi*r_ic(n_r)/r_cmb)
                END DO
            ELSE
                DO n_r=1,n_r_max
                    aj(l1m0,n_r)=aj(l1m0,n_r) + &
                                 b_tor*r(n_r)*dsin(pi*(r(n_r)-r_icb))
                END DO
            END IF
        END IF
         
    ELSE IF ( init_b1 == 3 ) THEN  ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
    ! toroidal field has again its maximum of amp_b1
    ! at mid-radius between r_icb and r_cmb for an
    ! insulating inner core and at r_cmb/2 for a
    ! conducting inner core
    ! The outer core poloidal field is defined by
    ! a homogeneous  current density, its maximum at
    ! the ICB is set to amp_b1.
    ! The inner core poloidal field is chosen accordingly.
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_tor=-4.d0/3.d0*amp_b1*dsqrt(pi/5.d0)
            IF ( l_cond_ic ) THEN
                b_pol=amp_b1*DSQRT(3.d0*pi)/(3.d0+r_cmb)
                DO n_r=1,n_r_max
                    b(l1m0,n_r)=b(l1m0,n_r) + &
                                b_pol*(r(n_r)**3 - 4.d0/3.d0*r_cmb*r(n_r)**2)
                END DO
                arg=pi*r_icb/r_cmb
                aj_ic1=(arg-2.d0*dsin(arg)*dcos(arg)) / &
                       (arg+dsin(arg)*dcos(arg))
                aj_ic2=(1.d0-aj_ic1)*r_icb*dsin(arg)/dcos(arg)
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * ( &
                                           r_ic(n_r)**2/2.d0/r_icb + &
                                                        r_icb/2.d0 - &
                                                   4.d0/3.d0*r_cmb )
                END DO
            ELSE
                b_pol=amp_b1*dsqrt(3.d0*pi)/4.d0
                DO n_r=1,n_r_max
                    b(l1m0,n_r)=b(l1m0,n_r)+ b_pol *     ( &
                                               r(n_r)**3 - &
                               4.d0/3.d0*r_cmb*r(n_r)**2 + &
                            1.d0/3.d0*r_icb**4/r(n_r)    )
                END DO
            END IF
        END IF

        IF ( lmStart <= l2m0 .AND. lmStop >= l2m0 ) THEN ! select processor
            b_tor=-4.d0/3.d0*amp_b1*dsqrt(pi/5.d0)
            IF ( l_cond_ic ) THEN
                b_pol=amp_b1*dsqrt(3.d0*pi)/(3.d0+r_cmb)
                DO n_r=1,n_r_max
                    aj(l2m0,n_r)=aj(l2m0,n_r) + &
                                 b_tor*r(n_r)*dsin(pi*(r(n_r)/r_cmb))
                END DO
                arg=pi*r_icb/r_cmb
                aj_ic1=(arg-2.d0*dsin(arg)*dcos(arg)) / &
                       (arg+dsin(arg)*dcos(arg))
                aj_ic2=(1.d0-aj_ic1)*r_icb*dsin(arg)/dcos(arg)
                DO n_r=1,n_r_ic_max
                    aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor*             ( &
                             aj_ic1*r_ic(n_r)*dsin(pi*r_ic(n_r)/r_cmb) + &
                                       aj_ic2*dcos(pi*r_ic(n_r)/r_cmb) )
                END DO
            ELSE
                b_pol=amp_b1*dsqrt(3.d0*pi)/4.d0
                DO n_r=1,n_r_max
                    aj(l2m0,n_r)=aj(l2m0,n_r) + b_tor * &
                                 r(n_r)*dsin(pi*(r(n_r)-r_icb))
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 4 .OR. imagcon == -1 ) THEN  ! l=1,m0 poloidal field
    ! with max field amplitude amp_b1 at r_icb
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_pol=-amp_b1*r_icb**3*dsqrt(pi/3.d0)
            DO n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol*or1(n_r)
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol/r_icb* &
                                   ( -3.d0/2.d0 + (r_ic(n_r)/r_icb)**2/2.d0 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 5 ) THEN  ! l=1,m0 poloidal field
    ! constant j density, defined max field value amp_v1 at r_cmb
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            IF ( l_cond_ic ) THEN
                b_pol=amp_b1*dsqrt(3.d0*pi)/r_cmb
                DO n_r=1,n_r_max
                    b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                                4.d0/3.d0*r_cmb * r(n_r)**2 )
                END DO
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * &
                       (-5.d0/6.d0*r_icb-4.d0/3.d0+r_ic(n_r)**2/r_icb/2.d0)
                END DO
            ELSE
                b_pol=amp_b1*dsqrt(3.d0*pi)/(r_cmb*(1.d0-radratio**4))
                DO n_r=1,n_r_max
                    b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                                     4.d0/3.d0*r_cmb * r(n_r)**2 + &
                                     1.d0/3.d0*r_icb**4 / r(n_r)    )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 6 ) THEN  ! l=1,m=0 poloidal field , constant in r !
    ! no potential at r_cmb but simple
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_pol=amp_b1
            DO n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol*r(n_r)**2
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 7 .OR. imagcon == -2 ) THEN  ! l=1,m0 poloidal field
    ! which is potential field at r_cmb
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_pol=amp_b1*5.d0/2.d0*dsqrt(pi/3.d0)*r_icb**2
            DO n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                            ( 1.d0 - 3.d0/5.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                   ( 1.d0 - 3.d0/5.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 8 ) THEN  ! l=1,m0 pol. field, l=2,m=0 toroidal field
    ! which is potential field at r_cmb
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_pol=amp_b1*5.d0/2.d0*dsqrt(pi/3.d0)*r_icb**2
            DO n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                            ( 1.d0 - 3.d0/5.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                   ( 1.d0 - 3.d0/5.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

        IF ( lmStart <= l2m0 .AND. lmStop >= l2m0 ) THEN ! select processor
            b_tor=amp_b1*3.d0/2.d0*dsqrt(pi/5.d0)*r_icb**2*radratio
            DO n_r=1,n_r_max
                aj(l2m0,n_r)=aj(l2m0,n_r)+b_tor*(r(n_r)/r_icb)**3 * &
                             ( 1.d0 - (r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor * &
                                    ( 1.d0 - (r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 9 ) THEN  ! l=2,m0 poloidal field
    ! which is potential field at r_cmb
        IF ( lmStart <= l2m0 .AND. lmStop >= l2m0 ) THEN ! select processor
            b_pol=amp_b1*7.d0/6.d0*dsqrt(pi/5.d0)*r_icb**2*radratio
            DO n_r=1,n_r_max
                b(l2m0,n_r)=b(l2m0,n_r)+b_pol*(r(n_r)/r_icb)**3 * &
                            ( 1.d0 - 5.d0/7.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l2m0,n_r)=b_ic(l2m0,n_r)+b_pol * &
                                   ( 1.d0 - 5.d0/7.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 10 ) THEN  ! only equatorial dipole

        IF ( lmStart <= l1m1 .AND. lmStop >= l1m1 ) THEN ! select processor
            IF ( l1m1 <= 0 ) THEN
                WRITE(*,*) '! Can not initialize l=1,m=1 !'
                STOP
            END IF
            b_pol=amp_b1*5.d0/2.d0*dsqrt(pi/3.d0)*r_icb**2
            DO n_r=1,n_r_max
                b(l1m1,n_r)=b(l1m1,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                            ( 1.d0 - 3.d0/5.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol * &
                                   ( 1.d0 - 3.d0/5.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 < 0 ) THEN  ! l,m mixture, random init

        bExp=IABS(init_b1)
        DO n_r=1,n_r_max
            b1(n_r)=(r(n_r)/r_cmb)**2 * &
                    ( 1.D0-3.D0/5.D0*(r(n_r)/r_cmb)**2 )
        END DO
        IF ( l_cond_ic ) THEN
            DO n_r=1,n_r_ic_max
                b1_ic(n_r)= ( 1.D0-3.D0/5.D0*(r_ic(n_r)/r_cmb)**2 )
            END DO
        END IF

    !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
        rr=random(1.D0)
        DO lm=lmStart00,lmStop
            bR=(-1.d0+2.d0*random(0.d0))*amp_b1/D_l(lm)**(bExp-1)
            bI=(-1.d0+2.d0*random(0.d0))*amp_b1/D_l(lm)**(bExp-1)
            IF ( lm2m(lm) == 0 ) bI=0.D0
            DO n_r=1,n_r_max
                b(lm,n_r)=b(lm,n_r) + CMPLX(bR*b1(n_r),bI*b1(n_r),KIND=KIND(0d0))
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(lm,n_r)=b_ic(lm,n_r) + &
                                 CMPLX(bR*b1_ic(n_r),bI*b1_ic(n_r),KIND=KIND(0d0))
                END DO
            END IF
        END DO

    ELSE IF ( init_b1 == 11 ) THEN  ! axial and equatorial dipole

        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            b_pol=amp_b1*5.d0/2.d0*dsqrt(pi/3.d0)*r_icb**2
            DO n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                            ( 1.d0 - 3.d0/5.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                   ( 1.d0 - 3.d0/5.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

        IF ( lmStart <= l1m1 .AND. lmStop >= l1m1 ) THEN ! select processor
            IF ( l1m1 <= 0 ) THEN
                WRITE(*,*) '! Cannot initialize l=1,m=1 !'
                STOP
            END IF
            b_pol=amp_b1*5.d0/2.d0*dsqrt(pi/3.d0)*r_icb**2
            DO n_r=1,n_r_max
                b(l1m1,n_r)=b(l1m1,n_r) +                   &
                            b_pol/10.D0*(r(n_r)/r_icb)**2 * &
                            ( 1.d0 - 3.d0/5.d0*(r(n_r)/r_cmb)**2 )
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol/5.D0 * &
                                   ( 1.d0 - 3.d0/5.d0*(r_ic(n_r)/r_cmb)**2 )
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 21 ) THEN ! toroidal field created by inner core rotation
    ! equatoriall symmetric
        IF ( lmStart <= l1m0 .AND. lmStop >= l1m0 ) THEN ! select processor
            DO n_r=1,n_r_max
                aj0(n_r)=amp_b1*(r_icb/r(n_r))**6
            END DO
            DO n_r=1,n_r_max             ! Diffusive toroidal field
                aj(l1m0,n_r)=aj(l1m0,n_r)+aj0(n_r)
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r)+aj0(n_r_icb)
                END DO
            END IF
        END IF

    ELSE IF ( init_b1 == 22 ) THEN ! toroidal field created by inner core rotation
    ! equatoriall asymmetric
        IF ( lmStart <= l2m0 .AND. lmStop >= l2m0 ) THEN ! select processor
            DO n_r=1,n_r_max
                aj0(n_r)=amp_b1*(r_icb/r(n_r))**6
            END DO
            DO n_r=1,n_r_max             ! Diffusive toroidal field
                aj(l2m0,n_r)=aj(l2m0,n_r)+aj0(n_r)
            END DO
            IF ( l_cond_ic ) THEN
                DO n_r=1,n_r_ic_max
                    aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+aj0(n_r_icb)
                END DO
            END IF
        END IF

    END IF

!-- Too lazy to calculate these:
    lorentz_torque_ic=0.d0
    lorentz_torque_ma=0.d0


    RETURN
    end SUBROUTINE initB

!-----------------------------------------------------------------------
