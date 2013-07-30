!$Id$
!***********************************************************************
    SUBROUTINE get_B_surface(b_r,b_t,b_p, &
                             bCMB,n_theta_start,n_theta_block)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Upward continuation of laplacian field to Earths surface.        |
!  |  Field is given by poloidal harmonic coefficients b at CMB.       |
!  |  Spherical harmonic transforms of upward continued field          |
!  |  to r/theta/phi vector components for all logitudes and           |
!  |  latitude are returned in br/bt/bp.                               |
!  |  Note that this routine given the real components of the magnetic |
!  |  fields while other transforms in the code provide only:          |
!  |   r**2 br, r**2 sin(theta) bt, r**2 sin(theta) bp                 |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif

    IMPLICIT NONE

!-- Input of constant parameters:
! include 'truncation.f'
! include 'c_blocking.f'
! include 'c_radial.f'    ! contains radial functions
! include 'c_horizontal.f'! contains horizontal functions

!-- Input of variables:
    INTEGER :: n_theta_start   ! No. of theta to start with
    INTEGER :: n_theta_block   ! Size of theta block

!-- Input of poloidal magnetic field scalar:
    COMPLEX(kind=8) :: bCMB(lm_max)

!-- Output:
    REAL(kind=8) :: b_r(nrp,*)       ! Radial magnetic field in (phi,theta)-space
    REAL(kind=8) :: b_t(nrp,*)       ! Latitudinal magnetic field
    REAL(kind=8) :: b_p(nrp,*)       ! Azimuthal magnetic field.

!-- local:
    INTEGER :: n_theta         ! No. of theta
    INTEGER :: n_theta_nhs     ! Counter for theta in northern hemisphere
    INTEGER :: l,m,lm,mc       ! degree/order,counter

    REAL(kind=8) :: r_ratio          ! r_cmb/r_surface
    REAL(kind=8) :: sign             ! Sign for southern hemisphere
    REAL(kind=8) :: r_dep(l_max)     ! Radial dependence
    REAL(kind=8) :: O_sint           ! 1/sin(theta)
    COMPLEX(kind=8) :: cs1(lm_max),cs2(lm_max) ! help arrays
    COMPLEX(kind=8) :: b_r_1,b_t_1,b_p_1
    COMPLEX(kind=8) :: b_r_n,b_t_n,b_p_n
    COMPLEX(kind=8) :: b_r_s,b_t_s,b_p_s

!-- end of declaration
!---------------------------------------------------------------------

!-- Radial dependence:
    r_ratio=r_cmb/r_surface
    r_dep(1)=r_ratio/(r_surface*r_surface)  ! l=1 term
    DO l=2,l_max
        r_dep(l)=r_dep(l-1)*r_ratio
    END DO
     
!-- Construct help arrays containing radial dependence
!   and l dependence: dLh=l*(l+1)
    cs1(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    cs2(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
    DO lm=2,lm_max
        cs1(lm) = bCMB(lm)*dLh(lm)*r_dep(lm2l(lm))
        cs2(lm)= -bCMB(lm)*D_l(lm)*r_dep(lm2l(lm))
    END DO
     

!-- Build field components:

!----- Loop over colatitudes:

    DO n_theta=1,n_theta_block,2

        n_theta_nhs=(n_theta_start+n_theta)/2
        O_sint     =osn1(n_theta_nhs)

    !------- Loop over degrees and orders:
        DO mc=1,n_m_max   ! Numbers ms
            m=(mc-1)*minc
            sign=-1.D0

            b_r_n=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            b_t_n=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            b_p_n=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            b_r_s=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            b_t_s=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            b_p_s=CMPLX(0.D0,0.D0,KIND=KIND(0d0))

            DO l=m,l_max
                lm=lm2(l,m)
                sign=-sign

                b_r_1=         cs1(lm)*Plm(lm,n_theta_nhs)
                b_t_1=         cs2(lm)*dPlm(lm,n_theta_nhs)
                b_p_1=dPhi(lm)*cs2(lm)*Plm(lm,n_theta_nhs)

            !-------- Northern hemisphere:
                b_r_n=b_r_n+b_r_1
                b_t_n=b_t_n+b_t_1
                b_p_n=b_p_n+b_p_1

            !-------- Southern hemisphere:
                b_r_s=b_r_s+sign*b_r_1
                b_t_s=b_t_s-sign*b_t_1
                b_p_s=b_p_s+sign*b_p_1

            END DO  ! Loop over order

            b_r(2*mc-1,n_theta)  =REAL(b_r_n)
            b_r(2*mc,n_theta)    =AIMAG(b_r_n)
            b_t(2*mc-1,n_theta)  =REAL(b_t_n)
            b_t(2*mc,n_theta)    =AIMAG(b_t_n)
            b_p(2*mc-1,n_theta)  =REAL(b_p_n)
            b_p(2*mc,n_theta)    =AIMAG(b_p_n)
            b_r(2*mc-1,n_theta+1)=REAL(b_r_s)
            b_r(2*mc,n_theta+1)  =AIMAG(b_r_s)
            b_t(2*mc-1,n_theta+1)=REAL(b_t_s)
            b_t(2*mc,n_theta+1)  =AIMAG(b_t_s)
            b_p(2*mc-1,n_theta+1)=REAL(b_p_s)
            b_p(2*mc,n_theta+1)  =AIMAG(b_p_s)

        END DO     ! Loop over degree

        DO mc=1,2*n_m_max
            b_t(mc,n_theta)  =O_sint*b_t(mc,n_theta)
            b_p(mc,n_theta)  =O_sint*b_p(mc,n_theta)
            b_t(mc,n_theta+1)=O_sint*b_t(mc,n_theta+1)
            b_p(mc,n_theta+1)=O_sint*b_p(mc,n_theta+1)
        END DO
        DO mc=2*n_m_max+1,nrp
            b_r(mc,n_theta)  =0.D0
            b_t(mc,n_theta)  =0.D0
            b_p(mc,n_theta)  =0.D0
            b_r(mc,n_theta+1)=0.D0
            b_t(mc,n_theta+1)=0.D0
            b_p(mc,n_theta+1)=0.D0
        END DO

    END DO        ! Loop over colatitudes

!-- Transform m 2 phi:
    CALL fft_thetab(b_r,1)
    CALL fft_thetab(b_t,1)
    CALL fft_thetab(b_p,1)


    RETURN
    end SUBROUTINE get_B_surface


!-------------------------------------------------------------------------
