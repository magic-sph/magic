!$Id$
!***********************************************************************
    subroutine get_Btor(Tlm,Bt,Bp, &
                        rT,n_theta_start,n_theta_block,lIC)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculated the components       |
!  |  Bt and Bp of the toroidal magnetic field  Tlm (in l,m space)     |
!  |  at the radial grid point r=rT and the block of theta grid points |
!  |  from n_theta=n_theta_start to                                    |
!  |       n_theta=n_theta_start+n_theta_block-1                       |
!  |  and for all phis.                                                |
!  |  For lIC=.true. the inner core field is calculated,               |
!  |  to get the IC field for a conducting inner core Plm has to be    |
!  |  the toroidal field at the ICB.                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif

    IMPLICIT NONE

!-- Input:

! include 'truncation.f'
! include 'c_blocking.f'
! include 'c_num_param.f'
! include 'c_radial.f'
! include 'c_horizontal.f'
! include 'c_logic.f'

    integer :: n_theta_start  ! first theta to be treated
    integer :: n_theta_block  ! last theta
    real(kind=8) ::  rT             ! radius
    logical :: lIC            ! true for inner core, special rDep !
    complex(kind=8) :: Tlm(lm_max) ! field in (l,m)-space for rT

!-- Output:
    real(kind=8) :: Bt(nrp,*),Bp(nrp,*)

!-- Local:
    integer :: lm,mc,m,l
    integer :: n_theta,n_theta_nhs
    real(kind=8) :: rRatio,rDep(0:l_max)
    real(kind=8) :: O_sint
    real(kind=8) :: sign
    complex(kind=8) :: cs1(lm_max)
    complex(kind=8) :: Bt_1,Bp_1
    complex(kind=8) :: Bt_n,Bp_n
    complex(kind=8) :: Bt_s,Bp_s

    complex(kind=8) :: ci


!-- end of declaration
!----------------------------------------------------------------------

    CI=CMPLX(0.D0,1.D0,KIND=KIND(0d0))

!-- Zero output field:
    do n_theta=1,n_theta_block
        do mc=1,nrp
            Bt(mc,n_theta)=0.d0
            Bp(mc,n_theta)=0.d0
        end do
    end do

    if ( lIC .AND. .NOT. l_cond_ic )    return

!-- Calculate radial dependencies: (r_ic(1)=r(n_r_max)=inner core boundary)
!   Radial dependence = (r/r_ICB)**(l-1) / r_ICB**2
    if ( lIC ) then
        rRatio=rT/r_icb
        rDep(0)=1.d0/r_icb
        do l=1,l_max
            rDep(l)=rDep(l-1)*rRatio
        end do
    else
        do l=0,l_max
            rDep(l)=1.d0/rT
        end do
    end if
!--- NOTE: mantle potential field may be included by adding special rDep
!          and using poloidal field at r_cmb


!-- Get coeffs with radial dependence:
    do m=0,m_max,minc
        do l=m,l_max
            lm=lm2(l,m)
            cs1(lm)=rDep(l)*Tlm(lm)
        end do
    end do

!-- Calculate northern and southern hemisphere components:

    do n_theta=1,n_theta_block,2 ! loop over thetas in northers HS

        n_theta_nhs=(n_theta_start+n_theta)/2
        O_sint     =osn1(n_theta_nhs)

        mc=0
        do m=0,m_max,minc   ! Numbers ms
            mc=mc+1
            sign=-1.d0

            Bt_n =CMPLX(0.d0,0.d0,KIND=KIND(0d0))
            Bp_n =CMPLX(0.d0,0.d0,KIND=KIND(0d0))
            Bt_s =CMPLX(0.d0,0.d0,KIND=KIND(0d0))
            Bp_s =CMPLX(0.d0,0.d0,KIND=KIND(0d0))

            do l=m,l_max
                lm=lm2(l,m)
                sign=-sign
                                 
                Bt_1=ci*dble(m)*cs1(lm)*Plm(lm,n_theta_nhs)
                Bp_1=          -cs1(lm)*dPlm(lm,n_theta_nhs)

            !-------- Northern hemisphere:
                Bt_n=Bt_n+Bt_1
                Bp_n=Bp_n+Bp_1

            !-------- Southern hemisphere:
                Bt_s=Bt_s+sign*Bt_1 ! Different equatoria symmetry
                Bp_s=Bp_s-sign*Bp_1

            end do  ! Loop over degree

        !-------- Divide theta and phi component by sin(theta):
            Bt_n =O_sint*Bt_n
            Bp_n =O_sint*Bp_n
            Bt_s =O_sint*Bt_s
            Bp_s =O_sint*Bp_s

        !-------- Copy to real array:
            Bt(2*mc-1,n_theta)   =REAL(Bt_n)
            Bt(2*mc  ,n_theta)   =AIMAG(Bt_n)
            Bp(2*mc-1,n_theta)   =REAL(Bp_n)
            Bp(2*mc  ,n_theta)   =AIMAG(Bp_n)
            Bt(2*mc-1,n_theta+1) =REAL(Bt_s)
            Bt(2*mc  ,n_theta+1) =AIMAG(Bt_s)
            Bp(2*mc-1,n_theta+1) =REAL(Bp_s)
            Bp(2*mc  ,n_theta+1) =AIMAG(Bp_s)

        end do     ! Loop over order

    !----- Zero remaining elements in array:
        do mc=2*n_m_max+1,nrp
            Bt(mc,n_theta)  =0.d0
            Bp(mc,n_theta)  =0.d0
            Bt(mc,n_theta+1)=0.d0
            Bp(mc,n_theta+1)=0.d0
        end do

    end do        ! Loop over colatitudes


!-- Transform m 2 phi:
    CALL fft_thetab(Bt,1)
    CALL fft_thetab(Bp,1)
            

    RETURN
    end subroutine get_Btor

!-------------------------------------------------------------------------------------
