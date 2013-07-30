!$Id$
!***********************************************************************
    SUBROUTINE get_br_v_bcs(br,vt,vp,omega,O_r_E_2,O_rho, &
                               n_theta_min,n_theta_block, &
                                       br_vt_lm,br_vp_lm)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculate the nonlinear term    |
!  |  of the magnetic boundary condition for a conducting mantle or    |
!  |  inner core in space (r,lm).                                      |
!  |  Calculation is performed for the theta block:                    |
!  |          n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        |
!  |  On input br, vt and vp are given on all phi points and           |
!  |  thetas in the specific block.                                    |
!  |  On output the contribution of these grid points to all           |
!  |  degree and orders is stored in br_vt_lm and br_vp_lm.            |
!  |  Output is [r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]            |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE horizontal_data
#if (FFTLIB==JW)
  USE fft_JW
#elif (FFTLIB==MKL)
  USE fft_MKL
#endif

    IMPLICIT NONE

!-- input:
! include 'truncation.f' ! contains dimensions
! include 'c_horizontal.f' ! contains O_sin_theta

    real(kind=8) :: br(nrp,*)       ! r**2 * B_r
    real(kind=8) :: vt(nrp,*)       ! r*sin(theta) U_theta
    real(kind=8) :: vp(nrp,*)       ! r*sin(theta) U_phi
    real(kind=8) :: omega           ! rotation rate of mantle or IC
    real(kind=8) :: O_r_E_2         ! 1/r**2
    real(kind=8) :: O_rho           ! 1/rho0 (anelastic)
    integer :: n_theta_min    ! start of theta block
    integer :: n_theta_block  ! size of theta_block

!-- output:
    complex(kind=8) :: br_vt_lm(lmP_max) ! br*vt/(sin(theta)**2*r**2)
    complex(kind=8) :: br_vp_lm(lmP_max) ! br*(vp/(sin(theta)**2*r**2)-omega_ma)

!-- local variables:
    integer :: n_theta     ! number of theta position
    integer :: n_theta_rel ! number of theta position in block
    integer :: n_phi       ! number of longitude
    real(kind=8) :: br_vt(nrp,n_theta_block)
    real(kind=8) :: br_vp(nrp,n_theta_block)
    real(kind=8) :: fac          ! 1/( r**2 sin(theta)**2 )

!-- end of declaration
!-----------------------------------------------------------------------


    n_theta=n_theta_min-1 ! n_theta needed for O_sin_theta_E_2

    do n_theta_rel=1,n_theta_block
        n_theta=n_theta+1           ! absolute number of theta

        fac=O_sin_theta(n_theta)*O_sin_theta(n_theta)*O_r_E_2*O_rho

        do n_phi=1,n_phi_max

            br_vt(n_phi,n_theta_rel)=  fac * &
                 br(n_phi,n_theta_rel)*vt(n_phi,n_theta_rel)

            br_vp(n_phi,n_theta_rel)=    &
                 br(n_phi,n_theta_rel) * &
                 ( fac*vp(n_phi,n_theta_rel) - omega )

        end do

    end do

!-- Fourier transform phi 2 m (real 2 complex!)
    CALL fft_thetab(br_vt,-1)
    CALL fft_thetab(br_vp,-1)

!-- Legendre transform contribution of thetas in block:
    call legTF2(n_theta_min,br_vt_lm,br_vp_lm,br_vt,br_vp)
            
    return
    end SUBROUTINE get_br_v_bcs

!----------------------------------------------------------------------------
