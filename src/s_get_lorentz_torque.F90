!$Id$
!***********************************************************************
    SUBROUTINE get_lorentz_torque(lorentz_torque, &
                          nThetaStart,sizeThetaB, &
                                        br,bp,nR)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this subroutine is to calculate the lorentz torque    |
!  |  on mantle or inner core respectively.                            |
!  |  Blocking in theta can be used to increased performance.          |
!  |  If no blocking required set n_theta_block=n_theta_max,           |
!  |  where n_theta_max is the absolut number of thetas used.          |
!  |  Note: lorentz_torque must be set to zero before loop over        |
!  |        theta blocks is started.                                   |
!  |  WARNING: subroutine returns -lorentz_torque if used at CMB       |
!  |        to calculate torque on mantle because if the inward        |
!  |        surface normal vector.
!  |  The Prandtl number is always the Prandtl number of the outer     |
!  |  core. This comes in via scaling of the magnetic field.           |
!  |  Theta alternates between northern and southern hemisphere in     |
!  |  br and bp but not in gauss. This has to be cared for, and we     |
!  |  use: gauss(latitude)=gauss(-latitude) here.                      |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE Grenoble
    USE horizontal_data

    IMPLICIT NONE

    REAL(kind=8) :: lorentz_torque  ! lorentz_torque for theta(1:n_theta)
    INTEGER :: nThetaStart    ! first number of theta in block
    INTEGER :: sizeThetaB     ! size of theta bloching
    REAL(kind=8) :: br(nrp,*)       ! array containing
! r**2 B_r(1:n_phi_max,1:n_theta_block)
    REAL(kind=8) :: bp(nrp,*)       ! array containing
! r*sin(theta) B_phi(1:n_phi_max,1:n_theta_block)
    INTEGER :: nR

!-- output:
!    lorentz_torque increased by the contribution of the theta block

!-- local variables:
    INTEGER :: nTheta,nPhi,nThetaNHS
    INTEGER :: nThetaB
    !REAL(kind=8) :: lorentz_torque_local
    REAL(kind=8) :: fac,b0r

!-- end of declaration
!-----------------------------------------------------------------------

    ! to avoid rounding errors for different theta blocking, we do not
    ! calculate sub sums with lorentz_torque_local, but keep on adding
    ! the contributions to the total lorentz_torque given as argument.

    IF ( nThetaStart == 1 ) THEN
       lorentz_torque=0.D0
    END IF

    !lorentz_torque_local=0.D0
    fac=8.D0*DATAN(1.D0)/DBLE(n_phi_max) ! 2 pi/n_phi_max

    nTheta=nThetaStart-1
    DO nThetaB=1,sizeThetaB
        nTheta=nTheta+1
        nThetaNHS=(nTheta+1)/2 ! northern hemisphere=odd n_theta
        IF ( lGrenoble ) THEN
            IF ( r(nR) == r_icb ) THEN
                b0r=2.D0*BIC*r_icb**2*cosTheta(nTheta)
            ELSE IF ( r(nR) == r_cmb ) THEN
                b0r=2.D0*BIC*r_icb**2*cosTheta(nTheta)*(r_icb/r_cmb)
            END IF
        ELSE
            b0r=0.D0
        END IF

        DO nPhi=1,n_phi_max
           !lorentz_torque_local=lorentz_torque_local + &
           !                          gauss(nThetaNHS) * &
           !       (br(nPhi,nThetaB)-b0r)*bp(nPhi,nThetaB)
           lorentz_torque=lorentz_torque + fac * gauss(nThetaNHS) * &
                  (br(nPhi,nThetaB)-b0r)*bp(nPhi,nThetaB)
        END DO
        !lorentz_torque_local = lorentz_torque_local + gauss(nThetaNHS)*phisum
    END DO

!-- normalisation of phi-integration and division by Pm:
    !lorentz_torque=lorentz_torque+fac*lorentz_torque_local
            
    RETURN
    end SUBROUTINE get_lorentz_torque
!-----------------------------------------------------------------------
