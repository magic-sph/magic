!$Id$
!***********************************************************************
SUBROUTINE get_viscous_torque(viscous_torque,z10,dz10,r)
  !***********************************************************************

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this subroutine is to calculate the viscous torque    |
  !  |  on mantle or inner core respectively.                            |
  !  |  NOTE: sign is wrong for torque on mantle!                        |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  IMPLICIT NONE

  !-- Input:
  COMPLEX(kind=8),INTENT(IN) :: z10,dz10    ! z10 coefficient and its radial deriv.
  ! at ICB or CMB
  REAL(kind=8),INTENT(IN) :: r               ! radius (ICB or CMB)

  !-- Output:
  REAL(kind=8),INTENT(OUT) :: viscous_torque

  !-- Local:
  REAL(kind=8) :: pi

  !-- end of declaration
  !-----------------------------------------------------------------------

  pi=4.D0*DATAN(1.D0)
  viscous_torque=-4.D0*DSQRT(pi/3.D0)*r *( 2.D0*REAL(z10) - r*REAL(dz10) )

  RETURN
end SUBROUTINE get_viscous_torque

!-----------------------------------------------------------------------
