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
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
            
    IMPLICIT NONE

!-- Input:
    COMPLEX(kind=8) :: z10,dz10    ! z10 coefficient and its radial deriv.
! at ICB or CMB
    REAL(kind=8) :: r               ! radius (ICB or CMB)

!-- Output:
    REAL(kind=8) :: viscous_torque

!-- Local:
    REAL(kind=8) :: pi

!-- end of declaration
!-----------------------------------------------------------------------

    pi=4.D0*DATAN(1.D0)
    viscous_torque=-4.D0*DSQRT(pi/3.D0)*r *( 2.D0*REAL(z10) - r*REAL(dz10) )
            
    RETURN
    end SUBROUTINE get_viscous_torque

!-----------------------------------------------------------------------
