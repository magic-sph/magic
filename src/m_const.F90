!$Id$
!*********************************************************************
!  module containing constants used for rotation of inner
!  core and mantle.
!  See subroutine s_get_new_constants for explanations.
!*********************************************************************

!    !------------ This is release 1 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

MODULE const

  IMPLICIT NONE

  REAL(kind=8) :: c_z10_omega_ic,c_z10_omega_ma
  REAL(kind=8) :: c_dt_z10_ic,c_dt_z10_ma
  REAL(kind=8) :: c_lorentz_ic,c_lorentz_ma
  REAL(kind=8) :: vol_ic,vol_oc,surf_cmb
  REAL(kind=8) :: pi
  REAL(kind=8) :: c_moi_ic,c_moi_ma,c_moi_oc
  REAL(kind=8) :: mass
  REAL(kind=8) :: y10_norm,y11_norm
  COMPLEX(kind=8),PARAMETER :: zero=(0.0D0,0.0D0)

  !---------------------------------------------------------------------
CONTAINS
  SUBROUTINE initialize_const
    
    pi=4.d0*DATAN(1.D0)
    
  END SUBROUTINE initialize_const
END MODULE const
