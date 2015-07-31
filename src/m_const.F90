!$Id$
module const
   !---------------------------------------------------------------
   !  module containing constants used for rotation of inner
   !  core and mantle.
   !  See subroutine s_get_new_constants for explanations.
   !---------------------------------------------------------------

   implicit none
 
   real(kind=8) :: c_z10_omega_ic,c_z10_omega_ma
   real(kind=8) :: c_dt_z10_ic,c_dt_z10_ma
   real(kind=8) :: c_lorentz_ic,c_lorentz_ma
   real(kind=8) :: vol_ic,vol_oc,surf_cmb
   real(kind=8), parameter :: pi=4.d0*datan(1.D0)
   real(kind=8) :: c_moi_ic,c_moi_ma,c_moi_oc
   real(kind=8) :: mass
   real(kind=8) :: y10_norm,y11_norm
   complex(kind=8),parameter :: zero=(0.0D0,0.0D0)

end module const
