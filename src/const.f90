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
   complex(kind=8), parameter :: zero=(0.0D0,0.0D0)

   real(kind=8), parameter :: sin36=dsin(36.d0*pi/180.d0)
   real(kind=8), parameter :: sin60=0.5d0*dsqrt(3.d0)
   real(kind=8), parameter :: sin72=dsin(72.d0*pi/180.d0)
   real(kind=8), parameter :: cos36=dcos(36.d0*pi/180.d0)
   real(kind=8), parameter :: cos72=dcos(72.d0*pi/180.d0)

end module const
