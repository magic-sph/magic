!$Id$
module const
   !---------------------------------------------------------------
   !  module containing constants used for rotation of inner
   !  core and mantle.
   !  See subroutine s_get_new_constants for explanations.
   !---------------------------------------------------------------

   use precision_mod

   implicit none
 
   real(cp) :: c_z10_omega_ic,c_z10_omega_ma
   real(cp) :: c_dt_z10_ic,c_dt_z10_ma
   real(cp) :: c_lorentz_ic,c_lorentz_ma
   real(cp) :: vol_ic,vol_oc,surf_cmb
   real(cp), parameter :: pi=4.0_cp*atan(1.0_cp)
   real(cp), parameter :: sq4pi=sqrt(4.0_cp*pi)
   real(cp), parameter :: osq4pi=1.0_cp/sq4pi
   real(cp) :: c_moi_ic,c_moi_ma,c_moi_oc
   real(cp) :: mass
   real(cp) :: y10_norm,y11_norm

   real(cp), parameter :: one  =1.0_cp
   real(cp), parameter :: two  =2.0_cp
   real(cp), parameter :: three=3.0_cp
   real(cp), parameter :: four =4.0_cp
   real(cp), parameter :: half =0.5_cp
   real(cp), parameter :: third=one/three
   complex(cp), parameter :: zero=(0.0_cp,0.0_cp)
   complex(cp), parameter :: ci=(0.0_cp,1.0_cp)

   real(cp), parameter :: sin36=sin(36.0_cp*pi/180.0_cp)
   real(cp), parameter :: sin60=0.5_cp*sqrt(3.0_cp)
   real(cp), parameter :: sin72=sin(72.0_cp*pi/180.0_cp)
   real(cp), parameter :: cos36=cos(36.0_cp*pi/180.0_cp)
   real(cp), parameter :: cos72=cos(72.0_cp*pi/180.0_cp)

   character(len=4), parameter :: codeVersion='5.1'

end module const
