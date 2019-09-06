module constants
   !
   ! module containing constants and parameters
   ! used in the code.
   !

   use precision_mod

   implicit none

   real(cp) :: c_z10_omega_ic,c_z10_omega_ma
   real(cp) :: c_dt_z10_ic,c_dt_z10_ma
   real(cp) :: c_lorentz_ic,c_lorentz_ma
   real(cp) :: vol_ic ! Volume of the inner core
   real(cp) :: vol_oc ! Volume of the outer core
   real(cp) :: surf_cmb ! Outer boundary surface
   real(cp), parameter :: pi=4.0_cp*atan(1.0_cp) ! pi=3.1415..
   real(cp), parameter :: sq4pi=sqrt(4.0_cp*pi)  ! :math:`\sqrt{4\pi}`
   real(cp), parameter :: osq4pi=1.0_cp/sq4pi    ! :math:`1/\sqrt{4\pi}`
   real(cp) :: c_moi_ic ! Moment of inertia of the inner core
   real(cp) :: c_moi_oc ! Moment of inertia of the outer core
   real(cp) :: c_moi_ma ! Moment of inertia of the mantle
   real(cp) :: mass ! Mass of the outer core
   real(cp) :: y10_norm,y11_norm

   real(cp), parameter :: one  =1.0_cp ! 1
   real(cp), parameter :: two  =2.0_cp ! 2
   real(cp), parameter :: three=3.0_cp ! 3
   real(cp), parameter :: four =4.0_cp ! 4
   real(cp), parameter :: half =0.5_cp ! 0.5
   real(cp), parameter :: third=one/three ! 1/3
   complex(cp), parameter :: zero=(0.0_cp,0.0_cp) ! cmplx(0.0, 0.0)
   complex(cp), parameter :: ci=(0.0_cp,1.0_cp) ! cmplx(0.0, 1.0)

   real(cp), parameter :: sin36=sin(36.0_cp*pi/180.0_cp) ! :math:`\sin{36\pi/180}`
   real(cp), parameter :: sin60=0.5_cp*sqrt(3.0_cp)      ! :math:`\sin{60\pi/180}`
   real(cp), parameter :: sin72=sin(72.0_cp*pi/180.0_cp) ! :math:`\sin{72\pi/180}`
   real(cp), parameter :: cos36=cos(36.0_cp*pi/180.0_cp) ! :math:`\cos{36\pi/180}`
   real(cp), parameter :: cos72=cos(72.0_cp*pi/180.0_cp) ! :math:`\cos{72\pi/180}`

   character(len=4), parameter :: codeVersion='5.8'

end module constants
