!$Id$
module num_param
   !--------------------------------------------------------------
   !  Module containing numerical and control parameters
   !--------------------------------------------------------------

   use truncation, only: n_r_max
   use precision_mod

   implicit none

   private
 
   !-- Time step control:
   integer, public :: n_time_steps
   real(cp), public :: alpha
   real(cp), public :: dtstart,dtMin,dtMax
   real(cp), public :: timeStart,tEND
 
   !-- Z-angular momentum at start of integration:
   real(cp), public :: AMstart
 
   !-- Courant criteria:
   integer, public :: n_cour_step ! step for controlling  Courant criteria
   real(cp), public :: courfac   ! input 
   real(cp), public :: alffac    ! input 
   real(cp), public :: intfac    ! input
   real(cp), public, allocatable :: delxr2(:),delxh2(:) ! ??
 
   !-- Hyperdiffusivity:
   integer, public :: ldif,ldifexp
   real(cp), public :: difeta,difnu,difkap
 
   !-- Scalings:
   real(cp), public :: tScale,lScale,vScale,pScale,eScale  ! scales
   real(cp), public :: enscale         ! (input) scale for energies !
   integer, public :: n_tScale       ! controlls time scale
   integer, public :: n_lScale       ! controlls length scale
 
   !-- Stop signal:
   integer, public :: istop
 
   !-- Controlling run time:
   integer, public :: runTimeLimit(4),runTime(4),runTimeStart(4)!,resetTime(4)

   public :: initialize_num_param

contains

   subroutine initialize_num_param

      allocate( delxr2(n_r_max),delxh2(n_r_max) )

   end subroutine initialize_num_param
!-------------------------------------------------------------------------------
end module num_param
