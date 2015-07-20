!$Id$
module num_param
!****************************************************************
!  Module containing numerical and control parameters
!****************************************************************
   implicit none
 
   !-- Time step control:
   integer :: n_time_steps
   real(kind=8) :: alpha
   real(kind=8) :: dtstart,dtMin,dtMax
   real(kind=8) :: timeStart,tEND
 
   !-- Z-angular momentum at start of integration:
   real(kind=8) :: AMstart
 
   !-- Courant criteria:
   integer :: n_cour_step ! step for controlling  Courant criteria
   real(kind=8) :: courfac   ! input 
   real(kind=8) :: alffac    ! input 
   real(kind=8) :: intfac    ! input
   real(kind=8), allocatable :: delxr2(:),delxh2(:) ! ??
 
   !-- Hyperdiffusivity:
   integer :: ldif,ldifexp
   real(kind=8) :: difeta,difnu,difkap
 
   !-- Scalings:
   real(kind=8) :: tScale,lScale,vScale,pScale,eScale  ! scales
   real(kind=8) :: enscale         ! (input) scale for energies !
   integer :: n_tScale       ! controlls time scale
   integer :: n_lScale       ! controlls length scale
 
   !-- Stop signal:
   integer :: istop
 
   !-- Controlling run time:
   integer :: runTimeLimit(4),runTime(4),runTimeStart(4)!,resetTime(4)

contains
   subroutine initialize_num_param

      use truncation, only: n_r_max

      allocate( delxr2(n_r_max),delxh2(n_r_max) )

   end subroutine initialize_num_param
end module num_param
