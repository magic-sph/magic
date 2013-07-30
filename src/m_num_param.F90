!$Id$
!****************************************************************
!  Module containing numerical and control parameters
!****************************************************************

!   !------------ This is release 1 level 1  --------------!
!   !------------ Created on 1/17/02  by JW. --------------!


MODULE num_param
  use truncation
  implicit none

  !-- Time step control:
  INTEGER :: n_time_steps
  REAL(kind=8) :: alpha
  REAL(kind=8) :: dtstart,dtMin,dtMax
  REAL(kind=8) :: timeStart,tEND

  !-- Z-angular momentum at start of integration:
  REAL(kind=8) :: AMstart

  !-- Courant criteria:
  INTEGER :: n_cour_step ! step for controlling  Courant criteria
  REAL(kind=8) :: courfac   ! input 
  REAL(kind=8) :: alffac    ! input 
  REAL(kind=8) :: intfac    ! input
  REAL(kind=8),ALLOCATABLE :: delxr2(:),delxh2(:) ! ??

  !-- Hyperdiffusivity:
  INTEGER :: ldif,ldifexp
  REAL(kind=8) :: difeta,difnu,difkap


  !-- Scalings:
  REAL(kind=8) :: tScale,lScale,vScale,pScale,eScale  ! scales
  REAL(kind=8) :: enscale         ! (input) scale for energies !
  INTEGER :: n_tScale       ! controlls time scale
  INTEGER :: n_lScale       ! controlls length scale

  !-- Stop signal:
  INTEGER :: istop

  !-- Controlling run time:
  INTEGER :: runTimeLimit(4),runTime(4),runTimeStart(4),resetTime(4)

CONTAINS
  SUBROUTINE initialize_num_param
    allocate( delxr2(n_r_max),delxh2(n_r_max) )

  END SUBROUTINE initialize_num_param
END MODULE num_param
