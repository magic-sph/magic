!$Id$
!***********************************************************************
SUBROUTINE get_dtB(dtB,dtBLM,DimB1,DimB2,n_r, &
     &             n_theta_start,n_theta_block,l_ic)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  USE truncation
  USE radial_functions
  USE blocking
  USE horizontal_data
  USE logic

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n_r             ! No. of radial grid point
  integer,INTENT(IN) :: n_theta_start   ! No. of theta to start with
  integer,INTENT(IN) :: n_theta_block   ! Size of theta block
  logical,INTENT(IN) :: l_ic            ! =true if inner core field

  integer,INTENT(IN) :: DimB1,DimB2
  complex(kind=8),INTENT(IN) :: dtBLM(DimB1,DimB2)

  !-- output:
  real(kind=8),INTENT(OUT) ::  dtB(*)           ! Result Field with dim>=n_theta_block

  !-- local:
  integer :: n_theta         ! No. of theta
  integer :: n_theta_nhs     ! Counter for thetas in north HS
  integer :: l,lm            ! Degree, counter for degree/order combinations

  real(kind=8) :: sign
  real(kind=8) :: r_ratio          ! r/r_ICB
  real(kind=8) :: O_r              ! 1/r
  real(kind=8) :: O_sint           ! 1/sin(theta)
  real(kind=8) :: r_dep(l_max)     ! (r/r_ICB)**l / r_ICB
  real(kind=8) :: fl_s,fl_n,fl_1


  !-- end of declaration
  !---------------------------------------------------------------------


  !-- Radial dependence:


  !-- Calculate radial dependencies:
  !     for IC: (r/r_ICB)**l / r_ICB
  !     for OC: 1/r
  if ( l_ic ) then
     r_ratio =r_ic(n_r)/r_icb
     r_dep(1)=r_ratio/r_icb
     do l=2,l_max
        r_dep(l)=r_dep(l-1)*r_ratio
     end do
  else
     O_r=or1(n_r)
  end if


  !----- Loop over colatitudes:

  do n_theta=1,n_theta_block,2

     n_theta_nhs=(n_theta_start+n_theta)/2

     !------- Loop over degrees and orders:

     sign=1.d0
     fl_n=0.d0
     fl_s=0.d0

     lm=1
     do l=1,l_max
        lm=lm+1
        sign=-sign

        if ( l_ic ) then ! Inner Core
           fl_1=r_dep(l)*REAL(dtBLM(lm,n_r))*dPlm(lm,n_theta_nhs)
        else             ! Outer Core
           fl_1=O_r*REAL(dtBLM(lm,n_r))*dPlm(lm,n_theta_nhs)
        end if

        !-------- Northern hemisphere:
        fl_n=fl_n+fl_1

        !-------- Southern hemisphere:
        fl_s=fl_s-sign*fl_1

     end do  ! Loop over order

     O_sint=osn1(n_theta_nhs)
     dtB(n_theta)  =-O_sint*fl_n
     dtB(n_theta+1)=-O_sint*fl_s

  end do        ! Loop over colatitudes


  return
end SUBROUTINE get_dtB
!-------------------------------------------------------------------------
