!$Id$
!***********************************************************************
SUBROUTINE get_fl(fl,n_r,n_theta_start,n_theta_block,l_ic)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. -----------

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  Return field fl whose contourlines are the fields lines          |
  !  |  of the axisymmetric poloidal mangetic field.                     |
  !  |    fl(r,theta)=d_theta b(r,theta,m=0)/r                           |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  !  This routine is called for l_ic=.true. only from rank 0 with full
  !  field b_ic in standard lm ordering available.
  !  The case l_ic=.false. is called from all ranks and uses b_Rloc.

  USE truncation
  USE radial_functions
  USE Grenoble
  USE blocking
  USE horizontal_data
  USE logic
  USE fields

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n_r             ! No. of radial grid point
  integer,INTENT(IN) :: n_theta_start   ! No. of theta to start with
  integer,INTENT(IN) :: n_theta_block   ! Size of theta block
  logical,INTENT(IN) :: l_ic            ! =true if inner core field

  !-- output:
  real(kind=8),INTENT(OUT) ::  fl(*)    ! Field for field lines

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
     r_ratio =r_ic(n_r)/r_ic(1)
     r_dep(1)=r_ratio/r_ic(1)
     do l=2,l_max
        r_dep(l)=r_dep(l-1)*r_ratio
     end do
  else
     O_r=or1(n_r)
  end if

  !----- Loop over colatitudes:

  do n_theta=1,n_theta_block,2

     n_theta_nhs=(n_theta_start+n_theta)/2
     O_sint=osn1(n_theta_nhs)

     !------- Loop over degrees and orders:

     sign=1.d0
     fl_n=0.d0
     fl_s=0.d0

     lm=1
     do l=1,l_max
        lm=lm+1
        sign=-sign

        IF ( l_ic ) THEN ! Inner Core
           IF ( l_cond_ic ) then
              fl_1=r_dep(l) * &
                   REAL(b_ic(lm,n_r))*dPlm(lm,n_theta_nhs)
           ELSE
              !                   IF ( lGrenoble .AND. l == 1 ) THEN
              !                       fl_1=r_dep(l) * dPlm(lm,n_theta_nhs) *
              !     *                       (REAL(b(lm,n_r_icb))+b0(n_r_icb))
              !                   ELSE
              fl_1=r_dep(l) * dPlm(lm,n_theta_nhs) * REAL(b(lm,n_r_icb))
              !                    END IF
           END IF
        ELSE             ! Outer Core
           !                 IF ( lGrenoble .AND. l == 1 ) THEN
           !                   fl_1=O_r*dPlm(lm,n_theta_nhs) *
           !     *                       (REAL(b(lm,n_r))+b0(n_r))
           !                ELSE
           fl_1=O_r*dPlm(lm,n_theta_nhs) * REAL(b_Rloc(lm,n_r))
           !                END IF
        END IF

        !-------- Northern hemisphere:
        fl_n=fl_n+fl_1

        !-------- Southern hemisphere:
        fl_s=fl_s-sign*fl_1

     end do  ! Loop over order

     !-- Devide by sin(theta):
     fl(n_theta)  =-O_sint*fl_n
     fl(n_theta+1)=-O_sint*fl_s

  end do        ! Loop over colatitudes


  return
end subroutine get_fl
!-------------------------------------------------------------------------
