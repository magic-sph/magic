!$Id$
!***********************************************************************
    subroutine get_sl(sl,n_r,n_theta_start,n_theta_block)
!***********************************************************************

!------------ This is release 2 level 1  --------------!
!------------ Created on 1/17/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  Return field sl whose contourlines are the stream lines          |
!  |  of the axisymmetric poloidal velocity field.                     |
!  |    sl(r,theta)=d_theta v(r,theta,m=0)/r                           |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE logic
    USE fields

    IMPLICIT NONE

    integer :: n_r             ! No. of radial grid point
    integer :: n_theta_start   ! No. of theta to start with
    integer :: n_theta_block   ! Size of theta block

!-- output:
    real(kind=8) ::  sl(*)           ! Field for field lines

!-- local:
    integer :: n_theta         ! No. of theta
    integer :: n_theta_nhs     ! Counter for thetas in north HS
    integer :: l,lm            ! Degree, counter for degree/order combinations

    real(kind=8) :: sign
    real(kind=8) :: O_r              ! 1/r
    real(kind=8) :: O_sint           ! 1/sin(theta)
    real(kind=8) :: sl_s,sl_n,sl_1


!-- end of declaration
!---------------------------------------------------------------------


!-- Radial dependence:


!-- Calculate radial dependencies:
    O_r=or1(n_r)

!----- Loop over colatitudes:

    do n_theta=1,n_theta_block,2

        n_theta_nhs=(n_theta_start+n_theta)/2
        O_sint=osn1(n_theta_nhs)

    !------- Loop over degrees and orders:

        sign=1.d0
        sl_n=0.d0
        sl_s=0.d0

        lm=1
        do l=1,l_max
            lm=lm+1
            sign=-sign
            sl_1=O_r*REAL(w(lm,n_r))*dPlm(lm,n_theta_nhs)
        !-------- Northern hemisphere:
            sl_n=sl_n+sl_1
        !-------- Southern hemisphere:
            sl_s=sl_s-sign*sl_1
        end do  ! Loop over order

    !-- Devide by sin(theta):
        sl(n_theta)  =O_sint*sl_n
        sl(n_theta+1)=O_sint*sl_s

    end do        ! Loop over colatitudes
     
    return
    end subroutine get_sl
!-------------------------------------------------------------------------
