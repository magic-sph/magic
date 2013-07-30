!$Id$
!***********************************************************************
    SUBROUTINE get_drNS_even(f,df,n_f_max,n_f_start,n_f_stop, &
                             n_r_max,n_cheb_max,dr_fac,work1, &
                                 i_costf1_init,d_costf1_init, &
                                   i_costf2_init,d_costf2_init)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Returns first rarial derivative df and second radial             |
!  |  derivative ddf of the input function f.                          |
!  |  Array f(n_f_max,*) may contain several functions numbered by     |
!  |  the first index. The subroutine calculates the derivaties of     |
!  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
!  |  to a Chebychev representation using n_r_max radial grid points.  |
!  |  The cheb transforms have to be initialized by calling            |
!  |   init_costf1 and init_costf2.                                    |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    IMPLICIT NONE

!-- input:
    INTEGER :: n_f_max           ! first dim of f
    REAL(kind=8) :: f(n_f_max,*)
    INTEGER :: n_f_start         ! first function to be treated
    INTEGER :: n_f_stop          ! last function to be treated
    INTEGER :: n_r_max           ! number of radial grid points
    INTEGER :: n_cheb_max        ! number of cheb modes
    REAL(kind=8) :: dr_fac             ! mapping factor
    REAL(kind=8) :: work1(n_f_max,*)   ! work array needed for costf
    INTEGER :: i_costf1_init(*)
    INTEGER :: i_costf2_init(*)
    REAL(kind=8) ::  d_costf1_init(*)
    REAL(kind=8) ::  d_costf2_init(*)

!-- output:
    REAL(kind=8) :: df(n_f_max,*)  ! first derivative of f

!-- end of declaration
!-----------------------------------------------------------------------


!-- Transform f to cheb space:
    CALL costf1(f,n_f_max,n_f_start,n_f_stop, &
                work1,i_costf1_init,d_costf1_init)

!-- Get derivatives:
    CALL get_dcheb_even(f,df,n_f_max,n_f_start,n_f_stop, &
                        n_r_max,n_cheb_max,dr_fac)

!-- Transform back, note the different transform used for df,
!   cause df is odd:
    CALL costf1(f,n_f_max,n_f_start,n_f_stop, &
                work1,i_costf1_init,d_costf1_init)
    CALL costf2(df,n_f_max,n_f_start,n_f_stop, &
                work1,i_costf2_init,d_costf2_init,1)


    RETURN
    end SUBROUTINE get_drNS_even

!-----------------------------------------------------------------------
