!$Id$
!***********************************************************************
    subroutine get_dcheb(f,df,n_f_max,n_f_start,n_f_stop, &
                         n_r_max,n_cheb_max,d_fac)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Returns chebychev coeffitiens of first derivative df and second  |
!  |  derivative ddf for a function whose cheb-coeff. are given as     |
!  |  columns in array f(n_f_max,n_r_max).                             |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+
      
    implicit none

!-- input:
    integer :: n_f_start  ! No of function to start with
    integer :: n_f_stop   ! No of function to stop with
    integer :: n_f_max    ! Max no of functions
    integer :: n_r_max    ! second dimension of f,df,ddf
    integer :: n_cheb_max ! Number of cheb modes
    real(kind=8) ::  f(n_f_max,*)
    real(kind=8) ::  d_fac    ! factor for intervall mapping

!-- output:
    real(kind=8) ::  df(n_f_max,*)

!-- local variables:
    integer :: n_f,n_cheb
    real(kind=8) :: fac_cheb

!-- end of declaration
!-----------------------------------------------------------------------

!----- initialize derivatives:
    DO n_cheb=n_cheb_max,n_r_max
        DO n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=0.d0
        END DO
    END DO
    n_cheb  =n_cheb_max-1
    fac_cheb=d_fac*DBLE(2*n_cheb)
    DO n_f=n_f_start,n_f_stop
        df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
    END DO

!----- recursion
    DO n_cheb=n_cheb_max-2,1,-1
        fac_cheb=d_fac*dble(2*n_cheb)
        DO n_f=n_f_start,n_f_stop
            df(n_f,n_cheb)=df(n_f,n_cheb+2) + &
                           fac_cheb*f(n_f,n_cheb+1)
        END DO
    END DO

    RETURN
    end subroutine get_dcheb

!-----------------------------------------------------------------------
