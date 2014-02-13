!$Id$
!***********************************************************************
SUBROUTINE get_ddcheb(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
     &                n_r_max,n_cheb_max,d_fac)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Returns chebychev coefficents of first derivative df and second  |
  !  |  derivative ddf for a function whose cheb-coeff. are given as     |
  !  |  columns in array f(n_c_tot,n_r_max).                             |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  implicit none

  !-- input:
  integer,INTENT(IN) :: n_f_start  ! No of column to start with
  integer,INTENT(IN) :: n_f_stop   ! No of column to stop with
  integer,INTENT(IN) :: n_f_max    ! First dimension of f,df,ddf
  integer,INTENT(IN) :: n_r_max    ! second dimension of f,df,ddf
  integer,INTENT(IN) :: n_cheb_max ! Number of cheb modes
  REAL(kind=8),INTENT(IN) ::  f(n_f_max,*)
  real(kind=8),INTENT(IN) ::  d_fac    ! factor for interval mapping

  !-- output:
  real(kind=8),INTENT(OUT) ::  df(n_f_max,*)
  real(kind=8),INTENT(OUT) ::  ddf(n_f_max,*)

  !-- local variables:
  integer :: n_f,n_cheb
  real(kind=8) :: fac_cheb

  !-- end of declaration
  !-----------------------------------------------------------------------

  !----- initialize derivatives:
  do n_cheb=n_cheb_max,n_r_max
     do n_f=n_f_start,n_f_stop
        df(n_f,n_cheb)=0.d0
        ddf(n_f,n_cheb)=0.d0
     end do
  end do
  n_cheb=n_cheb_max-1
  fac_cheb=d_fac*dble(2*n_cheb)
  do n_f=n_f_start,n_f_stop
     df(n_f,n_cheb)=fac_cheb*f(n_f,n_cheb+1)
     ddf(n_f,n_cheb)=0.d0
  end do

  !----- recursion
  do n_cheb=n_cheb_max-2,1,-1
     fac_cheb=d_fac*dble(2*n_cheb)
     do n_f=n_f_start,n_f_stop
        df(n_f,n_cheb)=df(n_f,n_cheb+2) + &
             fac_cheb*f(n_f,n_cheb+1)
        ddf(n_f,n_cheb)=ddf(n_f,n_cheb+2) + &
             fac_cheb*df(n_f,n_cheb+1)
     end do
  end do

  return
end subroutine get_ddcheb

!-----------------------------------------------------------------------
