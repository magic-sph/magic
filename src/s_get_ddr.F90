!$Id$
!***********************************************************************
SUBROUTINE get_ddr(f,df,ddf,n_f_max,n_f_start,n_f_stop, &
     n_r_max,n_cheb_max,work1,work2, &
     i_costf_init,d_costf_init,drx,ddrx)
  !***********************************************************************

  !------------ This is release 2 level 1  --------------!
  !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Returns first radial derivative df and second radial             |
  !  |  derivative ddf of the input function f.                          |
  !  |  Array f(n_f_max,*) may contain several functions numbered by     |
  !  |  the first index. The subroutine calculates the derivatives of    |
  !  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
  !  |  to a Chebychev representation using n_r_max radial grid points.  |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  IMPLICIT NONE

  !-- input:
  INTEGER,intent(IN) :: n_f_max         ! first dim of f
  REAL(kind=8),intent(IN) :: f(n_f_max,*)
  INTEGER,intent(IN) :: n_f_start       ! first function to be treated
  INTEGER,intent(IN) :: n_f_stop        ! last function to be treated
  INTEGER,intent(IN) :: n_r_max         ! number of radial grid points
  INTEGER,intent(IN) :: n_cheb_max      ! number of cheb modes
  REAL(kind=8),intent(OUT) :: work1(n_f_max,*)  ! work array needed for costf
  REAL(kind=8),intent(OUT) :: work2(n_f_max,*)  ! work array for f transfer
  REAL(kind=8),intent(IN) :: drx(*)           ! first derivatives of x(r)
  REAL(kind=8),intent(IN) :: ddrx(*)          ! second derivatives of x(r)

  INTEGER,intent(IN) :: i_costf_init(*) ! info for costf
  REAL(kind=8),intent(IN) :: d_costf_init(*)  ! info for costf

  !-- output:
  REAL(kind=8),intent(OUT) :: df(n_f_max,*)  ! first derivative of f
  REAL(kind=8),intent(OUT) :: ddf(n_f_max,*) ! second derivative of f

  !-- Local:
  INTEGER :: n_r,n_f

  !-- end of declaration
  !-----------------------------------------------------------------------


  !-- Copy input functions:
  DO n_r=1,n_r_max
     DO n_f=n_f_start,n_f_stop
        work2(n_f,n_r)=f(n_f,n_r)
     END DO
  END DO

  !-- Transform f to cheb space:
  CALL costf1(work2,n_f_max,n_f_start,n_f_stop, &
       work1,i_costf_init,d_costf_init)

  !-- Get derivatives:
  CALL get_ddcheb(work2,df,ddf, &
       n_f_max,n_f_start,n_f_stop, &
       n_r_max,n_cheb_max,1.D0)

  !-- Transform back:
  CALL costf1(df,n_f_max,n_f_start,n_f_stop, &
       work1,i_costf_init,d_costf_init)
  CALL costf1(ddf,n_f_max,n_f_start,n_f_stop, &
       work1,i_costf_init,d_costf_init)

  !-- New map:
  DO n_r=1,n_r_max
     DO n_f=n_f_start,n_f_stop
        ddf(n_f,n_r)=   ddrx(n_r)*df(n_f,n_r) + &
             drx(n_r)**2*ddf(n_f,n_r)
        df(n_f,n_r) =    drx(n_r)*df(n_f,n_r)
     END DO
  END DO


  RETURN
end SUBROUTINE get_ddr

!-----------------------------------------------------------------------
