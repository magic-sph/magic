!$Id$
!***********************************************************************
SUBROUTINE get_dr(f,df,n_f_max,n_f_start,n_f_stop, &
     &            n_r_max,n_cheb_max,work1,work2, &
     &            i_costf_init,d_costf_init,drx)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Returns first radial derivative df of the input function f.      |
  !  |  Array f(n_f_max,*) may contain several functions numbered by     |
  !  |  the first index. The subroutine calculates the derivaties of     |
  !  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
  !  |  to a Chebychev representation using n_r_max radial grid points . |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  !USE communications,ONLY: get_global_sum

  IMPLICIT NONE

  !-- input:
  INTEGER,INTENT(IN) :: n_f_max          ! first dim of f
  REAL(kind=8),INTENT(IN) :: f(n_f_max,*)
  INTEGER,INTENT(IN) :: n_f_start        ! first function to be treated
  INTEGER,INTENT(IN) :: n_f_stop         ! last function to be treated
  INTEGER,INTENT(IN) :: n_r_max          ! number of radial grid points
  INTEGER,INTENT(IN) :: n_cheb_max       ! max number of cheb modes
  REAL(kind=8),INTENT(OUT) :: work1(n_f_max,*)  ! work array needed for costf
  REAL(kind=8),INTENT(OUT) :: work2(n_f_max,n_r_max)  ! work array for f transfer
  REAL(kind=8),INTENT(IN) :: drx(*)            ! first derivatives of x(r)

  INTEGER,INTENT(IN) :: i_costf_init(*)  ! info for costf
  REAL(kind=8),INTENT(IN) :: d_costf_init(*)   ! info for costf

  !-- Output:
  REAL(kind=8),INTENT(OUT) :: df(n_f_max,*)  ! first derivative of f

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

  !WRITE(*,"(A,2ES22.10)") "work2=",get_global_sum(work2)
  !-- Transform f to cheb space:
  CALL costf1(work2,n_f_max,n_f_start,n_f_stop, &
       work1,i_costf_init,d_costf_init)
  !WRITE(*,"(A,2ES22.10)") "work2 in cheb=",get_global_sum(work2)

  !-- Get derivatives:
  CALL get_dcheb(work2,df,n_f_max,n_f_start,n_f_stop, &
       n_r_max,n_cheb_max,1.D0)
  !WRITE(*,"(A,2ES22.10)") "dcheb df=",get_global_sum(df)

  !-- Transform back:
  CALL costf1(df,n_f_max,n_f_start,n_f_stop, &
       work1,i_costf_init,d_costf_init)

  !-- New map:
  DO n_r=1,n_r_max
     DO n_f=n_f_start,n_f_stop
        df(n_f,n_r)=drx(n_r)*df(n_f,n_r)
     END DO
  END DO


  RETURN
end subroutine get_dr

!-----------------------------------------------------------------------
