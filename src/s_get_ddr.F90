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
!  |  Returns first rarial derivative df and second radial             |
!  |  derivative ddf of the input function f.                          |
!  |  Array f(n_f_max,*) may contain several functions numbered by     |
!  |  the first index. The subroutine calculates the derivaties of     |
!  |  the functions f(n_f_start,*) to f(n_f_stop) by transforming      |
!  |  to a Chebychev representation using n_r_max radial grid points.  |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


    IMPLICIT NONE

!-- input:
    INTEGER :: n_f_max         ! first dim of f
    REAL(kind=8) :: f(n_f_max,*)
    INTEGER :: n_f_start       ! first function to be treated
    INTEGER :: n_f_stop        ! last function to be treated
    INTEGER :: n_r_max         ! number of radial grid points
    INTEGER :: n_cheb_max      ! number of cheb modes
    REAL(kind=8) :: work1(n_f_max,*)  ! work array needed for costf
    REAL(kind=8) :: work2(n_f_max,*)  ! work array for f transfer
    REAL(kind=8) :: drx(*)           ! first derivatives of x(r)
    REAL(kind=8) :: ddrx(*)          ! second derivatives of x(r)

    INTEGER :: i_costf_init(*) ! info for costf
    REAL(kind=8) :: d_costf_init(*)  ! info for costf

!-- output:
    REAL(kind=8) :: df(n_f_max,*)  ! first derivative of f
    REAL(kind=8) :: ddf(n_f_max,*) ! second derivative of f

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
