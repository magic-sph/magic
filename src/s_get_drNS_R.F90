!$Id$
!*************************************************************************
    SUBROUTINE get_drNS_R(f,df,n_fMax,n_fStart,n_fStop, &
                      n_r_max,n_cheb_max,dr_fac2,work1, &
                               i_costf_init,d_costf_init)
!*************************************************************************
!   Same as function s_get_dr but for a radial dependent mapping function
!   dr_fac2.
!---------------------------------------------------------------------------

    IMPLICIT NONE

!--- Input:
    INTEGER :: n_fMax,n_fStart,n_fStop
    INTEGER :: n_r_max
    REAL(kind=8) :: f(n_fMax,*)
    INTEGER :: n_cheb_max
    REAL(kind=8) :: dr_fac2(*)
    REAL(kind=8) :: work1(n_fMax,*)

!--- Output:
    REAL(kind=8) :: df(n_fMax,*)

    INTEGER :: i_costf_init(*)   ! info for transform
    REAL(kind=8) ::  d_costf_init(*)   ! info for tranform

!-- Local:
    REAL(kind=8) :: drx(n_r_max)
    INTEGER :: nR,nF

!--- End of declaration
!-------------------------------------------------------------------------

!--- Get derivative based on the new transform:
    DO nR=1,n_r_max
        drx(nR)=1.D0
    END DO
    CALL get_drNS(f,df,n_fMax,n_fStart,n_fStop, &
                      n_r_max,n_cheb_max,work1, &
                  i_costf_init,d_costf_init,drx)

!--- Apply mapping function:
    DO nR=1,n_r_max
        DO nF=n_fStart,n_fStop
            df(nF,nR)=dr_fac2(nR)*df(nF,nR)
        END DO
    END DO

    RETURN
    end SUBROUTINE get_drNS_R
         
!---------------------------------------------------------------------------

