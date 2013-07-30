!$Id$
!*************************************************************************
    SUBROUTINE init_rNB(r,n_r_max,n_cheb_max,rCut,rDea, &
                               r2,n_r_max2,n_cheb_max2, &
                   nS,dr_fac2,i_costf_init2,nDi_costf1, &
                                d_costf_init2,nDd_costf1)
!*************************************************************************
! Prepares the usage of a cut back radial grid where nS points
! on both boundaries are discarded.
! The aim actually is to discard boundary effects, but just
! not considering the boundary grid points does not work when
! you also want radial derivatives and integrals. For these
! we use the Chebychev transform which needs are particular number
! of grid points so that the fast cosine transform can be
! applied. Therefor more than just 2 points have to be
! thrown away, which may make sense anyway.
!-------------------------------------------------------------------------

    IMPLICIT NONE

!--- Input:
    REAL(kind=8) :: r(*),rCut,rDea
    INTEGER :: n_r_max,n_cheb_max
    INTEGER :: nDi_costf1,nDd_costf1
    REAL(kind=8) :: drx(n_r_max)      ! first derivatives of x(r)

!--- Output:
    INTEGER :: nS,n_r_max2,n_cheb_max2
    REAL(kind=8) :: r2(*),dr_fac2(*)
    INTEGER :: i_costf_init2(nDi_costf1)   ! info for transform
    REAL(kind=8) ::  d_costf_init2(nDd_costf1)   ! info for tranfor

! Local stuff
    INTEGER, PARAMETER :: n_r_maxL=192
    REAL(kind=8) :: r2C(n_r_maxL)
    REAL(kind=8) :: r_cheb2(n_r_maxL)
    REAL(kind=8) :: dr2(n_r_maxL)
    REAL(kind=8) :: w1(n_r_maxL),w2(n_r_maxL)
    REAL(kind=8) :: r_icb2,r_cmb2,dr_fac
    INTEGER :: nRs(14)

    INTEGER :: nR,n

!-------------------------------------------------------------------------

!--- New radial grid:

!--- Find number of points to be cut away at either side:
    DO nS=1,(n_r_max-1)/2
        IF ( r(1)-r(nS) > rCut ) GOTO 10
    END DO
    WRITE(*,*) 'No nS found in init_rNB!'
    STOP
    10 CONTINUE
    n_r_max2=n_r_max-2*nS

! Allowed number of radial grid points:
    nRs(1) =25
    nRs(2) =33
    nRs(3) =37
    nRs(4) =41
    nRs(5) =49
    nRs(6) =61
    nRs(7) =65
    nRs(8) =73
    nRs(9) =81
    nRs(10)=97
    nRs(11)=101
    nRs(12)=109
    nRs(13)=121
    nRs(14)=161
    DO n=14,1,-1
        IF ( nRs(n) <= n_r_max2 ) GOTO 20
    END DO
    WRITE(*,*) 'No n_r_max2 found in init_rNB!'
    STOP
    20 CONTINUE
    n_r_max2=nRs(n)
    nS=(n_r_max-n_r_max2)/2
    IF ( n_r_max2 > n_r_maxL ) THEN
        WRITE(*,*) 'Increase n_r_maxL in subroutine dRMS.f!'
        WRITE(*,*) 'Should be at least:',n_r_max2
        STOP
    END IF
    n_cheb_max2=MIN0(INT((1.D0-rDea)*n_r_max2),n_cheb_max)

    DO nR=1,n_r_max2
        r2(nR)=r(nR+nS)
    END DO
    r_icb2=r2(n_r_max2)
    r_cmb2=r2(1)
    CALL cheb_x_map_e(r_icb2,r_cmb2,n_r_max2-1,r2C,r_cheb2)
    CALL init_costf1(n_r_max2,i_costf_init2,nDi_costf1, &
                     d_costf_init2,nDd_costf1)
    dr_fac=1.D0
    DO nR=1,n_r_max
        drx(nR)=1.D0
    END DO
    CALL get_dr(r2,dr2,1,1,1,n_r_max2,n_cheb_max2, &
                w1,w2,i_costf_init2,d_costf_init2,drx)
    DO nR=1,n_r_max2
        dr_fac2(nR)=1.D0/dr2(nR)
    END DO


    RETURN
    end SUBROUTINE init_rNB

!-------------------------------------------------------------------------

