!$Id$
MODULE integration

  IMPLICIT NONE

  CONTAINS
!*********************************************************************
    REAL(kind=8) FUNCTION rInt(f,nRmax,drFac,i_costf_init,d_costf_init)
!*********************************************************************
!  This function performs the radial integral over a
!  function f that is given on the appropriate nRmax
!  radial Chebychev grid points.
!  The arrays i_costf_init,d_costf_init are
!  defined by calling init_costf1.
!  Note: drFac maps radius to cheb space [-1,1]
!        drFac=2.D0/(rMax-rMin)
!---------------------------------------------------------------------

    IMPLICIT NONE

!-- Input:
    REAL(kind=8) :: f(*)
    INTEGER :: nRmax
    REAL(kind=8) :: drFac
    INTEGER :: i_costf_init(*)
    REAL(kind=8) :: d_costf_init(*)

!-- Internal:
    INTEGER :: nCheb,nR
    REAL(kind=8) :: WORK(nRmax)
    REAL(kind=8) :: f2(nRmax)
!---------------------------------------------------------------------

!-- Save input:
    DO nR=1,nRmax
        f2(nR)=f(nR)
    END DO

!-- Transform to cheb space:
    CALL costf1(f2,1,1,1,work,i_costf_init,d_costf_init)
    f2(1)    =0.5D0*f2(1)
    f2(nRmax)=0.5D0*f2(nRmax)

!-- Sum contribution:
    rInt=f2(1)           ! This is zero order contribution
    DO nCheb=3,nRmax,2  ! Only even chebs contribute
        rInt=rInt-1.D0/DBLE(nCheb*(nCheb-2))*f2(nCheb)
    END DO

!-- Remaining renormalisation:
    rInt=2.D0/drFac*DSQRT(2.D0/DBLE(nRmax-1))*rInt

    RETURN
    END FUNCTION rInt
!---------------------------------------------------------------------
    REAL(kind=8) FUNCTION rIntIC(f,nRmax,drFac, &
                                 i_costf_init,d_costf_init)
!  This function performs the radial integral over a
!  function f that is given on the appropriate nRmax
!  radial Chebychev grid points.
!  The arrays i_costf_init,d_costf_init are
!  defined by calling init_costf1.
!---------------------------------------------------------------------

    IMPLICIT NONE

!-- Input:
    REAL(kind=8) :: f(*)
    INTEGER :: nRmax
    REAL(kind=8) :: drFac
    INTEGER :: i_costf_init(*)
    REAL(kind=8) :: d_costf_init(*)

!-- Internal:
    INTEGER :: nCheb,nChebInt
    REAL(kind=8) :: weight
    REAL(kind=8) :: work(nRmax)

!---------------------------------------------------------------------

    CALL costf1(f,1,1,1,work,i_costf_init,d_costf_init)
    f(1)    =0.5D0*f(1)
    f(nRmax)=0.5D0*f(nRmax)

!-- Sum contribution:
    rIntIC=f(1)           ! This is zero order contribution
    DO nCheb=2,nRmax      ! Only even chebs for IC
        nChebInt=2*nCheb-1
        weight  =-1.D0/DBLE(nChebInt*(nChebInt-2))
        rIntIC  =rIntIC+weight*f(nCheb)
    END DO

!-- Remaining renormalisation:
    rIntIC=DSQRT(2.D0/DBLE(nRmax-1))*rIntIC/drFac

    RETURN
    END FUNCTION rIntIC
!---------------------------------------------------------------------
    REAL(kind=8) FUNCTION rInt_R(f,n_r_max,n_cheb_max,dr_fac, &
                                    i_costf_init,d_costf_init)
!   Same as function rInt but for a radial dependent mapping function
!   dr_fac2.
!---------------------------------------------------------------------------
    IMPLICIT NONE
               
    REAL(kind=8),INTENT(IN) :: f(*)
    INTEGER,INTENT(IN) :: n_r_max,n_cheb_max
    REAL(kind=8),INTENT(IN) :: dr_fac(*)
    REAL(kind=8),INTENT(IN) :: d_costf_init(*)
    INTEGER,INTENT(IN) :: i_costf_init(*)
            
    REAL(kind=8) :: f2(n_r_max)
    REAL(kind=8) :: work(n_r_max)
    INTEGER :: nR,nCheb
               
!---------------------------------------------------------------------------
            
!--- Integrals:
    DO nR=1,n_r_max
        f2(nR)=f(nR)/dr_fac(nR)
    END DO

!-- Transform to cheb space:
    CALL costf1(f2,1,1,1,work,i_costf_init,d_costf_init)
    f2(1)      =0.5D0*f2(1)
    f2(n_r_max)=0.5D0*f2(n_r_max)

!-- Sum contribution:
    rInt_R=f2(1)            ! This is zero order contribution
    DO nCheb=3,n_cheb_max,2 ! Only even chebs contribute
        rInt_R=rInt_R-1.D0/DBLE(nCheb*(nCheb-2))*f2(nCheb)
    END DO

!-- Remaining renormalisation:
    rInt_R=2.D0*DSQRT(2.D0/DBLE(n_r_max-1))*rInt_R

    END FUNCTION rInt_R
!---------------------------------------------------------------------------
END MODULE integration
