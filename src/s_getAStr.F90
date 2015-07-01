!$Id$
!****************************************************************************
    SUBROUTINE getAStr(fZ,flmn,nZmax,nZmaxA,lmMax,lMax, &
                       rMin,rMax,nChebMax,rZ,Plm)
!****************************************************************************
!  Calculates function value at radii rZ(nZmax) and
!  colatitudes for which Plm(theta) is given from
!  the spherical hamornic/Chebychev coefficients of an
!  axisymmetric function (order=0).
!----------------------------------------------------------------------------

    IMPLICIT NONE

!--- Output:
    REAL(kind=8) :: fZ(*)

!--- Input:
    INTEGER :: lmMax,lMax
    REAL(kind=8) :: flmn(lmMax,*)
    INTEGER :: nZmax,nZmaxA
    REAL(kind=8) :: rMin,rMax
    INTEGER :: nChebMax
    REAL(kind=8) :: rZ(nZmaxA/2+1)
    REAL(kind=8) :: Plm(lmMax,nZmaxA/2+1)

!--- Local:
    INTEGER :: nCheb
    REAL(kind=8) :: cheb(nChebMax)
    INTEGER ::   l,nZS,nZN
    REAL(kind=8) :: x,chebNorm,flmr

!----------------------------------------------------------------------------

    chebNorm=DSQRT(2.D0/DBLE(nChebMax-1))

    DO nZN=1,nZmax
        fZ(nZN)=0.D0
    END DO

    DO nZN=1,nZmax/2 ! Loop over all z-points in south HS
        nZS=nZmax+1-nZN

    !--- Map r to cheb intervall [-1,1]:
    !    and calculate the cheb polynomia:
    !    Note: the factor chebNorm is needed
    !    for renormalisation. Its not needed if one used
    !    costf1 for the back transform.
        x=2.D0*(rZ(nZN)-0.5D0*(rMin+rMax))/(rMax-rMin)
        cheb(1)=1.D0*chebNorm
        cheb(2)=x*chebNorm
        DO nCheb=3,nChebMax
            cheb(nCheb)=2.D0*x*cheb(nCheb-1)-cheb(nCheb-2)
        END DO
        cheb(1)       =0.5D0*cheb(1)
        cheb(nChebMax)=0.5D0*cheb(nChebMax)

    !--- Loop to add all contribution functions:
        DO l=0,lMax
            flmr=0.D0
            DO nCheb=1,nChebMax
                flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            END DO
            fZ(nZN)=fZ(nZN)+flmr*Plm(l+1,nZN)
            IF ( MOD(l,2) == 0 ) THEN ! Even contribution
                fZ(nZS)=fZ(nZS)+flmr*Plm(l+1,nZN)
            ELSE
                fZ(nZS)=fZ(nZS)-flmr*Plm(l+1,nZN)
            END IF
        END DO

    END DO

    IF ( MOD(nZmax,2) == 1 ) THEN ! Remaining equatorial point
        nZS=(nZmax-1)/2+1

        x=2.D0*(rZ(nZS)-0.5D0*(rMin+rMax))/(rMax-rMin)
        cheb(1)=1.D0*chebNorm
        cheb(2)=x*chebNorm
        DO nCheb=3,nChebMax
            cheb(nCheb)=2.D0*x*cheb(nCheb-1)-cheb(nCheb-2)
        END DO
        cheb(1)       =0.5D0*cheb(1)
        cheb(nChebMax)=0.5D0*cheb(nChebMax)

        DO l=0,lMax
            flmr=0.D0
            DO nCheb=1,nChebMax
                flmr=flmr+flmn(l+1,nCheb)*cheb(nCheb)
            END DO
            fZ(nZS)=fZ(nZS)+flmr*Plm(l+1,nZS)
        END DO

    END IF


    end SUBROUTINE getAStr

!----------------------------------------------------------------------
