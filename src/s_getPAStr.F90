!$Id$
!****************************************************************************
    SUBROUTINE getPAStr(fZ,flmn,nZmax,nZmaxA,lmMax,lMax, &
                           rMin,rMax,nChebMax,rZ,dPlm,OsinTS)
!****************************************************************************
!  Calculates axisymmetric phi component for a given toroidal
!  potential flmn in spherical harmonic/Cheb space.
!  This is calculated at radii rZ(nZmax) and matching
!  colatitutes theta(nZmax) for which dPlm(theta) and OsinTS(theta)
!  are provided for the south hemisphere.
!  Points in the northern hemipshere use Plm symmetries.
!----------------------------------------------------------------------------

    IMPLICIT NONE

!--- Output:
    REAL(kind=8) ::    fZ(*)

!--- Input:
    INTEGER ::   lmMax,lMax
    REAL(kind=8) ::    flmn(lmMax,*)
    INTEGER ::   nZmax,nZmaxA
    REAL(kind=8) ::    rMin,rMax
    INTEGER ::   nChebMax
    REAL(kind=8) ::    rZ(nZmaxA/2+1)
    REAL(kind=8) ::    dPlm(lmMax,nZmaxA/2+1)
    REAL(kind=8) ::    OsinTS(nZmaxA/2+1)

!--- Local:
    INTEGER ::   nCheb
    REAL(kind=8) ::    cheb(nChebMax)
    INTEGER ::   l,nZS,nZN!,nZ
    REAL(kind=8) ::    x,chebNorm,fac,flmr

!----------------------------------------------------------------------------

    chebNorm=DSQRT(2.D0/DBLE(nChebMax-1))
            
    DO nZN=1,nZmax
        fZ(nZN)=0.D0
    END DO

    DO nZN=1,nZmax/2 ! Loop over all z-points in northern HS
        nZS=nZmax+1-nZN  ! point in southern HS
        fac=-OsinTS(nZN)/rZ(nZN)

    !--- Map r to cheb intervall [-1,1]:
    !    and calculate the cheb polynomia:
    !    Note: the factor chebNorm is needed
    !    for renormalisation. Its not needed if one used
    !    costf1 for the back transform.
        x=2.D0*(rZ(nZN)-0.5D0*(rMin+rMax))/(rMax-rMin)
        cheb(1)=1.D0*chebNorm*fac
        cheb(2)=x*chebNorm*fac
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
            fZ(nZN)=fZ(nZN)+flmr*dPlm(l+1,nZN)
            IF ( MOD(l,2) == 0 ) THEN ! Odd contribution
                fZ(nZS)=fZ(nZS)-flmr*dPlm(l+1,nZN)
            ELSE
                fZ(nZS)=fZ(nZS)+flmr*dPlm(l+1,nZN)
            END IF
        END DO

    END DO

    IF ( MOD(nZmax,2) == 1 ) THEN ! Remaining equatorial point
        nZS=(nZmax-1)/2+1
        fac=-OsinTS(nZS)/rZ(nZS)

        x=2.D0*(rZ(nZS)-0.5D0*(rMin+rMax))/(rMax-rMin)
        cheb(1)=1.D0*chebNorm*fac
        cheb(2)=x*chebNorm*fac
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
            fZ(nZS)=fZ(nZS)+flmr*dPlm(l+1,nZS)
        END DO

    END IF


    end SUBROUTINE getPAStr

!----------------------------------------------------------------------
