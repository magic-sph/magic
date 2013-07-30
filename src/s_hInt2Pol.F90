!$Id$
!********************************************************************
    SUBROUTINE hInt2Pol(Pol,nR,lmStart,lmStop,PolLMr, &
                        Pol2hInt,PolAs2hInt)
!********************************************************************

!--------------------------------------------------------------------

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE const
    USE usefull, ONLY: cc2real

    IMPLICIT NONE

!-- Input:
    COMPLEX(kind=8) :: Pol(lm_max)   ! Poloidal field Potential
    INTEGER :: nR,lmStart,lmStop

!-- Output:
    COMPLEX(kind=8) :: PolLMr(lm_max,n_r_max)
    REAL(kind=8) :: Pol2hInt,PolAs2hInt

!-- Local:
    REAL(kind=8) :: help,rE2

    INTEGER :: lm,m

!-- end of declaration
!---------------------------------------------------------------------

!--- Set zero before calling for first lm-block
!        Pol2hInt   =0.D0
!        PolAs2hInt =0.D0

    rE2=r(nR)*r(nR)
    DO lm=lmStart,lmStop
        m=lm2m(lm)
        help=rE2*cc2real(Pol(lm),m)
        IF ( m == 0 ) PolAs2hInt=PolAs2hInt+help
        Pol2hInt=Pol2hInt+help
        PolLMr(lm,nR)=rE2/dLh(lm)*Pol(lm)
    END DO

    RETURN
    end SUBROUTINE hInt2Pol

!---------------------------------------------------------------------
