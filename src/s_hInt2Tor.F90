!$Id$
!********************************************************************
    SUBROUTINE hInt2Tor(Tor,nR,lmStart,lmStop, &
                        Tor2hInt,TorAs2hInt)
!********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!--------------------------------------------------------------------

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE usefull, ONLY: cc2real

    IMPLICIT NONE

!-- Input:
    COMPLEX(kind=8) :: Tor(lm_max)   ! Toroidal field Potential
    INTEGER :: nR
    INTEGER :: lmStart,lmStop

!-- Output:
    REAL(kind=8) :: Tor2hInt,TorAs2hInt

!-- Local:
    REAL(kind=8) :: help,rE4

    INTEGER :: lm,m

!-- end of declaration
!---------------------------------------------------------------------

!--- Set zero before calling for first lm-block
!        rPol2hInt  =0.D0
!        Tor2hInt   =0.D0
!        TorAs2hInt =0.D0

    rE4=r(nR)**4
    DO lm=lmStart,lmStop
        m=lm2m(lm)
        help=rE4/dLh(lm)*cc2real(Tor(lm),m)
        IF ( m == 0 ) TorAs2hInt=TorAs2hInt+help
        Tor2hInt=Tor2hInt+help
    END DO

    RETURN
    end SUBROUTINE hInt2Tor

!-----------------------------------------------------------------------------
