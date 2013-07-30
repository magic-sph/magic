!$Id$
!********************************************************************
    SUBROUTINE hInt2Tor(Tor,lb,ub,nR,lmStart,lmStop, &
                        Tor2hInt,TorAs2hInt,map)
!********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!
!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE usefull, ONLY: cc2real
    USE LMmapping,only:mappings
    IMPLICIT NONE

!-- Input:
    INTEGER,INTENT(IN) :: lb,ub
    COMPLEX(kind=8),intent(IN) :: Tor(lb:ub)   ! Toroidal field Potential
    INTEGER,intent(IN) :: nR
    INTEGER,intent(IN) :: lmStart,lmStop
    TYPE(mappings),intent(IN) :: map

!-- Output:
    REAL(kind=8),intent(INOUT) :: Tor2hInt,TorAs2hInt

!-- Local:
    REAL(kind=8) :: help,rE4

    INTEGER :: lm,l,m

!-- end of declaration
!---------------------------------------------------------------------

!--- Set zero before calling for first lm-block
!        rPol2hInt  =0.D0
!        Tor2hInt   =0.D0
!        TorAs2hInt =0.D0

    rE4=r(nR)**4
    DO lm=lmStart,lmStop
       l=map%lm2l(lm)
       m=map%lm2m(lm)
       help=rE4/dLh(st_map%lm2(l,m))*cc2real(Tor(lm),m)
        IF ( m == 0 ) TorAs2hInt=TorAs2hInt+help
        Tor2hInt=Tor2hInt+help
    END DO

    RETURN
    end SUBROUTINE hInt2Tor

!-----------------------------------------------------------------------------
