!$Id$
!********************************************************************
    SUBROUTINE hInt2Pol(Pol,lb,ub,nR,lmStart,lmStop,PolLMr, &
                        Pol2hInt,PolAs2hInt,map)
!********************************************************************

!--------------------------------------------------------------------

    USE truncation
    USE radial_functions
    USE blocking
    USE horizontal_data
    USE const
    USE usefull, ONLY: cc2real
    USE LMmapping,only: mappings
    IMPLICIT NONE

!-- Input:
    INTEGER,INTENT(IN) :: lb,ub
    COMPLEX(kind=8),intent(IN) :: Pol(lb:ub)   ! Poloidal field Potential
    INTEGER,intent(IN) :: nR,lmStart,lmStop
    TYPE(mappings),intent(IN) :: map

!-- Output:
    COMPLEX(kind=8),intent(OUT) :: PolLMr(lm_max,n_r_max)
    REAL(kind=8),intent(INOUT) :: Pol2hInt,PolAs2hInt

!-- Local:
    REAL(kind=8) :: help,rE2

    INTEGER :: lm,l,m

!-- end of declaration
!---------------------------------------------------------------------

    rE2=r(nR)*r(nR)
    DO lm=lmStart,lmStop
       l=map%lm2l(lm)
       m=map%lm2m(lm)
       help=rE2*cc2real(Pol(lm),m)
       IF ( m == 0 ) PolAs2hInt=PolAs2hInt+help
       Pol2hInt=Pol2hInt+help
       PolLMr(lm,nR)=rE2/dLh(st_map%lm2(l,m))*Pol(lm)
    END DO

    RETURN
    end SUBROUTINE hInt2Pol

!---------------------------------------------------------------------
