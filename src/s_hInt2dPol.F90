!$Id$
!********************************************************************
SUBROUTINE hInt2dPol(dPol,lmStart,lmStop,Pol2hInt,PolAs2hInt,map)
  !********************************************************************

  !------------ This is release 2 level 1  --------------!
  !------------ Created on 1/17/02  by JW. --------------!
  !--------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE blocking
  USE horizontal_data
  USE usefull, ONLY: cc2real
  USE LMmapping,only: mappings
  IMPLICIT NONE

  !-- Input:
  COMPLEX(kind=8),INTENT(IN) :: dPol(lm_max)   ! Toroidal field Potential
  INTEGER,INTENT(IN) :: lmStart,lmStop
  TYPE(mappings),INTENT(IN) :: map

  !-- Output:
  REAL(kind=8),INTENT(INOUT) :: Pol2hInt,PolAs2hInt

  !-- Local:
  REAL(kind=8) :: help

  INTEGER :: lm,l,m

  !-- end of declaration
  !---------------------------------------------------------------------

  !--- Set zero before calling for first lm-block
  !        Pol2hInt   =0.D0
  !        PolAs2hInt =0.D0

  DO lm=lmStart,lmStop
     l=map%lm2l(lm)
     m=map%lm2m(lm)
     help=dLh(st_map%lm2(l,m))*cc2real(dPol(lm),m)
     IF ( m == 0 ) PolAs2hInt=PolAs2hInt+help
     Pol2hInt=Pol2hInt+help
  END DO

  RETURN
end SUBROUTINE hInt2dPol

!-----------------------------------------------------------------------------
