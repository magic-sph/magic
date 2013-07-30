!$Id$
!***********************************************************************
SUBROUTINE get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
     BtVpLM,BpVtLM,BrVZLM,BtVZLM, &
     BtVpCotLM,BpVtCotLM,BtVZcotLM, &
     BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |  Purpose of this routine is to calculate theta and phi            |
  !  |  derivative related terms of the magnetic production and          |
  !  |  advection terms and store them.                                  |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data, ONLY: dPhi,dLh,dTheta1S,dTheta1A
  USE logic
  USE dtB_mod,ONLY: PstrLM_Rloc, PadvLM_Rloc,TstrLM_Rloc,TadvLM_Rloc,TomeLM_Rloc,&
       & TstrRLM_Rloc,TadvRLM_Rloc,TomeRLM_Rloc

  IMPLICIT NONE

  !-- input:
  INTEGER,intent(IN) :: nR

  !-- Nonlinear B/v combinations:
  COMPLEX(kind=8),intent(IN) :: BtVrLM(*),BpVrLM(*)
  COMPLEX(kind=8),intent(IN) :: BrVtLM(*),BrVpLM(*)
  COMPLEX(kind=8),intent(IN) :: BtVpLM(*),BpVtLM(*)
  COMPLEX(kind=8),intent(IN) :: BtVpCotLM(*),BpVtCotLM(*)
  COMPLEX(kind=8),intent(IN) :: BtVpSn2LM(*),BpVtSn2LM(*)
  COMPLEX(kind=8),intent(IN) :: BtVZcotLM(*),BtVZsn2LM(*)
  COMPLEX(kind=8),intent(IN) :: BrVZLM(*),BtVZLM(*)

  !-- output:
  !COMPLEX(kind=8),INTENT(OUT) :: TstrRLM(lm_max_dtB),TadvRLM(lm_max_dtB)
  !COMPLEX(kind=8),intent(OUT) :: TomeRLM(lm_max_dtB)
  !----- Output of PstrLM, PadvLM, TstrLM, TadvLM via module dtB_mod

  !-- local:
  INTEGER :: l,m,lm,lmP,lmPS,lmPA
  REAL(kind=8) :: fac


  !-- end of declaration
  !-------------------------------------------------------------------------

  PstrLM_Rloc(1,nR)=0.d0
  PadvLM_Rloc(1,nR)=0.D0
  DO lm=2,lm_max
     l   =lm2l(lm)
     m   =lm2m(lm)
     lmP =lm2lmP(lm)
     lmPS=lmP2lmPS(lmP)
     lmPA=lmP2lmPA(lmP)
     IF ( l > m ) THEN
        PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
             dTheta1S(lm)*BtVrLM(lmPS) - &
             dTheta1A(lm)*BtVrLM(lmPA) + &
             dPhi(lm)*BpVrLM(lmP)  )
        PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
             dTheta1S(lm)*BrVtLM(lmPS) - &
             dTheta1A(lm)*BrVtLM(lmPA) + &
             dPhi(lm)*BrVpLM(lmP)  )
     ELSE IF ( l == m ) THEN
        PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
             - dTheta1A(lm)*BtVrLM(lmPA) + &
             dPhi(lm)*BpVrLM(lmP)  )
        PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
             - dTheta1A(lm)*BrVtLM(lmPA) + &
             dPhi(lm)*BrVpLM(lmP) )
     END IF
  END DO

  !--- Poloidal advection and stretching term finished for radial level nR !

  TstrLM_Rloc(1,nR)=0.D0
  TstrRLM_Rloc(1,nR)  =0.D0
  DO lm=2,lm_max
     l   =lm2l(lm)
     m   =lm2m(lm)
     lmP =lm2lmP(lm)
     lmPS=lmP2lmPS(lmP)
     lmPA=lmP2lmPA(lmP)
     fac=or2(nR)/dLh(lm)
     IF ( l > m ) THEN
        TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
             fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
             BpVtSn2LM(lmP) )   + &
             fac * ( &
             dTheta1S(lm) * ( or1(nR)*BpVrLM(lmPS) + &
             BpVtCotLM(lmPS) + &
             BtVpCotLM(lmPS) ) - &
             dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
             BpVtCotLM(lmPA) + &
             BtVpCotLM(lmPA) ) ) - &
             fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
        TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
             dTheta1S(lm)*BrVpLM(lmPS) - &
             dTheta1A(lm)*BrVpLM(lmPA) - &
             dPhi(lm)*BrVtLM(lmP)  )
     ELSE IF ( l == m ) THEN
        TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
             fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
             BpVtSn2LM(lmP) )   + &
             fac * ( &
             - dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
             BpVtCotLM(lmPA) + &
             BtVpCotLM(lmPA) ) ) - &
             fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
        TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
             - dTheta1A(lm)*BrVpLM(lmPA) - &
             dPhi(lm)*BrVtLM(lmP)  )
     END IF
  END DO

  TadvLM_Rloc(1,nR)=0.D0
  TadvRLM_Rloc(1,nR)  =0.D0
  DO lm=2,lm_max
     l   =lm2l(lm)
     m   =lm2m(lm)
     lmP =lm2lmP(lm)
     lmPS=lmP2lmPS(lmP)
     lmPA=lmP2lmPA(lmP)
     fac=or2(nR)/dLh(lm)
     IF ( l > m ) THEN
        TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
             fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
             BtVpSn2LM(lmP) )   + &
             fac * ( &
             dTheta1S(lm) * ( or1(nR)*BrVpLM(lmPS) + &
             BtVpCotLM(lmPS) + &
             BpVtCotLM(lmPS) ) - &
             dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
             BtVpCotLM(lmPA) + &
             BpVtCotLM(lmPA) ) ) - &
             fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
        TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
             dTheta1S(lm)*BpVrLM(lmPS) - &
             dTheta1A(lm)*BpVrLM(lmPA) - &
             dPhi(lm)*BtVrLM(lmP)   )
     ELSE IF ( l == m ) THEN
        TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
             fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
             BtVpSn2LM(lmP) )   + &
             fac * ( &
             - dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
             BtVpCotLM(lmPA) + &
             BpVtCotLM(lmPA) ) ) - &
             fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
        TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
             - dTheta1A(lm)*BpVrLM(lmPA) - &
             dPhi(lm)*BtVrLM(lmP)   )
     END IF
  END DO

  !--- TomeLM same as TstrLM but where ever Vp appeared
  !    it is replaced by its axisymmetric contribution VZ:
  TomeLM_Rloc(1,nR)=0.D0
  TomeRLM_Rloc(1,nR)  =0.D0
  DO lm=2,lm_max
     l  =lm2l(lm)
     m  =lm2m(lm)
     lmP=lm2lmP(lm)
     lmPS=lmP2lmPS(lmP)
     lmPA=lmP2lmPA(lmP)
     fac=or2(nR)/dLh(lm)
     IF ( l > m ) THEN
        TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
             fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    + &
             fac*( dTheta1S(lm)*BtVZCotLM(lmPS) - &
             dTheta1A(lm)*BtVZCotLM(lmPA) )
        TomeRLM_Rloc(lm,nR)=          or2(nR)/dLh(lm) * ( &
             dTheta1S(lm)*BrVZLM(lmPS) - &
             dTheta1A(lm)*BrVZLM(lmPA) )
     ELSE IF ( l == m ) THEN
        TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
             fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    - &
             fac*dTheta1A(lm)*BtVZCotLM(lmPA)
        TomeRLM_Rloc(lm,nR)=        - or2(nR)/dLh(lm) * &
             dTheta1A(lm)*BrVZLM(lmPA)
     END IF
  END DO


  RETURN
end subroutine get_dH_dtBLM

!------------------------------------------------------------------------
