!$Id$
!***********************************************************************
    subroutine get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
                               BtVpLM,BpVtLM,BrVZLM,BtVZLM, &
                             BtVpCotLM,BpVtCotLM,BtVZcotLM, &
                             BtVpSn2LM,BpVtSn2LM,BtVZsn2LM, &
                                   TstrRLM,TadvRLM,TomeRLM)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this routine is to calculated theta and phi           |
!  |  derivative related terms of the magnetic production and          |
!  |  andvection terms and store them.                                 |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE blocking
    USE horizontal_data
    USE logic
    USE dtB_mod

    IMPLICIT NONE

!-- input:
! include 'truncation.f'
! include 'c_blocking.f'
! include 'c_horizontal.f'
! include 'c_phys_param.f'
! include 'c_num_param.f'
! include 'c_radial.f'
! include 'c_logic.f'

    INTEGER :: nR

!-- Nonlinear B/v combinations:
    COMPLEX(kind=8) :: BtVrLM(*),BpVrLM(*)
    COMPLEX(kind=8) :: BrVtLM(*),BrVpLM(*)
    COMPLEX(kind=8) :: BtVpLM(*),BpVtLM(*)
    COMPLEX(kind=8) :: BtVpCotLM(*),BpVtCotLM(*)
    COMPLEX(kind=8) :: BtVpSn2LM(*),BpVtSn2LM(*)
    COMPLEX(kind=8) :: BtVZcotLM(*),BtVZsn2LM(*)
    COMPLEX(kind=8) :: BrVZLM(*),BtVZLM(*)

!-- output:
    COMPLEX(kind=8) :: TstrRLM(lm_max_dtB),TadvRLM(lm_max_dtB)
    COMPLEX(kind=8) :: TomeRLM(lm_max_dtB)
!----- Output of PstrLM, PadvLM, TstrLM, TadvLM via
! include 'c_dtB.f'


!-- local:
    INTEGER :: l,m,lm,lmP,lmPS,lmPA
    REAL(kind=8) :: fac


!-- end of declaration
!-------------------------------------------------------------------------


    PstrLM(1,nR)=0.d0
    PadvLM(1,nR)=0.D0
    DO lm=2,lm_max
        l   =lm2l(lm)
        m   =lm2m(lm)
        lmP =lm2lmP(lm)
        lmPS=lmP2lmPS(lmP)
        lmPA=lmP2lmPA(lmP)
        IF ( l > m ) THEN
            PstrLM(lm,nR)=or2(nR)/dLh(lm) *   ( &
                    dTheta1S(lm)*BtVrLM(lmPS) - &
                    dTheta1A(lm)*BtVrLM(lmPA) + &
                         dPhi(lm)*BpVrLM(lmP)  )
            PadvLM(lm,nR)=or2(nR)/dLh(lm) *   ( &
                    dTheta1S(lm)*BrVtLM(lmPS) - &
                    dTheta1A(lm)*BrVtLM(lmPA) + &
                         dPhi(lm)*BrVpLM(lmP)  )
        ELSE IF ( l == m ) THEN
            PstrLM(lm,nR)=or2(nR)/dLh(lm) *   ( &
                  - dTheta1A(lm)*BtVrLM(lmPA) + &
                         dPhi(lm)*BpVrLM(lmP)  )
            PadvLM(lm,nR)=or2(nR)/dLh(lm) *   ( &
                  - dTheta1A(lm)*BrVtLM(lmPA) + &
                         dPhi(lm)*BrVpLM(lmP) )
        END IF
    END DO

!--- Poloidal advection and stretching term finished for radial level nR !

    TstrLM(1,nR)=0.D0
    TstrRLM(1)  =0.D0
    DO lm=2,lm_max
        l   =lm2l(lm)
        m   =lm2m(lm)
        lmP =lm2lmP(lm)
        lmPS=lmP2lmPS(lmP)
        lmPA=lmP2lmPA(lmP)
        fac=or2(nR)/dLh(lm)
        IF ( l > m ) THEN
            TstrLM(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
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
            TstrRLM(lm)=             or1(nR)/dLh(lm) * ( &
                             dTheta1S(lm)*BrVpLM(lmPS) - &
                             dTheta1A(lm)*BrVpLM(lmPA) - &
                                  dPhi(lm)*BrVtLM(lmP)  )
        ELSE IF ( l == m ) THEN
            TstrLM(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
                    fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
                                        BpVtSn2LM(lmP) )   + &
                                                     fac * ( &
                   - dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
                                           BpVtCotLM(lmPA) + &
                                       BtVpCotLM(lmPA) ) ) - &
                          fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
            TstrRLM(lm)=             or1(nR)/dLh(lm) * ( &
                           - dTheta1A(lm)*BrVpLM(lmPA) - &
                                   dPhi(lm)*BrVtLM(lmP)  )
        END IF
    END DO

    TadvLM(1,nR)=0.D0
    TadvRLM(1)  =0.D0
    DO lm=2,lm_max
        l   =lm2l(lm)
        m   =lm2m(lm)
        lmP =lm2lmP(lm)
        lmPS=lmP2lmPS(lmP)
        lmPA=lmP2lmPA(lmP)
        fac=or2(nR)/dLh(lm)
        IF ( l > m ) THEN
            TadvLM(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
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
            TadvRLM(lm)=or2(nR)/dLh(lm) * ( &
                dTheta1S(lm)*BpVrLM(lmPS) - &
                dTheta1A(lm)*BpVrLM(lmPA) - &
                     dPhi(lm)*BtVrLM(lmP)   )
        ELSE IF ( l == m ) THEN
            TadvLM(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
                   fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
                                       BtVpSn2LM(lmP) )   + &
                                                    fac * ( &
                  - dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
                                          BtVpCotLM(lmPA) + &
                                      BpVtCotLM(lmPA) ) ) - &
                         fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
            TadvRLM(lm)=or2(nR)/dLh(lm) * ( &
              - dTheta1A(lm)*BpVrLM(lmPA) - &
                     dPhi(lm)*BtVrLM(lmP)   )
        END IF
    END DO

!--- TomeLM same as TstrLM but where ever Vp appeared
!    it is replaced by its axisymmetric contribution VZ:
    TomeLM(1,nR)=0.D0
    TomeRLM(1)  =0.D0
    DO lm=2,lm_max
        l  =lm2l(lm)
        m  =lm2m(lm)
        lmP=lm2lmP(lm)
        lmPS=lmP2lmPS(lmP)
        lmPA=lmP2lmPA(lmP)
        fac=or2(nR)/dLh(lm)
        IF ( l > m ) THEN
            TomeLM(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
                 fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    + &
                      fac*( dTheta1S(lm)*BtVZCotLM(lmPS) - &
                            dTheta1A(lm)*BtVZCotLM(lmPA) )
            TomeRLM(lm)=          or2(nR)/dLh(lm) * ( &
                          dTheta1S(lm)*BrVZLM(lmPS) - &
                          dTheta1A(lm)*BrVZLM(lmPA) )
        ELSE IF ( l == m ) THEN
            TomeLM(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
                 fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    - &
                      fac*dTheta1A(lm)*BtVZCotLM(lmPA)
            TomeRLM(lm)=        - or2(nR)/dLh(lm) * &
                          dTheta1A(lm)*BrVZLM(lmPA)
        END IF
    END DO


    RETURN
    end subroutine get_dH_dtBLM

!------------------------------------------------------------------------
