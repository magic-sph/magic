!$Id$
!***********************************************************************
    SUBROUTINE get_td(nR,nBc,lRmsCalc,dVSrLM,dVxBhLM, &
                       dwdt,dzdt,dpdt,dsdt,dbdt,djdt, &
              AdvrLM,AdvtLM,AdvpLM,LFrLM,LFtLM,LFpLM, &
              VSrLM,VStLM,VSpLM,VxBrLM,VxBtLM,VxBpLM, &
                                ViscHeatLM,OhmLossLM, &
                     dLhw,dLhdw,dLhz,sR,preR,dpR)
!***********************************************************************

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!  +-------------------------------------------------------------------+
!  |                                                                   |
!  |  Purpose of this to calculate time derivatives                    |
!  |  dwdt,dzdt,dpdt,dsdt,dbdt,djdt                                    |
!  |  and auxiliary arrays dVSrLM and dVxBhLM                          |
!  |  from non-linear terms in spectral form,                          |
!  |  contained in flmw1-3,flms1-3, flmb1-3 (input)                    |
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
    USE RMS

    IMPLICIT NONE

!-- Input of variables:
    INTEGER,intent(IN) :: nR
    INTEGER,intent(IN) :: nBc ! signifies baundary conditions
    LOGICAL,intent(IN) :: lRmsCalc

!----- Nonlinear terms:
    COMPLEX(kind=8),intent(IN) :: AdvrLM(lmP_max),AdvtLM(lmP_max),AdvpLM(lmP_max)
    COMPLEX(kind=8),intent(IN) :: LFrLM(lmP_max),LFtLM(lmP_max),LFpLM(lmP_max)
    COMPLEX(kind=8),intent(IN) :: VSrLM(lmP_max),VStLM(lmP_max),VSpLM(lmP_max)
    COMPLEX(kind=8),intent(IN) :: VxBrLM(lmP_max),VxBtLM(lmP_max),VxBpLM(lmP_max)
    COMPLEX(kind=8),intent(IN) :: ViscHeatLM(lmP_max)
    COMPLEX(kind=8),intent(IN) :: OhmLossLM(lmP_max)

!----- Scalar fields in radial distribute space:
    COMPLEX(kind=8),intent(IN) :: dLhw(lm_max),dLhdw(lm_max),dLhz(lm_max)
    COMPLEX(kind=8),intent(IN) :: sR(lm_max),preR(lm_max),dpR(lm_max)
    !COMPLEX(kind=8) :: dsR(lm_max)

!-- Output:
    COMPLEX(kind=8),intent(OUT) :: dwdt(lm_max),dzdt(lm_max)
    COMPLEX(kind=8),intent(OUT) :: dpdt(lm_max),dsdt(lm_max)
    COMPLEX(kind=8),intent(OUT) :: dbdt(lm_maxMag),djdt(lm_maxMag)
    COMPLEX(kind=8),intent(OUT) :: dVxBhLM(lm_maxMag)
    COMPLEX(kind=8),intent(OUT) :: dVSrLM(lm_max)

!-- local:
    INTEGER :: l,m,lm,lmS,lmA,lmP,lmPS,lmPA
    COMPLEX(kind=8) :: CorPol(lm_max),CorTor(lm_max)
    COMPLEX(kind=8) :: AdvPol(lm_max),AdvTor(lm_max)
    COMPLEX(kind=8) :: LFPol(lm_max),LFTor(lm_max)
    COMPLEX(kind=8) :: zR(lm_max),wR(lm_max),dwR(lm_max)

!-- end of declaration
!-------------------------------------------------------------------------


    DO lm=1,lm_max
        IF ( dLh(lm) > 0 ) THEN
            wR(lm) =dLhw(lm)/dLh(lm)
            dwR(lm)=dLhdw(lm)/dLh(lm)
            zR(lm) =dLhz(lm)/dLh(lm)
        ELSE
            wR(lm) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dwR(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            zR(lm) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        ENDIF
    END DO

    IF ( nBc > 0 ) THEN   ! boundary !

        IF ( l_mag_nl .OR. l_mag_kin ) THEN

        !----- Stress free boundary, only nl mag. term for poloidal field needed.
        !      Because the radial derivative will be taken, this will contribute to
        !      the other radial grid points.
            dVxBhLM(1)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            dVSrLM(1) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            DO lm=2,lm_max
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)   ! l-1
                lmPA=lmP2lmPA(lmP)   ! l+1
                IF ( l > m ) THEN
                    dVxBhLM(lm)=r(nR)*r(nR)* (                &
                                  dTheta1S(lm)*VxBtLM(lmPS) - &
                                  dTheta1A(lm)*VxBtLM(lmPA) + &
                                      dPhi(lm)*VxBpLM(lmP)  )
                ELSE IF ( l == m ) THEN ! (l-1) not allowed !
                    dVxBhLM(lm)=r(nR)*r(nR)* (                &
                                - dTheta1A(lm)*VxBtLM(lmPA) + &
                                      dPhi(lm)*VxBpLM(lmP)  )
                END IF
                dVSrLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END DO
                       
        ELSE
            DO lm=1,lm_max
                IF ( l_mag ) dVxBhLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dVSrLM(lm) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END DO
        END IF
        IF ( l_heat ) THEN
            DO lm=1,lm_max
                dVSrLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END DO
        END IF

    ELSE
         
        IF ( l_conv ) THEN  ! Convection

            lm =1   ! This is l=0,m=0
            lmA=lm2lmA(lm)
            lmP=1
            lmPA=lmP2lmPA(lmP)
            IF ( l_conv_nl ) THEN
                AdvPol(lm)=      or2(nR)*AdvrLM(lm)
                AdvTor(lm)=-dTheta1A(lm)*AdvpLM(lmPA)
            ELSE
                AdvPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                AdvTor(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
            IF ( l_corr ) THEN
                CorPol(lm)=2.D0*CorFac*or1(nR) * dTheta2A(lm)* zR(lmA)
                CorTor(lm)= 2.d0*CorFac*or2(nR) * ( &
                            dTheta3A(lm)*dwR(lmA) + &
                    or1(nR)*dTheta4A(lm)* wR(lmA) )
            ELSE
                CorPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                CorTor(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            END IF
            dwdt(lm)=AdvPol(lm)+CorPol(lm)
            dzdt(lm)=AdvTor(lm)+CorTor(lm)
            IF ( lRmsCalc .AND. l_mag_LF ) THEN
                LFPol(lm) =      or2(nR)*LFrLM(lm)
                LFTor(lm) =-dTheta1A(lm)*LFpLM(lmPA)
                AdvPol(lm)=AdvPol(lm)-LFPol(lm)
                AdvTor(lm)=AdvTor(lm)-LFTor(lm)
            END IF

            DO lm=2,lm_max
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmS =lm2lmS(lm)
                lmA =lm2lmA(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)
                lmPA=lmP2lmPA(lmP)

                IF ( l_corr ) THEN
                    IF ( l < l_max .AND. l > m ) THEN
                        CorPol(lm) =2.D0*CorFac*or1(nR) * ( &
                                       dPhi0(lm)*dwR(lm) +  & ! phi-deriv of dw/dr
                                    dTheta2A(lm)*zR(lmA) -  & ! sin(theta) dtheta z
                                    dTheta2S(lm)*zR(lmS) )
                    ELSE IF ( l == l_max ) THEN
                        CorPol(lm)= 2.D0*CorFac*or1(nR) * ( &
                                    dPhi0(lm)*dwR(lm)  )
                    ELSE IF ( l == m ) THEN
                        CorPol(lm) =2.D0*CorFac*or1(nR) * ( &
                                       dPhi0(lm)*dwR(lm)  + &
                                     dTheta2A(lm)*zR(lmA) )
                    END IF
                ELSE
                    CorPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                IF ( l_conv_nl ) THEN
                    AdvPol(lm)=or2(nR)*AdvrLM(lmP)
                ELSE
                    AdvPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                dwdt(lm)=AdvPol(lm)+CorPol(lm)
                IF ( lRmsCalc .AND. l_mag_LF ) THEN
                !------ When RMS values are required, the Lorentz force is treated
                !       seperately:
                !--------- FLPol = ( curl(B)xB )_r
                    LFPol(lm) =or2(nR)*LFrLM(lmP)
                    AdvPol(lm)=AdvPol(lm)-LFPol(lm)
                END IF
                   
                IF ( l_corr ) THEN
                    IF ( l < l_max .AND. l > m ) THEN
                        CorTor(lm)=2.d0*CorFac*or2(nR) * ( &
                                      dPhi0(lm)*zR(lm)   + &
                                  dTheta3A(lm)*dwR(lmA)  + &
                          or1(nR)*dTheta4A(lm)* wR(lmA)  + &
                                  dTheta3S(lm)*dwR(lmS)  - &
                          or1(nR)*dTheta4S(lm)* wR(lmS)  )
                    ELSE IF ( l == l_max ) THEN
                        CorTor(lm)=2.d0*CorFac*or2(nR) * ( &
                                   dPhi0(lm)*zR(lm)   )
                    ELSE IF ( l == m ) THEN
                        CorTor(lm)=2.d0*CorFac*or2(nR) * ( &
                                      dPhi0(lm)*zR(lm)   + &
                                  dTheta3A(lm)*dwR(lmA)  + &
                          or1(nR)*dTheta4A(lm)* wR(lmA)  )
                    END IF
                ELSE
                    CorTor(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                IF ( l_conv_nl ) THEN
                    IF ( l > m ) THEN
                        AdvTor(lm)=   -dPhi(lm)*AdvtLM(lmP)  + &
                                   dTheta1S(lm)*AdvpLM(lmPS) - &
                                   dTheta1A(lm)*AdvpLM(lmPA)
                    ELSE IF ( l == m ) THEN
                        AdvTor(lm)=   -dPhi(lm)*AdvtLM(lmP)  - &
                                   dTheta1A(lm)*AdvpLM(lmPA)
                    END IF
                ELSE
                    AdvTor(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                dzdt(lm)=AdvTor(lm)+CorTor(lm)
                IF ( lRmsCalc .AND. l_mag_LF ) THEN
                !------ When RMS values are required, the Lorentz force is treated
                !       seperately:
                    IF ( l > m ) THEN
                    !------- LFTor= 1/(E*Pm) * curl( curl(B) x B )_r
                        LFTor(lm) =   -dPhi(lm)*LFtLM(lmP)  + &
                                   dTheta1S(lm)*LFpLM(lmPS) - &
                                   dTheta1A(lm)*LFpLM(lmPA)
                    ELSE IF ( l == m ) THEN
                        LFTor(lm) =   -dPhi(lm)*LFtLM(lmP)  - &
                                   dTheta1A(lm)*LFpLM(lmPA)
                    END IF
                    AdvTor(lm)=AdvTor(lm)-LFTor(lm)
                END IF
                    
            END DO

            IF ( lRmsCalc ) THEN

                IF ( l_conv_nl ) THEN
                    CALL hInt2Pol(AdvPol,nR,2,lm_max,AdvPolLMr, &
                                  AdvPol2hInt(nR),AdvPolAs2hInt(nR))
                    CALL hInt2Tor(AdvTor,nR,2,lm_max, &
                                  AdvTor2hInt(nR),AdvTorAs2hInt(nR))
                END IF
            ! rho* grad(p/rho) = grad(p) - beta*p
                CALL hInt2Pol(dpR-beta(nR)*preR,nR,2,lm_max,PreLMr, &
                              Pre2hInt(nR),PreAs2hInt(nR))
                IF ( ra /= 0.D0 ) &
                CALL hInt2Pol(sR,nR,2,lm_max,BuoLMr, &
                              Buo2hInt(nR),BuoAs2hInt(nR))
                IF ( l_corr ) THEN
                    CALL hInt2Pol(CorPol,nR,2,lm_max,CorPolLMr, &
                                  CorPol2hInt(nR),CorPolAs2hInt(nR))
                    CALL hInt2Tor(CorTor,nR,2,lm_max, &
                                  CorTor2hInt(nR),CorTorAs2hInt(nR))
                END IF
                IF ( l_mag_LF ) THEN
                    CALL hInt2Pol(LFPol,nR,2,lm_max,LFPolLMr, &
                                  LFPol2hInt(nR),LFPolAs2hInt(nR))
                    CALL hInt2Tor(LFTor,nR,2,lm_max, &
                                  LFTor2hInt(nR),LFTorAs2hInt(nR))
                END IF

            !----- Calculate balances: recycle AdvPol, AdvTor for this:
                DO lm=2,lm_max
                    IF ( l_RMStest ) THEN
                        AdvPol(lm)=AdvPol(lm)+CorPol(lm)+LFPol(lm)- &
                                         dpR(lm)+beta(nR)*preR(lm)+ &
                                         rho0(nR)*rgrav(nR)*sR(lm)
                        AdvTor(lm)=AdvTor(lm)+CorTor(lm)+LFTor(lm)
                        CorPol(lm)=dLh(lm)*or2(nR)*wR(lm)
                    ELSE
                        AdvPol(lm)=CorPol(lm)-dpR(lm)+beta(nR)*preR(lm)
                        AdvTor(lm)=AdvPol(lm)+LFPol(lm)
                        CorPol(lm)=AdvTor(lm)+rho0(nR)*rgrav(nR)*sR(lm)
                    END IF
                END DO
                CALL hInt2Pol(AdvPol,nR,2,lm_max,GeoLMr, &
                              Geo2hInt(nR),GeoAs2hInt(nR))
                IF ( l_RMStest ) THEN
                    CALL hInt2Tor(AdvTor,nR,2,lm_max, &
                                  Mag2hInt(nR),MagAs2hInt(nR))
                ELSE
                    CALL hInt2Pol(AdvTor,nR,2,lm_max,MagLMr, &
                                  Mag2hInt(nR),MagAs2hInt(nR))
                END IF
                CALL hInt2Pol(CorPol,nR,2,lm_max,ArcLMr, &
                              Arc2hInt(nR),ArcAs2hInt(nR))

            END IF

            DO lm=2,lm_max
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmS =lm2lmS(lm)
                lmA =lm2lmA(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)
                lmPA=lmP2lmPA(lmP)

            !------ Recycle CorPol and AdvPol:
                IF ( l_corr ) THEN
                    IF ( l < l_max .AND. l > m ) THEN
                        CorPol(lm)=  2.d0*CorFac*or2(nR) * ( &
                                  -dPhi0(lm) * ( dwR(lm)   + &
                                  or1(nR)*dLh(lm)*wR(lm) ) + &
                                     dTheta3A(lm)*zR(lmA)  + &
                                     dTheta3S(lm)*zR(lmS)  )
                    ELSE IF ( l == l_max ) THEN
                        CorPol(lm)=  2.d0*CorFac*or2(nR) * ( &
                                  -dPhi0(lm) * ( dwR(lm)   + &
                                  or1(nR)*dLh(lm)*wR(lm) ) )
                    ELSE IF ( l == m ) THEN
                        CorPol(lm)=  2.d0*CorFac*or2(nR) * ( &
                                  -dPhi0(lm) * ( dwR(lm)   + &
                                  or1(nR)*dLh(lm)*wR(lm) ) + &
                                     dTheta3A(lm)*zR(lmA)  )
                    END IF
                ELSE
                    CorPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                IF ( l_conv_nl ) THEN
                    IF ( l > m ) THEN
                        AdvPol(lm)= dTheta1S(lm)*AdvtLM(lmPS) - &
                                    dTheta1A(lm)*AdvtLM(lmPA) + &
                                        dPhi(lm)*AdvpLM(lmP)
                    ELSE IF ( l == m ) THEN
                        AdvPol(lm)=-dTheta1A(lm)*AdvtLM(lmPA) + &
                                        dPhi(lm)*AdvpLM(lmP)
                    END IF
                ELSE
                    AdvPol(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END IF
                dpdt(lm)=AdvPol(lm)+CorPol(lm)
                 
            END DO ! lm loop

        ELSE
            DO lm=2,lm_max
                dwdt(lm) =0.D0
                dzdt(lm) =0.D0
                dpdt(lm) =0.D0
            END DO
        END IF ! l_conv ?

        IF ( l_heat ) THEN
               
            dsdt(1)  =epsc*epscProf(nR)
            dVSrLM(1)=VSrLM(1)
            IF ( l_anel ) THEN
                IF ( l_mag_nl ) THEN
                    dsdt(1)=dsdt(1)+ViscHeatFac*hdif_V(1)*ViscHeatLM(1)+ &
                            OhmLossFac*hdif_B(1)*OhmLossLM(1)
                ELSE
                    dsdt(1)=dsdt(1)+ViscHeatFac*hdif_V(1)*ViscHeatLM(1)
                END IF
            END IF
            DO lm=2,lm_max
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)
                lmPA=lmP2lmPA(lmP)
            !------ This is horizontal heat advection:
                IF ( l > m ) THEN
                    dsdt(lm)= -dTheta1S(lm)*VStLM(lmPS) + &
                               dTheta1A(lm)*VStLM(lmPA) - &
                                   dPhi(lm)*VSpLM(lmP)
                ELSE IF ( l == m ) THEN
                    dsdt(lm)=  dTheta1A(lm)*VStLM(lmPA) - &
                                   dPhi(lm)*VSpLM(lmP)
                END IF

                IF ( l_anel ) THEN
                    IF ( l_mag_nl ) THEN
                        dsdt(lm)=dsdt(lm)+ &
                                 ViscHeatFac*hdif_V(lm)*ViscHeatLM(lmP)+ &
                                 OhmLossFac*hdif_B(lm)*OhmLossLM(lmP)
                    ELSE
                        dsdt(lm)=dsdt(lm)+ &
                                 ViscHeatFac*hdif_V(lm)*ViscHeatLM(lmP)
                    END IF
                END IF
            !-----   simplified form for linear onset !
            !        not ds not saved in the current program form!
            !                 dsdt(lm)=
            !                    -dLh(lm)*w(lm,nR)*or2(nR)*dsR(1)
                dVSrLM(lm)=VSrLM(lmP)
            END DO

        ELSE
            DO lm=2,lm_max
                dsdt(lm)  =0.D0
                dVSrLM(lm)=0.D0
            END DO
        END IF

        IF ( l_mag_nl .OR. l_mag_kin  ) THEN

            lm=1
            lmP=1
            lmPa=lmP2lmPA(lmP)
            dVxBhLM(lm)= -r(nR)*r(nR)* &
                          dTheta1A(lm)*VxBtLM(lmPA)
            dbdt(lm)   = -dTheta1A(lm)*VxBpLM(lmPA)
            djdt(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))

            DO lm=2,lm_max
                l   =lm2l(lm)
                m   =lm2m(lm)
                lmP =lm2lmP(lm)
                lmPS=lmP2lmPS(lmP)
                lmPA=lmP2lmPA(lmP)

            !------- This is the radial part of the dynamo terms \curl(VxB)
                IF ( l > m ) THEN
                    dbdt(lm)=  dTheta1S(lm)*VxBpLM(lmPS) - &
                               dTheta1A(lm)*VxBpLM(lmPA) - &
                                   dPhi(lm)*VxBtLM(lmP)
                ELSE IF ( l == m ) THEN
                    dbdt(lm)= -dTheta1A(lm)*VxBpLM(lmPA) - &
                                   dPhi(lm)*VxBtLM(lmP)
                END IF

            !------- Radial component of
            !           \curl\curl(UxB) = \grad\div(UxB) - \laplace(VxB)
                 
            !------- This is the radial part of \laplace (UxB)
                djdt(lm)=dLh(lm)*or4(nR)*VxBrLM(lmP)

            !------- This is r^2 * horizontal divergence of (UxB)
            !        Radial derivative performed in get_dr_td
                IF ( l > m ) THEN
                    dVxBhLM(lm)=          r(nR)*r(nR)* ( &
                             dTheta1S(lm)*VxBtLM(lmPS) - &
                             dTheta1A(lm)*VxBtLM(lmPA) + &
                                 dPhi(lm)*VxBpLM(lmP)  )
                ELSE IF ( l == m ) THEN
                    dVxBhLM(lm)=          r(nR)*r(nR)* ( &
                           - dTheta1A(lm)*VxBtLM(lmPA) + &
                                 dPhi(lm)*VxBpLM(lmP)  )
                END IF
            END DO

        ELSE
            IF ( l_mag ) THEN
                DO lm=1,lm_max
                   dbdt(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                   djdt(lm)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                   dVxBhLM(lm)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
               END DO
            END IF
        END IF

    END IF  ! boundary ? lvelo ?


    RETURN
    end SUBROUTINE get_td
!-----------------------------------------------------------------------------
