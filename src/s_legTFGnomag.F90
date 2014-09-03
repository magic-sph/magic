!$Id$
!********************************************************************************
SUBROUTINE legTFGnomag(nBc,lDeriv,nThetaStart, &
     &                 vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc, &
     &                 dvrdtc,dvrdpc,dvtdpc,dvpdpc,sc,drSc, &
     &                                      dsdtc,dsdpc,pc, &
     &                 dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG,dvhdrC,sR,dsR,pR)
  !********************************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !-----------------------------------------------------------------------------------

  !    Legendre transform from (nR,l,m) to (nR,nTheta,m) [spectral to grid]
  !    where nTheta numbers the colatitudes and l is the degree of
  !    the spherical harmonic representation.

  !    Transforms entropy, velocity and magnetic field components
  !    and terms involving spatial derivatives.
  !    The symmetry properties of the P_lm with respect to the equator
  !    are used. The equatorially anti-symmetric (EA) contribution is
  !    added to (subracted from ) the equatorially symmetric (ES) contribution
  !    in northern (southern) hemisphere.

  !     nBc            : (input) accounts for special conditions on radial boundaries
  !        nBc=2       : we are dealing with a no slip boundary, v_r and v_theta are
  !                      CMPLX(0.D0,0.D0,KIND=KIND(0d0)) and v_phi=r sin(theta) omega, where
  !                      omega is the rotation rate of the boundary (mantle of IC),
  !                      only magn. field terms are calculated, v is set later.
  !        nBc=1       : a free slip bounday: v_r is zero, derivatives of v and B are
  !                      not needed, only components of v,B and entropy are calculated
  !        nBc=0       : normal case, interior grid point
  !     lDeric=.TRUE.  : (input) calculate derivatives
  !     nThetaStart    : (input) transformation is done for the range of
  !                      points nThetaStart <= nTheta <= nThetaStart-1+sizeThetaB
  !     Plm            : associated Legendre polynomials
  !     dPlm           : sin(theta) d Plm / d theta
  !     osn2           : 1/sin(theta)^2
  !     vrc, ...., drSc: (output) components in (nTheta,m)-space
  !     dLhw,....,cbhC : (input) help arrays calculated in s_legPrep.f

  !-----------------------------------------------------------------------------------

  USE truncation
  USE blocking
  USE horizontal_data
  USE logic
  USE const
  IMPLICIT NONE

  !-- input:
  INTEGER,INTENT(IN) :: nBc
  LOGICAL,INTENT(IN) :: lDeriv
  INTEGER,INTENT(IN) :: nThetaStart

  !----- Stuff precomputed in legPrep:
  COMPLEX(kind=8),INTENT(IN) :: dLhw(lm_max),dLhdw(lm_max),dLhz(lm_max)
  COMPLEX(kind=8),INTENT(IN) :: vhG(lm_max),vhC(lm_max)
  COMPLEX(kind=8),INTENT(IN) :: dvhdrG(lm_max),dvhdrC(lm_max)
  COMPLEX(kind=8),INTENT(IN) :: sR(lm_max),dsR(lm_max)
  COMPLEX(kind=8),INTENT(IN) :: pR(lm_max)

  !------ Legendre Polynomials in c_horizontal.f
  !       REAL(kind=8) Plm(lm_max,n_theta_max/2)
  !       REAL(kind=8) dPlm(lm_max,n_theta_max/2)
  !       REAL(kind=8) osn2(*)
  REAL(kind=8) :: PlmG(lm_max)
  REAL(kind=8) :: PlmC(lm_max)

  !-- output: field on grid (theta,m) for the radial grid point nR
  !           and equatorially symmetric and antisymmetric contribution
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: vrc,vtc,vpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dvrdrc,dvtdrc,dvpdrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cvrc,sc,drSc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dsdtc,dsdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: pc

  !-- local:
  COMPLEX(kind=8) :: vrES,vrEA,dvrdrES,dvrdrEA,dvrdtES,dvrdtEA,cvrES,cvrEA
  complex(kind=8) :: sES,sEA,drsES,drsEA,pES,pEA
  COMPLEX(kind=8) :: dsdtES,dsdtEA

  INTEGER :: nThetaN,nThetaS,nThetaNHS
  INTEGER :: mc,lm,lmS
  REAL(kind=8) :: dm,dmT

  COMPLEX(kind=8) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
  COMPLEX(kind=8) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
  COMPLEX(kind=8) :: dvhdrN1M(n_m_max),dvhdrN2M(n_m_max),dvhdrN
  COMPLEX(kind=8) :: dvhdrS1M(n_m_max),dvhdrS2M(n_m_max),dvhdrS
  COMPLEX(kind=8) :: dvhdrN1,dvhdrN2,dvhdrS1,dvhdrS2

  !-- end of declaration
  !----------------------------------------------------------------------


  nThetaNHS=(nThetaStart-1)/2

  IF ( nBc == 0 .OR. lDeriv ) THEN ! not a boundary or derivs required

     DO nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
        nThetaS  =nThetaN+1      ! same theta but for southern HS
        nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point

        IF ( l_heat ) THEN
           DO mc=1,n_m_max
              lmS=lStop(mc)
              sES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
              sEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
              DO lm=lStart(mc),lmS-1,2
                 sES=sES+sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+sR(lmS)*Plm(lmS,nThetaNHS)
              sc(mc,nThetaN)=sES+sEA
              sc(mc,nThetaS)=sES-sEA
           END DO


           IF ( l_fluxProfs ) THEN
              DO mc=1,n_m_max
                 lmS=lStop(mc)
                 pES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
                 pEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
                 DO lm=lStart(mc),lmS-1,2
                    pES=pES+pR(lm)  *Plm(lm,nThetaNHS)
                    pEA=pEA+pR(lm+1)*Plm(lm+1,nThetaNHS)
                 END DO
                 IF ( lmOdd(mc) ) pES=pES+pR(lmS)*Plm(lmS,nThetaNHS)
                 pc(mc,nThetaN)=pES+pEA
                 pc(mc,nThetaS)=pES-pEA
              END DO
           END IF

           IF ( l_viscBcCalc ) THEN
              DO mc=1,n_m_max
                 dm =D_mc2m(mc)
                 lmS=lStop(mc)
                 dsdtES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 dsdtEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 DO lm=lStart(mc),lmS-1,2
                    dsdtEA =dsdtEA + sR(lm)*  dPlm(lm,nThetaNHS)
                    dsdtES =dsdtES + sR(lm+1)*dPlm(lm+1,nThetaNHS)
                 END DO
                 IF ( lmOdd(mc) ) THEN
                    dsdtEA =dsdtEA + sR(lmS)*dPlm(lmS,nThetaNHS)
                 END IF
                 dsdtc(mc,nThetaN)=dsdtES+dsdtEA
                 dsdtc(mc,nThetaS)=dsdtES-dsdtEA
              END DO

              DO mc=1,n_m_max
                 dm=D_mc2m(mc)
                 dsdpc(mc,nThetaN)= &
                  CMPLX(-dm*AIMAG(sc(mc,nThetaN)), &
                          dm*REAL(sc(mc,nThetaN)),KIND=KIND(0d0))
                 dsdpc(mc,nThetaS)= &
                  CMPLX(-dm*AIMAG(sc(mc,nThetaS)), &
                          dm*REAL(sc(mc,nThetaS)),KIND=KIND(0d0))
              END DO

           END IF ! thermal dissipation layer

        END IF

        !--- Loop over all oders m: (numbered by mc)
        DO mc=1,n_m_max
           lmS=lStop(mc)
           cvrES  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           cvrEA  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvrdrES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvrdrEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 4 add/mult
           DO lm=lStart(mc),lmS-1,2
              cvrES  =cvrES   +  dLhz(lm)  *Plm(lm,nThetaNHS)
              dvrdrES=dvrdrES + dLhdw(lm)  *Plm(lm,nThetaNHS)
              cvrEA  =cvrEA   +  dLhz(lm+1)*Plm(lm+1,nThetaNHS)
              dvrdrEA=dvrdrEA + dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              cvrES  =cvrES   +  dLhz(lmS)*Plm(lmS,nThetaNHS)
              dvrdrES=dvrdrES + dLhdw(lmS)*Plm(lmS,nThetaNHS)
           END IF
           cvrc(mc,nThetaN)  =cvrES  +cvrEA
           cvrc(mc,nThetaS)  =cvrES  -cvrEA
           dvrdrc(mc,nThetaN)=dvrdrES+dvrdrEA
           dvrdrc(mc,nThetaS)=dvrdrES-dvrdrEA
        END DO

        DO mc=1,n_m_max
           dm =D_mc2m(mc)
           lmS=lStop(mc)
           vrES   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vrEA   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvrdtES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvrdtEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 29 DBLE words
           DO lm=lStart(mc),lmS-1,2
              vrES    =vrES    + dLhw(lm)*   Plm(lm,nThetaNHS)
              dvrdtEA =dvrdtEA + dLhw(lm)*  dPlm(lm,nThetaNHS)
              PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
              vrEA    =vrEA    + dLhw(lm+1)* Plm(lm+1,nThetaNHS)
              dvrdtES =dvrdtES + dLhw(lm+1)*dPlm(lm+1,nThetaNHS)
              PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
              PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              vrES    =vrES    + dLhw(lmS)* Plm(lmS,nThetaNHS)
              dvrdtEA =dvrdtEA + dLhw(lmS)*dPlm(lmS,nThetaNHS)
              PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
              PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
           END IF
           vrc(mc,nThetaN)   =vrES   +vrEA
           vrc(mc,nThetaS)   =vrES   -vrEA
           dvrdtc(mc,nThetaN)=dvrdtES+dvrdtEA
           dvrdtc(mc,nThetaS)=dvrdtES-dvrdtEA
        END DO

        !--- Now the stuff using generalized harmonics:
        DO mc=1,n_m_max
           lmS=lStop(mc)
           vhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              vhN1=vhN1+vhG(lm)*PlmG(lm)+vhG(lm+1)*PlmG(lm+1)
              vhS1=vhS1-vhG(lm)*PlmC(lm)+vhG(lm+1)*PlmC(lm+1)
              vhN2=vhN2+vhC(lm)*PlmC(lm)+vhC(lm+1)*PlmC(lm+1)
              vhS2=vhS2-vhC(lm)*PlmG(lm)+vhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              vhN1=vhN1+vhG(lmS)*PlmG(lmS)
              vhS1=vhS1-vhG(lmS)*PlmC(lmS)
              vhN2=vhN2+vhC(lmS)*PlmC(lmS)
              vhS2=vhS2-vhC(lmS)*PlmG(lmS)
           END IF
           vhN1M(mc)=0.5D0*vhN1
           vhS1M(mc)=0.5D0*vhS1
           vhN2M(mc)=0.5D0*vhN2
           vhS2M(mc)=0.5D0*vhS2
        END DO

        DO mc=1,n_m_max
           lmS=lStop(mc)
           dvhdrN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              dvhdrN1=dvhdrN1+dvhdrG(lm)  *PlmG(lm) + &
                   dvhdrG(lm+1)*PlmG(lm+1)
              dvhdrS1=dvhdrS1-dvhdrG(lm)  *PlmC(lm) + &
                   dvhdrG(lm+1)*PlmC(lm+1)
              dvhdrN2=dvhdrN2+dvhdrC(lm)  *PlmC(lm) + &
                   dvhdrC(lm+1)*PlmC(lm+1)
              dvhdrS2=dvhdrS2-dvhdrC(lm)  *PlmG(lm) + &
                   dvhdrC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              dvhdrN1=dvhdrN1+dvhdrG(lmS)*PlmG(lmS)
              dvhdrS1=dvhdrS1-dvhdrG(lmS)*PlmC(lmS)
              dvhdrN2=dvhdrN2+dvhdrC(lmS)*PlmC(lmS)
              dvhdrS2=dvhdrS2-dvhdrC(lmS)*PlmG(lmS)
           END IF
           dvhdrN1M(mc)=0.5D0*dvhdrN1
           dvhdrS1M(mc)=0.5D0*dvhdrS1
           dvhdrN2M(mc)=0.5D0*dvhdrN2
           dvhdrS2M(mc)=0.5D0*dvhdrS2
        END DO

        !--- Unscramble:
        !--- 6 add/mult, 20 DBLE words
        DO mc=1,n_m_max
           vtc(mc,nThetaN)=vhN1M(mc)+vhN2M(mc)
           vhN            =vhN1M(mc)-vhN2M(mc)
           vtc(mc,nThetaS)=vhS1M(mc)+vhS2M(mc)
           vhS            =vhS1M(mc)-vhS2M(mc)
           vpc(mc,nThetaN)=CMPLX(AIMAG(vhN),-REAL(vhN),KIND=KIND(0d0))
           vpc(mc,nThetaS)=CMPLX(AIMAG(vhS),-REAL(vhS),KIND=KIND(0d0))
        END DO
        !--- 6 add/mult, 20 DBLE words
        DO mc=1,n_m_max
           dvtdrc(mc,nThetaN)=dvhdrN1M(mc)+dvhdrN2M(mc)
           dvhdrN            =dvhdrN1M(mc)-dvhdrN2M(mc)
           dvtdrc(mc,nThetaS)=dvhdrS1M(mc)+dvhdrS2M(mc)
           dvhdrS            =dvhdrS1M(mc)-dvhdrS2M(mc)
           dvpdrc(mc,nThetaN)= &
                CMPLX(AIMAG(dvhdrN),-REAL(dvhdrN),KIND=KIND(0d0))
           dvpdrc(mc,nThetaS)= &
                CMPLX(AIMAG(dvhdrS),-REAL(dvhdrS),KIND=KIND(0d0))
        END DO

        !--- Calculate phi derivatives:
        DO mc=1,n_m_max
           dm=D_mc2m(mc)
           dvrdpc(mc,nThetaN)= &
                CMPLX(-dm*AIMAG(vrc(mc,nThetaN)), &
                dm*REAL(vrc(mc,nThetaN)),KIND=KIND(0d0))
           dvrdpc(mc,nThetaS)= &
                CMPLX(-dm*AIMAG(vrc(mc,nThetaS)), &
                dm*REAL(vrc(mc,nThetaS)),KIND=KIND(0d0))
        END DO
        DO mc=1,n_m_max
           dmT=D_mc2m(mc)*osn2(nThetaNHS)
           dvtdpc(mc,nThetaN)= &
                CMPLX(-dmT*AIMAG(vtc(mc,nThetaN)), &
                dmT*REAL(vtc(mc,nThetaN)),KIND=KIND(0d0))
           dvtdpc(mc,nThetaS)= &
                CMPLX(-dmT*AIMAG(vtc(mc,nThetaS)), &
                dmT*REAL(vtc(mc,nThetaS)),KIND=KIND(0d0))
           dvpdpc(mc,nThetaN)= &
                CMPLX(-dmT*AIMAG(vpc(mc,nThetaN)), &
                dmT*REAL(vpc(mc,nThetaN)),KIND=KIND(0d0))
           dvpdpc(mc,nThetaS)= &
                CMPLX(-dmT*AIMAG(vpc(mc,nThetaS)), &
                dmT*REAL(vpc(mc,nThetaS)),KIND=KIND(0d0))
        END DO   ! End of loop over oder m numbered by mc

     END DO      ! End global loop over nTheta


     !-- Zero out terms with index mc > n_m_max:
     IF ( n_m_max < ncp ) THEN
        DO nThetaN=1,sizeThetaB
           DO mc=n_m_max+1,ncp
              sc(mc,nThetaN)    =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              IF ( l_fluxProfs) THEN
                 pc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              END IF
              IF ( l_viscBcCalc) THEN
                 dsdtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 dsdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              END IF
              vrc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vtc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vpc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              cvrc(mc,nThetaN)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvrdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvtdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvpdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvrdtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvrdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvtdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              dvpdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END DO
        END DO  ! loop over nThetaN (theta)
     END IF


  ELSE   ! boundary ?


     !-- Calculation for boundary r_cmb or r_icb:

     DO nThetaN=1,sizeThetaB,2
        nThetaS=nThetaN+1
        nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point

        IF ( l_heat ) THEN
           DO mc=1,n_m_max
              lmS=lStop(mc)
              sES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
              sEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
              DO lm=lStart(mc),lmS-1,2
                 sES=sES+sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+sR(lmS)*Plm(lmS,nThetaNHS)
              sc(mc,nThetaN)=sES+sEA
              sc(mc,nThetaS)=sES-sEA
           END DO
        END IF

        IF ( nBc == 1 ) THEN

           !--- Horizontal velocity components for nBc=1
           DO mc=1,n_m_max
              dm =D_mc2m(mc)
              lmS=lStop(mc)
              vhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              DO lm=lStart(mc),lmS-1,2
                 PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                 PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                 PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                 PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
                 vhN1=vhN1+vhG(lm)*PlmG(lm)+vhG(lm+1)*PlmG(lm+1)
                 vhS1=vhS1-vhG(lm)*PlmC(lm)+vhG(lm+1)*PlmC(lm+1)
                 vhN2=vhN2+vhC(lm)*PlmC(lm)+vhC(lm+1)*PlmC(lm+1)
                 vhS2=vhS2-vhC(lm)*PlmG(lm)+vhC(lm+1)*PlmG(lm+1)
              END DO
              IF ( lmOdd(mc) ) THEN
                 PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                 PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
                 vhN1=vhN1+vhG(lmS)*PlmG(lmS)
                 vhS1=vhS1-vhG(lmS)*PlmC(lmS)
                 vhN2=vhN2+vhC(lmS)*PlmC(lmS)
                 vhS2=vhS2-vhC(lmS)*PlmG(lmS)
              END IF
              vtc(mc,nThetaN)=0.5D0*vhN1+0.5D0*vhN2
              vtc(mc,nThetaS)=0.5D0*vhS1+0.5D0*vhS2
              vhN            =0.5D0*vhN1-0.5D0*vhN2
              vhS            =0.5D0*vhS1-0.5D0*vhS2
              vpc(mc,nThetaN)=CMPLX(AIMAG(vhN),-REAL(vhN),KIND=KIND(0d0))
              vpc(mc,nThetaS)=CMPLX(AIMAG(vhS),-REAL(vhS),KIND=KIND(0d0))
           END DO ! Loop over m

        END IF   ! nBc.eq.1 ? vrc and nBc=2 cared for later !

     END DO    ! End loop over nThetaN

     !-- Zero out terms with index mc > n_m_max :
     IF ( n_m_max < ncp ) THEN
        DO nThetaN=1,sizeThetaB
           DO mc=n_m_max+1,ncp
              sc(mc,nThetaN) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END DO
        END DO
        IF ( nBc == 1 ) THEN
           DO nThetaN=1,sizeThetaB
              DO mc=n_m_max+1,ncp
                 vtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 vpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              END DO
           END DO
        END IF
     END IF

  END IF  ! boundary ? nBc?


  IF ( l_HT .OR. l_viscBcCalc ) THEN    ! For movie output !
     nThetaNHS=(nThetaStart-1)/2

     !-- Caculate radial derivate of S for heatflux:
     DO nThetaN=1,sizeThetaB,2   ! Loop over thetas for one HS
        nThetaS  =nThetaN+1  ! same theta but at other HS
        nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
        DO mc=1,n_m_max
           lmS=lStop(mc)
           drsES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           drsEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           DO lm=lStart(mc),lmS-1,2
              drsES=drsES+dsR(lm)*Plm(lm,nThetaNHS)
              drsEA=drsEA+dsR(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) drsES=drsES+dsR(lmS)*Plm(lmS,nThetaNHS)
           drSc(mc,nThetaN)=drsES+drsEA
           drSc(mc,nThetaS)=drsES-drsEA
        END DO
     END DO
     !-- Zero out terms with index mc > n_m_max:
     IF ( n_m_max < ncp ) THEN
        DO nThetaN=1,sizeThetaB
           DO mc=n_m_max+1,ncp
              drSc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END DO
        END DO  ! loop over nThetaN (theta)
     END IF

  END IF


  RETURN
end SUBROUTINE legTFGnomag

!------------------------------------------------------------------------
