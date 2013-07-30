!$Id$
!*******************************************************************************
    SUBROUTINE legTFG(nBc,lDeriv,nThetaStart, &
       vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc, &
                 dvrdtc,dvrdpc,dvtdpc,dvpdpc, &
          brc,btc,bpc,cbrc,cbtc,cbpc,sc,drSc, &
              dLhw,dLhdw,dLhz,vhG,vhC,dvhdrG, &
      dvhdrC,dLhb,dLhj,bhG,bhC,cbhG,cbhC,sR,dsR)
!*******************************************************************************

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
    INTEGER,intent(IN) :: nBc
    LOGICAL,intent(IN) :: lDeriv
    INTEGER,intent(IN) :: nThetaStart


!----- Stuff precomputed in legPrep:
    COMPLEX(kind=8),intent(IN) :: dLhw(lm_max),dLhdw(lm_max),dLhz(lm_max)
    COMPLEX(kind=8),intent(IN) :: vhG(lm_max),vhC(lm_max)
    COMPLEX(kind=8),intent(IN) :: dvhdrG(lm_max),dvhdrC(lm_max)
    COMPLEX(kind=8),intent(IN) :: dLhb(lm_max),dLhj(lm_max)
    COMPLEX(kind=8),intent(IN) :: bhG(lm_max),bhC(lm_max)
    COMPLEX(kind=8),intent(IN) :: cbhG(lm_max),cbhC(lm_max)
    COMPLEX(kind=8),intent(IN) :: sR(lm_max),dsR(lm_max)

!------ Legendre Polynomials in c_horizontal.f
    REAL(kind=8) :: PlmG(lm_max)
    REAL(kind=8) :: PlmC(lm_max)

!-- output: field on grid (theta,m) for the radial grid point nR
!           and equatorially symmetric and antisymmetric contribution
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: vrc,vtc,vpc
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dvrdrc,dvtdrc,dvpdrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cvrc
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: brc,btc,bpc
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cbrc,cbtc,cbpc
    COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: sc,drSc
     
!-- local:
    COMPLEX(kind=8) :: vrES,vrEA,dvrdrES,dvrdrEA,dvrdtES,dvrdtEA,cvrES,cvrEA
    complex(kind=8) :: brES,brEA,cbrES,cbrEA,sES,sEA,drsES,drsEA
    INTEGER :: nThetaN,nThetaS,nThetaNHS
    INTEGER :: mc,lm,lmS
    REAL(kind=8) :: dm,dmT

    COMPLEX(kind=8) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
    COMPLEX(kind=8) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
    COMPLEX(kind=8) :: dvhdrN1M(n_m_max),dvhdrN2M(n_m_max),dvhdrN
    COMPLEX(kind=8) :: dvhdrS1M(n_m_max),dvhdrS2M(n_m_max),dvhdrS
    COMPLEX(kind=8) :: dvhdrN1,dvhdrN2,dvhdrS1,dvhdrS2
    COMPLEX(kind=8) :: bhN1M(n_m_max),bhN2M(n_m_max),bhN,bhN1,bhN2
    COMPLEX(kind=8) :: bhS1M(n_m_max),bhS2M(n_m_max),bhS,bhS1,bhS2
    COMPLEX(kind=8) :: cbhN1M(n_m_max),cbhN2M(n_m_max),cbhN,cbhN1,cbhN2
    COMPLEX(kind=8) :: cbhS1M(n_m_max),cbhS2M(n_m_max),cbhS,cbhS1,cbhS2

     
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
            END IF

        !--- Loop over all oders m: (numbered by mc)
            DO mc=1,n_m_max
                lmS=lStop(mc)
                cvrES  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cvrEA  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dvrdrES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dvrdrEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                brES   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                brEA   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            !--- 6 add/mult, 26 DBLE words
                DO lm=lStart(mc),lmS-1,2
                    cvrES  =cvrES   +  dLhz(lm)  *Plm(lm,nThetaNHS)
                    dvrdrES=dvrdrES + dLhdw(lm)  *Plm(lm,nThetaNHS)
                    brES   =brES    +  dLhb(lm)  *Plm(lm,nThetaNHS)
                    cvrEA  =cvrEA   +  dLhz(lm+1)*Plm(lm+1,nThetaNHS)
                    dvrdrEA=dvrdrEA + dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
                    brEA   =brEA    +  dLhb(lm+1)*Plm(lm+1,nThetaNHS)
                END DO
                IF ( lmOdd(mc) ) THEN
                    cvrES  =cvrES   +  dLhz(lmS)*Plm(lmS,nThetaNHS)
                    dvrdrES=dvrdrES + dLhdw(lmS)*Plm(lmS,nThetaNHS)
                    brES   =brES    +  dLhb(lmS)*Plm(lmS,nThetaNHS)
                END IF
                cvrc(mc,nThetaN)  =cvrES  +cvrEA
                cvrc(mc,nThetaS)  =cvrES  -cvrEA
                dvrdrc(mc,nThetaN)=dvrdrES+dvrdrEA
                dvrdrc(mc,nThetaS)=dvrdrES-dvrdrEA
                brc(mc,nThetaN)   =brES   +brEA
                brc(mc,nThetaS)   =brES   -brEA
            END DO

            DO mc=1,n_m_max
                dm =D_mc2m(mc)
                lmS=lStop(mc)
                vrES   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                vrEA   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dvrdtES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                dvrdtEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cbrES  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cbrEA  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            !--- 8 add/mult, 29 DBLE words
                DO lm=lStart(mc),lmS-1,2
                    vrES    =vrES    + dLhw(lm)*   Plm(lm,nThetaNHS)
                    dvrdtEA =dvrdtEA + dLhw(lm)*  dPlm(lm,nThetaNHS)
                    cbrES   =cbrES   + dLhj(lm)*   Plm(lm,nThetaNHS)
                    PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                    PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                    vrEA    =vrEA    + dLhw(lm+1)* Plm(lm+1,nThetaNHS)
                    dvrdtES =dvrdtES + dLhw(lm+1)*dPlm(lm+1,nThetaNHS)
                    cbrEA   =cbrEA   + dLhj(lm+1)* Plm(lm+1,nThetaNHS)
                    PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                    PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
                END DO
                IF ( lmOdd(mc) ) THEN
                    vrES    =vrES    + dLhw(lmS)* Plm(lmS,nThetaNHS)
                    dvrdtEA =dvrdtEA + dLhw(lmS)*dPlm(lmS,nThetaNHS)
                    cbrES   =cbrES   + dLhj(lmS)* Plm(lmS,nThetaNHS)
                    PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                    PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
                END IF
                vrc(mc,nThetaN)   =vrES   +vrEA
                vrc(mc,nThetaS)   =vrES   -vrEA
                dvrdtc(mc,nThetaN)=dvrdtES+dvrdtEA
                dvrdtc(mc,nThetaS)=dvrdtES-dvrdtEA
                cbrc(mc,nThetaN)  =cbrES  +cbrEA
                cbrc(mc,nThetaS)  =cbrES  -cbrEA
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

            DO mc=1,n_m_max
                lmS=lStop(mc)
                bhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            !--- 8 add/mult, 20 DBLE words
                DO lm=lStart(mc),lmS-1,2
                    bhN1=bhN1+bhG(lm)*PlmG(lm)+bhG(lm+1)*PlmG(lm+1)
                    bhS1=bhS1-bhG(lm)*PlmC(lm)+bhG(lm+1)*PlmC(lm+1)
                    bhN2=bhN2+bhC(lm)*PlmC(lm)+bhC(lm+1)*PlmC(lm+1)
                    bhS2=bhS2-bhC(lm)*PlmG(lm)+bhC(lm+1)*PlmG(lm+1)
                END DO
                IF ( lmOdd(mc) ) THEN
                    bhN1=bhN1+bhG(lmS)*PlmG(lmS)
                    bhS1=bhS1-bhG(lmS)*PlmC(lmS)
                    bhN2=bhN2+bhC(lmS)*PlmC(lmS)
                    bhS2=bhS2-bhC(lmS)*PlmG(lmS)
                END IF
                bhN1M(mc)=0.5D0*bhN1
                bhS1M(mc)=0.5D0*bhS1
                bhN2M(mc)=0.5D0*bhN2
                bhS2M(mc)=0.5D0*bhS2
            END DO

            DO mc=1,n_m_max
                lmS=lStop(mc)
                cbhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cbhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cbhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                cbhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
            !--- 8 add/mult, 20 DBLE words
                DO lm=lStart(mc),lmS-1,2
                    cbhN1=cbhN1+cbhG(lm)*PlmG(lm)+cbhG(lm+1)*PlmG(lm+1)
                    cbhS1=cbhS1-cbhG(lm)*PlmC(lm)+cbhG(lm+1)*PlmC(lm+1)
                    cbhN2=cbhN2+cbhC(lm)*PlmC(lm)+cbhC(lm+1)*PlmC(lm+1)
                    cbhS2=cbhS2-cbhC(lm)*PlmG(lm)+cbhC(lm+1)*PlmG(lm+1)
                END DO
                IF ( lmOdd(mc) ) THEN
                    cbhN1=cbhN1+cbhG(lmS)*PlmG(lmS)
                    cbhS1=cbhS1-cbhG(lmS)*PlmC(lmS)
                    cbhN2=cbhN2+cbhC(lmS)*PlmC(lmS)
                    cbhS2=cbhS2-cbhC(lmS)*PlmG(lmS)
                END IF
                cbhN1M(mc)=0.5D0*cbhN1
                cbhS1M(mc)=0.5D0*cbhS1
                cbhN2M(mc)=0.5D0*cbhN2
                cbhS2M(mc)=0.5D0*cbhS2
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
        !--- 6 add/mult, 20 DBLE words
            DO mc=1,n_m_max
                btc(mc,nThetaN)=bhN1M(mc)+bhN2M(mc)
                bhN            =bhN1M(mc)-bhN2M(mc)
                btc(mc,nThetaS)=bhS1M(mc)+bhS2M(mc)
                bhS            =bhS1M(mc)-bhS2M(mc)
                bpc(mc,nThetaN)=CMPLX(AIMAG(bhN),-REAL(bhN),KIND=KIND(0d0))
                bpc(mc,nThetaS)=CMPLX(AIMAG(bhS),-REAL(bhS),KIND=KIND(0d0))
            END DO
        !--- 6 add/mult, 20 DBLE words
            DO mc=1,n_m_max
                cbtc(mc,nThetaN)=cbhN1M(mc)+cbhN2M(mc)
                cbhN            =cbhN1M(mc)-cbhN2M(mc)
                cbtc(mc,nThetaS)=cbhS1M(mc)+cbhS2M(mc)
                cbhS            =cbhS1M(mc)-cbhS2M(mc)
                cbpc(mc,nThetaN)=CMPLX(AIMAG(cbhN),-REAL(cbhN),KIND=KIND(0d0))
                cbpc(mc,nThetaS)=CMPLX(AIMAG(cbhS),-REAL(cbhS),KIND=KIND(0d0))
            END DO ! Loop over order m


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
                    vrc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    vtc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    vpc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cvrc(mc,nThetaN)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvrdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvtdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvpdrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvrdtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvrdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                END DO
                DO mc=n_m_max+1,ncp
                    dvtdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    dvpdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    brc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    btc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    bpc(mc,nThetaN)   =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cbrc(mc,nThetaN)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cbtc(mc,nThetaN)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    cbpc(mc,nThetaN)  =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
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

            DO mc=1,n_m_max
                dm =D_mc2m(mc)
                lmS=lStop(mc)

            !------ br = r^2 B_r , bt = r sin(theta) B_theta , bp= r sin(theta) B_phi
                brES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                brEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                DO lm=lStart(mc),lmS-1,2
                    brES=brES + dLhb(lm)  *Plm(lm,nThetaNHS)
                    brEA=brEA + dLhb(lm+1)*Plm(lm+1,nThetaNHS)
                    PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                    PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                    PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                    PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
                END DO
                IF ( lmOdd(mc) ) THEN
                    brES=brES+dLhb(lm)*Plm(lm,nThetaNHS)
                    PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                    PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                END IF
                brc(mc,nThetaN)=brES+brEA
                brc(mc,nThetaS)=brES-brEA

                bhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                bhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                DO lm=lStart(mc),lmS-1,2
                    bhN1=bhN1+bhG(lm)*PlmG(lm)+bhG(lm+1)*PlmG(lm+1)
                    bhS1=bhS1-bhG(lm)*PlmC(lm)+bhG(lm+1)*PlmC(lm+1)
                    bhN2=bhN2+bhC(lm)*PlmC(lm)+bhC(lm+1)*PlmC(lm+1)
                    bhS2=bhS2-bhC(lm)*PlmG(lm)+bhC(lm+1)*PlmG(lm+1)
                END DO
                IF ( lmOdd(mc) ) THEN
                    bhN1=bhN1+bhG(lmS)*PlmG(lmS)
                    bhS1=bhS1-bhG(lmS)*PlmC(lmS)
                    bhN2=bhN2+bhC(lmS)*PlmC(lmS)
                    bhS2=bhS2-bhC(lmS)*PlmG(lmS)
                END IF
                btc(mc,nThetaN)=0.5D0*bhN1+0.5D0*bhN2
                btc(mc,nThetaS)=0.5D0*bhS1+0.5D0*bhS2
                bhN            =0.5D0*bhN1-0.5D0*bhN2
                bhS            =0.5D0*bhS1-0.5D0*bhS2
                bpc(mc,nThetaN)=CMPLX(AIMAG(bhN),-REAL(bhN),KIND=KIND(0d0))
                bpc(mc,nThetaS)=CMPLX(AIMAG(bhS),-REAL(bhS),KIND=KIND(0d0))

            END DO
             
            IF ( nBc == 1 ) THEN

            !--- Horizontal velocity components for nBc=1
                DO mc=1,n_m_max
                    lmS=lStop(mc)
                    vhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    vhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    vhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    vhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
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
                    brc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    btc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                    bpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
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


    IF ( l_HT ) THEN    ! For movie output !
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
    end SUBROUTINE legTFG

!------------------------------------------------------------------------
