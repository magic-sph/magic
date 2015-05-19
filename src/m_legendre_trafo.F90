!$Id$
!*******************************************************************************
#include "perflib_preproc.cpp"
MODULE legendre_trafo
  USE truncation, ONLY: ncp,lm_max,n_m_max
  USE blocking, ONLY: nfs,sizeThetaB,lm2mc
  USE horizontal_data, ONLY: Plm,dPlm,lStart,lStop,lmOdd,D_mc2m,osn2
  USE logic, ONLY: l_heat,l_ht
  USE parallel_mod,only: rank
#ifdef WITH_LIKWID
# include "likwid_f90.h"
#endif
  USE leg_helper_mod,only: leg_helper_t
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: legTFG,legTFGnomag
 CONTAINS
!*******************************************************************************
 SUBROUTINE legTFG(nBc,lDeriv,lViscBcCalc,lFluxProfCalc,nThetaStart, &
     &            vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc, &
     &            dvrdtc,dvrdpc,dvtdpc,dvpdpc, &
     &            brc,btc,bpc,cbrc,cbtc,cbpc,sc,&
     &            drSc,dsdtc,dsdpc,pc,    &
     &            leg_helper)
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
  

  !-- input:
  INTEGER,intent(IN) :: nBc
  LOGICAL,intent(IN) :: lDeriv,lFluxProfCalc,lViscBcCalc
  INTEGER,intent(IN) :: nThetaStart


  !----- Stuff precomputed in legPrep:
  type(leg_helper_t) :: leg_helper

  !-- output: field on grid (theta,m) for the radial grid point nR
  !           and equatorially symmetric and antisymmetric contribution
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: vrc,vtc,vpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dvrdrc,dvtdrc,dvpdrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cvrc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: brc,btc,bpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cbrc,cbtc,cbpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: sc,drSc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dsdtc,dsdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: pc

  !------ Legendre Polynomials in m_horizontal.F90
  REAL(kind=8) :: PlmG(lm_max)
  REAL(kind=8) :: PlmC(lm_max)

  !-- local:
  COMPLEX(kind=8) :: vrES,vrEA,dvrdrES,dvrdrEA,dvrdtES,dvrdtEA,cvrES,cvrEA
  COMPLEX(kind=8) :: brES,brEA,cbrES,cbrEA,sES,sEA,drsES,drsEA,pES,pEA
  COMPLEX(kind=8) :: dsdtES,dsdtEA
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

  !COMPLEX(kind=8) :: infield(lm_max,3),outfield(ncp,nfs,3)
  !-- end of declaration
  !----------------------------------------------------------------------
  !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  nThetaNHS=(nThetaStart-1)/2

  IF ( nBc == 0 .OR. lDeriv ) THEN ! not a boundary or derivs required
     !PERFON('TFG_inn')
     !PERFON('TFG_thl')
     DO nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
        nThetaS  =nThetaN+1      ! same theta but for southern HS
        nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point

        PERFON_I('TFG_1')
        IF ( l_heat ) THEN
           ! the original version has shown to be the fastest
           ! putting 
           DO mc=1,n_m_max
              lmS=lStop(mc)
              sES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
              sEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
              DO lm=lStart(mc),lmS-1,2
                 sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
              sc(mc,nThetaN)=sES+sEA
              sc(mc,nThetaS)=sES-sEA
           END DO

           IF ( lFluxProfCalc ) THEN
              DO mc=1,n_m_max
                  lmS=lStop(mc)
                  pES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
                  pEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
                  DO lm=lStart(mc),lmS-1,2
                     pES=pES+leg_helper%preR(lm)  *Plm(lm,nThetaNHS)
                     pEA=pEA+leg_helper%preR(lm+1)*Plm(lm+1,nThetaNHS)
                  END DO
                  IF ( lmOdd(mc) ) pES=pES+leg_helper%preR(lmS)*Plm(lmS,nThetaNHS)
                  pc(mc,nThetaN)=pES+pEA
                  pc(mc,nThetaS)=pES-pEA
               END DO
           END IF

           IF ( lViscBcCalc ) THEN
              DO mc=1,n_m_max
                 dm =D_mc2m(mc)
                 lmS=lStop(mc)
                 dsdtES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 dsdtEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 DO lm=lStart(mc),lmS-1,2
                    dsdtEA =dsdtEA + leg_helper%sR(lm)*  dPlm(lm,nThetaNHS)
                    dsdtES =dsdtES + leg_helper%sR(lm+1)*dPlm(lm+1,nThetaNHS)
                 END DO
                 IF ( lmOdd(mc) ) THEN
                    dsdtEA =dsdtEA + leg_helper%sR(lmS)*dPlm(lmS,nThetaNHS)
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
        PERFOFF_I
        !--- Loop over all orders m: (numbered by mc)
        PERFON_I('TFG_2')
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
              cvrES  =cvrES   +  leg_helper%dLhz(lm)  *Plm(lm,nThetaNHS)
              dvrdrES=dvrdrES + leg_helper%dLhdw(lm)  *Plm(lm,nThetaNHS)
              brES   =brES    +  leg_helper%dLhb(lm)  *Plm(lm,nThetaNHS)
              cvrEA  =cvrEA   +  leg_helper%dLhz(lm+1)*Plm(lm+1,nThetaNHS)
              dvrdrEA=dvrdrEA + leg_helper%dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
              brEA   =brEA    +  leg_helper%dLhb(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              cvrES  =cvrES   +  leg_helper%dLhz(lmS)*Plm(lmS,nThetaNHS)
              dvrdrES=dvrdrES + leg_helper%dLhdw(lmS)*Plm(lmS,nThetaNHS)
              brES   =brES    +  leg_helper%dLhb(lmS)*Plm(lmS,nThetaNHS)
           END IF
           cvrc(mc,nThetaN)  =cvrES  +cvrEA
           cvrc(mc,nThetaS)  =cvrES  -cvrEA
           dvrdrc(mc,nThetaN)=dvrdrES+dvrdrEA
           dvrdrc(mc,nThetaS)=dvrdrES-dvrdrEA
           brc(mc,nThetaN)   =brES   +brEA
           brc(mc,nThetaS)   =brES   -brEA
        END DO
        PERFOFF_I
        PERFON_I('TFG_3')
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
              vrES    =vrES    + leg_helper%dLhw(lm)*   Plm(lm,nThetaNHS)
              dvrdtEA =dvrdtEA + leg_helper%dLhw(lm)*  dPlm(lm,nThetaNHS)
              cbrES   =cbrES   + leg_helper%dLhj(lm)*   Plm(lm,nThetaNHS)
              PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
              vrEA    =vrEA    + leg_helper%dLhw(lm+1)* Plm(lm+1,nThetaNHS)
              dvrdtES =dvrdtES + leg_helper%dLhw(lm+1)*dPlm(lm+1,nThetaNHS)
              cbrEA   =cbrEA   + leg_helper%dLhj(lm+1)* Plm(lm+1,nThetaNHS)
              PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
              PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              vrES    =vrES    + leg_helper%dLhw(lmS)* Plm(lmS,nThetaNHS)
              dvrdtEA =dvrdtEA + leg_helper%dLhw(lmS)*dPlm(lmS,nThetaNHS)
              cbrES   =cbrES   + leg_helper%dLhj(lmS)* Plm(lmS,nThetaNHS)
              PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
              PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
           END IF
           vrc(mc,nThetaN)   =vrES   +vrEA
           vrc(mc,nThetaS)   =vrES   -vrEA
           !WRITE(*,"(2(A,I3),A,2ES20.12)") "vrc(",mc,",",nThetaN,") = ",vrc(mc,nThetaN)
           !WRITE(*,"(2(A,I3),A,2ES20.12)") "vrc(",mc,",",nThetaS,") = ",vrc(mc,nThetaS)
           dvrdtc(mc,nThetaN)=dvrdtES+dvrdtEA
           dvrdtc(mc,nThetaS)=dvrdtES-dvrdtEA
           cbrc(mc,nThetaN)  =cbrES  +cbrEA
           cbrc(mc,nThetaS)  =cbrES  -cbrEA
        END DO
        PERFOFF_I
        PERFON_I('TFG_4')

        !--- Now the stuff using generalized harmonics:
        DO mc=1,n_m_max
           lmS=lStop(mc)
           vhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              vhN1=vhN1+leg_helper%vhG(lm)*PlmG(lm)+leg_helper%vhG(lm+1)*PlmG(lm+1)
              vhS1=vhS1-leg_helper%vhG(lm)*PlmC(lm)+leg_helper%vhG(lm+1)*PlmC(lm+1)
              vhN2=vhN2+leg_helper%vhC(lm)*PlmC(lm)+leg_helper%vhC(lm+1)*PlmC(lm+1)
              vhS2=vhS2-leg_helper%vhC(lm)*PlmG(lm)+leg_helper%vhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
              vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
              vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
              vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
           END IF
           vhN1M(mc)=0.5D0*vhN1
           vhS1M(mc)=0.5D0*vhS1
           vhN2M(mc)=0.5D0*vhN2
           vhS2M(mc)=0.5D0*vhS2
        END DO
        PERFOFF_I
        PERFON_I('TFG_5')

        DO mc=1,n_m_max
           lmS=lStop(mc)
           dvhdrN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           dvhdrS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lm)  *PlmG(lm) + &
                   leg_helper%dvhdrG(lm+1)*PlmG(lm+1)
              dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lm)  *PlmC(lm) + &
                   leg_helper%dvhdrG(lm+1)*PlmC(lm+1)
              dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lm)  *PlmC(lm) + &
                   leg_helper%dvhdrC(lm+1)*PlmC(lm+1)
              dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lm)  *PlmG(lm) + &
                   leg_helper%dvhdrC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lmS)*PlmG(lmS)
              dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lmS)*PlmC(lmS)
              dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lmS)*PlmC(lmS)
              dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lmS)*PlmG(lmS)
           END IF
           dvhdrN1M(mc)=0.5D0*dvhdrN1
           dvhdrS1M(mc)=0.5D0*dvhdrS1
           dvhdrN2M(mc)=0.5D0*dvhdrN2
           dvhdrS2M(mc)=0.5D0*dvhdrS2
        END DO
        PERFOFF_I
        PERFON_I('TFG_6')

        DO mc=1,n_m_max
           lmS=lStop(mc)
           bhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           bhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           bhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           bhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              bhN1=bhN1+leg_helper%bhG(lm)*PlmG(lm)+leg_helper%bhG(lm+1)*PlmG(lm+1)
              bhS1=bhS1-leg_helper%bhG(lm)*PlmC(lm)+leg_helper%bhG(lm+1)*PlmC(lm+1)
              bhN2=bhN2+leg_helper%bhC(lm)*PlmC(lm)+leg_helper%bhC(lm+1)*PlmC(lm+1)
              bhS2=bhS2-leg_helper%bhC(lm)*PlmG(lm)+leg_helper%bhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              bhN1=bhN1+leg_helper%bhG(lmS)*PlmG(lmS)
              bhS1=bhS1-leg_helper%bhG(lmS)*PlmC(lmS)
              bhN2=bhN2+leg_helper%bhC(lmS)*PlmC(lmS)
              bhS2=bhS2-leg_helper%bhC(lmS)*PlmG(lmS)
           END IF
           bhN1M(mc)=0.5D0*bhN1
           bhS1M(mc)=0.5D0*bhS1
           bhN2M(mc)=0.5D0*bhN2
           bhS2M(mc)=0.5D0*bhS2
        END DO
        PERFOFF_I
        PERFON_I('TFG_7')

        DO mc=1,n_m_max
           lmS=lStop(mc)
           cbhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           cbhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           cbhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           cbhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           !--- 8 add/mult, 20 DBLE words
           DO lm=lStart(mc),lmS-1,2
              cbhN1=cbhN1+leg_helper%cbhG(lm)*PlmG(lm)+leg_helper%cbhG(lm+1)*PlmG(lm+1)
              cbhS1=cbhS1-leg_helper%cbhG(lm)*PlmC(lm)+leg_helper%cbhG(lm+1)*PlmC(lm+1)
              cbhN2=cbhN2+leg_helper%cbhC(lm)*PlmC(lm)+leg_helper%cbhC(lm+1)*PlmC(lm+1)
              cbhS2=cbhS2-leg_helper%cbhC(lm)*PlmG(lm)+leg_helper%cbhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              cbhN1=cbhN1+leg_helper%cbhG(lmS)*PlmG(lmS)
              cbhS1=cbhS1-leg_helper%cbhG(lmS)*PlmC(lmS)
              cbhN2=cbhN2+leg_helper%cbhC(lmS)*PlmC(lmS)
              cbhS2=cbhS2-leg_helper%cbhC(lmS)*PlmG(lmS)
           END IF
           cbhN1M(mc)=0.5D0*cbhN1
           cbhS1M(mc)=0.5D0*cbhS1
           cbhN2M(mc)=0.5D0*cbhN2
           cbhS2M(mc)=0.5D0*cbhS2
        END DO
        PERFOFF_I
        PERFON_I('TFG_8')

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
        PERFOFF_I
        PERFON_I('TFG_9')

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
        PERFOFF_I
     END DO      ! End global loop over nTheta

     !PERFOFF

     !-- Zero out terms with index mc > n_m_max:
     IF ( n_m_max < ncp ) THEN
        DO nThetaN=1,sizeThetaB
           DO mc=n_m_max+1,ncp
              sc(mc,nThetaN)    =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              IF ( lViscBcCalc) THEN
                 dsdtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 dsdpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              END IF
              IF ( lFluxProfCalc ) THEN
                 pc(mc,nThetaN)    =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
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
     !PERFOFF

  ELSE   ! boundary ?

     !PERFON('TFG_bnd')
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
                 sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
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
              brES=brES + leg_helper%dLhb(lm)  *Plm(lm,nThetaNHS)
              brEA=brEA + leg_helper%dLhb(lm+1)*Plm(lm+1,nThetaNHS)
              PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
              PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
              PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              brES=brES+leg_helper%dLhb(lm)*Plm(lm,nThetaNHS)
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
              bhN1=bhN1+leg_helper%bhG(lm)*PlmG(lm)+leg_helper%bhG(lm+1)*PlmG(lm+1)
              bhS1=bhS1-leg_helper%bhG(lm)*PlmC(lm)+leg_helper%bhG(lm+1)*PlmC(lm+1)
              bhN2=bhN2+leg_helper%bhC(lm)*PlmC(lm)+leg_helper%bhC(lm+1)*PlmC(lm+1)
              bhS2=bhS2-leg_helper%bhC(lm)*PlmG(lm)+leg_helper%bhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              bhN1=bhN1+leg_helper%bhG(lmS)*PlmG(lmS)
              bhS1=bhS1-leg_helper%bhG(lmS)*PlmC(lmS)
              bhN2=bhN2+leg_helper%bhC(lmS)*PlmC(lmS)
              bhS2=bhS2-leg_helper%bhC(lmS)*PlmG(lmS)
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
                 vhN1=vhN1+leg_helper%vhG(lm)*PlmG(lm)+leg_helper%vhG(lm+1)*PlmG(lm+1)
                 vhS1=vhS1-leg_helper%vhG(lm)*PlmC(lm)+leg_helper%vhG(lm+1)*PlmC(lm+1)
                 vhN2=vhN2+leg_helper%vhC(lm)*PlmC(lm)+leg_helper%vhC(lm+1)*PlmC(lm+1)
                 vhS2=vhS2-leg_helper%vhC(lm)*PlmG(lm)+leg_helper%vhC(lm+1)*PlmG(lm+1)
              END DO
              IF ( lmOdd(mc) ) THEN
                 vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                 vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                 vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                 vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
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
     !PERFOFF
  END IF  ! boundary ? nBc?


  IF ( l_HT .OR. lViscBcCalc ) THEN    ! For movie output !
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
              drsES=drsES+leg_helper%dsR(lm)*Plm(lm,nThetaNHS)
              drsEA=drsEA+leg_helper%dsR(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) drsES=drsES+leg_helper%dsR(lmS)*Plm(lmS,nThetaNHS)
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

 END SUBROUTINE legTFG
!*******************************************************************************
 SUBROUTINE legTFGnomag(nBc,lDeriv,lViscBcCalc,lFluxProfCalc, & 
     &                 nThetaStart,vrc,vtc,vpc,dvrdrc,dvtdrc, &
     &                 dvpdrc,cvrc,dvrdtc,dvrdpc,dvtdpc,      & 
     &                 dvpdpc,sc,drSc,dsdtc,dsdpc,pc,         &
     &                 leg_helper)

  !-- input:
  INTEGER,INTENT(IN) :: nBc
  LOGICAL,INTENT(IN) :: lDeriv,lViscBcCalc,lFluxProfCalc
  INTEGER,INTENT(IN) :: nThetaStart

  !----- Stuff precomputed in legPrep:
  type(leg_helper_t) :: leg_helper

  !-- output: field on grid (theta,m) for the radial grid point nR
  !           and equatorially symmetric and antisymmetric contribution
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: vrc,vtc,vpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dvrdrc,dvtdrc,dvpdrc,dvrdtc,dvrdpc,dvtdpc,dvpdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: cvrc,sc,drSc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: dsdtc,dsdpc
  COMPLEX(kind=8),DIMENSION(ncp,nfs),INTENT(OUT) :: pc

  !------ Legendre Polynomials in m_horizontal.F90
  REAL(kind=8) :: PlmG(lm_max)
  REAL(kind=8) :: PlmC(lm_max)

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
                 sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
              sc(mc,nThetaN)=sES+sEA
              sc(mc,nThetaS)=sES-sEA
           END DO


           IF ( lFluxProfCalc ) THEN
              DO mc=1,n_m_max
                 lmS=lStop(mc)
                 pES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! One equatorial symmetry
                 pEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))    ! The other equatorial symmetry
                 DO lm=lStart(mc),lmS-1,2
                    pES=pES+leg_helper%preR(lm)  *Plm(lm,nThetaNHS)
                    pEA=pEA+leg_helper%preR(lm+1)*Plm(lm+1,nThetaNHS)
                 END DO
                 IF ( lmOdd(mc) ) pES=pES+leg_helper%preR(lmS)*Plm(lmS,nThetaNHS)
                 pc(mc,nThetaN)=pES+pEA
                 pc(mc,nThetaS)=pES-pEA
              END DO
           END IF

           IF ( lViscBcCalc ) THEN
              DO mc=1,n_m_max
                 dm =D_mc2m(mc)
                 lmS=lStop(mc)
                 dsdtES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 dsdtEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
                 DO lm=lStart(mc),lmS-1,2
                    dsdtEA =dsdtEA + leg_helper%sR(lm)*  dPlm(lm,nThetaNHS)
                    dsdtES =dsdtES + leg_helper%sR(lm+1)*dPlm(lm+1,nThetaNHS)
                 END DO
                 IF ( lmOdd(mc) ) THEN
                    dsdtEA =dsdtEA + leg_helper%sR(lmS)*dPlm(lmS,nThetaNHS)
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
              cvrES  =cvrES   +  leg_helper%dLhz(lm)  *Plm(lm,nThetaNHS)
              dvrdrES=dvrdrES + leg_helper%dLhdw(lm)  *Plm(lm,nThetaNHS)
              cvrEA  =cvrEA   +  leg_helper%dLhz(lm+1)*Plm(lm+1,nThetaNHS)
              dvrdrEA=dvrdrEA + leg_helper%dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              cvrES  =cvrES   +  leg_helper%dLhz(lmS)*Plm(lmS,nThetaNHS)
              dvrdrES=dvrdrES + leg_helper%dLhdw(lmS)*Plm(lmS,nThetaNHS)
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
              vrES    =vrES    + leg_helper%dLhw(lm)*   Plm(lm,nThetaNHS)
              dvrdtEA =dvrdtEA + leg_helper%dLhw(lm)*  dPlm(lm,nThetaNHS)
              PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
              vrEA    =vrEA    + leg_helper%dLhw(lm+1)* Plm(lm+1,nThetaNHS)
              dvrdtES =dvrdtES + leg_helper%dLhw(lm+1)*dPlm(lm+1,nThetaNHS)
              PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
              PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) THEN
              vrES    =vrES    + leg_helper%dLhw(lmS)* Plm(lmS,nThetaNHS)
              dvrdtEA =dvrdtEA + leg_helper%dLhw(lmS)*dPlm(lmS,nThetaNHS)
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
              vhN1=vhN1+leg_helper%vhG(lm)*PlmG(lm)+leg_helper%vhG(lm+1)*PlmG(lm+1)
              vhS1=vhS1-leg_helper%vhG(lm)*PlmC(lm)+leg_helper%vhG(lm+1)*PlmC(lm+1)
              vhN2=vhN2+leg_helper%vhC(lm)*PlmC(lm)+leg_helper%vhC(lm+1)*PlmC(lm+1)
              vhS2=vhS2-leg_helper%vhC(lm)*PlmG(lm)+leg_helper%vhC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
              vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
              vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
              vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
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
              dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lm)  *PlmG(lm) + &
                   leg_helper%dvhdrG(lm+1)*PlmG(lm+1)
              dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lm)  *PlmC(lm) + &
                   leg_helper%dvhdrG(lm+1)*PlmC(lm+1)
              dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lm)  *PlmC(lm) + &
                   leg_helper%dvhdrC(lm+1)*PlmC(lm+1)
              dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lm)  *PlmG(lm) + &
                   leg_helper%dvhdrC(lm+1)*PlmG(lm+1)
           END DO
           IF ( lmOdd(mc) ) THEN
              dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lmS)*PlmG(lmS)
              dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lmS)*PlmC(lmS)
              dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lmS)*PlmC(lmS)
              dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lmS)*PlmG(lmS)
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
              IF ( lFluxProfCalc ) THEN
                 pc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              END IF
              IF ( lViscBcCalc) THEN
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
                 sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                 sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
              END DO
              IF ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
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
                 vhN1=vhN1+leg_helper%vhG(lm)*  PlmG(lm)+    &
                       &   leg_helper%vhG(lm+1)*PlmG(lm+1)
                 vhS1=vhS1-leg_helper%vhG(lm)*  PlmC(lm)+    &
                       &   leg_helper%vhG(lm+1)*PlmC(lm+1)
                 vhN2=vhN2+leg_helper%vhC(lm)  *PlmC(lm)+    &
                       &   leg_helper%vhC(lm+1)*PlmC(lm+1)
                 vhS2=vhS2-leg_helper%vhC(lm)  *PlmG(lm)+    &
                       &   leg_helper%vhC(lm+1)*PlmG(lm+1)
              END DO
              IF ( lmOdd(mc) ) THEN
                 PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                 PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
                 vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                 vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                 vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                 vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
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


  IF ( l_HT .OR. lViscBcCalc ) THEN    ! For movie output !
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
              drsES=drsES+leg_helper%dsR(lm)*Plm(lm,nThetaNHS)
              drsEA=drsEA+leg_helper%dsR(lm+1)*Plm(lm+1,nThetaNHS)
           END DO
           IF ( lmOdd(mc) ) drsES=drsES+leg_helper%dsR(lmS)*Plm(lmS,nThetaNHS)
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
 END SUBROUTINE legTFGnomag
!------------------------------------------------------------------------

#if 0
#if 1
  SUBROUTINE compute_l_sum_NS_3(nThetaN,nThetaS,infield,Plm_slice,outfield)
    INTEGER,PARAMETER :: nFields=3
    INTEGER,INTENT(IN) :: nThetaN,nThetaS
    COMPLEX(kind=8),INTENT(IN) :: infield(lm_max,nFields)
    REAL(kind=8),intent(IN) :: Plm_slice(lm_max)
    COMPLEX(kind=8),INTENT(INOUT),DIMENSION(ncp,nfs,nFields) :: outfield

    ! Local variables
    INTEGER :: mc,lmS,lm,nf
    COMPLEX(kind=8),DIMENSION(nFields) :: sum_sym,sum_asym

    DO mc=1,n_m_max
       lmS=lStop(mc)
       sum_sym  = CMPLX(0.D0,0.D0)    ! One equatorial symmetry
       sum_asym = CMPLX(0.D0,0.D0)    ! The other equatorial symmetry
       DO lm=lStart(mc),lmS-1,2
          DO nf=1,nFields
             sum_sym(nf) = sum_sym(nf) + infield(lm,nf) * Plm_slice(lm)
             sum_asym(nf)= sum_asym(nf)+ infield(lm+1,nf)*Plm_slice(lm+1)
          END DO
       END DO
       IF ( lmOdd(mc) ) THEN
          DO nf=1,nFields
             sum_sym(nf) = sum_sym(nf) + infield(lmS,nf) * Plm_slice(lmS)
          END DO
       END IF
       DO nf=1,nFields
          outfield(mc,nThetaN,nf) = sum_sym(nf)+sum_asym(nf)
          outfield(mc,nThetaS,nf) = sum_sym(nf)-sum_asym(nf)
       END DO
    END DO
  END SUBROUTINE compute_l_sum_NS_3

  SUBROUTINE compute_l_sum_NS_1(nThetaN,nThetaS,infield,Plm_slice,outfield)
    INTEGER,INTENT(IN) :: nThetaN,nThetaS
    COMPLEX(kind=8),INTENT(IN) :: infield(lm_max)
    REAL(kind=8),intent(IN) :: Plm_slice(lm_max)
    COMPLEX(kind=8),INTENT(INOUT),DIMENSION(ncp,nfs) :: outfield

    ! Local variables
    INTEGER :: mc,lmS,lm,nf
    COMPLEX(kind=8) :: sum_sym,sum_asym

    DO mc=1,n_m_max
       lmS=lStop(mc)
       sum_sym  = CMPLX(0.D0,0.D0)    ! One equatorial symmetry
       sum_asym = CMPLX(0.D0,0.D0)    ! The other equatorial symmetry
       DO lm=lStart(mc),lmS-1,2
          sum_sym = sum_sym + infield(lm) * Plm_slice(lm)
          sum_asym= sum_asym+ infield(lm+1)*Plm_slice(lm+1)
       END DO
       IF ( lmOdd(mc) ) THEN
          sum_sym = sum_sym + infield(lmS) * Plm_slice(lmS)
       END IF
       outfield(mc,nThetaN) = sum_sym+sum_asym
       outfield(mc,nThetaS) = sum_sym-sum_asym
    END DO
  END SUBROUTINE compute_l_sum_NS_1
#else
  SUBROUTINE compute_l_sum_NS(infield,Plm_slice,outN,outS)
    COMPLEX(kind=8),INTENT(IN) :: infield(lm_max)
    REAL(kind=8),intent(IN) :: Plm_slice(lm_max)
    COMPLEX(kind=8),INTENT(OUT),DIMENSION(ncp) :: outN,outS

    ! Local variables
    INTEGER :: mc,lmS,lm,mc_old
    !COMPLEX(kind=8) :: sES,sEA,temp
    COMPLEX(kind=8),dimension(lm_max) :: temp_field

    temp_field=infield*Plm_slice
    DO mc=1,n_m_max
       outN(mc)=SUM(temp_field(lStart(mc):lStop(mc)))
       outS(mc)=SUM(temp_field(lStart(mc):lStop(mc):2)) &
            & - SUM(temp_field(lStart(mc)+1:lStop(mc):2))
    END DO
#if 0
    sES=CMPLX(0.D0,0.D0)    ! One equatorial symmetry
    sEA=CMPLX(0.D0,0.D0)    ! The other equatorial symmetry
    mc_old=1
    DO lm=1,lm_max
       mc=lm2mc(lm)
       IF (mc.NE.mc_old) THEN
          ! next m
          outN(mc_old) = sES + sEA
          outS(mc_old) = sES - sEA
          sES=CMPLX(0.D0,0.D0)    ! One equatorial symmetry
          sEA=CMPLX(0.D0,0.D0)    ! The other equatorial symmetry
          mc_old = mc
       END IF
       temp = infield(lm)*Plm_slice(lm)
       IF (MODULO(lm-lStart(mc),2).EQ.0) THEN
          sES = sES + temp
       ELSE
          sEA = sEA + temp
       END IF
    END DO
    outN(mc_old) = sES + sEA
    outS(mc_old) = sES - sEA
#endif
  END SUBROUTINE compute_l_sum_NS
#endif
#endif
END MODULE legendre_trafo
  !------------------------------------------------------------------------
