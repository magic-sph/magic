!$Id$
!*******************************************************************************
SUBROUTINE legTF(dLhw,vhG,vhC,dLhz,cvhG,cvhC, &
     &           l_max,minc,nThetaStart,sizeThetaB, &
     &           Plm,dPlm,lm_max,ncp,lHor,lDeriv, &
     &           vrc,vtc,vpc,cvrc,cvtc,cvpc)
  !*******************************************************************************

  !-----------------------------------------------------------------------------------

  !    'Legendre transform' from (nR,l,m) to (nR,nTheta,m) [spectral to grid]
  !    where nTheta numbers the colatitudes and l and m are degree and
  !    order of the spherical harmonic representation.

  !    Calculates all three spherical components vrc,vtc,vpc of a field as
  !    well as its curl (cvrc,cvtc,cvpc) that is given a spherical harmonis poloidal
  !    toroidal decomposition. s_legPrep.f has to be called first and
  !    provides the input fields dLhW, .....
  !    The symmetry properties of the P_lm with respect to the equator
  !    are used. The equatorially anti-symmetric (EA) contribution is
  !    added to (subracted from ) the equatorially symmetric (ES) contribution
  !    in northern (southern) hemisphere.

  !    Output is given for all sizeThetaB colatitudes in a colatitude block
  !    that starts with colatitude nThetaStart. At output, each component
  !    in the northern hemisphere is followed by the component in the
  !    southern hemisphere.
  !    The Plms and dPlms=sin(theta) d Plm / d theta are only given
  !    for the colatitudes in the northern hemisphere.

  !     dLhw,....,cvhC : (input) arrays provided by s_legPrep.f
  !     l_max          : (input) maximum spherical harmonic degree
  !     minc           : (input) azimuthal symmetry
  !     nThetaStart    : (input) transformation is done for the range of
  !                      points nThetaStart <= nTheta <= nThetaStart-1+sizeThetaB
  !     sizeThetaB     : (input) size theta block
  !     Plm            : (input) associated Legendre polynomials
  !     dPlm           : (input) sin(theta) d Plm / d theta
  !     lm_max         : leading dimension of Plm and dPlm
  !     lHor=.TRUE.    : (input) calculate horizontal componenst
  !     lDeric=.TRUE.  : (input) calculate curl of field
  !     vrc, ....,cvpc : (output) components in (nTheta,m)-space
  !     ncp            : (input) leading dimension of complex vrc

  !-----------------------------------------------------------------------------------

  USE truncation, ONLY: n_m_max
  IMPLICIT NONE

  !-- INPUT:

  !----- Stuff precomputed in legPrep:
  COMPLEX(kind=8),INTENT(IN) :: dLhw(*),dLhz(*)
  COMPLEX(kind=8),INTENT(IN) :: vhG(*),vhC(*)
  COMPLEX(kind=8),INTENT(IN) :: cvhG(*),cvhC(*)

  !----- Defining theta block
  INTEGER,INTENT(IN) :: nThetaStart,sizeThetaB

  !------ Legendre Polynomials in c_horizontal.f
  INTEGER,INTENT(IN) :: l_max,minc,lm_max
  REAL(kind=8),INTENT(IN) :: Plm(lm_max,*)
  REAL(kind=8),INTENT(IN) :: dPlm(lm_max,*)

  !----- What should I calculate?
  LOGICAL,INTENT(IN) :: lHor
  LOGICAL,INTENT(IN) :: lDeriv

  INTEGER,INTENT(IN) :: ncp

  !-- Output: field on grid (theta,m) for the radial grid point nR
  !           and equatorially symmetric and antisymmetric contribution
  COMPLEX(kind=8),INTENT(OUT) :: vrc(ncp,*)
  COMPLEX(kind=8),INTENT(OUT) :: vtc(ncp,*)
  COMPLEX(kind=8),INTENT(OUT) :: vpc(ncp,*)
  COMPLEX(kind=8),INTENT(OUT) :: cvrc(ncp,*)
  COMPLEX(kind=8),INTENT(OUT) :: cvtc(ncp,*)
  COMPLEX(kind=8),INTENT(OUT) :: cvpc(ncp,*)

  !-- Local:
  COMPLEX(kind=8) :: vrES,vrEA,cvrES,cvrEA
  INTEGER :: nThetaN,nThetaS,nThetaNHS
  INTEGER :: m,l,mc,lm
  REAL(kind=8) :: PlmG(lm_max)
  REAL(kind=8) :: PlmC(lm_max)

  COMPLEX(kind=8) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
  COMPLEX(kind=8) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
  COMPLEX(kind=8) :: cvhN1M(n_m_max),cvhN2M(n_m_max)
  COMPLEX(kind=8) :: cvhS1M(n_m_max),cvhS2M(n_m_max)

  !-- end of declaration
  !----------------------------------------------------------------------

  !-- Theta blocking possible here.
  !   Note that the theta order mixed, i.e. each
  !   value for a theta in the northern hemisphere is
  !   followed by its counterpart in the southern hemisphere.
  nThetaNHS=(nThetaStart-1)/2
  DO nThetaN=1,sizeThetaB,2 ! Loop over thetas for north HS
     nThetaS=nThetaN+1
     nThetaNHS=nThetaNHS+1

     !--- Loop over all oders m: (numbered by mc)
     lm=0
     mc=0
     DO m=0,l_max,minc
        mc=mc+1
        IF ( mc > ncp ) THEN
           WRITE(*,*) 'ncp too small in s_legTF!'
           WRITE(*,*) 'Increase ncp in calling routine!'
           STOP
        END IF
        vrES =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        vrEA =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvrES=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        cvrEA=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
        DO l=m,l_max-1,2
           lm=lm+2
           vrES=vrES+dLhw(lm-1)*Plm(lm-1,nThetaNHS)
           IF ( lHor ) THEN
              PlmG(lm-1)=dPlm(lm-1,nThetaNHS) - &
                   m*Plm(lm-1,nThetaNHS)
              PlmC(lm-1)=dPlm(lm-1,nThetaNHS) + &
                   m*Plm(lm-1,nThetaNHS)
           END IF
           IF ( lDeriv ) cvrES=cvrES+dLhz(lm-1)*Plm(lm-1,nThetaNHS)
           vrEA=vrEA+dLhw(lm)*Plm(lm,nThetaNHS)
           IF ( lHor ) THEN
              PlmG(lm)=dPlm(lm,nThetaNHS)-m*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+m*Plm(lm,nThetaNHS)
           END IF
           IF ( lDeriv ) cvrEA=cvrEA+dLhz(lm)*Plm(lm,nThetaNHS)
        END DO
        IF ( MOD(l_max-m+1,2) == 1 ) THEN
           lm=lm+1
           vrES=vrES+dLhw(lm)*Plm(lm,nThetaNHS)
           IF ( lHor ) THEN
              PlmG(lm)=dPlm(lm,nThetaNHS)-m*Plm(lm,nThetaNHS)
              PlmC(lm)=dPlm(lm,nThetaNHS)+m*Plm(lm,nThetaNHS)
           END IF
           IF ( lDeriv ) cvrES=cvrES+dLhz(lm)*Plm(lm,nThetaNHS)
        END IF
        vrc(mc,nThetaN)=vrES+vrEA
        vrc(mc,nThetaS)=vrES-vrEA
        IF ( lDeriv ) THEN
           cvrc(mc,nThetaN)=cvrES+cvrEA
           cvrc(mc,nThetaS)=cvrES-cvrEA
        END IF
     END DO
     n_m_max=mc

     !--- Now the stuff using generalized harmonics:
     IF ( lHor ) THEN
        lm=0
        mc=0
        DO m=0,l_max,minc
           mc=mc+1
           vhN1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS1=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhN2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS2=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           DO l=m,l_max-1,2
              lm=lm+2
              vhN1=vhN1+vhG(lm-1)*PlmG(lm-1)+vhG(lm)*PlmG(lm)
              vhS1=vhS1-vhG(lm-1)*PlmC(lm-1)+vhG(lm)*PlmC(lm)
              vhN2=vhN2+vhC(lm-1)*PlmC(lm-1)+vhC(lm)*PlmC(lm)
              vhS2=vhS2-vhC(lm-1)*PlmG(lm-1)+vhC(lm)*PlmG(lm)
           END DO
           IF ( MOD(l_max-m+1,2) == 1 ) THEN
              lm=lm+1
              vhN1=vhN1+vhG(lm)*PlmG(lm)
              vhS1=vhS1-vhG(lm)*PlmC(lm)
              vhN2=vhN2+vhC(lm)*PlmC(lm)
              vhS2=vhS2-vhC(lm)*PlmG(lm)
           END IF
           vhN1M(mc)=0.5D0*vhN1
           vhS1M(mc)=0.5D0*vhS1
           vhN2M(mc)=0.5D0*vhN2
           vhS2M(mc)=0.5D0*vhS2
        END DO
        !--- Unscramble:
        n_m_max=mc
        DO mc=1,n_m_max
           vtc(mc,nThetaN)=vhN1M(mc)+vhN2M(mc)
           vhN            =vhN1M(mc)-vhN2M(mc)
           vtc(mc,nThetaS)=vhS1M(mc)+vhS2M(mc)
           vhS            =vhS1M(mc)-vhS2M(mc)
           vpc(mc,nThetaN)=CMPLX(AIMAG(vhN),-REAL(vhN),KIND=KIND(0d0))
           vpc(mc,nThetaS)=CMPLX(AIMAG(vhS),-REAL(vhS),KIND=KIND(0d0))
        END DO
     END IF

     IF ( lDeriv ) THEN
        lm=0
        mc=0
        DO m=0,l_max,minc
           mc=mc+1
           vhN1 =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS1 =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhN2 =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           vhS2 =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           DO l=m,l_max-1,2
              lm=lm+2
              vhN1=vhN1+cvhG(lm-1)*PlmG(lm-1)+cvhG(lm)*PlmG(lm)
              vhS1=vhS1-cvhG(lm-1)*PlmC(lm-1)+cvhG(lm)*PlmC(lm)
              vhN2=vhN2+cvhC(lm-1)*PlmC(lm-1)+cvhC(lm)*PlmC(lm)
              vhS2=vhS2-cvhC(lm-1)*PlmG(lm-1)+cvhC(lm)*PlmG(lm)
           END DO
           IF ( MOD(l_max-m+1,2) == 1 ) THEN
              lm=lm+1
              vhN1=vhN1+cvhG(lm)*PlmG(lm)
              vhS1=vhS1-cvhG(lm)*PlmC(lm)
              vhN2=vhN2+cvhC(lm)*PlmC(lm)
              vhS2=vhS2-cvhC(lm)*PlmG(lm)
           END IF
           cvhN1M(mc)=0.5D0*vhN1
           cvhS1M(mc)=0.5D0*vhS1
           cvhN2M(mc)=0.5D0*vhN2
           cvhS2M(mc)=0.5D0*vhS2
        END DO
        n_m_max=mc
        DO mc=1,n_m_max
           cvtc(mc,nThetaN)=cvhN1M(mc)+cvhN2M(mc)
           vhN             =cvhN1M(mc)-cvhN2M(mc)
           cvtc(mc,nThetaS)=cvhS1M(mc)+cvhS2M(mc)
           vhS             =cvhS1M(mc)-cvhS2M(mc)
           cvpc(mc,nThetaN)=CMPLX(AIMAG(vhN),-REAL(vhN),KIND=KIND(0d0))
           cvpc(mc,nThetaS)=CMPLX(AIMAG(vhS),-REAL(vhS),KIND=KIND(0d0))
        END DO
     END IF

  END DO

  !-- Zero out terms with index mc > n_m_max:
  IF ( n_m_max < ncp ) THEN
     DO nThetaN=1,sizeThetaB
        DO mc=n_m_max+1,ncp
           vrc(mc,nThetaN) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           IF ( lHor ) THEN
              vtc(mc,nThetaN) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              vpc(mc,nThetaN) =CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
           IF ( lDeriv ) THEN
              cvrc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              cvtc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
              cvpc(mc,nThetaN)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
           END IF
        END DO
     END DO  ! loop over nThetaN (theta)
  END IF


  RETURN
end SUBROUTINE legTF

!------------------------------------------------------------------------
