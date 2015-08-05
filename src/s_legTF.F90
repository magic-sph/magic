!$Id$
subroutine legTF(dLhw,vhG,vhC,dLhz,cvhG,cvhC,l_max,minc, &
           &     nThetaStart,sizeThetaB,Plm,dPlm,lHor,   &
           &     lDeriv,vrc,vtc,vpc,cvrc,cvtc,cvpc)
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
   !     lHor=.true.    : (input) calculate horizontal componenst
   !     lDeric=.true.  : (input) calculate curl of field
   !     vrc, ....,cvpc : (output) components in (nTheta,m)-space
 
   !-----------------------------------------------------------------------------------
 
   use truncation, only: n_m_max, nrp, lm_max
 
   implicit none
 
   !-- Input:
 
   !----- Stuff precomputed in legPrep:
   complex(kind=8), intent(in) :: dLhw(*),dLhz(*)
   complex(kind=8), intent(in) :: vhG(*),vhC(*)
   complex(kind=8), intent(in) :: cvhG(*),cvhC(*)
 
   !----- Defining theta block
   integer,         intent(in) :: nThetaStart,sizeThetaB
 
   !------ Legendre Polynomials in c_horizontal.f
   integer,         intent(in) :: l_max,minc
   real(kind=8),    intent(in) :: Plm(lm_max,*)
   real(kind=8),    intent(in) :: dPlm(lm_max,*)
 
   !----- What should I calculate?
   logical,         intent(in) :: lHor
   logical,         intent(in) :: lDeriv
 
   !-- Output: field on grid (theta,m) for the radial grid point nR
   !           and equatorially symmetric and antisymmetric contribution
   real(kind=8), intent(out) :: vrc(nrp,*)
   real(kind=8), intent(out) :: vtc(nrp,*)
   real(kind=8), intent(out) :: vpc(nrp,*)
   real(kind=8), intent(out) :: cvrc(nrp,*)
   real(kind=8), intent(out) :: cvtc(nrp,*)
   real(kind=8), intent(out) :: cvpc(nrp,*)
 
   !-- Local variables:
   complex(kind=8) :: vrES,vrEA,cvrES,cvrEA
   integer :: nThetaN,nThetaS,nThetaNHS
   integer :: m,l,mc,lm
   real(kind=8) :: PlmG(lm_max)
   real(kind=8) :: PlmC(lm_max)
 
   complex(kind=8) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
   complex(kind=8) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
   complex(kind=8) :: cvhN1M(n_m_max),cvhN2M(n_m_max)
   complex(kind=8) :: cvhS1M(n_m_max),cvhS2M(n_m_max)
 
   !-- Theta blocking possible here.
   !   Note that the theta order mixed, i.e. each
   !   value for a theta in the northern hemisphere is
   !   followed by its counterpart in the southern hemisphere.
   nThetaNHS=(nThetaStart-1)/2
   do nThetaN=1,sizeThetaB,2 ! Loop over thetas for north HS
      nThetaS=nThetaN+1
      nThetaNHS=nThetaNHS+1
 
      !--- Loop over all oders m: (numbered by mc)
      lm=0
      mc=0
      do m=0,l_max,minc
         mc=mc+1
         if ( mc > nrp/2 ) then
            write(*,*) 'nrp too small in s_legTF!'
            write(*,*) 'Increase nrp in calling routine!'
            stop
         end if
         vrES =cmplx(0.D0,0.D0,kind=kind(0d0))
         vrEA =cmplx(0.D0,0.D0,kind=kind(0d0))
         cvrES=cmplx(0.D0,0.D0,kind=kind(0d0))
         cvrEA=cmplx(0.D0,0.D0,kind=kind(0d0))
         do l=m,l_max-1,2
            lm=lm+2
            vrES=vrES+dLhw(lm-1)*Plm(lm-1,nThetaNHS)
            if ( lHor ) then
               PlmG(lm-1)=dPlm(lm-1,nThetaNHS) - m*Plm(lm-1,nThetaNHS)
               PlmC(lm-1)=dPlm(lm-1,nThetaNHS) + m*Plm(lm-1,nThetaNHS)
            end if
            if ( lDeriv ) cvrES=cvrES+dLhz(lm-1)*Plm(lm-1,nThetaNHS)
            vrEA=vrEA+dLhw(lm)*Plm(lm,nThetaNHS)
            if ( lHor ) then
               PlmG(lm)=dPlm(lm,nThetaNHS)-m*Plm(lm,nThetaNHS)
               PlmC(lm)=dPlm(lm,nThetaNHS)+m*Plm(lm,nThetaNHS)
            end if
            if ( lDeriv ) cvrEA=cvrEA+dLhz(lm)*Plm(lm,nThetaNHS)
         end do
         if ( mod(l_max-m+1,2) == 1 ) then
            lm=lm+1
            vrES=vrES+dLhw(lm)*Plm(lm,nThetaNHS)
            if ( lHor ) then
               PlmG(lm)=dPlm(lm,nThetaNHS)-m*Plm(lm,nThetaNHS)
               PlmC(lm)=dPlm(lm,nThetaNHS)+m*Plm(lm,nThetaNHS)
            end if
            if ( lDeriv ) cvrES=cvrES+dLhz(lm)*Plm(lm,nThetaNHS)
         end if
         vrc(2*mc-1,nThetaN)= real(vrES+vrEA)
         vrc(2*mc  ,nThetaN)=aimag(vrES+vrEA)
         vrc(2*mc-1,nThetaS)= real(vrES-vrEA)
         vrc(2*mc  ,nThetaS)=aimag(vrES-vrEA)
         if ( lDeriv ) then
            cvrc(2*mc-1,nThetaN)= real(cvrES+cvrEA)
            cvrc(2*mc  ,nThetaN)=aimag(cvrES+cvrEA)
            cvrc(2*mc-1,nThetaS)= real(cvrES-cvrEA)
            cvrc(2*mc  ,nThetaS)=aimag(cvrES-cvrEA)
         end if
      end do
      n_m_max=mc
 
      !--- Now the stuff using generalized harmonics:
      if ( lHor ) then
         lm=0
         mc=0
         do m=0,l_max,minc
            mc=mc+1
            vhN1=cmplx(0.D0,0.D0,kind=kind(0d0))
            vhS1=cmplx(0.D0,0.D0,kind=kind(0d0))
            vhN2=cmplx(0.D0,0.D0,kind=kind(0d0))
            vhS2=cmplx(0.D0,0.D0,kind=kind(0d0))
            do l=m,l_max-1,2
               lm=lm+2
               vhN1=vhN1+vhG(lm-1)*PlmG(lm-1)+vhG(lm)*PlmG(lm)
               vhS1=vhS1-vhG(lm-1)*PlmC(lm-1)+vhG(lm)*PlmC(lm)
               vhN2=vhN2+vhC(lm-1)*PlmC(lm-1)+vhC(lm)*PlmC(lm)
               vhS2=vhS2-vhC(lm-1)*PlmG(lm-1)+vhC(lm)*PlmG(lm)
            end do
            if ( mod(l_max-m+1,2) == 1 ) then
               lm=lm+1
               vhN1=vhN1+vhG(lm)*PlmG(lm)
               vhS1=vhS1-vhG(lm)*PlmC(lm)
               vhN2=vhN2+vhC(lm)*PlmC(lm)
               vhS2=vhS2-vhC(lm)*PlmG(lm)
            end if
            vhN1M(mc)=0.5D0*vhN1
            vhS1M(mc)=0.5D0*vhS1
            vhN2M(mc)=0.5D0*vhN2
            vhS2M(mc)=0.5D0*vhS2
         end do
         !--- Unscramble:
         n_m_max=mc
         do mc=1,n_m_max
            vtc(2*mc-1,nThetaN)= real(vhN1M(mc)+vhN2M(mc))
            vtc(2*mc  ,nThetaN)=aimag(vhN1M(mc)+vhN2M(mc))
            vhN                =vhN1M(mc)-vhN2M(mc)
            vtc(2*mc-1,nThetaS)= real(vhS1M(mc)+vhS2M(mc))
            vtc(2*mc  ,nThetaS)=aimag(vhS1M(mc)+vhS2M(mc))
            vhS                =vhS1M(mc)-vhS2M(mc)
            vpc(2*mc-1,nThetaN)=aimag(vhN)
            vpc(2*mc  ,nThetaN)=-real(vhN)
            vpc(2*mc-1,nThetaS)=aimag(vhS)
            vpc(2*mc  ,nThetaS)=-real(vhS)
         end do
      end if
 
      if ( lDeriv ) then
         lm=0
         mc=0
         do m=0,l_max,minc
            mc=mc+1
            vhN1 =cmplx(0.D0,0.D0,kind=kind(0d0))
            vhS1 =cmplx(0.D0,0.D0,kind=kind(0d0))
            vhN2 =cmplx(0.D0,0.D0,kind=kind(0d0))
            vhS2 =cmplx(0.D0,0.D0,kind=kind(0d0))
            do l=m,l_max-1,2
               lm=lm+2
               vhN1=vhN1+cvhG(lm-1)*PlmG(lm-1)+cvhG(lm)*PlmG(lm)
               vhS1=vhS1-cvhG(lm-1)*PlmC(lm-1)+cvhG(lm)*PlmC(lm)
               vhN2=vhN2+cvhC(lm-1)*PlmC(lm-1)+cvhC(lm)*PlmC(lm)
               vhS2=vhS2-cvhC(lm-1)*PlmG(lm-1)+cvhC(lm)*PlmG(lm)
            end do
            if ( mod(l_max-m+1,2) == 1 ) then
               lm=lm+1
               vhN1=vhN1+cvhG(lm)*PlmG(lm)
               vhS1=vhS1-cvhG(lm)*PlmC(lm)
               vhN2=vhN2+cvhC(lm)*PlmC(lm)
               vhS2=vhS2-cvhC(lm)*PlmG(lm)
            end if
            cvhN1M(mc)=0.5D0*vhN1
            cvhS1M(mc)=0.5D0*vhS1
            cvhN2M(mc)=0.5D0*vhN2
            cvhS2M(mc)=0.5D0*vhS2
         end do
         n_m_max=mc
         do mc=1,n_m_max
            cvtc(2*mc-1,nThetaN)= real(cvhN1M(mc)+cvhN2M(mc))
            cvtc(2*mc  ,nThetaN)=aimag(cvhN1M(mc)+cvhN2M(mc))
            vhN                 =cvhN1M(mc)-cvhN2M(mc)
            cvtc(2*mc-1,nThetaS)= real(cvhS1M(mc)+cvhS2M(mc))
            cvtc(2*mc  ,nThetaS)=aimag(cvhS1M(mc)+cvhS2M(mc))
            vhS                 =cvhS1M(mc)-cvhS2M(mc)
            cvpc(2*mc-1,nThetaN)=aimag(vhN)
            cvpc(2*mc  ,nThetaN)=-real(vhN)
            cvpc(2*mc-1,nThetaS)=aimag(vhS)
            cvpc(2*mc  ,nThetaS)=-real(vhS)
         end do
      end if
 
   end do
 
   !-- Zero out terms with index mc > n_m_max:
   if ( n_m_max < nrp/2 ) then
      do nThetaN=1,sizeThetaB
         do mc=2*n_m_max+1,nrp
            vrc(mc,nThetaN) =0.d0
            if ( lHor ) then
               vtc(mc,nThetaN) =0.d0
               vpc(mc,nThetaN) =0.d0
            end if
            if ( lDeriv ) then
               cvrc(mc,nThetaN)=0.d0
               cvtc(mc,nThetaN)=0.d0
               cvpc(mc,nThetaN)=0.d0
            end if
         end do
      end do  ! loop over nThetaN (theta)
   end if

end subroutine legTF
!------------------------------------------------------------------------
