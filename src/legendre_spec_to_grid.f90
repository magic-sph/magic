!$Id$
#include "perflib_preproc.cpp"
module legendre_spec_to_grid

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use precision_mod
   use truncation, only: lm_max, n_m_max, nrp, l_max
   use blocking, only: nfs, sizeThetaB, lm2mc, lm2
   use horizontal_data, only: Plm, dPlm, lStart, lStop, lmOdd, D_mc2m, &
                              osn2, dPlmt, osn1
   use logic, only: l_heat, l_ht
   use const, only: zero, half, one
   use parallel_mod, only: rank
   use leg_helper_mod, only: leg_helper_t

   implicit none
 
   private

   public :: legTFG, legTFGnomag, legTF, lmAS2pt

contains

   subroutine legTFG(nBc,lDeriv,lViscBcCalc,lFluxProfCalc,nThetaStart, &
     &               vrc,vtc,vpc,dvrdrc,dvtdrc,dvpdrc,cvrc,            &
     &               dvrdtc,dvrdpc,dvtdpc,dvpdpc,                      &
     &               brc,btc,bpc,cbrc,cbtc,cbpc,sc,                    &
     &               drSc,dsdtc,dsdpc,pc,leg_helper)
      !
      !    Legendre transform from (nR,l,m) to (nR,nTheta,m) [spectral to grid]
      !    where nTheta numbers the colatitudes and l is the degree of
      !    the spherical harmonic representation.
      !
      !    Transforms entropy, velocity and magnetic field components
      !    and terms involving spatial derivatives.
      !    The symmetry properties of the P_lm with respect to the equator
      !    are used. The equatorially anti-symmetric (EA) contribution is
      !    added to (subracted from ) the equatorially symmetric (ES) contribution
      !    in northern (southern) hemisphere.
      !
      !      * nBc            : (input) accounts for special conditions on radial boundaries
      !
      !         -nBc=2       : we are dealing with a no slip boundary, v_r and v_theta are
      !          zero and v_phi=r sin(theta) omega, where omega is the rotation rate of the 
      !          boundary (mantle of IC), only magn. field terms are calculated, v is 
      !          set later.
      !
      !         -nBc=1       : a free slip bounday: v_r is zero, derivatives of v and B 
      !                      are not needed, only components of v,B and entropy 
      !                      are calculated
      !
      !         -nBc=0       : normal case, interior grid point
      !
      !      * lDeriv=.true.  : (input) calculate derivatives
      !      * nThetaStart    : (input) transformation is done for the range of
      !        points nThetaStart <= nTheta <= nThetaStart-1+sizeThetaB
      !      * Plm            : associated Legendre polynomials
      !      * dPlm           : sin(theta) d Plm / d theta
      !      * osn2           : 1/sin(theta)^2
      !      * vrc, ...., drSc: (output) components in (nTheta,m)-space
      !      * dLhw,....,cbhC : (input) help arrays calculated in s_legPrep.f
      !
      
      !-- Input variables:
      integer, intent(in) :: nBc
      logical, intent(in) :: lDeriv,lFluxProfCalc,lViscBcCalc
      integer, intent(in) :: nThetaStart
    
      !----- Stuff precomputed in legPrep:
      type(leg_helper_t) :: leg_helper
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: vrc(nrp,nfs), vtc(nrp,nfs), vpc(nrp,nfs)
      real(cp), intent(out) :: dvrdrc(nrp,nfs), dvtdrc(nrp,nfs), dvpdrc(nrp,nfs)
      real(cp), intent(out) :: dvrdtc(nrp,nfs), dvrdpc(nrp,nfs)
      real(cp), intent(out) :: dvtdpc(nrp,nfs), dvpdpc(nrp, nfs)
      real(cp), intent(out) :: cvrc(nrp,nfs)
      real(cp), intent(out) :: brc(nrp,nfs), btc(nrp,nfs), bpc(nrp,nfs)
      real(cp), intent(out) :: cbrc(nrp,nfs), cbtc(nrp,nfs), cbpc(nrp,nfs)
      real(cp), intent(out) :: sc(nrp,nfs), drSc(nrp,nfs), pc(nrp,nfs)
      real(cp), intent(out) :: dsdtc(nrp,nfs), dsdpc(nrp,nfs)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max)
      real(cp) :: PlmC(lm_max)
    
      !-- Local variables:
      complex(cp) :: vrES,vrEA,dvrdrES,dvrdrEA,dvrdtES,dvrdtEA,cvrES,cvrEA
      complex(cp) :: brES,brEA,cbrES,cbrEA,sES,sEA,drsES,drsEA,pES,pEA
      complex(cp) :: dsdtES,dsdtEA
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS
      real(cp) :: dm
    
      complex(cp) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
      complex(cp) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
      complex(cp) :: dvhdrN1M(n_m_max),dvhdrN2M(n_m_max),dvhdrN
      complex(cp) :: dvhdrS1M(n_m_max),dvhdrS2M(n_m_max),dvhdrS
      complex(cp) :: dvhdrN1,dvhdrN2,dvhdrS1,dvhdrS2
      complex(cp) :: bhN1M(n_m_max),bhN2M(n_m_max),bhN,bhN1,bhN2
      complex(cp) :: bhS1M(n_m_max),bhS2M(n_m_max),bhS,bhS1,bhS2
      complex(cp) :: cbhN1M(n_m_max),cbhN2M(n_m_max),cbhN,cbhN1,cbhN2
      complex(cp) :: cbhS1M(n_m_max),cbhS2M(n_m_max),cbhS,cbhS1,cbhS2
    
      !call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
      nThetaNHS=(nThetaStart-1)/2
    
      if ( nBc == 0 .or. lDeriv ) then ! not a boundary or derivs required
         !PERFON('TFG_inn')
         !PERFON('TFG_thl')
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
    
            PERFON_I('TFG_1')
            if ( l_heat ) then
               ! the original version has shown to be the fastest
               ! putting 
               do mc=1,n_m_max
                  lmS=lStop(mc)
                  sES=zero  ! One equatorial symmetry
                  sEA=zero  ! The other equatorial symmetry
                  do lm=lStart(mc),lmS-1,2
                     sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                     sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
                  end do
                  if ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
                  sc(2*mc-1,nThetaN)= real(sES+sEA)
                  sc(2*mc  ,nThetaN)=aimag(sES+sEA)
                  sc(2*mc-1,nThetaS)= real(sES-sEA)
                  sc(2*mc  ,nThetaS)=aimag(sES-sEA)
               end do
    
               if ( lFluxProfCalc ) then
                  do mc=1,n_m_max
                      lmS=lStop(mc)
                      pES=zero ! One equatorial symmetry
                      pEA=zero ! The other equatorial symmetry
                      do lm=lStart(mc),lmS-1,2
                         pES=pES+leg_helper%preR(lm)  *Plm(lm,nThetaNHS)
                         pEA=pEA+leg_helper%preR(lm+1)*Plm(lm+1,nThetaNHS)
                      end do
                      if ( lmOdd(mc) ) pES=pES+leg_helper%preR(lmS)*Plm(lmS,nThetaNHS)
                      pc(2*mc-1,nThetaN)= real(pES+pEA)
                      pc(2*mc  ,nThetaN)=aimag(pES+pEA)
                      pc(2*mc-1,nThetaS)= real(pES-pEA)
                      pc(2*mc  ,nThetaS)=aimag(pES-pEA)
                   end do
               end if
    
               if ( lViscBcCalc ) then
                  do mc=1,n_m_max
                     dm =D_mc2m(mc)*osn1(nThetaNHS)
                     lmS=lStop(mc)
                     dsdtES=zero
                     dsdtEA=zero
                     do lm=lStart(mc),lmS-1,2
                        dsdtEA =dsdtEA + leg_helper%sR(lm)*  dPlmt(lm,nThetaNHS)
                        dsdtES =dsdtES + leg_helper%sR(lm+1)*dPlmt(lm+1,nThetaNHS)
                     end do
                     if ( lmOdd(mc) ) then
                        dsdtEA =dsdtEA + leg_helper%sR(lmS)*dPlmt(lmS,nThetaNHS)
                     end if
                     dsdtc(2*mc-1,nThetaN)= real(dsdtES+dsdtEA)
                     dsdtc(2*mc  ,nThetaN)=aimag(dsdtES+dsdtEA)
                     dsdtc(2*mc-1,nThetaS)= real(dsdtES-dsdtEA)
                     dsdtc(2*mc  ,nThetaS)=aimag(dsdtES-dsdtEA)
                  end do
    
                  do mc=1,n_m_max
                     dm=D_mc2m(mc)
                     dsdpc(2*mc-1,nThetaN)=-dm*sc(2*mc  ,nThetaN)
                     dsdpc(2*mc  ,nThetaN)= dm*sc(2*mc-1,nThetaN)
                     dsdpc(2*mc-1,nThetaS)=-dm*sc(2*mc  ,nThetaS)
                     dsdpc(2*mc  ,nThetaS)= dm*sc(2*mc-1,nThetaS)
                  end do
    
               end if ! thermal dissipation layer
    
            end if
            PERFOFF_I
            !--- Loop over all orders m: (numbered by mc)
            PERFON_I('TFG_2')
            do mc=1,n_m_max
               lmS=lStop(mc)
               cvrES  =zero
               cvrEA  =zero
               dvrdrES=zero
               dvrdrEA=zero
               brES   =zero
               brEA   =zero
               !--- 6 add/mult, 26 dble words
               do lm=lStart(mc),lmS-1,2
                  cvrES  =cvrES   +  leg_helper%dLhz(lm)  *Plm(lm,nThetaNHS)
                  dvrdrES=dvrdrES + leg_helper%dLhdw(lm)  *Plm(lm,nThetaNHS)
                  brES   =brES    +  leg_helper%dLhb(lm)  *Plm(lm,nThetaNHS)
                  cvrEA  =cvrEA   +  leg_helper%dLhz(lm+1)*Plm(lm+1,nThetaNHS)
                  dvrdrEA=dvrdrEA + leg_helper%dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
                  brEA   =brEA    +  leg_helper%dLhb(lm+1)*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  cvrES  =cvrES   +  leg_helper%dLhz(lmS)*Plm(lmS,nThetaNHS)
                  dvrdrES=dvrdrES + leg_helper%dLhdw(lmS)*Plm(lmS,nThetaNHS)
                  brES   =brES    +  leg_helper%dLhb(lmS)*Plm(lmS,nThetaNHS)
               end if
               cvrc(2*mc-1,nThetaN)  = real(cvrES  +cvrEA)
               cvrc(2*mc  ,nThetaN)  =aimag(cvrES  +cvrEA)
               cvrc(2*mc-1,nThetaS)  = real(cvrES  -cvrEA)
               cvrc(2*mc  ,nThetaS)  =aimag(cvrES  -cvrEA)
               dvrdrc(2*mc-1,nThetaN)= real(dvrdrES+dvrdrEA)
               dvrdrc(2*mc  ,nThetaN)=aimag(dvrdrES+dvrdrEA)
               dvrdrc(2*mc-1,nThetaS)= real(dvrdrES-dvrdrEA)
               dvrdrc(2*mc  ,nThetaS)=aimag(dvrdrES-dvrdrEA)
               brc(2*mc-1,nThetaN)   = real(brES   +brEA)
               brc(2*mc  ,nThetaN)   =aimag(brES   +brEA)
               brc(2*mc-1,nThetaS)   = real(brES   -brEA)
               brc(2*mc  ,nThetaS)   =aimag(brES   -brEA)
            end do
            PERFOFF_I
            PERFON_I('TFG_3')
            do mc=1,n_m_max
               dm =D_mc2m(mc)*osn1(nThetaNHS)
               lmS=lStop(mc)
               vrES   =zero
               vrEA   =zero
               dvrdtES=zero
               dvrdtEA=zero
               cbrES  =zero
               cbrEA  =zero
               !--- 8 add/mult, 29 dble words
               do lm=lStart(mc),lmS-1,2
                  vrES    =vrES    + leg_helper%dLhw(lm)*   Plm(lm,nThetaNHS)
                  dvrdtEA =dvrdtEA + leg_helper%dLhw(lm)*  dPlmt(lm,nThetaNHS)
                  cbrES   =cbrES   + leg_helper%dLhj(lm)*   Plm(lm,nThetaNHS)
                  PlmG(lm)=dPlmt(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlmt(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                  vrEA    =vrEA    + leg_helper%dLhw(lm+1)* Plm(lm+1,nThetaNHS)
                  dvrdtES =dvrdtES + leg_helper%dLhw(lm+1)*dPlmt(lm+1,nThetaNHS)
                  cbrEA   =cbrEA   + leg_helper%dLhj(lm+1)* Plm(lm+1,nThetaNHS)
                  PlmG(lm+1)=dPlmt(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                  PlmC(lm+1)=dPlmt(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  vrES    =vrES    + leg_helper%dLhw(lmS)* Plm(lmS,nThetaNHS)
                  dvrdtEA =dvrdtEA + leg_helper%dLhw(lmS)*dPlmt(lmS,nThetaNHS)
                  cbrES   =cbrES   + leg_helper%dLhj(lmS)* Plm(lmS,nThetaNHS)
                  PlmG(lmS)=dPlmt(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                  PlmC(lmS)=dPlmt(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
               end if
               vrc(2*mc-1,nThetaN)   = real(vrES+vrEA)
               vrc(2*mc  ,nThetaN)   =aimag(vrES+vrEA)
               vrc(2*mc-1,nThetaS)   = real(vrES-vrEA)
               vrc(2*mc,nThetaS)     =aimag(vrES-vrEA)
               dvrdtc(2*mc-1,nThetaN)= real(dvrdtES+dvrdtEA)
               dvrdtc(2*mc  ,nThetaN)=aimag(dvrdtES+dvrdtEA)
               dvrdtc(2*mc-1,nThetaS)= real(dvrdtES-dvrdtEA)
               dvrdtc(2*mc  ,nThetaS)=aimag(dvrdtES-dvrdtEA)
               cbrc(2*mc-1,nThetaN)  = real(cbrES  +cbrEA)
               cbrc(2*mc  ,nThetaN)  =aimag(cbrES  +cbrEA)
               cbrc(2*mc-1,nThetaS)  = real(cbrES  -cbrEA)
               cbrc(2*mc  ,nThetaS)  =aimag(cbrES  -cbrEA)
            end do
            PERFOFF_I
            PERFON_I('TFG_4')
    
            !--- Now the stuff using generalized harmonics:
            do mc=1,n_m_max
               lmS=lStop(mc)
               vhN1=zero
               vhS1=zero
               vhN2=zero
               vhS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  vhN1=vhN1+leg_helper%vhG(lm)*PlmG(lm)+leg_helper%vhG(lm+1)*PlmG(lm+1)
                  vhS1=vhS1-leg_helper%vhG(lm)*PlmC(lm)+leg_helper%vhG(lm+1)*PlmC(lm+1)
                  vhN2=vhN2+leg_helper%vhC(lm)*PlmC(lm)+leg_helper%vhC(lm+1)*PlmC(lm+1)
                  vhS2=vhS2-leg_helper%vhC(lm)*PlmG(lm)+leg_helper%vhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                  vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                  vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                  vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
               end if
               vhN1M(mc)=half*vhN1
               vhS1M(mc)=half*vhS1
               vhN2M(mc)=half*vhN2
               vhS2M(mc)=half*vhS2
            end do
            PERFOFF_I
            PERFON_I('TFG_5')
    
            do mc=1,n_m_max
               lmS=lStop(mc)
               dvhdrN1=zero
               dvhdrS1=zero
               dvhdrN2=zero
               dvhdrS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lm)  *PlmG(lm) + &
                       leg_helper%dvhdrG(lm+1)*PlmG(lm+1)
                  dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lm)  *PlmC(lm) + &
                       leg_helper%dvhdrG(lm+1)*PlmC(lm+1)
                  dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lm)  *PlmC(lm) + &
                       leg_helper%dvhdrC(lm+1)*PlmC(lm+1)
                  dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lm)  *PlmG(lm) + &
                       leg_helper%dvhdrC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lmS)*PlmG(lmS)
                  dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lmS)*PlmC(lmS)
                  dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lmS)*PlmC(lmS)
                  dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lmS)*PlmG(lmS)
               end if
               dvhdrN1M(mc)=half*dvhdrN1
               dvhdrS1M(mc)=half*dvhdrS1
               dvhdrN2M(mc)=half*dvhdrN2
               dvhdrS2M(mc)=half*dvhdrS2
            end do
            PERFOFF_I
            PERFON_I('TFG_6')
    
            do mc=1,n_m_max
               lmS=lStop(mc)
               bhN1=zero
               bhS1=zero
               bhN2=zero
               bhS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  bhN1=bhN1+leg_helper%bhG(lm)*PlmG(lm)+leg_helper%bhG(lm+1)*PlmG(lm+1)
                  bhS1=bhS1-leg_helper%bhG(lm)*PlmC(lm)+leg_helper%bhG(lm+1)*PlmC(lm+1)
                  bhN2=bhN2+leg_helper%bhC(lm)*PlmC(lm)+leg_helper%bhC(lm+1)*PlmC(lm+1)
                  bhS2=bhS2-leg_helper%bhC(lm)*PlmG(lm)+leg_helper%bhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  bhN1=bhN1+leg_helper%bhG(lmS)*PlmG(lmS)
                  bhS1=bhS1-leg_helper%bhG(lmS)*PlmC(lmS)
                  bhN2=bhN2+leg_helper%bhC(lmS)*PlmC(lmS)
                  bhS2=bhS2-leg_helper%bhC(lmS)*PlmG(lmS)
               end if
               bhN1M(mc)=half*bhN1
               bhS1M(mc)=half*bhS1
               bhN2M(mc)=half*bhN2
               bhS2M(mc)=half*bhS2
            end do
            PERFOFF_I
            PERFON_I('TFG_7')
    
            do mc=1,n_m_max
               lmS=lStop(mc)
               cbhN1=zero
               cbhS1=zero
               cbhN2=zero
               cbhS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  cbhN1=cbhN1+leg_helper%cbhG(lm)  *PlmG(lm)+ &
                              leg_helper%cbhG(lm+1)*PlmG(lm+1)
                  cbhS1=cbhS1-leg_helper%cbhG(lm)  *PlmC(lm)+ &
                              leg_helper%cbhG(lm+1)*PlmC(lm+1)
                  cbhN2=cbhN2+leg_helper%cbhC(lm)  *PlmC(lm)+ &
                              leg_helper%cbhC(lm+1)*PlmC(lm+1)
                  cbhS2=cbhS2-leg_helper%cbhC(lm)  *PlmG(lm)+ &
                              leg_helper%cbhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  cbhN1=cbhN1+leg_helper%cbhG(lmS)*PlmG(lmS)
                  cbhS1=cbhS1-leg_helper%cbhG(lmS)*PlmC(lmS)
                  cbhN2=cbhN2+leg_helper%cbhC(lmS)*PlmC(lmS)
                  cbhS2=cbhS2-leg_helper%cbhC(lmS)*PlmG(lmS)
               end if
               cbhN1M(mc)=half*cbhN1
               cbhS1M(mc)=half*cbhS1
               cbhN2M(mc)=half*cbhN2
               cbhS2M(mc)=half*cbhS2
            end do
            PERFOFF_I
            PERFON_I('TFG_8')
    
            !--- Unscramble:
            !--- 6 add/mult, 20 dble words
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
            !--- 6 add/mult, 20 dble words
            do mc=1,n_m_max
               dvtdrc(2*mc-1,nThetaN)= real(dvhdrN1M(mc)+dvhdrN2M(mc))
               dvtdrc(2*mc  ,nThetaN)=aimag(dvhdrN1M(mc)+dvhdrN2M(mc))
               dvhdrN                =dvhdrN1M(mc)-dvhdrN2M(mc)
               dvtdrc(2*mc-1,nThetaS)= real(dvhdrS1M(mc)+dvhdrS2M(mc))
               dvtdrc(2*mc  ,nThetaS)=aimag(dvhdrS1M(mc)+dvhdrS2M(mc))
               dvhdrS                =dvhdrS1M(mc)-dvhdrS2M(mc)
               dvpdrc(2*mc-1,nThetaN)=aimag(dvhdrN)
               dvpdrc(2*mc  ,nThetaN)=-real(dvhdrN)
               dvpdrc(2*mc-1,nThetaS)=aimag(dvhdrS)
               dvpdrc(2*mc  ,nThetaS)=-real(dvhdrS)
            end do
            !--- 6 add/mult, 20 dble words
            do mc=1,n_m_max
               btc(2*mc-1,nThetaN)= real(bhN1M(mc)+bhN2M(mc))
               btc(2*mc  ,nThetaN)=aimag(bhN1M(mc)+bhN2M(mc))
               bhN                =bhN1M(mc)-bhN2M(mc)
               btc(2*mc-1,nThetaS)= real(bhS1M(mc)+bhS2M(mc))
               btc(2*mc  ,nThetaS)=aimag(bhS1M(mc)+bhS2M(mc))
               bhS                =bhS1M(mc)-bhS2M(mc)
               bpc(2*mc-1,nThetaN)=aimag(bhN)
               bpc(2*mc  ,nThetaN)=-real(bhN)
               bpc(2*mc-1,nThetaS)=aimag(bhS)
               bpc(2*mc  ,nThetaS)=-real(bhS)
            end do
            !--- 6 add/mult, 20 dble words
            do mc=1,n_m_max
               cbtc(2*mc-1,nThetaN)= real(cbhN1M(mc)+cbhN2M(mc))
               cbtc(2*mc  ,nThetaN)=aimag(cbhN1M(mc)+cbhN2M(mc))
               cbhN                =cbhN1M(mc)-cbhN2M(mc)
               cbtc(2*mc-1,nThetaS)= real(cbhS1M(mc)+cbhS2M(mc))
               cbtc(2*mc  ,nThetaS)=aimag(cbhS1M(mc)+cbhS2M(mc))
               cbhS                =cbhS1M(mc)-cbhS2M(mc)
               cbpc(2*mc-1,nThetaN)=aimag(cbhN)
               cbpc(2*mc  ,nThetaN)=-real(cbhN)
               cbpc(2*mc-1,nThetaS)=aimag(cbhS)
               cbpc(2*mc  ,nThetaS)=-real(cbhS)
            end do ! Loop over order m
            PERFOFF_I
            PERFON_I('TFG_9')
    
            !--- Calculate phi derivatives:
            do mc=1,n_m_max
               dm=D_mc2m(mc)
               dvrdpc(2*mc-1,nThetaN)=-dm*vrc(2*mc  ,nThetaN)
               dvrdpc(2*mc  ,nThetaN)= dm*vrc(2*mc-1,nThetaN)
               dvrdpc(2*mc-1,nThetaS)=-dm*vrc(2*mc  ,nThetaS)
               dvrdpc(2*mc  ,nThetaS)= dm*vrc(2*mc-1,nThetaS)
            end do
            do mc=1,n_m_max
               dm=D_mc2m(mc)
               dvtdpc(2*mc-1,nThetaN)=-dm*vtc(2*mc  ,nThetaN)
               dvtdpc(2*mc  ,nThetaN)= dm*vtc(2*mc-1,nThetaN)
               dvtdpc(2*mc-1,nThetaS)=-dm*vtc(2*mc  ,nThetaS)
               dvtdpc(2*mc  ,nThetaS)= dm*vtc(2*mc-1,nThetaS)
               dvpdpc(2*mc-1,nThetaN)=-dm*vpc(2*mc  ,nThetaN)
               dvpdpc(2*mc  ,nThetaN)= dm*vpc(2*mc-1,nThetaN)
               dvpdpc(2*mc-1,nThetaS)=-dm*vpc(2*mc  ,nThetaS)
               dvpdpc(2*mc  ,nThetaS)= dm*vpc(2*mc-1,nThetaS)
            end do   ! End of loop over oder m numbered by mc
            PERFOFF_I
         end do      ! End global loop over nTheta
    
         !PERFOFF
    
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  sc(mc,nThetaN)=0.0_cp
                  if ( lViscBcCalc) then
                     dsdtc(mc,nThetaN)=0.0_cp
                     dsdpc(mc,nThetaN)=0.0_cp
                  end if
                  if ( lFluxProfCalc ) then
                     pc(mc,nThetaN)=0.0_cp
                  end if
                  vrc(mc,nThetaN)   =0.0_cp
                  vtc(mc,nThetaN)   =0.0_cp
                  vpc(mc,nThetaN)   =0.0_cp
                  cvrc(mc,nThetaN)  =0.0_cp
                  dvrdrc(mc,nThetaN)=0.0_cp
                  dvtdrc(mc,nThetaN)=0.0_cp
                  dvpdrc(mc,nThetaN)=0.0_cp
                  dvrdtc(mc,nThetaN)=0.0_cp
                  dvrdpc(mc,nThetaN)=0.0_cp
               end do
               do mc=2*n_m_max+1,nrp
                  dvtdpc(mc,nThetaN)=0.0_cp
                  dvpdpc(mc,nThetaN)=0.0_cp
                  brc(mc,nThetaN)   =0.0_cp
                  btc(mc,nThetaN)   =0.0_cp
                  bpc(mc,nThetaN)   =0.0_cp
                  cbrc(mc,nThetaN)  =0.0_cp
                  cbtc(mc,nThetaN)  =0.0_cp
                  cbpc(mc,nThetaN)  =0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
         !PERFOFF
    
      else   ! boundary ?
    
         !PERFON('TFG_bnd')
         !-- Calculation for boundary r_cmb or r_icb:
    
         do nThetaN=1,sizeThetaB,2
            nThetaS=nThetaN+1
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
    
            if ( l_heat ) then
               do mc=1,n_m_max
                  lmS=lStop(mc)
                  sES=zero    ! One equatorial symmetry
                  sEA=zero    ! The other equatorial symmetry
                  do lm=lStart(mc),lmS-1,2
                     sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                     sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
                  end do
                  if ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
                  sc(2*mc-1,nThetaN)= real(sES+sEA)
                  sc(2*mc  ,nThetaN)=aimag(sES+sEA)
                  sc(2*mc-1,nThetaS)= real(sES-sEA)
                  sc(2*mc  ,nThetaS)=aimag(sES-sEA)
               end do
            end if
    
            do mc=1,n_m_max
               dm =D_mc2m(mc)*osn1(nThetaNHS)
               lmS=lStop(mc)
    
               !------ br = r^2 B_r , bt = r sin(theta) B_theta , bp= r sin(theta) B_phi
               brES=zero
               brEA=zero
               do lm=lStart(mc),lmS-1,2
                  brES=brES + leg_helper%dLhb(lm)  *Plm(lm,nThetaNHS)
                  brEA=brEA + leg_helper%dLhb(lm+1)*Plm(lm+1,nThetaNHS)
                  PlmG(lm)=dPlmt(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlmt(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                  PlmG(lm+1)=dPlmt(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                  PlmC(lm+1)=dPlmt(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  brES=brES+leg_helper%dLhb(lm)*Plm(lm,nThetaNHS)
                  PlmG(lm)=dPlmt(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlmt(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
               end if
               brc(2*mc-1,nThetaN)=real(brES+brEA)
               brc(2*mc  ,nThetaN)=aimag(brES+brEA)
               brc(2*mc-1,nThetaS)=real(brES-brEA)
               brc(2*mc  ,nThetaS)=aimag(brES-brEA)
    
               bhN1=zero
               bhS1=zero
               bhN2=zero
               bhS2=zero
               do lm=lStart(mc),lmS-1,2
                  bhN1=bhN1+leg_helper%bhG(lm)*PlmG(lm)+leg_helper%bhG(lm+1)*PlmG(lm+1)
                  bhS1=bhS1-leg_helper%bhG(lm)*PlmC(lm)+leg_helper%bhG(lm+1)*PlmC(lm+1)
                  bhN2=bhN2+leg_helper%bhC(lm)*PlmC(lm)+leg_helper%bhC(lm+1)*PlmC(lm+1)
                  bhS2=bhS2-leg_helper%bhC(lm)*PlmG(lm)+leg_helper%bhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  bhN1=bhN1+leg_helper%bhG(lmS)*PlmG(lmS)
                  bhS1=bhS1-leg_helper%bhG(lmS)*PlmC(lmS)
                  bhN2=bhN2+leg_helper%bhC(lmS)*PlmC(lmS)
                  bhS2=bhS2-leg_helper%bhC(lmS)*PlmG(lmS)
               end if
               btc(2*mc-1,nThetaN)=real(half*bhN1+half*bhN2)
               btc(2*mc  ,nThetaN)=aimag(half*bhN1+half*bhN2)
               btc(2*mc-1,nThetaS)=real(half*bhS1+half*bhS2)
               btc(2*mc  ,nThetaS)=aimag(half*bhS1+half*bhS2)
               bhN                =half*bhN1-half*bhN2
               bhS                =half*bhS1-half*bhS2
               bpc(2*mc-1,nThetaN)=aimag(bhN)
               bpc(2*mc  ,nThetaN)=-real(bhN)
               bpc(2*mc-1,nThetaS)=aimag(bhS)
               bpc(2*mc  ,nThetaS)=-real(bhS)
    
            end do
    
            if ( nBc == 1 ) then
    
               !--- Horizontal velocity components for nBc=1
               do mc=1,n_m_max
                  lmS=lStop(mc)
                  vhN1=zero
                  vhS1=zero
                  vhN2=zero
                  vhS2=zero
                  do lm=lStart(mc),lmS-1,2
                     vhN1=vhN1+leg_helper%vhG(lm)  *PlmG(lm)+ &
                               leg_helper%vhG(lm+1)*PlmG(lm+1)
                     vhS1=vhS1-leg_helper%vhG(lm)  *PlmC(lm)+ &
                               leg_helper%vhG(lm+1)*PlmC(lm+1)
                     vhN2=vhN2+leg_helper%vhC(lm)  *PlmC(lm)+ &
                               leg_helper%vhC(lm+1)*PlmC(lm+1)
                     vhS2=vhS2-leg_helper%vhC(lm)  *PlmG(lm)+ &
                               leg_helper%vhC(lm+1)*PlmG(lm+1)
                  end do
                  if ( lmOdd(mc) ) then
                     vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                     vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                     vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                     vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
                  end if
                  vtc(2*mc-1,nThetaN)=real(half*vhN1+half*vhN2)
                  vtc(2*mc  ,nThetaN)=aimag(half*vhN1+half*vhN2)
                  vtc(2*mc-1,nThetaS)=real(half*vhS1+half*vhS2)
                  vtc(2*mc  ,nThetaS)=aimag(half*vhS1+half*vhS2)
                  vhN            =half*vhN1-half*vhN2
                  vhS            =half*vhS1-half*vhS2
                  vpc(2*mc-1,nThetaN)=aimag(vhN)
                  vpc(2*mc  ,nThetaN)=-real(vhN)
                  vpc(2*mc-1,nThetaS)=aimag(vhS)
                  vpc(2*mc  ,nThetaS)=-real(vhS)
               end do ! Loop over m
    
            end if   ! nBc == 1 ? vrc and nBc=2 cared for later !
    
         end do    ! End loop over nThetaN
    
         !-- Zero out terms with index mc > n_m_max :
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  sc(mc,nThetaN) =0.0_cp
                  brc(mc,nThetaN)=0.0_cp
                  btc(mc,nThetaN)=0.0_cp
                  bpc(mc,nThetaN)=0.0_cp
               end do
            end do
            if ( nBc == 1 ) then
               do nThetaN=1,sizeThetaB
                  do mc=2*n_m_max+1,nrp
                     vtc(mc,nThetaN)=0.0_cp
                     vpc(mc,nThetaN)=0.0_cp
                  end do
               end do
            end if
         end if
         !PERFOFF
      end if  ! boundary ? nBc?
    
    
      if ( l_HT .or. lViscBcCalc ) then    ! For movie output !
         nThetaNHS=(nThetaStart-1)/2
    
         !-- Caculate radial derivate of S for heatflux:
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for one HS
            nThetaS  =nThetaN+1  ! same theta but at other HS
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
            do mc=1,n_m_max
               lmS=lStop(mc)
               drsES=zero
               drsEA=zero
               do lm=lStart(mc),lmS-1,2
                  drsES=drsES+leg_helper%dsR(lm)*Plm(lm,nThetaNHS)
                  drsEA=drsEA+leg_helper%dsR(lm+1)*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) drsES=drsES+leg_helper%dsR(lmS)*Plm(lmS,nThetaNHS)
               drSc(2*mc-1,nThetaN)= real(drsES+drsEA)
               drSc(2*mc  ,nThetaN)=aimag(drsES+drsEA)
               drSc(2*mc-1,nThetaS)= real(drsES-drsEA)
               drSc(2*mc  ,nThetaS)=aimag(drsES-drsEA)
            end do
         end do
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  drSc(mc,nThetaN)=0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
    
      end if

   end subroutine legTFG
!------------------------------------------------------------------------------
   subroutine legTFGnomag(nBc,lDeriv,lViscBcCalc,lFluxProfCalc, & 
       &                 nThetaStart,vrc,vtc,vpc,dvrdrc,dvtdrc, &
       &                 dvpdrc,cvrc,dvrdtc,dvrdpc,dvtdpc,      & 
       &                 dvpdpc,sc,drSc,dsdtc,dsdpc,pc,         &
       &                 leg_helper)
      !
      ! Same as legTFG for non-magnetic cases
      !

      !-- Input:
      integer, intent(in) :: nBc
      logical, intent(in) :: lDeriv,lViscBcCalc,lFluxProfCalc
      integer, intent(in) :: nThetaStart
    
      !----- Stuff precomputed in legPrep:
      type(leg_helper_t) :: leg_helper
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: vrc(nrp,nfs), vtc(nrp,nfs), vpc(nrp,nfs)
      real(cp), intent(out) :: dvrdrc(nrp,nfs), dvtdrc(nrp,nfs), dvpdrc(nrp,nfs)
      real(cp), intent(out) :: dvrdtc(nrp,nfs), dvrdpc(nrp, nfs)
      real(cp), intent(out) :: dvtdpc(nrp,nfs), dvpdpc(nrp, nfs)
      real(cp), intent(out) :: cvrc(nrp,nfs), pc(nrp,nfs)
      real(cp), intent(out) :: sc(nrp,nfs), drSc(nrp, nfs)
      real(cp), intent(out) :: dsdtc(nrp,nfs), dsdpc(nrp,nfs)
    
      !------ Legendre Polynomials
      real(cp) :: PlmG(lm_max)
      real(cp) :: PlmC(lm_max)
    
      !-- Local:
      complex(cp) :: vrES,vrEA,dvrdrES,dvrdrEA,dvrdtES,dvrdtEA,cvrES,cvrEA
      complex(cp) :: sES,sEA,drsES,drsEA,pES,pEA
      complex(cp) :: dsdtES,dsdtEA
    
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS
      real(cp) :: dm
    
      complex(cp) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
      complex(cp) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
      complex(cp) :: dvhdrN1M(n_m_max),dvhdrN2M(n_m_max),dvhdrN
      complex(cp) :: dvhdrS1M(n_m_max),dvhdrS2M(n_m_max),dvhdrS
      complex(cp) :: dvhdrN1,dvhdrN2,dvhdrS1,dvhdrS2
    
    
      nThetaNHS=(nThetaStart-1)/2
    
      if ( nBc == 0 .or. lDeriv ) then ! not a boundary or derivs required
    
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
    
            if ( l_heat ) then
               do mc=1,n_m_max
                  lmS=lStop(mc)
                  sES=zero  ! One equatorial symmetry
                  sEA=zero  ! The other equatorial symmetry
                  do lm=lStart(mc),lmS-1,2
                     sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                     sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
                  end do
                  if ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
                  sc(2*mc-1,nThetaN)= real(sES+sEA)
                  sc(2*mc  ,nThetaN)=aimag(sES+sEA)
                  sc(2*mc-1,nThetaS)= real(sES-sEA)
                  sc(2*mc  ,nThetaS)=aimag(sES-sEA)
               end do
    
    
               if ( lFluxProfCalc ) then
                  do mc=1,n_m_max
                     lmS=lStop(mc)
                     pES=zero ! One equatorial symmetry
                     pEA=zero ! The other equatorial symmetry
                     do lm=lStart(mc),lmS-1,2
                        pES=pES+leg_helper%preR(lm)  *Plm(lm,nThetaNHS)
                        pEA=pEA+leg_helper%preR(lm+1)*Plm(lm+1,nThetaNHS)
                     end do
                     if ( lmOdd(mc) ) pES=pES+leg_helper%preR(lmS)*Plm(lmS,nThetaNHS)
                     pc(2*mc-1,nThetaN)= real(pES+pEA)
                     pc(2*mc  ,nThetaN)=aimag(pES+pEA)
                     pc(2*mc-1,nThetaS)= real(pES-pEA)
                     pc(2*mc  ,nThetaS)=aimag(pES-pEA)
                  end do
               end if
    
               if ( lViscBcCalc ) then
                  do mc=1,n_m_max
                     dm =D_mc2m(mc)*osn1(nThetaNHS)
                     lmS=lStop(mc)
                     dsdtES=zero
                     dsdtEA=zero
                     do lm=lStart(mc),lmS-1,2
                        dsdtEA =dsdtEA + leg_helper%sR(lm)*  dPlm(lm,nThetaNHS)
                        dsdtES =dsdtES + leg_helper%sR(lm+1)*dPlmt(lm+1,nThetaNHS)
                     end do
                     if ( lmOdd(mc) ) then
                        dsdtEA =dsdtEA + leg_helper%sR(lmS)*dPlmt(lmS,nThetaNHS)
                     end if
                     dsdtc(2*mc-1,nThetaN)= real(dsdtES+dsdtEA)
                     dsdtc(2*mc  ,nThetaN)=aimag(dsdtES+dsdtEA)
                     dsdtc(2*mc-1,nThetaS)= real(dsdtES-dsdtEA)
                     dsdtc(2*mc  ,nThetaS)=aimag(dsdtES-dsdtEA)
                  end do
    
                  do mc=1,n_m_max
                     dm=D_mc2m(mc)
                     dsdpc(2*mc-1,nThetaN)=-dm*sc(2*mc  ,nThetaN)
                     dsdpc(2*mc  ,nThetaN)= dm*sc(2*mc-1,nThetaN)
                     dsdpc(2*mc-1,nThetaS)=-dm*sc(2*mc  ,nThetaS)
                     dsdpc(2*mc  ,nThetaS)= dm*sc(2*mc-1,nThetaS)
                  end do
    
               end if ! thermal dissipation layer
    
            end if
    
            !--- Loop over all oders m: (numbered by mc)
            do mc=1,n_m_max
               lmS=lStop(mc)
               cvrES  =zero
               cvrEA  =zero
               dvrdrES=zero
               dvrdrEA=zero
               !--- 4 add/mult
               do lm=lStart(mc),lmS-1,2
                  cvrES  =cvrES   +  leg_helper%dLhz(lm)  *Plm(lm,nThetaNHS)
                  dvrdrES=dvrdrES + leg_helper%dLhdw(lm)  *Plm(lm,nThetaNHS)
                  cvrEA  =cvrEA   +  leg_helper%dLhz(lm+1)*Plm(lm+1,nThetaNHS)
                  dvrdrEA=dvrdrEA + leg_helper%dLhdw(lm+1)*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  cvrES  =cvrES   +  leg_helper%dLhz(lmS)*Plm(lmS,nThetaNHS)
                  dvrdrES=dvrdrES + leg_helper%dLhdw(lmS)*Plm(lmS,nThetaNHS)
               end if
               cvrc(2*mc-1,nThetaN)  = real(cvrES  +cvrEA)
               cvrc(2*mc  ,nThetaN)  =aimag(cvrES  +cvrEA)
               cvrc(2*mc-1,nThetaS)  = real(cvrES  -cvrEA)
               cvrc(2*mc  ,nThetaS)  =aimag(cvrES  -cvrEA)
               dvrdrc(2*mc-1,nThetaN)= real(dvrdrES+dvrdrEA)
               dvrdrc(2*mc  ,nThetaN)=aimag(dvrdrES+dvrdrEA)
               dvrdrc(2*mc-1,nThetaS)= real(dvrdrES-dvrdrEA)
               dvrdrc(2*mc  ,nThetaS)=aimag(dvrdrES-dvrdrEA)
            end do
    
            do mc=1,n_m_max
               dm =D_mc2m(mc)*osn1(nThetaNHS)
               lmS=lStop(mc)
               vrES   =zero
               vrEA   =zero
               dvrdtES=zero
               dvrdtEA=zero
               !--- 8 add/mult, 29 dble words
               do lm=lStart(mc),lmS-1,2
                  vrES    =vrES    + leg_helper%dLhw(lm)*   Plm(lm,nThetaNHS)
                  dvrdtEA =dvrdtEA + leg_helper%dLhw(lm)*  dPlmt(lm,nThetaNHS)
                  PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                  vrEA    =vrEA    + leg_helper%dLhw(lm+1)* Plm(lm+1,nThetaNHS)
                  dvrdtES =dvrdtES + leg_helper%dLhw(lm+1)*dPlmt(lm+1,nThetaNHS)
                  PlmG(lm+1)=dPlmt(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                  PlmC(lm+1)=dPlmt(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  vrES    =vrES    + leg_helper%dLhw(lmS)* Plm(lmS,nThetaNHS)
                  dvrdtEA =dvrdtEA + leg_helper%dLhw(lmS)*dPlmt(lmS,nThetaNHS)
                  PlmG(lmS)=dPlmt(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                  PlmC(lmS)=dPlmt(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
               end if
               vrc(2*mc-1,nThetaN)   = real(vrES+vrEA)
               vrc(2*mc  ,nThetaN)   =aimag(vrES+vrEA)
               vrc(2*mc-1,nThetaS)   = real(vrES-vrEA)
               vrc(2*mc,nThetaS)     =aimag(vrES-vrEA)
               dvrdtc(2*mc-1,nThetaN)= real(dvrdtES+dvrdtEA)
               dvrdtc(2*mc  ,nThetaN)=aimag(dvrdtES+dvrdtEA)
               dvrdtc(2*mc-1,nThetaS)= real(dvrdtES-dvrdtEA)
               dvrdtc(2*mc  ,nThetaS)=aimag(dvrdtES-dvrdtEA)
            end do
    
            !--- Now the stuff using generalized harmonics:
            do mc=1,n_m_max
               lmS=lStop(mc)
               vhN1=zero
               vhS1=zero
               vhN2=zero
               vhS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  vhN1=vhN1+leg_helper%vhG(lm)*PlmG(lm)+leg_helper%vhG(lm+1)*PlmG(lm+1)
                  vhS1=vhS1-leg_helper%vhG(lm)*PlmC(lm)+leg_helper%vhG(lm+1)*PlmC(lm+1)
                  vhN2=vhN2+leg_helper%vhC(lm)*PlmC(lm)+leg_helper%vhC(lm+1)*PlmC(lm+1)
                  vhS2=vhS2-leg_helper%vhC(lm)*PlmG(lm)+leg_helper%vhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                  vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                  vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                  vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
               end if
               vhN1M(mc)=half*vhN1
               vhS1M(mc)=half*vhS1
               vhN2M(mc)=half*vhN2
               vhS2M(mc)=half*vhS2
            end do
    
            do mc=1,n_m_max
               lmS=lStop(mc)
               dvhdrN1=zero
               dvhdrS1=zero
               dvhdrN2=zero
               dvhdrS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lm)  *PlmG(lm) + &
                       leg_helper%dvhdrG(lm+1)*PlmG(lm+1)
                  dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lm)  *PlmC(lm) + &
                       leg_helper%dvhdrG(lm+1)*PlmC(lm+1)
                  dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lm)  *PlmC(lm) + &
                       leg_helper%dvhdrC(lm+1)*PlmC(lm+1)
                  dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lm)  *PlmG(lm) + &
                       leg_helper%dvhdrC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  dvhdrN1=dvhdrN1+leg_helper%dvhdrG(lmS)*PlmG(lmS)
                  dvhdrS1=dvhdrS1-leg_helper%dvhdrG(lmS)*PlmC(lmS)
                  dvhdrN2=dvhdrN2+leg_helper%dvhdrC(lmS)*PlmC(lmS)
                  dvhdrS2=dvhdrS2-leg_helper%dvhdrC(lmS)*PlmG(lmS)
               end if
               dvhdrN1M(mc)=half*dvhdrN1
               dvhdrS1M(mc)=half*dvhdrS1
               dvhdrN2M(mc)=half*dvhdrN2
               dvhdrS2M(mc)=half*dvhdrS2
            end do
    
            !--- Unscramble:
            !--- 6 add/mult, 20 dble words
            do mc=1,n_m_max
               vtc(2*mc-1,nThetaN)= real(vhN1M(mc)+vhN2M(mc))
               vtc(2*mc  ,nThetaN)=aimag(vhN1M(mc)+vhN2M(mc))
               vhN                =vhN1M(mc)-vhN2M(mc)
               vtc(2*mc-1,nThetaS)= real(vhS1M(mc)+vhS2M(mc))
               vtc(2*mc,nThetaS)  =aimag(vhS1M(mc)+vhS2M(mc))
               vhS                =vhS1M(mc)-vhS2M(mc)
               vpc(2*mc-1,nThetaN)=aimag(vhN)
               vpc(2*mc  ,nThetaN)=-real(vhN)
               vpc(2*mc-1,nThetaS)=aimag(vhS)
               vpc(2*mc  ,nThetaS)=-real(vhS)
            end do
            !--- 6 add/mult, 20 dble words
            do mc=1,n_m_max
               dvtdrc(2*mc-1,nThetaN)= real(dvhdrN1M(mc)+dvhdrN2M(mc))
               dvtdrc(2*mc  ,nThetaN)=aimag(dvhdrN1M(mc)+dvhdrN2M(mc))
               dvhdrN                =dvhdrN1M(mc)-dvhdrN2M(mc)
               dvtdrc(2*mc-1,nThetaS)= real(dvhdrS1M(mc)+dvhdrS2M(mc))
               dvtdrc(2*mc  ,nThetaS)=aimag(dvhdrS1M(mc)+dvhdrS2M(mc))
               dvhdrS                =dvhdrS1M(mc)-dvhdrS2M(mc)
               dvpdrc(2*mc-1,nThetaN)=aimag(dvhdrN)
               dvpdrc(2*mc  ,nThetaN)=-real(dvhdrN)
               dvpdrc(2*mc-1,nThetaS)=aimag(dvhdrS)
               dvpdrc(2*mc  ,nThetaS)=-real(dvhdrS)
            end do
    
            !--- Calculate phi derivatives:
            do mc=1,n_m_max
               dm=D_mc2m(mc)
               dvrdpc(2*mc-1,nThetaN)=-dm*vrc(2*mc  ,nThetaN)
               dvrdpc(2*mc  ,nThetaN)= dm*vrc(2*mc-1,nThetaN)
               dvrdpc(2*mc-1,nThetaS)=-dm*vrc(2*mc  ,nThetaS)
               dvrdpc(2*mc  ,nThetaS)= dm*vrc(2*mc-1,nThetaS)
            end do
            do mc=1,n_m_max
               dm=D_mc2m(mc)
               dvtdpc(2*mc-1,nThetaN)=-dm*vtc(2*mc  ,nThetaN)
               dvtdpc(2*mc  ,nThetaN)= dm*vtc(2*mc-1,nThetaN)
               dvtdpc(2*mc-1,nThetaS)=-dm*vtc(2*mc  ,nThetaS)
               dvtdpc(2*mc  ,nThetaS)= dm*vtc(2*mc-1,nThetaS)
               dvpdpc(2*mc-1,nThetaN)=-dm*vpc(2*mc  ,nThetaN)
               dvpdpc(2*mc  ,nThetaN)= dm*vpc(2*mc-1,nThetaN)
               dvpdpc(2*mc-1,nThetaS)=-dm*vpc(2*mc  ,nThetaS)
               dvpdpc(2*mc  ,nThetaS)= dm*vpc(2*mc-1,nThetaS)
            end do   ! End of loop over oder m numbered by mc
    
         end do      ! End global loop over nTheta
    
    
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  sc(mc,nThetaN)=0.0_cp
                  if ( lFluxProfCalc ) then
                     pc(mc,nThetaN)=0.0_cp
                  end if
                  if ( lViscBcCalc) then
                     dsdtc(mc,nThetaN)=0.0_cp
                     dsdpc(mc,nThetaN)=0.0_cp
                  end if
                  vrc(mc,nThetaN)   =0.0_cp
                  vtc(mc,nThetaN)   =0.0_cp
                  vpc(mc,nThetaN)   =0.0_cp
                  cvrc(mc,nThetaN)  =0.0_cp
                  dvrdrc(mc,nThetaN)=0.0_cp
                  dvtdrc(mc,nThetaN)=0.0_cp
                  dvpdrc(mc,nThetaN)=0.0_cp
                  dvrdtc(mc,nThetaN)=0.0_cp
                  dvrdpc(mc,nThetaN)=0.0_cp
                  dvtdpc(mc,nThetaN)=0.0_cp
                  dvpdpc(mc,nThetaN)=0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
    
    
      else   ! boundary ?
    
    
         !-- Calculation for boundary r_cmb or r_icb:
    
         do nThetaN=1,sizeThetaB,2
            nThetaS=nThetaN+1
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
    
            if ( l_heat ) then
               do mc=1,n_m_max
                  lmS=lStop(mc)
                  sES=zero    ! One equatorial symmetry
                  sEA=zero    ! The other equatorial symmetry
                  do lm=lStart(mc),lmS-1,2
                     sES=sES+leg_helper%sR(lm)  *Plm(lm,nThetaNHS)
                     sEA=sEA+leg_helper%sR(lm+1)*Plm(lm+1,nThetaNHS)
                  end do
                  if ( lmOdd(mc) ) sES=sES+leg_helper%sR(lmS)*Plm(lmS,nThetaNHS)
                  sc(2*mc-1,nThetaN)= real(sES+sEA)
                  sc(2*mc  ,nThetaN)=aimag(sES+sEA)
                  sc(2*mc-1,nThetaS)= real(sES-sEA)
                  sc(2*mc  ,nThetaS)=aimag(sES-sEA)
               end do
            end if
    
            if ( nBc == 1 ) then
    
               !--- Horizontal velocity components for nBc=1
               do mc=1,n_m_max
                  dm =D_mc2m(mc)*osn1(nThetaNHS)
                  lmS=lStop(mc)
                  vhN1=zero
                  vhS1=zero
                  vhN2=zero
                  vhS2=zero
                  do lm=lStart(mc),lmS-1,2
                     PlmG(lm)=dPlmt(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                     PlmC(lm)=dPlmt(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                     PlmG(lm+1)=dPlmt(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                     PlmC(lm+1)=dPlmt(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
                     vhN1=vhN1+leg_helper%vhG(lm)*  PlmG(lm)+    &
                           &   leg_helper%vhG(lm+1)*PlmG(lm+1)
                     vhS1=vhS1-leg_helper%vhG(lm)*  PlmC(lm)+    &
                           &   leg_helper%vhG(lm+1)*PlmC(lm+1)
                     vhN2=vhN2+leg_helper%vhC(lm)  *PlmC(lm)+    &
                           &   leg_helper%vhC(lm+1)*PlmC(lm+1)
                     vhS2=vhS2-leg_helper%vhC(lm)  *PlmG(lm)+    &
                           &   leg_helper%vhC(lm+1)*PlmG(lm+1)
                  end do
                  if ( lmOdd(mc) ) then
                     PlmG(lmS)=dPlmt(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                     PlmC(lmS)=dPlmt(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
                     vhN1=vhN1+leg_helper%vhG(lmS)*PlmG(lmS)
                     vhS1=vhS1-leg_helper%vhG(lmS)*PlmC(lmS)
                     vhN2=vhN2+leg_helper%vhC(lmS)*PlmC(lmS)
                     vhS2=vhS2-leg_helper%vhC(lmS)*PlmG(lmS)
                  end if
                  vtc(2*mc-1,nThetaN)=real(half*vhN1+half*vhN2)
                  vtc(2*mc,nThetaN)  =aimag(half*vhN1+half*vhN2)
                  vtc(2*mc-1,nThetaS)=real(half*vhS1+half*vhS2)
                  vtc(2*mc,nThetaS)  =aimag(half*vhS1+half*vhS2)
                  vhN            =half*vhN1-half*vhN2
                  vhS            =half*vhS1-half*vhS2
                  vpc(2*mc-1,nThetaN)=aimag(vhN)
                  vpc(2*mc  ,nThetaN)=-real(vhN)
                  vpc(2*mc-1,nThetaS)=aimag(vhS)
                  vpc(2*mc  ,nThetaS)=-real(vhS)
               end do ! Loop over m
    
            end if   ! nBc == 1 ? vrc and nBc=2 cared for later !
    
         end do    ! End loop over nThetaN
    
         !-- Zero out terms with index mc > n_m_max :
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  sc(mc,nThetaN)=0.0_cp
               end do
            end do
            if ( nBc == 1 ) then
               do nThetaN=1,sizeThetaB
                  do mc=2*n_m_max+1,nrp
                     vtc(mc,nThetaN)=0.0_cp
                     vpc(mc,nThetaN)=0.0_cp
                  end do
               end do
            end if
         end if
    
      end if  ! boundary ? nBc?
    
    
      if ( l_HT .or. lViscBcCalc ) then    ! For movie output !
         nThetaNHS=(nThetaStart-1)/2
    
         !-- Caculate radial derivate of S for heatflux:
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for one HS
            nThetaS  =nThetaN+1  ! same theta but at other HS
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
            do mc=1,n_m_max
               lmS=lStop(mc)
               drsES=zero
               drsEA=zero
               do lm=lStart(mc),lmS-1,2
                  drsES=drsES+leg_helper%dsR(lm)*Plm(lm,nThetaNHS)
                  drsEA=drsEA+leg_helper%dsR(lm+1)*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) drsES=drsES+leg_helper%dsR(lmS)*Plm(lmS,nThetaNHS)
               drSc(2*mc-1,nThetaN)= real(drsES+drsEA)
               drSc(2*mc  ,nThetaN)=aimag(drsES+drsEA)
               drSc(2*mc-1,nThetaS)= real(drsES-drsEA)
               drSc(2*mc  ,nThetaS)=aimag(drsES-drsEA)
            end do
         end do
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  drSc(mc,nThetaN)=0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
    
      end if

   end subroutine legTFGnomag
!------------------------------------------------------------------------------
   subroutine legTF(dLhw,vhG,vhC,dLhz,cvhG,cvhC,l_max,minc, &
              &     nThetaStart,sizeThetaB,Plm,dPlm,lHor,   &
              &     lDeriv,vrc,vtc,vpc,cvrc,cvtc,cvpc)
      !
      !    'Legendre transform' from (nR,l,m) to (nR,nTheta,m) [spectral to grid]
      !    where nTheta numbers the colatitudes and l and m are degree and
      !    order of the spherical harmonic representation.
      !
      !    Calculates all three spherical components vrc,vtc,vpc of a field as
      !    well as its curl (cvrc,cvtc,cvpc) that is given a spherical harmonis poloidal
      !    toroidal decomposition. s_legPrep.f has to be called first and
      !    provides the input fields dLhW, .....
      !    The symmetry properties of the P_lm with respect to the equator
      !    are used. The equatorially anti-symmetric (EA) contribution is
      !    added to (subracted from ) the equatorially symmetric (ES) contribution
      !    in northern (southern) hemisphere.
      !
      !    Output is given for all sizeThetaB colatitudes in a colatitude block
      !    that starts with colatitude nThetaStart. At output, each component
      !    in the northern hemisphere is followed by the component in the
      !    southern hemisphere.
      !    The Plms and dPlms=sin(theta) d Plm / d theta are only given
      !    for the colatitudes in the northern hemisphere.
      !
      !      * dLhw,....,cvhC : (input) arrays provided by s_legPrep.f
      !      * l_max          : (input) maximum spherical harmonic degree
      !      * minc           : (input) azimuthal symmetry
      !      * nThetaStart    : (input) transformation is done for the range of
      !        points nThetaStart <= nTheta <= nThetaStart-1+sizeThetaB
      !      * sizeThetaB     : (input) size theta block
      !      * Plm            : (input) associated Legendre polynomials
      !      * dPlm           : (input) sin(theta) d Plm / d theta
      !      * lHor=.true.    : (input) calculate horizontal componenst
      !      * lDeriv=.true.  : (input) calculate curl of field
      !      * vrc, ....,cvpc : (output) components in (nTheta,m)-space
      !
    
    
      !-- Input:
      !----- Stuff precomputed in legPrep:
      complex(cp), intent(in) :: dLhw(*),dLhz(*)
      complex(cp), intent(in) :: vhG(*),vhC(*)
      complex(cp), intent(in) :: cvhG(*),cvhC(*)
    
      !----- Defining theta block
      integer,     intent(in) :: nThetaStart,sizeThetaB
    
      !------ Legendre Polynomials 
      integer,     intent(in) :: l_max,minc
      real(cp),    intent(in) :: Plm(lm_max,*)
      real(cp),    intent(in) :: dPlm(lm_max,*)
    
      !----- What should I calculate?
      logical,     intent(in) :: lHor
      logical,     intent(in) :: lDeriv
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: vrc(nrp,*)
      real(cp), intent(out) :: vtc(nrp,*)
      real(cp), intent(out) :: vpc(nrp,*)
      real(cp), intent(out) :: cvrc(nrp,*)
      real(cp), intent(out) :: cvtc(nrp,*)
      real(cp), intent(out) :: cvpc(nrp,*)
    
      !-- Local variables:
      complex(cp) :: vrES,vrEA,cvrES,cvrEA
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: m,l,mc,lm
      real(cp) :: PlmG(lm_max)
      real(cp) :: PlmC(lm_max)
    
      complex(cp) :: vhN1M(n_m_max),vhN2M(n_m_max),vhN1,vhN2,vhN
      complex(cp) :: vhS1M(n_m_max),vhS2M(n_m_max),vhS1,vhS2,vhS
      complex(cp) :: cvhN1M(n_m_max),cvhN2M(n_m_max)
      complex(cp) :: cvhS1M(n_m_max),cvhS2M(n_m_max)
    
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
            vrES =zero
            vrEA =zero
            cvrES=zero
            cvrEA=zero
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
               vhN1=zero
               vhS1=zero
               vhN2=zero
               vhS2=zero
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
               vhN1M(mc)=half*vhN1
               vhS1M(mc)=half*vhS1
               vhN2M(mc)=half*vhN2
               vhS2M(mc)=half*vhS2
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
               vhN1 =zero
               vhS1 =zero
               vhN2 =zero
               vhS2 =zero
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
               cvhN1M(mc)=half*vhN1
               cvhS1M(mc)=half*vhS1
               cvhN2M(mc)=half*vhN2
               cvhS2M(mc)=half*vhS2
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
               vrc(mc,nThetaN) =0.0_cp
               if ( lHor ) then
                  vtc(mc,nThetaN) =0.0_cp
                  vpc(mc,nThetaN) =0.0_cp
               end if
               if ( lDeriv ) then
                  cvrc(mc,nThetaN)=0.0_cp
                  cvtc(mc,nThetaN)=0.0_cp
                  cvpc(mc,nThetaN)=0.0_cp
               end if
            end do
         end do  ! loop over nThetaN (theta)
      end if

   end subroutine legTF
!------------------------------------------------------------------------------
   subroutine lmAS2pt(alm,aij,nThetaStart,nThetaBlockSize)
      !
      ! Spherical harmonic transform from alm(l) to aij(theta)            
      ! Done within the range [nThetaStart,n_thetaStart+nThetaBlockSize-1]
      ! only considering axisymmetric contributions.                      
      ! alm contains only m=0 coefficients                                
      !

      !-- Input variables:
      integer,  intent(in) :: nThetaStart     ! first theta to be treated
      integer,  intent(in) :: nThetaBlockSize !
      real(cp), intent(in) :: alm(*)      ! field in (l,m)-space

      !-- Output variables:
      real(cp), intent(out) :: aij(*)  ! field in (theta,phi)-space

      !-- Local variables:
      integer :: nTheta        ! last theta treated
      integer :: nThetaNHS     ! counter for theta in one HS
      integer :: nThetaN       ! counter for theta in NHS
      integer :: nThetaS       ! counter for theta in SHS
      integer :: nThetaBlock   ! counter for theta in block
      integer :: l,lm          ! degree/order
      real(cp) ::  sign

      nTheta=nThetaStart-1  ! last theta
       
      !-- Zero output field:
      do nThetaBlock=1,nThetaBlockSize
         aij(nThetaBlock)=0.0_cp
      end do

      !-- Transform l 2 theta:
      nThetaNHS=nTheta/2  ! last theta in one HS

      do nThetaN=1,nThetaBlockSize,2
         nThetaS  =nThetaN+1
         nThetaNHS=nThetaNHS+1

         sign=-one
         do l=0,l_max
            lm=lm2(l,0)
            sign=-sign
            ! Northern hemisphere
            aij(nThetaN)=aij(nThetaN) +      alm(l+1)*Plm(lm,nThetaNHS)
            ! Southern hemisphere
            aij(nThetaS)=aij(nThetaS) + sign*alm(l+1)*Plm(lm,nThetaNHS)
         end do

      end do

   end subroutine lmAS2pt
!------------------------------------------------------------------------------
#if 0
#if 1
   subroutine compute_l_sum_NS_3(nThetaN,nThetaS,infield,Plm_slice,outfield)

      integer, parameter :: nFields=3

      ! Input variables:
      integer,     intent(in) :: nThetaN,nThetaS
      complex(cp), intent(in) :: infield(lm_max,nFields)
      real(cp),    intent(in) :: Plm_slice(lm_max)
      real(cp),    intent(inout) :: outfield(nrp,nfs,nFields)

      ! Local variables:
      integer :: mc,lmS,lm,nf
      complex(cp) :: sum_sym(nFields),sum_asym(nFields)

      do mc=1,n_m_max
         lmS=lStop(mc)
         sum_sym  = zero    ! One equatorial symmetry
         sum_asym = zero    ! The other equatorial symmetry
         do lm=lStart(mc),lmS-1,2
            do nf=1,nFields
               sum_sym(nf) = sum_sym(nf) + infield(lm,nf) * Plm_slice(lm)
               sum_asym(nf)= sum_asym(nf)+ infield(lm+1,nf)*Plm_slice(lm+1)
            end do
         end do
         if ( lmOdd(mc) ) then
            do nf=1,nFields
               sum_sym(nf) = sum_sym(nf) + infield(lmS,nf) * Plm_slice(lmS)
            end do
         end if
         do nf=1,nFields
            outfield(2*mc-1,nThetaN,nf) = real(sum_sym(nf)+sum_asym(nf))
            outfield(2*mc  ,nThetaN,nf) =aimag(sum_sym(nf)+sum_asym(nf))
            outfield(2*mc-1,nThetaS,nf) = real(sum_sym(nf)-sum_asym(nf))
            outfield(2*mc  ,nThetaS,nf) =aimag(sum_sym(nf)-sum_asym(nf))
         end do
      end do

   end subroutine compute_l_sum_NS_3
!------------------------------------------------------------------------------
   subroutine compute_l_sum_NS_1(nThetaN,nThetaS,infield,Plm_slice,outfield)

      !-- Input variables:
      integer,     intent(in) :: nThetaN,nThetaS
      complex(cp), intent(in) :: infield(lm_max)
      real(cp),    intent(in) :: Plm_slice(lm_max)
      real(cp),    intent(inout) :: outfield(nrp,nfs)

      !-- Local variables:
      integer :: mc,lmS,lm,nf
      complex(cp) :: sum_sym,sum_asym

      do mc=1,n_m_max
         lmS=lStop(mc)
         sum_sym  = zero    ! One equatorial symmetry
         sum_asym = zero    ! The other equatorial symmetry
         do lm=lStart(mc),lmS-1,2
            sum_sym = sum_sym + infield(lm) * Plm_slice(lm)
            sum_asym= sum_asym+ infield(lm+1)*Plm_slice(lm+1)
         end do
         if ( lmOdd(mc) ) then
            sum_sym = sum_sym + infield(lmS) * Plm_slice(lmS)
         end if
         outfield(2*mc-1,nThetaN) = real(sum_sym+sum_asym)
         outfield(2*mc  ,nThetaN) =aimag(sum_sym+sum_asym)
         outfield(2*mc-1,nThetaS) = real(sum_sym-sum_asym)
         outfield(2*mc  ,nThetaS) =aimag(sum_sym-sum_asym)
      end do

   end subroutine compute_l_sum_NS_1
!------------------------------------------------------------------------------
#else
   subroutine compute_l_sum_NS(infield,Plm_slice,outN,outS)

      !-- Input variables:
      complex(cp), intent(in) :: infield(lm_max)
      real(cp),    intent(in) :: Plm_slice(lm_max)
      
      !-- Output variables:
      real(cp), intent(out) :: outN(nrp), outS(nrp)

      !-- Local variables:
      integer :: mc,lmS,lm,mc_old
      complex(cp) :: temp_field(lm_max)

      temp_field=infield*Plm_slice
      do mc=1,n_m_max
         outN(mc)=sum(temp_field(lStart(mc):lStop(mc)))
         outS(mc)=sum(temp_field(lStart(mc):lStop(mc):2)) &
              & - sum(temp_field(lStart(mc)+1:lStop(mc):2))
      end do
#if 0
      sES=zero    ! One equatorial symmetry
      sEA=zero    ! The other equatorial symmetry
      mc_old=1
      do lm=1,lm_max
         mc=lm2mc(lm)
         if ( mc /= mc_old ) then
            ! next m
            outN(mc_old) = sES + sEA
            outS(mc_old) = sES - sEA
            sES=zero    ! One equatorial symmetry
            sEA=zero    ! The other equatorial symmetry
            mc_old = mc
         end if
         temp = infield(lm)*Plm_slice(lm)
         if ( modulo(lm-lStart(mc),2) == 0 ) then
            sES = sES + temp
         else
            sEA = sEA + temp
         end if
      end do
      outN(2*mc_old-1) = real(sES + sEA)
      outN(2*mc_old  ) =aimag(sES + sEA)
      outS(2*mc_old-1) = real(sES - sEA)
      outS(2*mc_old  ) =aimag(sES - sEA)
#endif
   end subroutine compute_l_sum_NS
#endif
#endif
!------------------------------------------------------------------------------
end module legendre_spec_to_grid
