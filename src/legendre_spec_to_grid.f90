#include "perflib_preproc.cpp"
module legendre_spec_to_grid

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use precision_mod
   use truncation, only: lm_max, n_m_max, nrp, l_max, l_axi, n_theta_max, minc
   use LMmapping, only: map_dist_st
   use blocking, only: nfs, sizeThetaB, lm2mc, lm2
   use horizontal_data, only: Plm, dPlm, lStart, lStop, lmOdd, D_mc2m, &
       &                      osn2, Plm_loc, dPlm_loc
   use constants, only: zero, half, one, ci
   use useful, only: abortRun

   implicit none
 
   private

   public :: lmAS2pt, leg_scal_to_grad_spat, leg_scal_to_spat, sh_to_spat_ml, &
   &         leg_polsphtor_to_spat, leg_pol_to_grad_spat, leg_dphi_vec

contains

   subroutine leg_polsphtor_to_spat(l_calc, nThetaStart, Qlm, bhG, bhC, brc, &
              &                     btc, bpc, n_thetas)
      !
      ! Take Q,S,T and transform them to vectors
      !
      
      !-- Input variables:
      logical,     intent(in) :: l_calc
      integer,     intent(in) :: nThetaStart
      complex(cp), intent(in) :: Qlm(lm_max) ! Poloidal
      complex(cp), intent(in) :: bhG(lm_max)
      complex(cp), intent(in) :: bhC(lm_max)
      integer, optional, intent(in) :: n_thetas
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: brc(nrp,nfs), btc(nrp,nfs), bpc(nrp,nfs)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max)
      real(cp) :: PlmC(lm_max)
    
      !-- Local variables:
      complex(cp) :: brES,brEA
      integer :: nThetaN,nThetaS,nThetaNHS,sThetaB
      integer :: mc,lm,lmS
      real(cp) :: dm
    
      complex(cp) :: bhN1M(n_m_max),bhN2M(n_m_max),bhN,bhN1,bhN2
      complex(cp) :: bhS1M(n_m_max),bhS2M(n_m_max),bhS,bhS1,bhS2
    
      nThetaNHS=(nThetaStart-1)/2

      if ( present(n_thetas) ) then
         sThetaB = n_thetas
      else
         sThetaB = sizeThetaB
      end if
    
      if ( l_calc ) then ! not a boundary or derivs required
         do nThetaN=1,sThetaB,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
    
            !--- Loop over all orders m: (numbered by mc)
            PERFON_I('TFG_2')
            do mc=1,n_m_max
               dm = D_mc2m(mc)
               lmS=lStop(mc)
               brES   =zero
               brEA   =zero
               !--- 6 add/mult, 26 dble words
               do lm=lStart(mc),lmS-1,2
                  brES   =brES + Qlm(lm)  *Plm(lm,nThetaNHS)
                  brEA   =brEA + Qlm(lm+1)*Plm(lm+1,nThetaNHS)
                  PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                  PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                  PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  brES   =brES + Qlm(lmS)*Plm(lmS,nThetaNHS)
                  PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                  PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
               end if
               brc(2*mc-1,nThetaN)   = real(brES+brEA)
               brc(2*mc  ,nThetaN)   =aimag(brES+brEA)
               brc(2*mc-1,nThetaS)   = real(brES-brEA)
               brc(2*mc  ,nThetaS)   =aimag(brES-brEA)
            end do
    
            do mc=1,n_m_max
               lmS=lStop(mc)
               bhN1=zero
               bhS1=zero
               bhN2=zero
               bhS2=zero
               !--- 8 add/mult, 20 dble words
               do lm=lStart(mc),lmS-1,2
                  bhN1=bhN1+bhG(lm)*PlmG(lm)+bhG(lm+1)*PlmG(lm+1)
                  bhS1=bhS1-bhG(lm)*PlmC(lm)+bhG(lm+1)*PlmC(lm+1)
                  bhN2=bhN2+bhC(lm)*PlmC(lm)+bhC(lm+1)*PlmC(lm+1)
                  bhS2=bhS2-bhC(lm)*PlmG(lm)+bhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  bhN1=bhN1+bhG(lmS)*PlmG(lmS)
                  bhS1=bhS1-bhG(lmS)*PlmC(lmS)
                  bhN2=bhN2+bhC(lmS)*PlmC(lmS)
                  bhS2=bhS2-bhC(lmS)*PlmG(lmS)
               end if
               bhN1M(mc)=half*bhN1
               bhS1M(mc)=half*bhS1
               bhN2M(mc)=half*bhN2
               bhS2M(mc)=half*bhS2
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
    
         end do      ! End global loop over nTheta
    
         !PERFOFF
    
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sThetaB
               do mc=2*n_m_max+1,nrp
                  brc(mc,nThetaN)   =0.0_cp
                  btc(mc,nThetaN)   =0.0_cp
                  bpc(mc,nThetaN)   =0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
         !PERFOFF
    
      else   ! boundary ?
    
         !-- Calculation for boundary r_cmb or r_icb:
         do nThetaN=1,sThetaB,2
            nThetaS=nThetaN+1
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
    
            do mc=1,n_m_max
               dm =D_mc2m(mc)
               lmS=lStop(mc)
    
               !------ br = r^2 B_r , bt = r sin(theta) B_theta , bp= r sin(theta) B_phi
               brES=zero
               brEA=zero
               do lm=lStart(mc),lmS-1,2
                  brES=brES + Qlm(lm)  *Plm(lm,nThetaNHS)
                  brEA=brEA + Qlm(lm+1)*Plm(lm+1,nThetaNHS)
                  PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
                  PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
                  PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  brES=brES+Qlm(lm)*Plm(lm,nThetaNHS)
                  PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
                  PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
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
                  bhN1=bhN1+bhG(lm)*PlmG(lm)+bhG(lm+1)*PlmG(lm+1)
                  bhS1=bhS1-bhG(lm)*PlmC(lm)+bhG(lm+1)*PlmC(lm+1)
                  bhN2=bhN2+bhC(lm)*PlmC(lm)+bhC(lm+1)*PlmC(lm+1)
                  bhS2=bhS2-bhC(lm)*PlmG(lm)+bhC(lm+1)*PlmG(lm+1)
               end do
               if ( lmOdd(mc) ) then
                  bhN1=bhN1+bhG(lmS)*PlmG(lmS)
                  bhS1=bhS1-bhG(lmS)*PlmC(lmS)
                  bhN2=bhN2+bhC(lmS)*PlmC(lmS)
                  bhS2=bhS2-bhC(lmS)*PlmG(lmS)
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
    
         end do    ! End loop over nThetaN
    
         !-- Zero out terms with index mc > n_m_max :
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sThetaB
               do mc=2*n_m_max+1,nrp
                  brc(mc,nThetaN)=0.0_cp
                  btc(mc,nThetaN)=0.0_cp
                  bpc(mc,nThetaN)=0.0_cp
               end do
            end do
         end if

      end if  ! boundary ? nBc?
    
   end subroutine leg_polsphtor_to_spat
!------------------------------------------------------------------------------
   subroutine leg_pol_to_grad_spat(l_calc,nThetaStart,Qlm,dvrdtc)
      !
      ! Poloidal to d/dth
      !

      !-- Input:
      logical,     intent(in) :: l_calc
      integer,     intent(in) :: nThetaStart
      complex(cp), intent(in) :: Qlm(lm_max)
    
      !-- Output fields
      real(cp), intent(out) :: dvrdtc(nrp,nfs)
    
      !------ Legendre Polynomials
      real(cp) :: PlmG(lm_max)
      real(cp) :: PlmC(lm_max)
    
      !-- Local:
      complex(cp) :: dvrdtEA,dvrdtES
    
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS
      real(cp) :: dm
    
    
      nThetaNHS=(nThetaStart-1)/2
    
      if ( l_calc ) then ! not a boundary
    
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
    
            do mc=1,n_m_max
               dm =D_mc2m(mc)
               lmS=lStop(mc)
               dvrdtES=zero
               dvrdtEA=zero
               !--- 8 add/mult, 29 dble words
               do lm=lStart(mc),lmS-1,2
                  dvrdtEA =dvrdtEA + Qlm(lm)*  dPlm(lm,nThetaNHS)
                  dvrdtES =dvrdtES + Qlm(lm+1)*dPlm(lm+1,nThetaNHS)
               end do
               if ( lmOdd(mc) ) then
                  dvrdtEA =dvrdtEA + Qlm(lmS)*dPlm(lmS,nThetaNHS)
                  PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
                  PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
               end if
               dvrdtc(2*mc-1,nThetaN)= real(dvrdtES+dvrdtEA)
               dvrdtc(2*mc  ,nThetaN)=aimag(dvrdtES+dvrdtEA)
               dvrdtc(2*mc-1,nThetaS)= real(dvrdtES-dvrdtEA)
               dvrdtc(2*mc  ,nThetaS)=aimag(dvrdtES-dvrdtEA)
            end do

         end do
    
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  dvrdtc(mc,nThetaN)=0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
    
      end if  ! boundary ? nBc?
    
   end subroutine leg_pol_to_grad_spat
!------------------------------------------------------------------------------
   subroutine leg_dphi_vec(l_calc,nThetaStart,vrc,vtc,vpc,dvrdpc,dvtdpc,dvpdpc)
      !
      ! Take the phi derivative of an input vector
      !

      !-- Input:
      logical, intent(in) :: l_calc
      integer, intent(in) :: nThetaStart
      real(cp), intent(in) :: vrc(nrp,nfs), vtc(nrp,nfs), vpc(nrp,nfs)
    
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: dvrdpc(nrp, nfs), dvtdpc(nrp,nfs), dvpdpc(nrp, nfs)
    
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc
      real(cp) :: dm,dmT
    
      nThetaNHS=(nThetaStart-1)/2
    
      if ( l_calc  ) then ! not a boundary or derivs required
    
         do nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
    
            !--- Calculate phi derivatives:
            do mc=1,n_m_max
               dm=D_mc2m(mc)
               dvrdpc(2*mc-1,nThetaN)=-dm*vrc(2*mc  ,nThetaN)
               dvrdpc(2*mc  ,nThetaN)= dm*vrc(2*mc-1,nThetaN)
               dvrdpc(2*mc-1,nThetaS)=-dm*vrc(2*mc  ,nThetaS)
               dvrdpc(2*mc  ,nThetaS)= dm*vrc(2*mc-1,nThetaS)
            end do
            do mc=1,n_m_max
               dmT=D_mc2m(mc)*osn2(nThetaNHS)
               dvtdpc(2*mc-1,nThetaN)=-dmT*vtc(2*mc  ,nThetaN)
               dvtdpc(2*mc  ,nThetaN)= dmT*vtc(2*mc-1,nThetaN)
               dvtdpc(2*mc-1,nThetaS)=-dmT*vtc(2*mc  ,nThetaS)
               dvtdpc(2*mc  ,nThetaS)= dmT*vtc(2*mc-1,nThetaS)
               dvpdpc(2*mc-1,nThetaN)=-dmT*vpc(2*mc  ,nThetaN)
               dvpdpc(2*mc  ,nThetaN)= dmT*vpc(2*mc-1,nThetaN)
               dvpdpc(2*mc-1,nThetaS)=-dmT*vpc(2*mc  ,nThetaS)
               dvpdpc(2*mc  ,nThetaS)= dmT*vpc(2*mc-1,nThetaS)
            end do   ! End of loop over oder m numbered by mc
    
         end do      ! End global loop over nTheta
    
    
         !-- Zero out terms with index mc > n_m_max:
         if ( n_m_max < nrp/2 ) then
            do nThetaN=1,sizeThetaB
               do mc=2*n_m_max+1,nrp
                  dvrdpc(mc,nThetaN)=0.0_cp
                  dvtdpc(mc,nThetaN)=0.0_cp
                  dvpdpc(mc,nThetaN)=0.0_cp
               end do
            end do  ! loop over nThetaN (theta)
         end if
    
      end if  ! boundary ? nBc?
    
   end subroutine leg_dphi_vec
!------------------------------------------------------------------------------
   subroutine qst_to_spat_ml(m, Qlm, Slm, Tlm, brc, btc, bpc, lcut)
      !
      ! Take Q,S,T and transform them to a vector
      !
      
      !-- Input variables:
      integer,     intent(in) :: m
      integer,     intent(in) :: lcut
      complex(cp), intent(in) :: Qlm(:) ! Poloidal
      complex(cp), intent(in) :: Slm(:)
      complex(cp), intent(in) :: Tlm(:)
    
      complex(cp), intent(out) :: brc(n_theta_max), btc(n_theta_max), bpc(n_theta_max)
    
      !------ Legendre Polynomials helpers
      complex(cp) :: bhG(size(Slm)), bhC(size(Slm))
      real(cp) :: PlmG(size(Slm)), PlmC(size(Slm))
    
      !-- Local variables:
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lj,lb,lu
      real(cp) :: dm
      complex(cp) :: brES,brEA,bhN,bhN1,bhN2,bhS,bhS1,bhS2
    
      nThetaNHS=0
      lb = map_dist_st%lm2(m, m)
      lu = map_dist_st%lm2(l_max, m)
      mc = m/minc+1
      dm = real(m,cp)
      
      do lj=1,lu-lb+1
         bhG(lj)=Slm(lj)-ci*Tlm(lj)
         bhC(lj)=Slm(lj)+ci*Tlm(lj)
      end do

      do nThetaN=1,n_theta_max,2   ! Loop over thetas for north HS
         nThetaS  =nThetaN+1      ! same theta but for southern HS
         nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
         !--- Loop over all orders m: (numbered by mc)
         brES=zero
         brEA=zero
         !--- 6 add/mult, 26 dble words
         do lj=1,lu-lb,2
            brES   =brES + Qlm(lj)  *Plm_loc(lb+lj-1,nThetaNHS)
            brEA   =brEA + Qlm(lj+1)*Plm_loc(lb+lj,nThetaNHS)
            PlmG(lj)=dPlm_loc(lb+lj-1,nThetaNHS)-dm*Plm_loc(lb+lj-1,nThetaNHS)
            PlmC(lj)=dPlm_loc(lb+lj-1,nThetaNHS)+dm*Plm_loc(lb+lj-1,nThetaNHS)
            PlmG(lj+1)=dPlm_loc(lb+lj,nThetaNHS)-dm*Plm_loc(lb+lj,nThetaNHS)
            PlmC(lj+1)=dPlm_loc(lb+lj,nThetaNHS)+dm*Plm_loc(lb+lj,nThetaNHS)
         end do
         if ( mod(lu-lb,2) == 0 ) then
            brES   =brES + Qlm(lu-lb+1)*Plm(lu,nThetaNHS)
            PlmG(lu-lb+1)=dPlm(lu,nThetaNHS)-dm*Plm(lu,nThetaNHS)
            PlmC(lu-lb+1)=dPlm(lu,nThetaNHS)+dm*Plm(lu,nThetaNHS)
         end if
         brc(nThetaN)=brES+brEA
         brc(nThetaS)=brES-brEA
 
         bhN1=zero
         bhS1=zero
         bhN2=zero
         bhS2=zero
         !--- 8 add/mult, 20 dble words
         do lj=1,lu-lb,2
            bhN1=bhN1+bhG(lj)*PlmG(lj)+bhG(lj+1)*PlmG(lj+1)
            bhS1=bhS1-bhG(lj)*PlmC(lj)+bhG(lj+1)*PlmC(lj+1)
            bhN2=bhN2+bhC(lj)*PlmC(lj)+bhC(lj+1)*PlmC(lj+1)
            bhS2=bhS2-bhC(lj)*PlmG(lj)+bhC(lj+1)*PlmG(lj+1)
         end do
         if ( mod(lu-lb,2) == 0 ) then
            bhN1=bhN1+bhG(lu-lb+1)*PlmG(lu-lb+1)
            bhS1=bhS1-bhG(lu-lb+1)*PlmC(lu-lb+1)
            bhN2=bhN2+bhC(lu-lb+1)*PlmC(lu-lb+1)
            bhS2=bhS2-bhC(lu-lb+1)*PlmG(lu-lb+1)
         end if
         btc(nThetaN)=half*(bhN1+bhN2)
         btc(nThetaS)=half*(bhS1+bhS2)
         bhN         =half*(bhN1-bhN2)
         bhS         =half*(bhS1-bhS2)
         bpc(nThetaN)=cmplx(aimag(bhN),-real(bhN),cp)
         bpc(nThetaS)=cmplx(aimag(bhS),-real(bhS),cp)
 
      end do      ! End global loop over nTheta
 
   end subroutine qst_to_spat_ml
!------------------------------------------------------------------------------
   subroutine sh_to_spat_ml(m, Slm, sc, lcut)
      !
      ! Legendre transform from (l) to (theta) for a scalar input field
      !


      !-- Input variable
      integer,     intent(in) :: m
      integer,     intent(in) :: lcut
      complex(cp), intent(in) :: Slm(:)

      !-- Output variables
      complex(cp), intent(out) :: sc(n_theta_max)

      !-- Local variables
      integer :: mc, lb, lu, lj
      integer :: nThetaN, nThetaS, nThetaNHS
      complex(cp) :: sES, sEA

      lb = map_dist_st%lm2(m, m)
      lu = map_dist_st%lm2(l_max, m)
      mc = m/minc+1
      nThetaNHS=0
      do nThetaN=1,n_theta_max,2   ! Loop over thetas for one HS
         nThetaS=nThetaN+1  ! same theta but at other HS
         nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
         sES=zero ! One equatorial symmetry
         sEA=zero ! The other equatorial symmetry
         do lj=1,lu-lb,2
            sES=sES+Slm(lj)  *Plm_loc(lb+lj-1,nThetaNHS)
            sEA=sEA+Slm(lj+1)*Plm_loc(lb+lj,nThetaNHS)
         end do
         if ( mod(lu-lb,2) == 0 ) sES=sES+Slm(lu-lb+1)*Plm_loc(lu,nThetaNHS)
         sc(nThetaN)=sES+sEA
         sc(nThetaS)=sES-SEA
      end do

   end subroutine sh_to_spat_ml
!------------------------------------------------------------------------------
   subroutine leg_scal_to_spat(nThetaStart, Slm, sc)
      !
      ! Legendre transform from (l) to (theta) for a scalar input field
      !


      !-- Input variable
      integer,     intent(in) :: nThetaStart
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: sc(nrp,nfs)

      !-- Local variables
      integer :: mc, lm, lmS
      integer :: nThetaN, nThetaS, nThetaNHS
      complex(cp) :: sES, sEA

      nThetaNHS=(nThetaStart-1)/2
      do nThetaN=1,sizeThetaB,2   ! Loop over thetas for one HS
         nThetaS  =nThetaN+1  ! same theta but at other HS
         nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
         do mc=1,n_m_max
            lmS=lStop(mc)
            sES=zero ! One equatorial symmetry
            sEA=zero ! The other equatorial symmetry
            do lm=lStart(mc),lmS-1,2
               sES=sES+Slm(lm)  *Plm(lm,nThetaNHS)
               sEA=sEA+Slm(lm+1)*Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(mc) ) sES=sES+Slm(lmS)*Plm(lmS,nThetaNHS)
            sc(2*mc-1,nThetaN)= real(sES+sEA)
            sc(2*mc  ,nThetaN)=aimag(sES+sEA)
            sc(2*mc-1,nThetaS)= real(sES-sEA)
            sc(2*mc  ,nThetaS)=aimag(sES-sEA)
         end do
      end do

      if ( n_m_max < nrp/2 ) then
         do nThetaN=1,sizeThetaB
            do mc=2*n_m_max+1,nrp
               sc(mc,nThetaN)=0.0_cp
            end do
         end do  ! loop over nThetaN (theta)
      end if

   end subroutine leg_scal_to_spat
!------------------------------------------------------------------------------
   subroutine leg_scal_to_grad_spat(nThetaStart, slm, gradtc, gradpc)
      !
      ! Transform s(l) into dsdt(theta) and dsdp(theta)
      !

      !-- Input variable
      integer,     intent(in) :: nThetaStart
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: gradtc(nrp,nfs), gradpc(nrp,nfs)

      !-- Local variables
      integer :: mc, lm, lmS
      integer :: nThetaN, nThetaS, nThetaNHS
      real(cp) :: dm
      complex(cp) :: gradtcES, gradtcEA, sES, sEA

      nThetaNHS=(nThetaStart-1)/2

      do nThetaN=1,sizeThetaB,2   ! Loop over thetas for north HS
         nThetaS  =nThetaN+1      ! same theta but for southern HS
         nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
         do mc=1,n_m_max
            dm =D_mc2m(mc)
            lmS=lStop(mc)
            sES=zero  ! One equatorial symmetry
            sEA=zero  ! The other equatorial symmetry
            gradtcES=zero
            gradtcEA=zero
            do lm=lStart(mc),lmS-1,2
               gradtcEA=gradtcEA + Slm(lm)*  dPlm(lm  ,nThetaNHS)
               gradtcES=gradtcES + Slm(lm+1)*dPlm(lm+1,nThetaNHS)
               sES     =sES      + Slm(lm)  * Plm(lm  ,nThetaNHS)
               sEA     =sEA      + Slm(lm+1)* Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(mc) ) then
               gradtcEA=gradtcEA + Slm(lmS)*dPlm(lmS,nThetaNHS)
               sES     =sES      + Slm(lmS)* Plm(lmS,nThetaNHS)
            end if
            gradtc(2*mc-1,nThetaN)= real(gradtcES+gradtcEA)
            gradtc(2*mc  ,nThetaN)=aimag(gradtcES+gradtcEA)
            gradtc(2*mc-1,nThetaS)= real(gradtcES-gradtcEA)
            gradtc(2*mc  ,nThetaS)=aimag(gradtcES-gradtcEA)

            gradpc(2*mc-1,nThetaN)=-dm*aimag(sES+sEA)
            gradpc(2*mc  ,nThetaN)= dm* real(sES+sEA)
            gradpc(2*mc-1,nThetaS)=-dm*aimag(sES-sEA)
            gradpc(2*mc  ,nThetaS)= dm* real(sES-sEA)
         end do
    
      end do

      if ( n_m_max < nrp/2 ) then
         do nThetaN=1,sizeThetaB
            do mc=2*n_m_max+1,nrp
               gradtc(mc,nThetaN)=0.0_cp
               gradpc(mc,nThetaN)=0.0_cp
            end do
         end do  ! loop over nThetaN (theta)
      end if

   end subroutine leg_scal_to_grad_spat
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
end module legendre_spec_to_grid
