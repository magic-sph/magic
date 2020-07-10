module legendre_spec_to_grid

   use precision_mod
   use truncation, only: lm_max, n_m_max, l_max, l_axi, n_theta_max, minc, n_phi_max
   use LMmapping, only: map_dist_st, map_glbl_st
   use horizontal_data, only: Plm, dPlm, lStart, lStop, lmOdd, D_mc2m, &
       &                      Plm_loc, dPlm_loc
   use constants, only: zero, half, one, ci
   use fft, only: ifft_many

   implicit none
 
   private

   public :: lmAS2pt, qst_to_spat_ml, sh_to_spat_ml, qst_to_spat, sphtor_to_spat,&
   &         sh_to_spat, sh_to_grad_spat

contains

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
   subroutine qst_to_spat(Qlm, Slm, Tlm, brc, btc, bpc)
      !
      ! Vector spherical harmonic transform: take Q,S,T and transform them 
      ! to a vector field
      !
      
      !-- Input variables:
      complex(cp), intent(in) :: Qlm(lm_max) ! Poloidal
      complex(cp), intent(in) :: Slm(lm_max) ! Spheroidal
      complex(cp), intent(in) :: Tlm(lm_max) ! Toroidal
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: brc(n_phi_max,n_theta_max)
      real(cp), intent(out) :: btc(n_phi_max,n_theta_max)
      real(cp), intent(out) :: bpc(n_phi_max,n_theta_max)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max),PlmC(lm_max)
      complex(cp) :: bhG(lm_max),bhC(lm_max)
      complex(cp) :: tmpr(n_phi_max/2+1,n_theta_max)
      complex(cp) :: tmpt(n_phi_max/2+1,n_theta_max)
      complex(cp) :: tmpp(n_phi_max/2+1,n_theta_max)
    
      !-- Local variables:
      complex(cp) :: brES,brEA
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS
      real(cp) :: dm
    
      complex(cp) :: bhN1M(n_m_max),bhN2M(n_m_max),bhN,bhN1,bhN2
      complex(cp) :: bhS1M(n_m_max),bhS2M(n_m_max),bhS,bhS1,bhS2
    
      bhG(1)=zero
      bhC(1)=zero
      do lm=2,lm_max
         bhG(lm)=Slm(lm)-ci*Tlm(lm)
         bhC(lm)=Slm(lm)+ci*Tlm(lm)
      end do

      nThetaNHS=0
      do nThetaN=1,n_theta_max,2   ! Loop over thetas for north HS
         nThetaS  =nThetaN+1      ! same theta but for southern HS
         nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
         !--- Loop over all orders m: (numbered by mc)
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
            tmpr(mc,nThetaN)=brES+brEA
            tmpr(mc,nThetaS)=brES-brEA
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
            tmpt(mc,nThetaN)=bhN1M(mc)+bhN2M(mc)
            bhN             =bhN1M(mc)-bhN2M(mc)
            tmpt(mc,nThetaS)=bhS1M(mc)+bhS2M(mc)
            bhS             =bhS1M(mc)-bhS2M(mc)
            tmpp(mc,nThetaN)=-ci*bhN
            tmpp(mc,nThetaS)=-ci*bhS
         end do
 
      end do      ! End global loop over nTheta
 
      !-- Zero out terms with index mc > n_m_max:
      if ( n_m_max < n_phi_max/2+1 ) then
         do nThetaN=1,n_theta_max
            tmpr(n_m_max+1:,nThetaN)=zero
            tmpt(n_m_max+1:,nThetaN)=zero
            tmpp(n_m_max+1:,nThetaN)=zero
         end do  ! loop over nThetaN (theta)
      end if

      if ( .not. l_axi ) then
         call ifft_many(tmpr,brc)
         call ifft_many(tmpt,btc)
         call ifft_many(tmpp,bpc)
      end if
    
   end subroutine qst_to_spat
!------------------------------------------------------------------------------
   subroutine sphtor_to_spat(Slm, Tlm, btc, bpc)
      !
      ! Use spheroidal and toroidal potentials to transform them to angular
      ! vector components btheta and bphi
      !
      
      !-- Input variables:
      complex(cp), intent(in) :: Slm(lm_max)
      complex(cp), intent(in) :: Tlm(lm_max)
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: btc(n_phi_max,n_theta_max)
      real(cp), intent(out) :: bpc(n_phi_max,n_theta_max)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max),PlmC(lm_max)
      complex(cp) :: bhG(lm_max),bhC(lm_max)
      complex(cp) :: tmpt(n_phi_max/2+1,n_theta_max)
      complex(cp) :: tmpp(n_phi_max/2+1,n_theta_max)
    
      !-- Local variables:
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS
      real(cp) :: dm
    
      complex(cp) :: bhN1M(n_m_max),bhN2M(n_m_max),bhN,bhN1,bhN2
      complex(cp) :: bhS1M(n_m_max),bhS2M(n_m_max),bhS,bhS1,bhS2
    
      nThetaNHS=0

      do lm=1,lm_max
         bhG(lm)=Slm(lm)-ci*Tlm(lm)
         bhC(lm)=Slm(lm)+ci*Tlm(lm)
      end do

      do nThetaN=1,n_theta_max,2   ! Loop over thetas for north HS
         nThetaS  =nThetaN+1      ! same theta but for southern HS
         nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
         !--- Loop over all orders m: (numbered by mc)
         do mc=1,n_m_max
            dm = D_mc2m(mc)
            lmS=lStop(mc)
            !--- 6 add/mult, 26 dble words
            do lm=lStart(mc),lmS-1,2
               PlmG(lm)=dPlm(lm,nThetaNHS)-dm*Plm(lm,nThetaNHS)
               PlmC(lm)=dPlm(lm,nThetaNHS)+dm*Plm(lm,nThetaNHS)
               PlmG(lm+1)=dPlm(lm+1,nThetaNHS)-dm*Plm(lm+1,nThetaNHS)
               PlmC(lm+1)=dPlm(lm+1,nThetaNHS)+dm*Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(mc) ) then
               PlmG(lmS)=dPlm(lmS,nThetaNHS)-dm*Plm(lmS,nThetaNHS)
               PlmC(lmS)=dPlm(lmS,nThetaNHS)+dm*Plm(lmS,nThetaNHS)
            end if
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
            tmpt(mc,nThetaN)=bhN1M(mc)+bhN2M(mc)
            bhN             =bhN1M(mc)-bhN2M(mc)
            tmpt(mc,nThetaS)=bhS1M(mc)+bhS2M(mc)
            bhS             =bhS1M(mc)-bhS2M(mc)
            tmpp(mc,nThetaN)=-ci*bhN
            tmpp(mc,nThetaS)=-ci*bhS
         end do
 
      end do      ! End global loop over nTheta
 
      !-- Zero out terms with index mc > n_m_max:
      if ( n_m_max < n_phi_max/2+1 ) then
         do nThetaN=1,n_theta_max
            tmpt(n_m_max+1:,nThetaN)=zero
            tmpp(n_m_max+1:,nThetaN)=zero
         end do  ! loop over nThetaN (theta)
      end if

      if ( .not. l_axi ) then
         call ifft_many(tmpt,btc)
         call ifft_many(tmpp,bpc)
      end if
    
   end subroutine sphtor_to_spat
!------------------------------------------------------------------------------
   subroutine sh_to_spat(Slm, sc)
      !
      ! Spherical Harmonic Transform for a scalar input field
      !

      !-- Input variable
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: sc(n_phi_max,n_theta_max)

      !-- Local variables
      complex(cp) :: tmp(n_phi_max/2+1,n_theta_max)
      integer :: mc, lm, lmS
      integer :: nThetaN, nThetaS, nThetaNHS
      complex(cp) :: sES, sEA

      nThetaNHS=0
      do nThetaN=1,n_theta_max,2   ! Loop over thetas for one HS
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
            tmp(mc,nThetaN)=sES+sEA
            tmp(mc,nThetaS)=sES-sEA
         end do
      end do

      if ( n_m_max < n_phi_max/2+1 ) then
         do nThetaN=1,n_theta_max
            tmp(n_m_max+1:,nThetaN)=zero
         end do  ! loop over nThetaN (theta)
      end if

      if ( .not. l_axi ) call ifft_many(tmp,sc)

   end subroutine sh_to_spat
!------------------------------------------------------------------------------
   subroutine sh_to_grad_spat(slm, gradtc, gradpc)
      !
      ! Transform s(l) into dsdt(theta) and dsdp(theta)
      !

      !-- Input variable
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_phi_max,n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max,n_theta_max)

      !-- Local variables
      complex(cp) :: tmpt(n_phi_max/2+1,n_theta_max),tmpp(n_phi_max/2+1,n_theta_max)
      integer :: mc, lm, lmS, nThetaN, nThetaS, nThetaNHS
      real(cp) :: dm
      complex(cp) :: gradtcES, gradtcEA, sES, sEA

      nThetaNHS=0
      do nThetaN=1,n_theta_max,2   ! Loop over thetas for north HS
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
            tmpt(mc,nThetaN)=gradtcES+gradtcEA
            tmpt(mc,nThetaS)=gradtcES-gradtcEA
            tmpp(mc,nThetaN)=dm*ci*(sES+sEA)
            tmpp(mc,nThetaS)=dm*ci*(sES-sEA)
         end do
      end do

      if ( n_m_max < n_phi_max/2+1 ) then
         do nThetaN=1,n_theta_max
            tmpt(n_m_max+1:,nThetaN)=zero
            tmpp(n_m_max+1:,nThetaN)=zero
         end do  ! loop over nThetaN (theta)
      end if

      if ( .not. l_axi ) then
         call ifft_many(tmpt,gradtc)
         call ifft_many(tmpp,gradpc)
      end if

   end subroutine sh_to_grad_spat
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
            lm=map_glbl_st%lm2(l,0)
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
