module shtransforms

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_m_max, l_max, l_axi, n_theta_max, minc, &
       &                 n_phi_max, lmP_max, m_max
   use blocking, only: lmP2l, lmP2lm
   use horizontal_data, only: gauleg, O_sin_theta_E2
   use plms_theta, only: plm_theta
   use constants, only: zero, half, one, ci, pi, two
   use fft, only: fft_many, ifft_many
   use parallel_mod, only: get_openmp_blocks

   implicit none

   !-- Legendres:
   real(cp), allocatable :: Plm(:,:)
   real(cp), allocatable :: wPlm(:,:)
   real(cp), allocatable :: wdPlm(:,:)
   real(cp), allocatable :: dPlm(:,:)

   !-- Limiting l for a given m, used in legtf
   real(cp), public, allocatable :: D_mc2m(:)
   integer, allocatable :: lStart(:),lStop(:)
   integer, allocatable :: lStartP(:),lStopP(:)
   logical, allocatable :: lmOdd(:),lmOddP(:)

   public :: initialize_transforms, finalize_transforms, native_qst_to_spat,   &
   &         native_sphtor_to_spat, native_sph_to_spat, native_spat_to_sph,    &
   &         native_spat_to_sph_tor, native_sph_to_grad_spat,                  &
   &         native_toraxi_to_spat, native_spat_to_SH_axi

contains

   subroutine initialize_transforms

      !-- Local variables
      integer :: n_theta, lmP, norm, l, lm, mc, m
      real(cp) :: colat,theta_ord(n_theta_max), gauss(n_theta_max)
      real(cp) :: plma(lmP_max), dtheta_plma(lmP_max)

      allocate( Plm(lm_max,n_theta_max/2) )
      allocate( wPlm(lmP_max,n_theta_max/2) )
      allocate( wdPlm(lmP_max,n_theta_max/2) )
      allocate( dPlm(lm_max,n_theta_max/2) )
      bytes_allocated = bytes_allocated+(lm_max*n_theta_max+ &
      &                 lmP_max*n_theta_max)*SIZEOF_DEF_REAL

      allocate( D_mc2m(n_m_max) )
      bytes_allocated = bytes_allocated+n_m_max*SIZEOF_DEF_REAL

      !-- Limiting l for a given m, used in legtf
      allocate( lStart(n_m_max),lStop(n_m_max) )
      allocate( lStartP(n_m_max),lStopP(n_m_max) )
      allocate( lmOdd(n_m_max),lmOddP(n_m_max) )
      bytes_allocated = bytes_allocated+6*n_m_max*SIZEOF_INTEGER

      norm=2 ! norm chosen so that a surface integral over
             ! any ylm**2 is 1.

      !-- Calculate grid points and weights for the
      !--   Gauss-Legendre integration of the plms:
      call gauleg(-one,one,theta_ord,gauss,n_theta_max)

      do n_theta=1,n_theta_max/2  ! Loop over colat in NHS
         colat=theta_ord(n_theta)
         !----- plmtheta calculates plms and their derivatives
         !      up to degree and order l_max+1 and m_max at
         !      the points cos(theta_ord(n_theta)):
         call plm_theta(colat,l_max+1,m_max,minc,plma,dtheta_plma,lmP_max,norm)
         do lmP=1,lmP_max
            l=lmP2l(lmP)
            if ( l <= l_max ) then
               lm=lmP2lm(lmP)
               Plm(lm,n_theta) =plma(lmP)
               dPlm(lm,n_theta)=dtheta_plma(lmP)
            end if
            wPlm(lmP,n_theta) =two*pi*gauss(n_theta)*plma(lmP)
            wdPlm(lmP,n_theta)=two*pi*gauss(n_theta)*dtheta_plma(lmP)
         end do
      end do

      !-- Build auxiliary index arrays for Legendre transform:
      !   lStartP, lStopP give start and end positions in lmP-block.
      !   lStart, lStop give start and end positions in lm-block.
      lStartP(1)=1
      lStopP(1) =l_max+2
      lStart(1) =1
      lStop(1)  =l_max+1
      D_mc2m(1)=0
      if ( mod(l_max,2) == 0 ) then
         lmOdd(1) =.true.
         lmOddP(1)=.false.
      else
         lmOdd(1) =.false.
         lmOddP(1)=.true.
      end if
      do mc=2,n_m_max
         m=(mc-1)*minc
         D_mc2m(mc) =real(m,cp)
         lStartP(mc)=lStopP(mc-1)+1
         lStopP(mc) =lStartP(mc) +l_max-m+1
         lStart(mc) =lStop(mc-1) +1
         lStop(mc)  =lStart(mc)  +l_max-m
         if ( mod(lStop(mc)-lStart(mc),2) == 0 ) then
            lmOdd(mc) =.true.
            lmOddP(mc)=.false.
         else
            lmOdd(mc) =.false.
            lmOddP(mc)=.true.
         end if
      end do

   end subroutine initialize_transforms
!------------------------------------------------------------------------------
   subroutine finalize_transforms

      deallocate( Plm, wPlm, wdPlm, dPlm)
      deallocate( D_mc2m, lStart, lStop, lStartP, lStopP, lmOdd, lmOddP)

   end subroutine finalize_transforms
!------------------------------------------------------------------------------
   subroutine native_qst_to_spat(Qlm, Slm, Tlm, brc, btc, bpc)
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
      real(cp), intent(out) :: brc(n_theta_max,n_phi_max)
      real(cp), intent(out) :: btc(n_theta_max,n_phi_max)
      real(cp), intent(out) :: bpc(n_theta_max,n_phi_max)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max),PlmC(lm_max)
      complex(cp) :: bhG(lm_max),bhC(lm_max)
      complex(cp) :: tmpr(n_theta_max,n_phi_max/2+1)
      complex(cp) :: tmpt(n_theta_max,n_phi_max/2+1)
      complex(cp) :: tmpp(n_theta_max,n_phi_max/2+1)
    
      !-- Local variables:
      complex(cp) :: brES,brEA
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS, nThStart,nThStop
      real(cp) :: dm
    
      complex(cp) :: bhN1M,bhN2M,bhN,bhN1,bhN2
      complex(cp) :: bhS1M,bhS2M,bhS,bhS1,bhS2
    
      bhG(1)=zero
      bhC(1)=zero
      do lm=2,lm_max
         bhG(lm)=Slm(lm)-ci*Tlm(lm)
         bhC(lm)=Slm(lm)+ci*Tlm(lm)
      end do

      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc) &
      !!$omp private(dm, lm, lmS, bhN, bhN1, bhN2, bhS, bhS1, bhS2)      &
      !!$omp private(brES, brEA, bhN1M, bhN2M, bhS1M, bhS2M)
      nThStart=1; nThStop=n_theta_max/2
      !call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThStart-1 ; nThStop=2*nThStop
      !--- Loop over all orders m: (numbered by mc)
      do mc=1,n_m_max
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
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
            tmpr(nThetaN,mc)=brES+brEA
            tmpr(nThetaS,mc)=brES-brEA
 
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
            bhN1M=half*bhN1
            bhS1M=half*bhS1
            bhN2M=half*bhN2
            bhS2M=half*bhS2
 
            !--- 6 add/mult, 20 dble words
            tmpt(nThetaN,mc)=bhN1M+bhN2M
            bhN             =bhN1M-bhN2M
            tmpt(nThetaS,mc)=bhS1M+bhS2M
            bhS             =bhS1M-bhS2M
            tmpp(nThetaN,mc)=-ci*bhN
            tmpp(nThetaS,mc)=-ci*bhS
         end do      ! End global loop over nTheta
      end do
 
      !-- Zero out terms with index mc > n_m_max:
      if ( n_m_max < n_phi_max/2+1 ) then
         do mc=n_m_max+1,n_phi_max/2+1
            do nThetaN=nThStart,nThStop
               tmpr(nThetaN,mc)=zero
               tmpt(nThetaN,mc)=zero
               tmpp(nThetaN,mc)=zero
            end do  ! loop over nThetaN (theta)
         end do
      end if
      !!$omp end parallel

      if ( .not. l_axi ) then
         call ifft_many(tmpr,brc)
         call ifft_many(tmpt,btc)
         call ifft_many(tmpp,bpc)
      end if
    
   end subroutine native_qst_to_spat
!------------------------------------------------------------------------------
   subroutine native_sphtor_to_spat(Slm, Tlm, btc, bpc)
      !
      ! Use spheroidal and toroidal potentials to transform them to angular
      ! vector components btheta and bphi
      !
      
      !-- Input variables:
      complex(cp), intent(in) :: Slm(lm_max)
      complex(cp), intent(in) :: Tlm(lm_max)
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: btc(n_theta_max,n_phi_max)
      real(cp), intent(out) :: bpc(n_theta_max,n_phi_max)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(lm_max),PlmC(lm_max)
      complex(cp) :: bhG(lm_max),bhC(lm_max)
      complex(cp) :: tmpt(n_theta_max,n_phi_max/2+1)
      complex(cp) :: tmpp(n_theta_max,n_phi_max/2+1)
    
      !-- Local variables:
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: mc,lm,lmS,nThStart,nThStop
      real(cp) :: dm
    
      complex(cp) :: bhN1M,bhN2M,bhN,bhN1,bhN2
      complex(cp) :: bhS1M,bhS2M,bhS,bhS1,bhS2
    
      do lm=1,lm_max
         bhG(lm)=Slm(lm)-ci*Tlm(lm)
         bhC(lm)=Slm(lm)+ci*Tlm(lm)
      end do

      !--- Loop over all orders m: (numbered by mc)
      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc) &
      !!$omp private(dm, lm, lmS, bhN, bhN1, bhN2, bhS, bhS1, bhS2)      &
      !!$omp private(bhN1M, bhN2M, bhS1M, bhS2M)
      nThStart=1; nThStop=n_theta_max/2
      !call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThstart-1 ; nThStop=2*nThStop
      do mc=1,n_m_max
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
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
            bhN1M=half*bhN1
            bhS1M=half*bhS1
            bhN2M=half*bhN2
            bhS2M=half*bhS2
 
            !--- 6 add/mult, 20 dble words
            tmpt(nThetaN,mc)=bhN1M+bhN2M
            bhN             =bhN1M-bhN2M
            tmpt(nThetaS,mc)=bhS1M+bhS2M
            bhS             =bhS1M-bhS2M
            tmpp(nThetaN,mc)=-ci*bhN
            tmpp(nThetaS,mc)=-ci*bhS
         end do
      end do      ! End global loop over mc
 
      !-- Zero out terms with index mc > n_m_max:
      if ( n_m_max < n_phi_max/2+1 ) then
         do mc=n_m_max+1,n_phi_max/2+1
            do nThetaN=nThStart,nThStop
               tmpt(nThetaN,mc)=zero
               tmpp(nThetaN,mc)=zero
            end do  ! loop over nThetaN (theta)
         end do
      end if
      !!$omp end parallel

      if ( .not. l_axi ) then
         call ifft_many(tmpt,btc)
         call ifft_many(tmpp,bpc)
      end if
    
   end subroutine native_sphtor_to_spat
!------------------------------------------------------------------------------
   subroutine native_toraxi_to_spat(Tlm, btc, bpc)
      !
      ! Use spheroidal and toroidal potentials to transform them to angular
      ! vector components btheta and bphi
      !
      
      !-- Input variables:
      complex(cp), intent(in) :: Tlm(l_max+1)
    
      !-- Output: field on grid (theta,m) for the radial grid point nR
      !           and equatorially symmetric and antisymmetric contribution
      real(cp), intent(out) :: btc(n_theta_max)
      real(cp), intent(out) :: bpc(n_theta_max)
    
      !------ Legendre Polynomials 
      real(cp) :: PlmG(l_max+1),PlmC(l_max+1)
      complex(cp) :: bhG(l_max+1),bhC(l_max+1)
    
      !-- Local variables:
      integer :: nThetaN,nThetaS,nThetaNHS
      integer :: lm,lmS,l
    
      complex(cp) :: bhN1M,bhN2M,bhN,bhN1,bhN2
      complex(cp) :: bhS1M,bhS2M,bhS,bhS1,bhS2
    
      do l=1,l_max+1
         bhG(l)=-ci*Tlm(l)
         bhC(l)= ci*Tlm(l)
      end do

      nThetaNHS=0
      do nThetaN=1,n_theta_max,2   ! Loop over thetas for north HS
         nThetaS  =nThetaN+1      ! same theta but for southern HS
         nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
 
         lmS=lStop(1)
         !--- 6 add/mult, 26 dble words
         do lm=lStart(1),lmS-1,2
            PlmG(lm)=dPlm(lm,nThetaNHS)
            PlmC(lm)=dPlm(lm,nThetaNHS)
            PlmG(lm+1)=dPlm(lm+1,nThetaNHS)
            PlmC(lm+1)=dPlm(lm+1,nThetaNHS)
         end do
         if ( lmOdd(1) ) then
            PlmG(lmS)=dPlm(lmS,nThetaNHS)
            PlmC(lmS)=dPlm(lmS,nThetaNHS)
         end if
 
         lmS=lStop(1)
         bhN1=zero
         bhS1=zero
         bhN2=zero
         bhS2=zero
         !--- 8 add/mult, 20 dble words
         do lm=lStart(1),lmS-1,2
            bhN1=bhN1+bhG(lm)*PlmG(lm)+bhG(lm+1)*PlmG(lm+1)
            bhS1=bhS1-bhG(lm)*PlmC(lm)+bhG(lm+1)*PlmC(lm+1)
            bhN2=bhN2+bhC(lm)*PlmC(lm)+bhC(lm+1)*PlmC(lm+1)
            bhS2=bhS2-bhC(lm)*PlmG(lm)+bhC(lm+1)*PlmG(lm+1)
         end do
         if ( lmOdd(1) ) then
            bhN1=bhN1+bhG(lmS)*PlmG(lmS)
            bhS1=bhS1-bhG(lmS)*PlmC(lmS)
            bhN2=bhN2+bhC(lmS)*PlmC(lmS)
            bhS2=bhS2-bhC(lmS)*PlmG(lmS)
         end if
         bhN1M=half*bhN1
         bhS1M=half*bhS1
         bhN2M=half*bhN2
         bhS2M=half*bhS2

         !--- 6 add/mult, 20 dble words
         btc(nThetaN)=real(bhN1M+bhN2M)
         bhN         =bhN1M-bhN2M
         btc(nThetaS)=real(bhS1M+bhS2M)
         bhS         =bhS1M-bhS2M
         bpc(nThetaN)=real(-ci*bhN)
         bpc(nThetaS)=real(-ci*bhS)
      end do
 
   end subroutine native_toraxi_to_spat
!------------------------------------------------------------------------------
   subroutine native_sph_to_spat(Slm, sc)
      !
      ! Spherical Harmonic Transform for a scalar input field
      !

      !-- Input variable
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: sc(n_theta_max,n_phi_max)

      !-- Local variables
      complex(cp) :: tmp(n_theta_max,n_phi_max/2+1)
      integer :: mc, lm, lmS
      integer :: nThetaN, nThetaS, nThetaNHS, nThStart, nThStop
      complex(cp) :: sES, sEA

      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc) &
      !!$omp private(sES, sEA, lm, lmS)
      nThStart=1; nThStop=n_theta_max/2
      call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThstart-1 ; nThStop=2*nThStop
      do mc=1,n_m_max
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2   ! Loop over thetas for one HS
            nThetaS  =nThetaN+1  ! same theta but at other HS
            nThetaNHS=nThetaNHS+1  ! ic-index of northern hemisph. point
            lmS=lStop(mc)
            sES=zero ! One equatorial symmetry
            sEA=zero ! The other equatorial symmetry
            do lm=lStart(mc),lmS-1,2
               sES=sES+Slm(lm)  *Plm(lm,nThetaNHS)
               sEA=sEA+Slm(lm+1)*Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(mc) ) sES=sES+Slm(lmS)*Plm(lmS,nThetaNHS)
            tmp(nThetaN,mc)=sES+sEA
            tmp(nThetaS,mc)=sES-sEA
         end do
      end do

      if ( n_m_max < n_phi_max/2+1 ) then
         do mc=n_m_max+1,n_phi_max/2+1
            do nThetaN=nThstart,nThStop
               tmp(nThetaN,mc)=zero
            end do ! loop over nThetaN (theta)
         end do  
      end if
      !!$omp end parallel

      if ( .not. l_axi ) call ifft_many(tmp,sc)

   end subroutine native_sph_to_spat
!------------------------------------------------------------------------------
   subroutine native_sph_to_grad_spat(slm, gradtc, gradpc)
      !
      ! Transform s(l) into dsdt(theta) and dsdp(theta)
      !

      !-- Input variable
      complex(cp), intent(in) :: Slm(lm_max)

      !-- Output variables
      real(cp), intent(out) :: gradtc(n_theta_max,n_phi_max)
      real(cp), intent(out) :: gradpc(n_theta_max,n_phi_max)

      !-- Local variables
      complex(cp) :: tmpt(n_theta_max,n_phi_max/2+1),tmpp(n_theta_max,n_phi_max/2+1)
      integer :: mc, lm, lmS, nThetaN, nThetaS, nThetaNHS, nThStart, nThStop
      real(cp) :: dm
      complex(cp) :: gradtcES, gradtcEA, sES, sEA

      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc) &
      !!$omp private(sES, sEA, gradtcES, gradtcEA, dm, lm, lmS)
      nThStart=1; nThStop=n_theta_max/2
      call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThstart-1 ; nThStop=2*nThStop
      do mc=1,n_m_max
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2   ! Loop over thetas for north HS
            nThetaS  =nThetaN+1      ! same theta but for southern HS
            nThetaNHS=nThetaNHS+1    ! theta-index of northern hemisph. point
            dm =D_mc2m(mc)
            lmS=lStop(mc)
            sES=zero  ! One equatorial symmetry
            sEA=zero  ! The other equatorial symmetry
            gradtcES=zero
            gradtcEA=zero
            do lm=lStart(mc),lmS-1,2
               gradtcEA=gradtcEA + Slm(lm)*  dPlm(lm,nThetaNHS)
               gradtcES=gradtcES + Slm(lm+1)*dPlm(lm+1,nThetaNHS)
               sES     =sES      + Slm(lm)  * Plm(lm,nThetaNHS)
               sEA     =sEA      + Slm(lm+1)* Plm(lm+1,nThetaNHS)
            end do
            if ( lmOdd(mc) ) then
               gradtcEA=gradtcEA + Slm(lmS)*dPlm(lmS,nThetaNHS)
               sES     =sES      + Slm(lmS)* Plm(lmS,nThetaNHS)
            end if
            tmpt(nThetaN,mc)=gradtcES+gradtcEA
            tmpt(nThetaS,mc)=gradtcES-gradtcEA
            tmpp(nThetaN,mc)=dm*ci*(sES+sEA)
            tmpp(nThetaS,mc)=dm*ci*(sES-sEA)
         end do
      end do

      if ( n_m_max < n_phi_max/2+1 ) then
         do mc=n_m_max+1,n_phi_max/2+1
            do nThetaN=nThStart,nThStop
               tmpt(nThetaN,mc)=zero
               tmpp(nThetaN,mc)=zero
            end do  ! loop over nThetaN (theta)
         end do
      end if
      !!$omp end parallel

      if ( .not. l_axi ) then
         call ifft_many(tmpt,gradtc)
         call ifft_many(tmpp,gradpc)
      end if

   end subroutine native_sph_to_grad_spat
!------------------------------------------------------------------------------
   subroutine native_spat_to_sph(scal,f1LM)
      !
      !  Legendre transform (n_r,n_theta,m) to (n_r,l,m)
      !  [grid to spectral] for 2 arrays
      !  f1TM (input) to f1LM (output)
      !  One call to this routine does part of the transform
      !  by summation over theta points in one theta block:
      !  nThetaStart,..,nThetaStart+n_theta_block-1
      !

      !-- Input variables:
      real(cp), intent(inout) :: scal(:,:)

      !-- Output variable:
      complex(cp), intent(out) :: f1LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)

      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode
      integer :: nThStart,nThStop

      complex(cp) :: f1TM(n_theta_max,n_phi_max/2+1)
      complex(cp) :: f1ES(n_theta_max/2,n_m_max),f1ES1,f1ES2
      complex(cp) :: f1EA(n_theta_max/2,n_m_max),f1EA1,f1EA2

      if ( .not. l_axi ) call fft_many(scal,f1TM)

      f1LM(:)=zero

      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc) &
      !!$omp private(lm, lmS, nTheta1, nTheta2, f1ES1, f1ES2, f1EA1, f1EA2)
      nThStart=1; nThStop=n_theta_max/2
      call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThstart-1 ; nThStop=2*nThStop
      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      do mc=1,n_m_max        ! counts spherical harmonic orders
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2 ! thetas in NHS
            nThetaS=nThetaN+1      ! thetas in SHS
            nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
            f1ES(nThetaNHS,mc)=f1TM(nThetaN,mc)+f1TM(nThetaS,mc)
            f1EA(nThetaNHS,mc)=f1TM(nThetaN,mc)-f1TM(nThetaS,mc)
         end do

         !-- Loop over half of the thetas with step 2 unrolling:
         do nTheta1=(nThStart+1)/2,nThStop/2,2
            nTheta2=nTheta1+1
            lmS=lStopP(mc)
            f1ES1=f1ES(nTheta1,mc)
            f1ES2=f1ES(nTheta2,mc)
            f1EA1=f1EA(nTheta1,mc)
            f1EA2=f1EA(nTheta2,mc)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)  +f1ES1*wPlm(lm,nTheta1)  +f1ES2*wPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1)+f1EA1*wPlm(lm+1,nTheta1)+f1EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + f1ES2*wPlm(lmS,nTheta2)
            end if
         end do
      end do 
      !!$omp end parallel

   end subroutine native_spat_to_sph
!------------------------------------------------------------------------------
   subroutine native_spat_to_sph_tor(vt,vp,f1LM,f2LM)
      !
      !  Vector Legendre transform
      !  vt(n_r,n_theta,m), vp(n_r,n_theta,m) to Spheroidal(n_r,l,m)
      !  and Toroidal(n_r,l,m)

      !-- Input variables:
      real(cp), intent(inout) :: vt(:,:)
      real(cp), intent(inout) :: vp(:,:)

      !-- Output variables:
      complex(cp), intent(out) :: f1LM(lmP_max),f2LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)

      integer :: l, nThStart, nThStop
      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode

      complex(cp) :: f1TM(n_theta_max,n_phi_max/2+1)
      complex(cp) :: f2TM(n_theta_max,n_phi_max/2+1)
      complex(cp) :: f1ES(n_theta_max/2,n_m_max),f1ES1,f1ES2
      complex(cp) :: f1EA(n_theta_max/2,n_m_max),f1EA1,f1EA2
      complex(cp) :: f2ES(n_theta_max/2,n_m_max),f2ES1,f2ES2
      complex(cp) :: f2EA(n_theta_max/2,n_m_max),f2EA1,f2EA2

      real(cp) :: dm

      !-- Not sure why there is a crossing between f2TM and f1TM here
      if ( .not. l_axi ) then
         call fft_many(vt,f2TM)
         call fft_many(vp,f1TM)
      end if

      f1LM(:)=zero
      f2LM(:)=zero

      !!$omp parallel default(shared) &
      !!$omp private(nThStart, nThStop, nThetaNHS, nThetaN, nThetaS, mc)    &
      !!$omp private(lm, lmS, nTheta1, nTheta2, f1ES1, f1ES2, f1EA1, f1EA2) &
      !!$omp private(f2ES1, f2ES2, f2EA1, f2EA2)
      nThStart=1; nThStop=n_theta_max/2
      call get_openmp_blocks(nThStart,nThStop)
      nThStart=2*nThstart-1 ; nThStop=2*nThStop

      !-- SHTns MagIC layout adopts the following convention
      !-- One has to multiplty by 1/sin(theta)**2 here
      do mc=1,n_phi_max/2+1
         do nThetaN=nThStart,nThStop
            f1TM(nThetaN,mc)=f1TM(nThetaN,mc)*O_sin_theta_E2(nThetaN)
            f2TM(nThetaN,mc)=f2TM(nThetaN,mc)*O_sin_theta_E2(nThetaN)
         end do
      end do

      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      do mc=1,n_m_max         ! counts spherical harmonic orders
         nThetaNHS=(nThStart-1)/2
         do nThetaN=nThStart,nThStop,2 ! thetas in NHS
            nThetaS=nThetaN+1       ! thetas in SHS
            nThetaNHS=nThetaNHS+1   ! thetas counted in NHS only
            f1ES(nThetaNHS,mc)=f1TM(nThetaN,mc)+f1TM(nThetaS,mc)
            f1EA(nThetaNHS,mc)=f1TM(nThetaN,mc)-f1TM(nThetaS,mc)
            f2ES(nThetaNHS,mc)=f2TM(nThetaN,mc)+f2TM(nThetaS,mc)
            f2EA(nThetaNHS,mc)=f2TM(nThetaN,mc)-f2TM(nThetaS,mc)
         end do

         !-- Loop over half of the thetas with step 2 unrolling:
         do nTheta1=(nThStart+1)/2,nThStop/2,2
            nTheta2 =nTheta1+1

            !dm = real((mc-1)*minc,cp)
            dm = D_mc2m(mc)
            lmS=lStopP(mc)
            f1ES1=f1ES(nTheta1,mc)
            f1ES2=f1ES(nTheta2,mc)
            f1EA1=f1EA(nTheta1,mc)
            f1EA2=f1EA(nTheta2,mc)
            f2ES1=f2ES(nTheta1,mc)
            f2ES2=f2ES(nTheta2,mc)
            f2EA1=f2EA(nTheta1,mc)
            f2EA2=f2EA(nTheta2,mc)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)-ci*dm*f1ES1* wPlm(lm,nTheta1)      &
               &                  -ci*dm*f1ES2* wPlm(lm,nTheta2)      &
               &                  +      f2EA1*wdPlm(lm,nTheta1)      &
               &                  +      f2EA2*wdPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1)-ci*dm*f1EA1* wPlm(lm+1,nTheta1)  &
               &                    -ci*dm*f1EA2* wPlm(lm+1,nTheta2)  &
               &                    +      f2ES1*wdPlm(lm+1,nTheta1)  &
               &                    +      f2ES2*wdPlm(lm+1,nTheta2)
               f2LM(lm)  =f2LM(lm)-      f1EA1*wdPlm(lm,nTheta1) &
               &                  -      f1EA2*wdPlm(lm,nTheta2) &
               &                  -ci*dm*f2ES1* wPlm(lm,nTheta1) &
               &                  -ci*dm*f2ES2* wPlm(lm,nTheta2)
               f2LM(lm+1)=f2LM(lm+1)-      f1ES1*wdPlm(lm+1,nTheta1) &
               &                    -      f1ES2*wdPlm(lm+1,nTheta2) &
               &                    -ci*dm*f2EA1* wPlm(lm+1,nTheta1) &
               &                    -ci*dm*f2EA2* wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS)-ci*dm*f1ES1* wPlm(lmS,nTheta1) &
               &                  -ci*dm*f1ES2* wPlm(lmS,nTheta2) &
               &                  +      f2EA1*wdPlm(lmS,nTheta1) &
               &                  +      f2EA2*wdPlm(lmS,nTheta2)
               f2LM(lmS)=f2LM(lmS)-      f1EA1*wdPlm(lmS,nTheta1) &
               &                  -      f1EA2*wdPlm(lmS,nTheta2) &
               &                  -ci*dm*f2ES1* wPlm(lmS,nTheta1) &
               &                  -ci*dm*f2ES2* wPlm(lmS,nTheta2)
            end if
         end do !  loop over theta in block

      end do 
      !!$omp end parallel

      !-- Division by l(l+1) except for (l=0,m=0)
      do lm=2,lmP_max
         l=lmP2l(lm)
         f1LM(lm)=f1LM(lm)/real(l*(l+1),cp)
         f2LM(lm)=f2LM(lm)/real(l*(l+1),cp)
      end do

   end subroutine native_spat_to_sph_tor
!------------------------------------------------------------------------------
   subroutine native_spat_to_SH_axi(ft1,flm1,lmMax)
      !
      !  Legendre transform (n_r,n_theta,m) to (n_r,l,m)
      !  [grid to spectral] for 2 arrays
      !  ancl1/2 (input) to flm1/2 (output)
      !  One call to this routine does part of the transform
      !  by summation over theta points in on theta block:
      !  n_theta_min,..,n_theta_min+n_theta_block-1
      !

      !-- Input variables:
      integer,  intent(in) :: lmMax          ! Number of modes to be processed
      real(cp), intent(in) :: ft1(*)

      !-- Output: transformed arrays anlc1,anlc2
      real(cp), intent(out) :: flm1(*)

      !-- Local variables:
      integer :: nThetaN,nThetaS,nThetaNHS,nTheta1,nTheta2,lm1,lm2
      real(cp) :: f1p(n_theta_max/2),f1m(n_theta_max/2)
      
      flm1(1:lmMax)=0.0_cp

      !-- Prepare arrays of sums and differences:
      nThetaNHS=0
      do nThetaN=1,n_theta_max,2 ! thetas in NHS
         nThetaS=nThetaN+1         ! thetas in SHS
         nThetaNHS=nThetaNHS+1      ! thetas in one HS
         f1p(nThetaNHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
         f1m(nThetaNHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
      end do

      do nTheta1=1,n_theta_max/2,2
         nTheta2=nTheta1+1

         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nTheta1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nTheta2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nTheta1) + &
            &                     f1m(nTheta2)*wPlm(lm2,nTheta2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nTheta1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nTheta2)
         end if
      end do

   end subroutine native_spat_to_SH_axi
!------------------------------------------------------------------------------
end module shtransforms