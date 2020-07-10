module legendre_grid_to_spec

   use precision_mod
   use truncation, only: n_m_max, nrp, lmP_max, n_theta_max, n_phi_max, l_axi, minc
   use truncation, only: l_max, m_max
   use LMmapping, only: map_glbl_st
   use blocking, only: nfs, sizeThetaB
   use horizontal_data, only: lStartP, wPlm, lmOddP, lStopP, wdPlm, D_mc2m, &
       &                      O_sin_theta_E2
   use constants, only: ci, zero
   use fft, only: fft_many

   implicit none

   private

   public :: spat_to_sph, legTFAS, legTFAS2, spat_to_sph_tor

contains

   subroutine spat_to_sph(scal,f1LM)
      !
      !  Legendre transform (n_r,n_theta,m) to (n_r,l,m)
      !  [grid to spectral] for 2 arrays
      !  f1TM (input) to f1LM (output)
      !  One call to this routine does part of the transform
      !  by summation over theta points in one theta block:
      !  nThetaStart,..,nThetaStart+n_theta_block-1
      !

      !-- Input variables:
      real(cp), intent(inout) :: scal(n_phi_max,n_theta_max)

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

      complex(cp) :: f1TM(n_phi_max/2+1,n_theta_max)
      complex(cp) :: f1ES(n_m_max,n_theta_max/2),f1ES1,f1ES2
      complex(cp) :: f1EA(n_m_max,n_theta_max/2),f1EA1,f1EA2

      if ( .not. l_axi ) call fft_many(scal,f1TM)

      f1LM(:)=zero

      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      nThetaNHS=0
      do nThetaN=1,n_theta_max,2 ! thetas in NHS
         nThetaS=nThetaN+1      ! thetas in SHS
         nThetaNHS=nThetaNHS+1  ! thetas counted in NHS only
         do mc=1,n_m_max        ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=f1TM(mc,nThetaN)+f1TM(mc,nThetaS)
            f1EA(mc,nThetaNHS)=f1TM(mc,nThetaN)-f1TM(mc,nThetaS)
         end do
      end do

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=1,n_theta_max/2,2
         nTheta2=nTheta1+1
         do mc=1,n_m_max
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nTheta1)
            f1ES2=f1ES(mc,nTheta2)
            f1EA1=f1EA(mc,nTheta1)
            f1EA2=f1EA(mc,nTheta2)
            do lm=lStartP(mc),lmS-1,2
               f1LM(lm)  =f1LM(lm)  +f1ES1*wPlm(lm,nTheta1)  +f1ES2*wPlm(lm,nTheta2)
               f1LM(lm+1)=f1LM(lm+1)+f1EA1*wPlm(lm+1,nTheta1)+f1EA2*wPlm(lm+1,nTheta2)
            end do
            if ( lmOddP(mc) ) then
               f1LM(lmS)=f1LM(lmS) + f1ES1*wPlm(lmS,nTheta1) + f1ES2*wPlm(lmS,nTheta2)
            end if
        end do
      end do 

   end subroutine spat_to_sph
!------------------------------------------------------------------------
   subroutine spat_to_sph_tor(vt,vp,f1LM,f2LM)
      !
      !  Vector Legendre transform
      !  vt(n_r,n_theta,m), vp(n_r,n_theta,m) to Spheroidal(n_r,l,m)
      !  and Toroidal(n_r,l,m)
      !

      !-- Input variables:
      real(cp), intent(inout) :: vt(n_phi_max,n_theta_max)
      real(cp), intent(inout) :: vp(n_phi_max,n_theta_max)

      !-- Output variables:
      complex(cp), intent(out) :: f1LM(lmP_max),f2LM(lmP_max)

      !-- Local variables:
      integer :: nThetaN     ! No. of theta in NHS
      integer :: nThetaS     ! No. of theta in SHS
      integer :: nThetaNHS   ! No. of thetas in one HS only
      integer :: nTheta1     ! No. of theta (in one HS)
      integer :: nTheta2     ! No. of theta (in one HS)

      integer :: l, m
      integer :: mc          ! counter of spherical order
      integer :: lmS,lm      ! counter of spherical mode

      complex(cp) :: f1TM(n_phi_max/2+1,n_theta_max)
      complex(cp) :: f2TM(n_phi_max/2+1,n_theta_max)
      complex(cp) :: f1ES(n_m_max,n_theta_max/2),f1ES1,f1ES2
      complex(cp) :: f1EA(n_m_max,n_theta_max/2),f1EA1,f1EA2
      complex(cp) :: f2ES(n_m_max,n_theta_max/2),f2ES1,f2ES2
      complex(cp) :: f2EA(n_m_max,n_theta_max/2),f2EA1,f2EA2

      real(cp) :: dm

      if ( .not. l_axi ) then
         call fft_many(vt,f2TM)
         call fft_many(vp,f1TM)
      end if

      do nThetaN=1,n_theta_max
         f1TM(:,nThetaN)=f1TM(:,nThetaN)*O_sin_theta_E2(nThetaN)
         f2TM(:,nThetaN)=f2TM(:,nThetaN)*O_sin_theta_E2(nThetaN)
      end do

      f1LM(:)=zero
      f2LM(:)=zero

      !-- Unscrambles equatorially symmetric and antisymmetric contributions:
      nThetaNHS=0
      do nThetaN=1,n_theta_max,2 ! thetas in NHS
         nThetaS=nThetaN+1       ! thetas in SHS
         nThetaNHS=nThetaNHS+1   ! thetas counted in NHS only
         do mc=1,n_m_max         ! counts spherical harmonic orders
            f1ES(mc,nThetaNHS)=(f1TM(mc,nThetaN)+f1TM(mc,nThetaS))
            f1EA(mc,nThetaNHS)=(f1TM(mc,nThetaN)-f1TM(mc,nThetaS))
            f2ES(mc,nThetaNHS)=(f2TM(mc,nThetaN)+f2TM(mc,nThetaS))
            f2EA(mc,nThetaNHS)=(f2TM(mc,nThetaN)-f2TM(mc,nThetaS))
         end do
      end do

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=1,n_theta_max/2,2
         nTheta2 =nTheta1+1

         do mc=1,n_m_max
            !dm = real((mc-1)*minc,cp)
            dm = D_mc2m(mc)
            lmS=lStopP(mc)
            f1ES1=f1ES(mc,nTheta1)
            f1ES2=f1ES(mc,nTheta2)
            f1EA1=f1EA(mc,nTheta1)
            f1EA2=f1EA(mc,nTheta2)
            f2ES1=f2ES(mc,nTheta1)
            f2ES2=f2ES(mc,nTheta2)
            f2EA1=f2EA(mc,nTheta1)
            f2EA2=f2EA(mc,nTheta2)
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
         end do

      end do  !  loop over theta in block

      do lm=2,lmP_max
         l=map_glbl_st%lmP2l(lm)
         f1LM(lm)=f1LM(lm)/real(l*(l+1),cp)
         f2LM(lm)=f2LM(lm)/real(l*(l+1),cp)
      end do

   end subroutine spat_to_sph_tor
!------------------------------------------------------------------------
   subroutine legTFAS(flm1,ft1,lmMax,nThetaStart,sizeThetaB)
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
      integer,  intent(in) :: nThetaStart    ! First no of theta on block
      integer,  intent(in) :: sizeThetaB     ! Size of theta block
      real(cp), intent(in) :: ft1(*)

      !-- Output: transformed arrays anlc1,anlc2
      real(cp), intent(out) :: flm1(*)

      !-- Local variables:
      integer :: nThetaN,nThetaS
      integer :: nThetaHS
      integer :: nTheta1,nTheta2
      integer :: nThetaC1,nThetaC2
      integer :: nThetaMin
      integer :: lm1,lm2

      real(cp) :: f1p(n_theta_max/2),f1m(n_theta_max/2)

      !-- Prepare arrays of sums and differences:
      nThetaHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS =nThetaN+1         ! thetas in SHS
         nThetaHS=nThetaHS+1        ! thetas in one HS
         f1p(nThetaHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
         f1m(nThetaHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
      end do

      !-- Start with first two thetas in first theta block:
      !   This initalizes the flm1 and flm2
      if ( nThetaStart == 1 ) then
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm1(lm2)=f1m(1)*wPlm(lm2,1)+f1m(2)*wPlm(lm2,2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
         end if

         if ( sizeThetaB <= 4 ) return

      end if

      !-- Following values of n_theta=3,4,5... when n_theta_min=1
      !   or all values of n_theta when n_theta_min > 0

      nThetaMin=1
      If ( nThetaStart == 1 ) nThetaMin=3 ! 2 already done

      !-- Calculating position for all thetas:
      !     (nThetaStart-1)/2 are previous blocks,
      !     (nThetaMin-1) is what has been done already
      nThetaC1=(nThetaStart-1)/2+nThetaMin-2

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=nThetaMin,sizeThetaB/2,2
         nTheta2=nTheta1+1
         nThetaC1=nThetaC1+2
         nThetaC2=nThetaC1+1

         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nThetaC1) + &
            &                     f1m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nThetaC2)
         end if

      end do


   end subroutine legTFAS
!------------------------------------------------------------------------
   subroutine legTFAS2(flm1,flm2,ft1,ft2,lmMax,nThetaStart,sizeThetaB)
      !
      !  Legendre transform (n_r,n_theta,m) to (n_r,l,m)
      !  [grid to spectral] for 2 arrays
      !  ancl1/2 (input) to flm1/2 (output)
      !  One call to this routine does part of the transform
      !  by summation over theta points in on theta block:
      !  n_theta_min,..,n_theta_min+n_theta_block-1
      !

      !-- Input variables
      integer,  intent(in) :: lmMax          ! Number of modes to be processed
      integer,  intent(in) :: nThetaStart    ! First no of theta on block
      integer,  intent(in) :: sizeThetaB     ! Size of theta block
      real(cp), intent(in) :: ft1(*),ft2(*)

      !-- Output: transformed arrays anlc1,anlc2
      real(cp), intent(out) :: flm1(*),flm2(*)

      !-- Local variables:
      integer :: nThetaN,nThetaS
      integer :: nThetaHS
      integer :: nTheta1,nTheta2
      integer :: nThetaC1,nThetaC2
      integer :: nThetaMin
      integer :: lm1,lm2

      real(cp) :: f1p(n_theta_max/2),f1m(n_theta_max/2)
      real(cp) :: f2p(n_theta_max/2),f2m(n_theta_max/2)

      !-- Prepare arrays of sums and differences:
      nThetaHS=0
      do nThetaN=1,sizeThetaB,2 ! thetas in NHS
         nThetaS =nThetaN+1         ! thetas in SHS
         nThetaHS=nThetaHS+1        ! thetas in one HS
         f1p(nThetaHS)=ft1(nThetaN)+ft1(nThetaS) ! Symm
         f1m(nThetaHS)=ft1(nThetaN)-ft1(nThetaS) ! ASymm
         f2p(nThetaHS)=ft2(nThetaN)+ft2(nThetaS)
         f2m(nThetaHS)=ft2(nThetaN)-ft2(nThetaS)
      end do

      !-- Start with first two thetas in first theta block:
      !   This initalizes the flm1 and flm2
      if ( nThetaStart == 1 ) then
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm1(lm2)=f1m(1)*wPlm(lm2,1)+f1m(2)*wPlm(lm2,2)
         end do
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
            flm2(lm2)=f2m(1)*wPlm(lm2,1)+f2m(2)*wPlm(lm2,2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=f1p(1)*wPlm(lm1,1)+f1p(2)*wPlm(lm1,2)
            flm2(lm1)=f2p(1)*wPlm(lm1,1)+f2p(2)*wPlm(lm1,2)
         end if

         if ( sizeThetaB <= 4 ) return

      end if

      !-- Following values of n_theta=3,4,5... when n_theta_min=1
      !   or all values of n_theta when n_theta_min > 0
      nThetaMin=1
      if ( nThetaStart == 1 ) nThetaMin=3 ! 2 already done

      !-- Calculating position for all thetas:
      !     (nThetaStart-1)/2 are previous blocks,
      !     (nThetaMin-1) is what has been done already
      nThetaC1=(nThetaStart-1)/2+nThetaMin-2

      !-- Loop over half of the thetas with step 2 unrolling:
      do nTheta1=nThetaMin,sizeThetaB/2,2
         nTheta2=nTheta1+1
         nThetaC1=nThetaC1+2
         nThetaC2=nThetaC1+1

         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm1(lm2)=flm1(lm2) + f1m(nTheta1)*wPlm(lm2,nThetaC1) + &
            &                     f1m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         do lm1=1,lmMax-1,2
            lm2=lm1+1
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f2p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm2)=flm2(lm2) + f2m(nTheta1)*wPlm(lm2,nThetaC1) + &
            &                     f2m(nTheta2)*wPlm(lm2,nThetaC2)
         end do
         if ( lm2 < lmMax ) then
            lm1=lmMax
            flm1(lm1)=flm1(lm1) + f1p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f1p(nTheta2)*wPlm(lm1,nThetaC2)
            flm2(lm1)=flm2(lm1) + f2p(nTheta1)*wPlm(lm1,nThetaC1) + &
            &                     f2p(nTheta2)*wPlm(lm1,nThetaC2)
         end if

      end do

   end subroutine legTFAS2
!------------------------------------------------------------------------
end module legendre_grid_to_spec
