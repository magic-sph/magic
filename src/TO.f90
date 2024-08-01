module torsional_oscillations
   !
   !  This module contains information for TO calculation and output
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_phi_max, n_r_max, l_max, n_theta_max, lm_max, &
       &                 nlat_padded
   use grid_blocking, only: radlatlon2spat
   use radial_data, only: nRstart, nRstop
   use radial_functions, only: or1, or2, or3, or4, beta, orho1, dbeta
   use physical_parameters, only: CorFac, kbotv, ktopv, epsPhase, penaltyFac
   use blocking, only: lm2
   use horizontal_data, only: sinTheta, cosTheta, hdif_V, dLh, &
       &                      n_theta_cal2ord, O_sin_theta
   use constants, only: one, two
   use logic, only: lVerbose, l_mag, l_parallel_solve, l_phase_field
   use sht, only: toraxi_to_spat

   implicit none

   private

   real(cp), allocatable :: zASL(:), dzASL(:)
   real(cp), public, allocatable :: ddzASL(:,:)
   real(cp), allocatable :: dzddVpLMr(:,:), dzdVpLMr(:,:)

   real(cp), public, allocatable :: dzStrAS_Rloc(:,:),dzRstrAS_Rloc(:,:)
   real(cp), public, allocatable :: dzAStrAS_Rloc(:,:),dzCorAS_Rloc(:,:)
   real(cp), public, allocatable :: dzLFAS_Rloc(:,:),dzdVpAS_Rloc(:,:)
   real(cp), public, allocatable :: dzddVpAS_Rloc(:,:), VAS_Rloc(:,:)
   real(cp), public, allocatable :: V2AS_Rloc(:,:),Bs2AS_Rloc(:,:)
   real(cp), public, allocatable :: BszAS_Rloc(:,:),BspAS_Rloc(:,:)
   real(cp), public, allocatable :: BpzAS_Rloc(:,:),BspdAS_Rloc(:,:)
   real(cp), public, allocatable :: BpsdAS_Rloc(:,:),BzpdAS_Rloc(:,:)
   real(cp), public, allocatable :: BpzdAS_Rloc(:,:),dzPenAS_Rloc(:,:)

   real(cp), allocatable :: BsLast(:,:,:), BpLast(:,:,:), BzLast(:,:,:)

   public :: initialize_TO, prep_TO_axi, getTO, getTOnext, getTOfinish, finalize_TO

contains

   subroutine initialize_TO
      !
      ! Allocate the memory needed
      !

      allocate( dzStrAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( dzRstrAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( dzAStrAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( dzCorAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( dzdVpAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( dzddVpAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( VAS_Rloc(n_theta_max,nRstart:nRstop) )
      allocate( V2AS_Rloc(n_theta_max,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+8*(nRstop-nRstart+1)*n_theta_max* &
      &                 SIZEOF_DEF_REAL
      if ( l_phase_field ) then
         allocate( dzPenAS_Rloc(n_theta_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+(nRstop-nRstart+1)*n_theta_max* &
         &                 SIZEOF_DEF_REAL
      end if

      allocate( dzdVpLMr(l_max+1,nRstart:nRstop), dzddVpLMr(l_max+1,nRstart:nRstop) )
      dzdVpLMr(:,:) =0.0_cp
      dzddVpLMr(:,:)=0.0_cp
      bytes_allocated = bytes_allocated+2*(l_max+1)*(nRstop-nRstart+1)*SIZEOF_DEF_REAL

      if ( l_parallel_solve ) then
         allocate( ddzASL(l_max+1,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+(l_max+1)*(nRstop-nRstart+1)*SIZEOF_DEF_REAL
      else
         allocate( ddzASL(l_max+1,n_r_max) )
         bytes_allocated = bytes_allocated+(l_max+1)*n_r_max*SIZEOF_DEF_REAL
      end if
      allocate( zASL(l_max+1), dzASL(l_max+1) )
      bytes_allocated = bytes_allocated+2*(l_max+1)*SIZEOF_DEF_REAL
      ddzASL(:,:)=0.0_cp
      dzASL(:)   =0.0_cp
      zASL(:)    =0.0_cp

      if ( l_mag ) then
         allocate( dzLFAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( Bs2AS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BszAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BspAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BpzAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BspdAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BpsdAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BzpdAS_Rloc(n_theta_max,nRstart:nRstop) )
         allocate( BpzdAS_Rloc(n_theta_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+9*(nRstop-nRstart+1)*n_theta_max* &
         &                 SIZEOF_DEF_REAL

         allocate( BsLast(n_theta_max,n_phi_max,nRstart:nRstop) )
         allocate( BpLast(n_theta_max,n_phi_max,nRstart:nRstop) )
         allocate( BzLast(n_theta_max,n_phi_max,nRstart:nRstop) )
         bytes_allocated = bytes_allocated+3*n_phi_max*n_theta_max* &
         &                 (nRstop-nRstart+1)*SIZEOF_DEF_REAL
      end if

   end subroutine initialize_TO
!-----------------------------------------------------------------------------
   subroutine finalize_TO
      !
      ! Deallocate the memory
      !

      deallocate( dzddVpAS_Rloc, dzdVpAS_Rloc, dzCorAS_Rloc )
      if ( l_phase_field ) deallocate(dzPenAS_Rloc)
      if ( l_mag ) then
         deallocate( BsLast, BpLast, BzLast, dzLFAS_Rloc )
         deallocate( BpzdAS_Rloc, BzpdAS_Rloc, BpsdAS_Rloc, BspdAS_Rloc )
         deallocate( BpzAS_Rloc, BspAS_Rloc, BszAS_Rloc, Bs2AS_Rloc )
      end if
      deallocate( dzAStrAS_Rloc, dzRstrAS_Rloc, dzStrAS_Rloc )
      deallocate( zASL, dzASL, ddzASL, VAS_Rloc, V2AS_Rloc )
      deallocate( dzddVpLMr, dzdVpLMr )

   end subroutine finalize_TO
!-----------------------------------------------------------------------------
   subroutine prep_TO_axi(z,dz)

      !-- Input variables
      complex(cp), intent(in) :: z(:)
      complex(cp), intent(in) :: dz(:)

      !-- Local variables
      integer :: l, lm

      do l=0,l_max
         lm=lm2(l,0)
         zASL(l+1) =real(z(lm))
         dzASL(l+1)=real(dz(lm))
      end do

   end subroutine prep_TO_axi
!-----------------------------------------------------------------------------
   subroutine getTO(vr,vt,vp,cvr,dvpdr,br,bt,bp,cbr,cbt,phase,dtLast,nR)
      !
      !  This program calculates various axisymmetric linear
      !  and nonlinear variables for a radial grid point nR and
      !  a theta-block.
      !  Input are the fields vr,vt,vp,cvr,dvpdr
      !  Output are linear azimuthally averaged  field VpAS (flow phi component),
      !  VpAS2 (square of flow phi component), V2AS (V*V),
      !  and Coriolis force Cor. These are give in (r,theta)-space.
      !  Also in (r,theta)-space are azimuthally averaged correlations of
      !  non-axisymmetric flow components and the respective squares:
      !  Vsp=Vs*Vp,Vzp,Vsz,VspC,VzpC,VszC. These are used to calulcate
      !  the respective correlations and Reynolds stress.
      !

      !-- Input of variables
      real(cp), intent(in) :: dtLast      ! last time step
      integer,  intent(in) :: nR          ! radial grid point
      real(cp), intent(in) :: vr(*),vt(*),vp(*)
      real(cp), intent(in) :: cvr(*),dvpdr(*)
      real(cp), intent(in) :: br(*),bt(*),bp(*)
      real(cp), intent(in) :: cbr(*),cbt(*),phase(*)

      !-- Local variables:
      integer :: nTheta,nPhi,nTheta1,nelem
      real(cp) :: VrMean,VtMean,VpMean,Vr2Mean,Vt2Mean,Vp2Mean
      real(cp) :: LFmean,cvrMean,dvpdrMean,VrdVpdrMean,VtcVrMean
      real(cp) :: Bs2Mean,BszMean,BspMean,BpzMean,BspdMean,BpsdMean
      real(cp) :: BzpdMean,BpzdMean,VpPhiMean
      real(cp) :: sinT,Osin,Osin2,cosT,phiNorm
      real(cp) :: BsL,BzL,BpL,Bs2F1,Bs2F2,Bs2F3,BspF2
      real(cp) :: BpzF1,BpzF2,BszF1,BszF2,BszF3

      if ( lVerbose ) write(output_unit,*) '! Starting getTO!'

      phiNorm=one/real(n_phi_max, kind=cp)

      !-- Set values to zero before filling it
      dzCorAS_Rloc(:,nR) =0.0_cp
      dzRstrAS_Rloc(:,nR)=0.0_cp
      dzAstrAS_Rloc(:,nR)=0.0_cp
      if ( l_mag ) dzLFAS_Rloc(:,nR)=0.0_cp
      if ( l_phase_field ) dzPenAS_Rloc(:,nR)=0.0_cp

      !-- Big loop over thetas in block:
      do nTheta=1,n_theta_max
         nTheta1=n_theta_cal2ord(nTheta)
         sinT =sinTheta(nTheta)
         cosT =cosTheta(nTheta)
         Osin =one/sinT
         Osin2=Osin*Osin
         if ( l_mag ) then
            Bs2F1=sinT*sinT*or4(nR)
            Bs2F2=cosT*cosT*Osin2*or2(nR)
            Bs2F3=two*cosT*or3(nR)
            BspF2=cosT*Osin2*or2(nR)
            BpzF1=cosT*Osin*or3(nR)
            BpzF2=or2(nR)*Osin
            BszF1=sinT*cosT*or4(nR)
            BszF2=(two*cosT*cosT-one)*Osin*or3(nR)
            BszF3=cosT*Osin*or2(nR)
         end if

         !--- Get zonal means of velocity and derivatives:
         VrMean     =0.0_cp
         VtMean     =0.0_cp
         VpMean     =0.0_cp
         Vr2Mean    =0.0_cp
         Vt2Mean    =0.0_cp
         Vp2Mean    =0.0_cp
         dVpdrMean  =0.0_cp
         cVrMean    =0.0_cp
         VrdVpdrMean=0.0_cp
         VtcVrMean  =0.0_cp
         VpPhiMean  =0.0_cp
         if ( l_mag ) then
            LFmean     =0.0_cp
            Bs2Mean    =0.0_cp
            BspMean    =0.0_cp
            BpzMean    =0.0_cp
            BszMean    =0.0_cp
            BspdMean   =0.0_cp
            BpsdMean   =0.0_cp
            BzpdMean   =0.0_cp
            BpzdMean   =0.0_cp
         end if
         do nPhi=1,n_phi_max
            nelem = radlatlon2spat(nTheta,nPhi,nR)
            VrMean =VrMean +vr(nelem)
            VtMean =VtMean +vt(nelem)
            VpMean =VpMean +vp(nelem)
            Vr2Mean=Vr2Mean+vr(nelem)*vr(nelem)
            Vt2Mean=Vt2Mean+vt(nelem)*vt(nelem)
            Vp2Mean=Vp2Mean+vp(nelem)*vp(nelem)
            dVpdrMean=dVpdrMean + orho1(nR)*( dvpdr(nelem) - beta(nR)*vp(nelem) )
            cVrMean  =cVrMean  +cvr(nelem)
            VrdVpdrMean=VrdVpdrMean  + orho1(nR)*           & ! rho * vr * dvp/dr
            &           vr(nelem)*(dvpdr(nelem) &
            &                -beta(nR)*vp(nelem))
            VtcVrMean=VtcVrMean    + orho1(nR)*  &  ! rho * vt * cvr
            &         vt(nelem)*cvr(nelem)
            if ( l_phase_field ) VpPhiMean=VpPhiMean+phase(nelem)*vp(nelem)
            if ( l_mag ) then
               LFmean=LFmean + cbr(nelem)*bt(nelem) - &
               &               cbt(nelem)*br(nelem)
               Bs2Mean=Bs2Mean + Bs2F1*br(nelem)*br(nelem) + &
               &                 Bs2F2*bt(nelem)*bt(nelem) + &
               &                 Bs2F3*br(nelem)*bt(nelem)
               BspMean=BspMean + or3(nR)*br(nelem)*bp(nelem) + &
               &                 BspF2  *bt(nelem)*bp(nelem)
               BpzMean=BpzMean + BpzF1*br(nelem)*bp(nelem) - &
               &                 BpzF2*bt(nelem)*bp(nelem)
               BszMean=BszMean + BszF1*br(nelem)*br(nelem) + &
               &                 BszF2*br(nelem)*bt(nelem) - &
               &                 BszF3*bt(nelem)*bt(nelem)
               BsL=sinT*or2(nR)*br(nelem) + cosT*Osin*or1(nR)*bt(nelem)
               BpL=Osin*or1(nR)*bp(nelem)
               BzL=cosT*or2(nR)*br(nelem) - or1(nR)*bt(nelem)
               BspdMean=BspdMean+BsL*(BpL-BpLast(nTheta,nPhi,nR))
               BpsdMean=BpsdMean+BpL*(BsL-BsLast(nTheta,nPhi,nR))
               BzpdMean=BzpdMean+BzL*(BpL-BpLast(nTheta,nPhi,nR))
               BpzdMean=BpzdMean+BpL*(BzL-BzLast(nTheta,nPhi,nR))
            end if
         end do

         Vr2Mean=phiNorm*or4(nR)*Vr2Mean
         Vt2Mean=phiNorm*or2(nR)*Osin2*Vt2Mean
         Vp2Mean=phiNorm*or2(nR)*Osin2*Vp2Mean

         V2AS_Rloc(nTheta1,nR)=Vr2Mean+Vt2Mean+Vp2Mean
         VpMean =phiNorm*or1(nR)*Osin*VpMean
         VAS_Rloc(nTheta1,nR)=orho1(nR)*VpMean
         !--- This is Coriolis force = 2\Omega u_s
         dzCorAS_Rloc(nTheta1,nR)=-phiNorm*two*CorFac *  &
         &                 (or2(nR)*sinT*VrMean+or1(nR)*cosT*Osin*VtMean)
         if ( l_phase_field ) then
            dzPenAS_Rloc(nTheta1,nR)=-phiNorm*VpPhiMean*or1(nR)*Osin/ &
            &                         epsPhase**2/penaltyFac**2
         end if
         if ( l_mag ) then
            !--- This is the Lorentz force
            dzLFAS_Rloc(nTheta1,nR)=phiNorm*or3(nR)*Osin*LFmean
            Bs2AS_Rloc(nTheta1,nR) =phiNorm*Bs2Mean
            BspAS_Rloc(nTheta1,nR) =phiNorm*BspMean
            BpzAS_Rloc(nTheta1,nR) =phiNorm*BpzMean
            BszAS_Rloc(nTheta1,nR) =phiNorm*BszMean
            BspdAS_Rloc(nTheta1,nR)=phiNorm*(BspdMean/dtLast)
            BpsdAS_Rloc(nTheta1,nR)=phiNorm*(BpsdMean/dtLast)
            BzpdAS_Rloc(nTheta1,nR)=phiNorm*(BzpdMean/dtLast)
            BpzdAS_Rloc(nTheta1,nR)=phiNorm*(BpzdMean/dtLast)
         end if

         ! dVpdrMean, VtcVrMean and VrdVpdrMean are already divided by rho
         dzRstrAS_Rloc(nTheta1,nR)=                           &
         &           -phiNorm * or3(nR)*Osin*  (VrdVpdrMean - &
         &                         phiNorm*VrMean*dVpdrMean + &
         &                                        VtcVrMean - &
         &                  orho1(nR)*phiNorm*VtMean*cVrMean )
         dzAstrAS_Rloc(nTheta1,nR)=                           &
         &  -phiNorm*or3(nR)*Osin*phiNorm*(VrMean*dVpdrMean + &
         &                         orho1(nR)*VtMean*cVrMean )

      end do ! Loop over thetas

      if ( lVerbose ) write(output_unit,*) '! End of getTO!'

   end subroutine getTO
!-----------------------------------------------------------------------------
   subroutine getTOnext(br,bt,bp,lTONext,lTONext2,dt,dtLast,nR)
      !
      !  Preparing TO calculation by storing flow and magnetic field
      !  contribution to build time derivative.
      !

      !-- Input of variables:
      real(cp), intent(in) :: dt,dtLast
      integer,  intent(in) :: nR
      logical,  intent(in) :: lTONext,lTONext2
      real(cp), intent(in) :: br(*),bt(*),bp(*)

      !-- Local variables:
      integer :: nPhi,nTh,nelem

      if ( lVerbose ) write(output_unit,*) '! Starting getTOnext!',dtLast

      if ( lTONext2 ) then

         dzddVpLMr(1,nR) =0.0_cp
         dzddVpLMr(1:,nR)=zASL(1:)

      else if ( lTOnext ) then

         if ( l_mag ) then
            !$omp parallel do default(shared)    &
            !$omp& private(nPhi,nTh)
            do nPhi=1,n_phi_max
               do nTh=1,n_theta_max
                  nelem = radlatlon2spat(nTh,nPhi,nR)
                  BsLast(nTh,nPhi,nR)=or2(nR)*sinTheta(nTh)*br(nelem) + &
                  &                 or1(nR)*cosTheta(nTh)*O_sin_theta(nTh)*bt(nelem)
                  BpLast(nTh,nPhi,nR)=or1(nR)*bp(nelem)*O_sin_theta(nTh)
                  BzLast(nTh,nPhi,nR)=cosTheta(nTh)*or2(nR)*br(nelem)-or1(nR)*bt(nelem)
               end do
            end do
            !$omp end parallel do
         end if

         dzdVpLMr(1,nR)  =0.0_cp
         dzdVpLMr(1:,nR) =zASL(1:)
         dzddVpLMr(1,nR) =0.0_cp
         dzddVpLMr(1:,nR)=( dzddVpLMr(1:,nR)-((dtLast+dt)/dt)*zASL(1:) )/dtLast
      end if

      if ( lVerbose ) write(output_unit,*) '! End of getTOnext!'

   end subroutine getTOnext
!-----------------------------------------------------------------------------
   subroutine getTOfinish(nR,dtLast)
      !
      ! This handles the computation of the axisymmetric viscous stress
      ! in spectral space using z and its radial derivatives, and then
      ! transform it to grid space using SHTs.
      !

      !-- Input of variables:
      integer,  intent(in) :: nR
      real(cp), intent(in) :: dtLast

      !-- Local variables:
      real(cp) :: dzStrLMr(l_max+1)
      integer :: l,lm

      !------ When all thetas are done calculate viscous stress in LM space:
      dzStrLMr(1)    =0.0_cp
      dzdVpLMr(1,nR) =0.0_cp
      dzddVpLMr(1,nR)=0.0_cp
      do l=1,l_max
         lm=lm2(l,0)
         dzStrLMr(l+1)= hdif_V(l) * (  ddzASL(l+1,nR) -  beta(nR)* dzASL(l+1) - &
         &         (dLh(lm)*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR))*zASL(l+1) )
         dzdVpLMr(l+1,nR) =(zASL(l+1)-dzdVpLMr(l+1,nR))/dtLast
         dzddVpLMr(l+1,nR)=(zASL(l+1)/dtLast+dzddVpLMr(l+1,nR))/dtLast
      end do

      !-- Bring back quantity to physical-space (also unscramble theta)
      call get_PAS(dzStrLMr(:), dzStrAS_Rloc(:,nR), nR)
      call get_PAS(dzdVpLMr(:,nR), dzdVpAS_Rloc(:,nR), nR)
      call get_PAS(dzddVpLMr(:,nR), dzddVpAS_Rloc(:,nR), nR)

   end subroutine getTOfinish
!-----------------------------------------------------------------------------
   subroutine get_PAS(Tlm,Bp,nR)
      !
      !  Purpose of this subroutine is to calculate the axisymmetric
      !  component Bp of an axisymmetric toroidal field Tlm
      !  given in spherical harmonic space (1:lmax+1).
      !  Unscrambling of theta is also ensured

      !-- Input variables
      integer,  intent(in) :: nR
      real(cp), intent(in) :: Tlm(:)    ! field in (l,m)-space for rT

      !-- Output variables:
      real(cp), intent(out) :: Bp(:)

      !-- Local variables:
      integer :: lm,l
      integer :: nTheta,nTheta1
      complex(cp) :: Tl_AX(1:l_max+1)
      real(cp) :: tmpt(nlat_padded), tmpp(nlat_padded)

      do l=0,l_max
         lm=lm2(l,0)
         Tl_AX(l+1)=cmplx(Tlm(lm),0.0_cp,kind=cp)
      end do

      call toraxi_to_spat(Tl_AX(1:l_max+1), tmpt(:), tmpp(:))

      !-- Unscramble theta
      do nTheta=1,n_theta_max
         nTheta1=n_theta_cal2ord(nTheta)
         Bp(nTheta1)=O_sin_theta(nTheta)*tmpp(nTheta)*or1(nR)
      end do

   end subroutine get_PAS
!-----------------------------------------------------------------------------
end module torsional_oscillations
