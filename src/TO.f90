module torsional_oscillations
   !
   !  This module contains information for TO calculation and output
   !

   use iso_fortran_env, only: output_unit
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_phi_maxStr, n_r_maxStr, l_max, n_theta_maxStr, lm_max, &
       &                 nlat_padded
   use radial_data, only: n_r_cmb, nRstart, nRstop
   use radial_functions, only: r, or1, or2, or3, or4, beta, orho1, dbeta
   use physical_parameters, only: CorFac, kbotv, ktopv
   use blocking, only: lm2, llmMag, ulmMag, lo_map, lm2l, lm2m
   use horizontal_data, only: sinTheta, cosTheta, hdif_V, dTheta1A, dTheta1S, dLh, &
       &                      n_theta_cal2ord, O_sin_theta
   use constants, only: one, two
   use logic, only: lVerbose, l_mag
   use sht, only: spat_to_SH_axi, toraxi_to_spat

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
   real(cp), public, allocatable :: BpzdAS_Rloc(:,:)

   real(cp), allocatable :: BsLast(:,:,:), BpLast(:,:,:), BzLast(:,:,:)

   public :: initialize_TO, prep_TO_axi, getTO, getTOnext, getTOfinish, finalize_TO

contains

   subroutine initialize_TO
      !
      ! Allocate the memory needed
      !

      allocate( dzStrAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzRstrAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzAStrAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzCorAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzLFAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzdVpAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( dzddVpAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+7*(nRstop-nRstart+1)*n_theta_maxStr* &
      &                 SIZEOF_DEF_REAL
      allocate( dzdVpLMr(l_max+1,nRstart:nRstop), dzddVpLMr(l_max+1,nRstart:nRstop) )
      dzdVpLMr(:,:) =0.0_cp
      dzddVpLMr(:,:)=0.0_cp
      bytes_allocated = bytes_allocated+2*(l_max+1)*(nRstop-nRstart+1)*SIZEOF_DEF_REAL

      allocate( VAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( V2AS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( Bs2AS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BszAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BspAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BpzAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BspdAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BpsdAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BzpdAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      allocate( BpzdAS_Rloc(n_theta_maxStr,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+10*(nRstop-nRstart+1)*n_theta_maxStr* &
      &                 SIZEOF_DEF_REAL

      allocate( zASL(l_max+1), dzASL(l_max+1), ddzASL(l_max+1,n_r_maxStr) )
      bytes_allocated = bytes_allocated+2*(l_max+1)+(l_max+1)*n_r_maxStr*&
      &                 SIZEOF_DEF_REAL
      ddzASL(:,:)=0.0_cp
      dzASL(:)   =0.0_cp
      zASL(:)    =0.0_cp

      allocate( BsLast(n_theta_maxStr,n_phi_maxStr,nRstart:nRstop) )
      allocate( BpLast(n_theta_maxStr,n_phi_maxStr,nRstart:nRstop) )
      allocate( BzLast(n_theta_maxStr,n_phi_maxStr,nRstart:nRstop) )
      bytes_allocated = bytes_allocated+3*n_phi_maxStr*n_theta_maxStr* &
      &                 (nRstop-nRstart+1)*SIZEOF_DEF_REAL

   end subroutine initialize_TO
!-----------------------------------------------------------------------------
   subroutine finalize_TO
      !
      ! Deallocate the memory
      !

      deallocate( ddzASL, BpzdAS_Rloc, BzpdAS_Rloc, BpsdAS_Rloc, BspdAS_Rloc )
      deallocate( BpzAS_Rloc, BspAS_Rloc, BszAS_Rloc, Bs2AS_Rloc, V2AS_Rloc )
      deallocate( dzddVpAS_Rloc, dzdVpAS_Rloc, dzLFAS_Rloc, dzCorAS_Rloc )
      deallocate( dzAStrAS_Rloc, dzRstrAS_Rloc, dzStrAS_Rloc )
      deallocate( zASL, dzASL, BsLast, BpLast, BzLast, VAS_Rloc )
      deallocate( dzddVpLMr, dzdVpLMr )

   end subroutine finalize_TO
!-----------------------------------------------------------------------------
   subroutine prep_TO_axi(z,dz)

      !-- Input variables
      complex(cp), intent(in) :: z(:)
      complex(cp), intent(in) :: dz(:)

      !-- Local variables
      integer :: l, m, lm

      do lm=1,lm_max
         l=lm2l(lm)
         m=lm2m(lm)
         if ( l <= l_max .and. m == 0 ) then
            zASL(l+1)  =real(z(lm))
            dzASL(l+1) =real(dz(lm))
         end if
      end do

   end subroutine prep_TO_axi
!-----------------------------------------------------------------------------
   subroutine getTO(vr,vt,vp,cvr,dvpdr,br,bt,bp,cbr,cbt, &
              &     dzRstrLM,dzAstrLM,dzCorLM,dzLFLM,dtLast,nR)
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
      !  In addition three output field are given in (lm,r) space:
      !  dzRstrLMr,dzAstrLMr,dzCorLM,dzLFLM.
      !
      !  These are used to calculate the total Reynolds stress,
      !  advection and viscous stress later. Their calculation
      !  retraces the calculations done in the time-stepping part
      !  of the code.
      !

      !-- Input of variables
      real(cp), intent(in) :: dtLast              ! last time step
      integer,  intent(in) :: nR                 ! radial grid point
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: cvr(:,:),dvpdr(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: cbr(:,:),cbt(:,:)

      !-- Output of arrays needing further treatment in s_getTOfinish.f:
      real(cp), intent(out) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(cp), intent(out) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      integer :: nTheta,nPhi,nTheta1
      real(cp) :: VrMean,VtMean,VpMean,Vr2Mean,Vt2Mean,Vp2Mean
      real(cp) :: LFmean,cvrMean,dvpdrMean,VrdVpdrMean,VtcVrMean
      real(cp) :: Bs2Mean,BszMean,BspMean,BpzMean,BspdMean,BpsdMean
      real(cp) :: BzpdMean,BpzdMean
      real(cp) :: Rmean(nlat_padded),Amean(nlat_padded)
      real(cp) :: dzCorMean(nlat_padded),dzLFmean(nlat_padded)
      real(cp) :: sinT,Osin,Osin2,cosT,cosOsin2,phiNorm
      real(cp) :: BsL,BzL,BpL,Bs2F1,Bs2F2,Bs2F3,BspF1,BspF2
      real(cp) :: BpzF1,BpzF2,BszF1,BszF2,BszF3
      real(cp) :: BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(output_unit,*) '! Starting getTO!'

      phiNorm=one/real(n_phi_maxStr, kind=cp)

      !-- Big loop over thetas in block:
      do nTheta=1,n_theta_maxStr
         nTheta1=n_theta_cal2ord(nTheta)
         sinT =sinTheta(nTheta)
         cosT =cosTheta(nTheta)
         Osin =one/sinT
         Osin2=Osin*Osin
         cosOsin2=cosT*Osin2
         Bs2F1=sinT*sinT*or4(nR)
         Bs2F2=cosT*cosT*Osin2*or2(nR)
         Bs2F3=two*cosT*or3(nR)
         BspF1=or3(nR)
         BspF2=cosT*Osin2*or2(nR)
         BpzF1=cosT*Osin*or3(nR)
         BpzF2=or2(nR)*Osin
         BszF1=sinT*cosT*or4(nR)
         BszF2=(two*cosT*cosT-one)*Osin*or3(nR)
         BszF3=cosT*Osin*or2(nR)
         BsF1 =sinT*or2(nR)
         BsF2 =cosT*Osin*or1(nR)
         BpF1 =Osin*or1(nR)
         BzF1 =cosT*or2(nR)
         BzF2 =or1(nR)

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
         LFmean     =0.0_cp
         Bs2Mean    =0.0_cp
         BspMean    =0.0_cp
         BpzMean    =0.0_cp
         BszMean    =0.0_cp
         BspdMean   =0.0_cp
         BpsdMean   =0.0_cp
         BzpdMean   =0.0_cp
         BpzdMean   =0.0_cp
         do nPhi=1,n_phi_maxStr
            VrMean =VrMean +vr(nTheta,nPhi)
            VtMean =VtMean +vt(nTheta,nPhi)
            VpMean =VpMean +vp(nTheta,nPhi)
            Vr2Mean=Vr2Mean+vr(nTheta,nPhi)*vr(nTheta,nPhi)
            Vt2Mean=Vt2Mean+vt(nTheta,nPhi)*vt(nTheta,nPhi)
            Vp2Mean=Vp2Mean+vp(nTheta,nPhi)*vp(nTheta,nPhi)
            dVpdrMean=dVpdrMean    + orho1(nR)*( dvpdr(nTheta,nPhi) - &
            &                              beta(nR)*vp(nTheta,nPhi) )
            cVrMean  =cVrMean  +cvr(nTheta,nPhi)
            VrdVpdrMean=VrdVpdrMean  + orho1(nR)*           & ! rho * vr * dvp/dr
            & vr(nTheta,nPhi)*(dvpdr(nTheta,nPhi) &
            &                -beta(nR)*vp(nTheta,nPhi))
            VtcVrMean=VtcVrMean    + orho1(nR)*  &  ! rho * vt * cvr
            &         vt(nTheta,nPhi)*cvr(nTheta,nPhi)
            if ( l_mag ) then
               LFmean=LFmean + cbr(nTheta,nPhi)*bt(nTheta,nPhi) - &
               &               cbt(nTheta,nPhi)*br(nTheta,nPhi)
               Bs2Mean=Bs2Mean + Bs2F1*br(nTheta,nPhi)*br(nTheta,nPhi) + &
               &                 Bs2F2*bt(nTheta,nPhi)*bt(nTheta,nPhi) + &
               &                 Bs2F3*br(nTheta,nPhi)*bt(nTheta,nPhi)
               BspMean=BspMean + BspF1*br(nTheta,nPhi)*bp(nTheta,nPhi) + &
               &                 BspF2*bt(nTheta,nPhi)*bp(nTheta,nPhi)
               BpzMean=BpzMean + BpzF1*br(nTheta,nPhi)*bp(nTheta,nPhi) - &
               &                 BpzF2*bt(nTheta,nPhi)*bp(nTheta,nPhi)
               BszMean=BszMean + BszF1*br(nTheta,nPhi)*br(nTheta,nPhi) + &
               &                 BszF2*br(nTheta,nPhi)*bt(nTheta,nPhi) - &
               &                 BszF3*bt(nTheta,nPhi)*bt(nTheta,nPhi)
               BsL=BsF1*br(nTheta,nPhi) + BsF2*bt(nTheta,nPhi)
               BpL=BpF1*bp(nTheta,nPhi)
               BzL=BzF1*br(nTheta,nPhi) - BzF2*bt(nTheta,nPhi)
               BspdMean=BspdMean+BsL*(BpL-BpLast(nTheta,nPhi,nR))
               BpsdMean=BpsdMean+BpL*(BsL-BsLast(nTheta,nPhi,nR))
               BzpdMean=BzpdMean+BzL*(BpL-BpLast(nTheta,nPhi,nR))
               BpzdMean=BpzdMean+BpL*(BzL-BzLast(nTheta,nPhi,nR))
            end if
         end do

         Vr2Mean=phiNorm*or4(nR)*Vr2Mean
         Vt2Mean=phiNorm*or2(nR)*Osin2*Vt2Mean
         Vp2Mean=phiNorm*or2(nR)*Osin2*Vp2Mean
         if ( nR == n_r_CMB ) then
            VrMean =0.0_cp
            Vr2Mean=0.0_cp
            if ( ktopv == 2 ) then
               VtMean =0.0_cp
               Vt2Mean=0.0_cp
               VpMean =0.0_cp
               Vp2Mean=0.0_cp
            end if
         end if
         if ( nR == n_r_CMB ) then
            VrMean =0.0_cp
            Vr2Mean=0.0_cp
            if ( kbotv == 2 ) then
               VtMean =0.0_cp
               Vt2Mean=0.0_cp
               VpMean =0.0_cp
               Vp2Mean=0.0_cp
            end if
         end if
         V2AS_Rloc(nTheta1,nR)=Vr2Mean+Vt2Mean+Vp2Mean
         VpMean =phiNorm*or1(nR)*Osin*VpMean
         VAS_Rloc(nTheta1,nR)=orho1(nR)*VpMean
         !--- This is Coriolis force / r*sin(theta)
         dzCorMean(nTheta)=phiNorm*two*CorFac * &
         &                 (or3(nR)*VrMean+or2(nR)*cosOsin2*VtMean)
         if ( l_mag ) then
            !--- This is Lorentz force/ r*sin(theta)
            dzLFmean(nTheta)      =phiNorm*or4(nR)*Osin2*LFmean
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
         Rmean(nTheta)=                                  &
         &            phiNorm * or4(nR)*Osin2* (VrdVpdrMean - &
         &                         phiNorm*VrMean*dVpdrMean + &
         &                                        VtcVrMean - &
         &                  orho1(nR)*phiNorm*VtMean*cVrMean )
         Amean(nTheta)=                                  &
         &  phiNorm*or4(nR)*Osin2*phiNorm*(VrMean*dVpdrMean + &
         &                         orho1(nR)*VtMean*cVrMean )

      end do ! Loop over thetas

      !--- Transform and Add to LM-Space:
      !------ Add contribution from thetas in block:
      !       note legtfAS2 returns modes l=0 -- l=l_max+1
      call spat_to_SH_axi(Rmean,dzRstrLM)
      call spat_to_SH_axi(Amean,dzAstrLM)
      call spat_to_SH_axi(dzCorMean,dzCorLM)
      call spat_to_SH_axi(dZLFmean,dzLFLM)

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
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)

      !-- Local variables:
      integer :: l,lm,nTheta,nPhi

      real(cp) :: sinT,cosT
      real(cp) :: BsF1,BsF2,BpF1,BzF1,BzF2

      if ( lVerbose ) write(output_unit,*) '! Starting getTOnext!',dtLast

      if ( lTONext2 ) then

         dzddVpLMr(1,nR)=0.0_cp
         do l=1,l_max
            lm=lm2(l,0)
            dzddVpLMr(l+1,nR)=zASL(l+1)
         end do

      else if ( lTOnext ) then

         !$omp parallel do default(shared)                     &
         !$omp& private(nTheta,nPhi,sinT,cosT,BsF1,BsF2,BpF1,BzF1,BzF2)
         do nPhi=1,n_phi_maxStr
            do nTheta=1,n_theta_maxStr
               sinT =sinTheta(nTheta)
               cosT =cosTheta(nTheta)
               BsF1 =sinT*or2(nR)
               BsF2 =cosT/sinT*or1(nR)
               BpF1 =or1(nR)/sinT
               BzF1 =cosT*or2(nR)
               BzF2 =or1(nR)

               BsLast(nTheta,nPhi,nR)=BsF1*br(nTheta,nPhi) + BsF2*bt(nTheta,nPhi)
               BpLast(nTheta,nPhi,nR)=BpF1*bp(nTheta,nPhi)
               BzLast(nTheta,nPhi,nR)=BzF1*br(nTheta,nPhi) - BzF2*bt(nTheta,nPhi)
            end do
         end do 
         !$omp end parallel do

         dzdVpLMr(1,nR) =0.0_cp
         dzddVpLMr(1,nR)=0.0_cp
         do l=1,l_max
            lm=lm2(l,0)
            dzdVpLMr(l+1,nR) = zASL(l+1)
            dzddVpLMr(l+1,nR)= ( dzddVpLMr(l+1,nR)-((dtLast+dt)/dt)*zASL(l+1) )/dtLast
         end do
      end if

      if ( lVerbose ) write(output_unit,*) '! End of getTOnext!'

   end subroutine getTOnext
!-----------------------------------------------------------------------------
   subroutine getTOfinish(nR,dtLast,dzRstrLM,dzAstrLM,dzCorLM,dzLFLM)
      !
      !  This program was previously part of getTO(...)
      !  It has now been separated to get it out of the theta-block loop.
      !

      !-- Input of variables:
      integer,  intent(in) :: nR
      real(cp), intent(in) :: dtLast
      real(cp), intent(in) :: dzRstrLM(l_max+2),dzAstrLM(l_max+2)
      real(cp), intent(in) :: dzCorLM(l_max+2),dzLFLM(l_max+2)

      !-- Local variables:
      real(cp) :: dzStrLMr(l_max+1), dzRstrLMr(l_max+1), dzAstrLMr(l_max+1)
      real(cp) :: dzCorLMr(l_max+1), dzLFLMr(l_max+1)
      integer :: l,lS,lA,lm

      !------ When all thetas are done calculate viscous stress in LM space:
      dzStrLMr(1) =0.0_cp
      dzRstrLMr(1)=0.0_cp
      dzAstrLMr(1)=0.0_cp
      dzCorLMr(1) =0.0_cp
      dzLFLMr(1)  =0.0_cp
      dzdVpLMr(1,nR) =0.0_cp
      dzddVpLMr(1,nR)=0.0_cp
      do l=1,l_max
         lS=(l-1)+1
         lA=(l+1)+1
         lm=lm2(l,0)
         dzStrLMr(l+1)= hdif_V(l) * (  ddzASL(l+1,nR) -  beta(nR)* dzASL(l+1) - &
         &         (dLh(lm)*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR))*zASL(l+1) )
         !---- -r**2/(l(l+1)) 1/sin(theta) dtheta sin(theta)**2
         !     minus sign to bring stuff on the RHS of NS equation !
         dzRstrLMr(l+1)=-r(nR)*r(nR)/dLh(lm) * ( dTheta1S(lm)*dzRstrLM(lS) - &
         &                                       dTheta1A(lm)*dzRstrLM(lA) )
         dzAstrLMr(l+1)=-r(nR)*r(nR)/dLh(lm) * ( dTheta1S(lm)*dzAstrLM(lS) - &
         &                                       dTheta1A(lm)*dzAstrLM(lA) )
         dzCorLMr(l+1) =-r(nR)*r(nR)/dLh(lm) * ( dTheta1S(lm)*dzCorLM(lS) - &
         &                                       dTheta1A(lm)*dzCorLM(lA) )
         dzLFLMr(l+1)  = r(nR)*r(nR)/dLh(lm) * ( dTheta1S(lm)*dzLFLM(lS) - &
         &                                       dTheta1A(lm)*dzLFLM(lA) )
         dzdVpLMr(l+1,nR) =(zASL(l+1)-dzdVpLMr(l+1,nR))/dtLast
         dzddVpLMr(l+1,nR)=(zASL(l+1)/dtLast+dzddVpLMr(l+1,nR))/dtLast
      end do

      !-- Bring back quantity to physical-space (also unscramble theta)
      call get_PAS(dzStrLMr(:), dzStrAS_Rloc(:,nR), r(nR))
      call get_PAS(dzRstrLMr(:), dzRstrAS_Rloc(:,nR), r(nR))
      call get_PAS(dzAstrLMr(:), dzAstrAS_Rloc(:,nR), r(nR))
      call get_PAS(dzCorLMr(:), dzCorAS_Rloc(:,nR), r(nR))
      call get_PAS(dzLFLMr(:), dzLFAS_Rloc(:,nR), r(nR))
      call get_PAS(dzdVpLMr(:,nR), dzdVpAS_Rloc(:,nR), r(nR))
      call get_PAS(dzddVpLMr(:,nR), dzddVpAS_Rloc(:,nR), r(nR))

   end subroutine getTOfinish
!-----------------------------------------------------------------------------
   subroutine get_PAS(Tlm,Bp,rT)
      !
      !  Purpose of this subroutine is to calculate the axisymmetric
      !  component Bp of an axisymmetric toroidal field Tlm
      !  given in spherical harmonic space (1:lmax+1).
      !  Unscrambling of theta is also ensured

      !-- Input variables
      real(cp), intent(in) :: rT             ! radius
      real(cp), intent(in) :: Tlm(:)         ! field in (l,m)-space for rT

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
      do nTheta=1,n_theta_maxStr
         nTheta1=n_theta_cal2ord(nTheta)
         Bp(nTheta1)=O_sin_theta(nTheta)*tmpp(nTheta)/rT
      end do

   end subroutine get_PAS
!-----------------------------------------------------------------------------
end module torsional_oscillations
