module mri
   !
   !  This module contains information for MRI calculation and output
   !

  use iso_fortran_env, only: output_unit
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_phi_maxStr, n_r_maxStr, l_max, n_theta_maxStr, lm_max, &
        &                 nlat_padded
   use radial_data, only: n_r_cmb, nRstart, nRstop
   use radial_functions, only: r, or1, or2, or3, or4, beta, orho1, &
       &                       dbeta,visc
   use physical_parameters, only: CorFac, kbotv, ktopv
   use blocking, only: lm2,llmMag, ulmMag, lo_map
   use horizontal_data, only: sinTheta, cosTheta, hdif_V, dTheta1A, dTheta1S, &
       &                      dLh,O_sin_theta_E2 !osn2
   use constants, only: one, two, third
   use logic, only: lVerbose, l_mag, l_MRICalc
   use sht, only: spat_to_SH_axi

   implicit none

   private

   real(cp), public, allocatable :: VspAS_Rloc(:,:) !VspAS = (Vs*Vp).mean(axis='phi')
   real(cp), public, allocatable :: BspAS_Rloc(:,:) !BspAS = (Bs*Bp).mean(axis='phi')
   !real(cp), public, allocatable :: EkinAS_Rloc(:,:)
   !real(cp), public, allocatable :: EmagAS_Rloc(:,:)
   real(cp), public, allocatable :: VsAS_Rloc(:,:)  !Vs
   real(cp), public, allocatable :: VpAS_Rloc(:,:)  !Vp
   real(cp), public, allocatable :: VzAS_Rloc(:,:)  !Vz
   real(cp), public, allocatable :: BsAS_Rloc(:,:)  !Bs
   real(cp), public, allocatable :: BpAS_Rloc(:,:)  !Bp
   real(cp), public, allocatable :: BzAS_Rloc(:,:)  !Bz
   real(cp), public, allocatable :: Vs2AS_Rloc(:,:)  !Vs
   real(cp), public, allocatable :: Vp2AS_Rloc(:,:)  !Vp
   real(cp), public, allocatable :: Vz2AS_Rloc(:,:)  !Vz
   real(cp), public, allocatable :: Bs2AS_Rloc(:,:)  !Bs
   real(cp), public, allocatable :: Bp2AS_Rloc(:,:)  !Bp
   real(cp), public, allocatable :: Bz2AS_Rloc(:,:)  !Bz
   real(cp), public, allocatable :: fvisc_nr_cmbAS(:),fvisc_nr_1AS(:)


   public :: initialize_MRI, get_MRI, finalize_MRI

contains
!-----------------------------------------------------------------------------
  subroutine initialize_MRI
      !
      ! Allocate the memory needed
      !


    allocate( VspAS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( VsAS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( VpAS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( VzAS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( Vs2AS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( Vp2AS_Rloc(nlat_padded,nRstart:nRstop) )
    allocate( Vz2AS_Rloc(nlat_padded,nRstart:nRstop) )
    if (l_mag) then
       allocate( BspAS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( BsAS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( BpAS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( BzAS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( Bs2AS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( Bp2AS_Rloc(nlat_padded,nRstart:nRstop) )
       allocate( Bz2AS_Rloc(nlat_padded,nRstart:nRstop) )
    end if
    allocate (fvisc_nr_cmbAS(nlat_padded))
    allocate (fvisc_nr_1AS(nlat_padded))
    !allocate( EkinAS_Rloc(nlat_padded,nRstart:nRstop) )
    !allocate( EmagAS_Rloc(nlat_padded,nRstart:nRstop) )
    bytes_allocated = bytes_allocated+ &
         &                 10*(nRstop-nRstart+1)*nlat_padded*SIZEOF_DEF_REAL

  end subroutine initialize_MRI
!-----------------------------------------------------------------------------
  subroutine finalize_MRI
      !
      ! Deallocate the memory
      !

    deallocate(VspAS_Rloc)
    deallocate(VsAS_Rloc)
    deallocate(VpAS_Rloc)
    deallocate(VzAS_Rloc)
    if (l_mag) then
       deallocate(BspAS_Rloc)
       deallocate(BsAS_Rloc)
       deallocate(BpAS_Rloc)
       deallocate(BzAS_Rloc)
    end if
    deallocate(fvisc_nr_cmbAS)
    deallocate(fvisc_nr_1AS)
    !deallocate(EkinAS_Rloc)
    !deallocate(EmagAS_Rloc)

  end subroutine finalize_MRI
!-----------------------------------------------------------------------------
  subroutine get_MRI(vr,vt,vp,br,bt,bp,dvrdr,dvtdr,dvpdr,dvrdt,dvrdp,nR)
     !  This program calculates various axisymmetric linear
      !  and nonlinear variables for a radial grid point nR and
      !  a theta-block.
      !  Input are the fields vr,vt,vp
      !  Output are non-linear azimuthally averaged  field VpsAS (flow phi component).
      !  These are give in (r,theta)-space.
      !  Also in (r,theta)-space are azimuthally averaged correlations of
      !  non-axisymmetric magnetic field component :
      !  Bsp=Bs*Bp. These are used to calulcate
      !  the respective Maxwell and Reynolds stress.
      !  Their calculation
      !  retraces the calculations done in the time-stepping part
      !  of the code.
      !

      !-- Input of variables
      integer,  intent(in) :: nR                 ! radial grid point
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)
      real(cp), intent(in) :: dvrdr(:,:),dvtdr(:,:),dvpdr(:,:)
      real(cp), intent(in) :: dvrdt(:,:),dvrdp(:,:)

      !-- Local variables:
      integer :: nTheta,nPhi, nThetaNHS
      real(cp) :: BspMean,VspMean
      real(cp) :: Vp2Mean,Vs2Mean,Vz2Mean
      real(cp) :: Bp2Mean,Bs2Mean,Bz2Mean
      real(cp) :: VsMean,VpMean,VzMean
      real(cp) :: BsMean,BpMean,BzMean
      real(cp) :: sinT,Osin,Osin2,cosT,cosOsin2
      real(cp) :: phiNorm
      real(cp) :: BspF1,BspF2
      real(cp) :: BsF1,BsF2,Bs2F1,Bs2F2,Bs2F3
      real(cp) :: BzF1,BzF2,BpF1,Bz2F1,Bz2F2,Bz2F3
      ! for visc flux at bcs
      real(cp) :: fvisc

      if ( lVerbose ) write(*,*) '! Starting getMRI!'

      phiNorm=one/real(n_phi_maxStr, kind=cp)

      if (nR == n_r_cmb) fvisc_nr_cmbAS=0.0_cp
      if (nR == n_r_cmb+1) fvisc_nr_1AS=0.0_cp

      !-- Big loop over thetas in block:
      do nTheta=1,n_theta_maxStr
         nThetaNHS = (nTheta+1)/2
         sinT =sinTheta(nTheta)
         cosT =cosTheta(nTheta)
         Osin =one/sinT
         Osin2=Osin*Osin
         BspF1=or3(nR)
         BspF2=cosT*Osin2*or2(nR)
         BsF1 =sinT*or2(nR)
         BsF2 =cosT*Osin*or1(nR)
         BzF1 =cosT*or2(nR)
         BzF2 =or1(nR)
         BpF1 =or1(nR)*Osin
         Bs2F1=sinT*sinT*or4(nR)
         Bs2F2= cosT*cosT*Osin2*or2(nR)
         Bs2F3=two*cosT*or3(nR)
         Bz2F1=cosT*cosT*or4(nR)
         Bz2F2=or2(nR)
         Bz2F3=two*cosT*or3(nR)
         !--- Get zonal means of velocity and derivatives:
         BspMean    =0.0_cp
         VspMean    =0.0_cp
         Vpmean     =0.0_cp
         Vp2Mean   =0.0_cp
         Vs2Mean   =0.0_cp
         Vz2Mean   =0.0_cp
         Bp2Mean   =0.0_cp
         Bs2Mean   =0.0_cp
         Bz2Mean   =0.0_cp
         VsMean   =0.0_cp
         VpMean   =0.0_cp
         VzMean   =0.0_cp
         BsMean   =0.0_cp
         BpMean   =0.0_cp
         BzMean   =0.0_cp
         fvisc = 0
         do nPhi=1,n_phi_maxStr
            Vs2Mean=Vs2Mean   +Bs2F1*vr(nTheta,nPhi)*vr(nTheta,nPhi) &
                   &           +Bs2F2*vt(nTheta,nPhi)*vt(nTheta,nPhi)&
                   &           +Bs2F3*vr(nTheta,nPhi)*vt(nTheta,nPhi)
            Vz2Mean=Vz2Mean   +Bz2F1*vr(nTheta,nPhi)*vr(nTheta,nPhi) &
                    &          +Bz2F2*vt(nTheta,nPhi)*vt(nTheta,nPhi)&
                    &          +Bz2F3*vr(nTheta,nPhi)*vt(nTheta,nPhi)
            VsMean = VsMean + BsF1*vr(nTheta,nPhi) + BsF2*vt(nTheta,nPhi)
            VpMean = VpMean + vp(nTheta,nPhi)
            VzMean = VzMean + BzF1*vr(nTheta,nPhi) - BzF2*vt(nTheta,nPhi)
            if (nR == n_r_cmb) then
               fvisc = fvisc - visc(nR)*orho1(nR)*vp(nTheta,nPhi)*         &
                    &                              O_sin_theta_E2(nThetaNHS)* (           &
                    &                              dvpdr(nTheta,nPhi)          &
                    &       -(two*or1(nR)+beta(nR))*vp(nTheta,nPhi) )

            else if (nR== n_r_cmb+1) then
               fvisc = fvisc -two*visc(nR)*orho1(nR)*vr(nTheta,nPhi)*or2(nR)* (&
                    &                             dvrdr(nTheta,nPhi)                 &
                    & -(two*or1(nR)+two*third*beta(nR))*vr(nTheta,nPhi) )-       &
                    &                       visc(nR)*orho1(nR)*vt(nTheta,nPhi)*  &
                    &                            O_sin_theta_E2(nThetaNHS)* (                   &
                    &                       or2(nR)*dvrdt(nTheta,nPhi)             &
                    &                              +dvtdr(nTheta,nPhi)               &
                    &       -(two*or1(nR)+beta(nR))*vt(nTheta,nPhi) )  -         &
                    &       visc(nR)*orho1(nR)*vp(nTheta,nPhi)*                  &
                    &                               O_sin_theta_E2(nThetaNHS)* (                &
                    &                       or2(nR)*dvrdp(nTheta,nPhi)             &
                    &                              +dvpdr(nTheta,nPhi)               &
                    &       -(two*or1(nR)+beta(nR))*vp(nTheta,nPhi) )

            end if

            if (l_mag) then
               Bs2Mean=Bs2Mean+Bs2F1*br(nTheta,nPhi)*br(nTheta,nPhi) &
                     &         +Bs2F2*bt(nTheta,nPhi)*bt(nTheta,nPhi)&
                     &         +Bs2F3*br(nTheta,nPhi)*bt(nTheta,nPhi)
               Bz2Mean=Bz2Mean+Bz2F1*br(nTheta,nPhi)*br(nTheta,nPhi) &
                     &         +Bz2F2*bt(nTheta,nPhi)*bt(nTheta,nPhi)&
                     &         +Bz2F3*br(nTheta,nPhi)*bt(nTheta,nPhi)
               Bp2Mean=Bp2Mean+bp(nTheta,nPhi)*bp(nTheta,nPhi)
               BsMean = BsMean + BsF1*br(nTheta,nPhi) + BsF2*bt(nTheta,nPhi)
               BpMean = BpMean + bp(nTheta,nPhi)
               BzMean = BzMean + BzF1*br(nTheta,nPhi) - BzF2*bt(nTheta,nPhi)
            end if
         end do
         Vpmean = phiNorm*Vpmean
         do nPhi=1,n_phi_maxStr
            Vp2Mean=Vp2Mean+(vp(nTheta,nPhi)-Vpmean)*(vp(nTheta,nPhi)-Vpmean)
            VspMean=VspMean + BspF1*vr(nTheta,nPhi)*(vp(nTheta,nPhi)-Vpmean) + &
                 &                 BspF2*vt(nTheta,nPhi)*(vp(nTheta,nPhi)-Vpmean)
            if ( l_mag ) then
               BspMean=BspMean + BspF1*br(nTheta,nPhi)*bp(nTheta,nPhi) + &
                    &                 BspF2*bt(nTheta,nPhi)*bp(nTheta,nPhi)
            end if

         end do
         VspAS_Rloc(nTheta,nR) =phiNorm*VspMean
         !EkinAS_Rloc(nTheta,nR) = Vr2Mean+Vt2Mean +Vp2Mean
         !EmagAS_Rloc(nTheta,nR) = Br2Mean+Bt2Mean +Bp2Mean
         VsAS_Rloc(nTheta,nR) = phinorm*VsMean
         VpAs_Rloc(nTheta,nR) = BpF1*VpMean
         VzAS_Rloc(nTheta,nR) = phinorm*VzMean
         Vs2AS_Rloc(nTheta,nR) = phiNorm*Vs2Mean
         Vp2As_Rloc(nTheta,nR) = phiNorm*or2(nR)*Osin2*Vp2Mean
         Vz2AS_Rloc(nTheta,nR) = phiNorm*Vz2Mean
         if (l_mag) then
            BspAS_Rloc(nTheta,nR) =phiNorm*BspMean
            BsAs_Rloc(nTheta,nR) = phinorm*BsMean
            BpAS_Rloc(nTheta,nR) = phinorm*BpF1*BpMean
            BzAs_Rloc(nTheta,nR) = phinorm*BzMean
            Bs2As_Rloc(nTheta,nR) = phinorm*Bs2Mean
            Bp2AS_Rloc(nTheta,nR) = phinorm*or2(nR)*Osin2*Bp2Mean
            Bz2As_Rloc(nTheta,nR) = phinorm*Bz2Mean
         end if
         if (nR == n_r_cmb) then
            fvisc_nr_cmbAS(nTheta) = phinorm*fvisc
         else if (nR == n_r_cmb+1) then
            fvisc_nr_1AS(nTheta)   = phinorm*fvisc
         end if
      end do

      if ( lVerbose ) write(*,*) '! End of getMRI!'

  end subroutine get_MRI
!----------------------------------------------------------------------------

!  subroutine get_viscflux(vr,vt,vp,dvrdr,dvtdr,dvpdr,        &
!       &          dvrdt,dvrdp, fvisc_glob, &
!       &          nR,nThetaStart)
!    !
!    !   Calculates the flux:
!    !
!    !     * Viscous flux: :math:`F_= -(u \cdot S )_r`)
!
!    !-- Input of variables
!    integer,  intent(in) :: nR
!    integer,  intent(in) :: nThetaStart
!    real(cp), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
!    real(cp), intent(in) :: dvrdr(nrp,nfs),dvtdr(nrp,nfs),dvpdr(nrp,nfs)
!    real(cp), intent(in) :: dvrdt(nrp,nfs),dvrdp(nrp,nfs)
!
!    !-- Output variables:
!    real(cp), intent(out) :: fvisc_glob
!
!    !-- Local variables:
!      integer :: nTheta,nThetaB,nThetaNHS
!      integer :: nPhi
!
!      phiNorm=two*pi/real(n_phi_max,cp)
!
!      nTheta=nThetaStart-1
!#ifdef WITH_SHTNS
!      !$OMP PARALLEL DO default(shared)         &
!      !$OMP& private(nThetaB, nTheta, nPhi)     &
!      !$OMP& private(fkin, fconv, fvisc)
!#endif
!      do nThetaB=1,sizeThetaB
!         nTheta=nThetaStart+nThetaB-1
!         nThetaNHS=(nTheta+1)/2
!         fkinAS(nThetaB) =0.0_cp
!         fconvAS(nThetaB)=0.0_cp
!         fviscAS(nThetaB)=0.0_cp
!         fkin=0.0_cp
!         fconv=0.0_cp
!         fvisc=0.0_cp
!         do nPhi=1,n_phi_max
!
!            if (nR == n_r_cmb):
!               fvisc =  - visc(nR+1)*orho1(nR+1)*vp(nPhi,nThetaB)*         &
!               &                               O_sin_theta_E2(nThetaNHS)* (           &
!               &                              dvpdr(nPhi,nThetaB)           &
!               &       -(two*or1(nR+1)+beta(nR+1))*vp(nPhi,nThetaB) )
!           else:
!               fvisc=-two*visc(nR+1)*orho1(nR+1)*vr(nPhi,nThetaB)*or2(nR+1)* (     &
!               &                             dvrdr(nPhi,nThetaB)             &
!               & -(two*or1(nR+1)+two*third*beta(nR+1))*vr(nPhi,nThetaB) )-       &
!               &                       visc(nR+1)*orho1(nR+1)*vt(nPhi,nThetaB)*  &
!               &                            O_sin_theta_E2(nThetaNHS)* (               &
!               &                       or2(nR+1)*dvrdt(nPhi,nThetaB)           &
!               &                              +dvtdr(nPhi,nThetaB)           &
!               &       -(two*or1(nR+1)+beta(nR+1))*vt(nPhi,nThetaB) )  -         &
!               &       visc(nR+1)*orho1(nR+1)*vp(nPhi,nThetaB)*                  &
!               &                               O_sin_theta_E2(nThetaNHS)* (            &
!               &                       or2(nR+1)*dvrdp(nPhi,nThetaB)           &
!               &                              +dvpdr(nPhi,nThetaB)           &
!               &       -(two*or1(nR+1)+beta(nR+1))*vp(nPhi,nThetaB) )
!
!           fviscAS(nThetaB)=fviscAS(nThetaB)+fvisc
!         end do
!         fviscAS(nThetaB)=phiNorm*fviscAS(nThetaB)
!         f_visc_glob = f_visc_glob + gauss(nThetaNHS)*fviscAS
!      end do
!#ifdef WITH_SHTNS
!      !$OMP END PARALLEL DO
!#endif
!
!    end subroutine get_viscflux

end module mri
