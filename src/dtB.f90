module dtB_mod
   !
   !  This module contains magnetic field stretching and advection terms
   !  plus a separate omega-effect.
   !  It is used for movie output.
   !
   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: nrp, n_r_maxMag, n_r_ic_maxMag, n_r_max, lm_max_dtB, &
       &                 n_r_max_dtB, n_r_ic_max_dtB, lm_max, n_cheb_max,     &
       &                 n_r_ic_max, l_max, n_phi_max, ldtBmem, l_axi,        &
       &                 nRstart, nRstop, nThetaStart, nThetaStop
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC
   use mpi_transp, only: type_mpitransp
   use mpi_ptop_mod, only: type_mpiptop
   use physical_parameters, only: opm,O_sr
   use radial_functions, only: O_r_ic, lambda, or2, dLlambda, rscheme_oc, &
       &                       or1, orho1, l_R
   use horizontal_data, only: dPhi, dLh, hdif_B, osn2, cosn2, osn1, &
       &                      dTheta1S, dTheta1A
   use logic, only: l_cond_ic, l_DTrMagSpec, l_dtBmovie
   use blocking, only: lo_map, st_map, l2lmAS, lm2l, lm2m, lmP2lmPS, lmP2lmPA, &
                       lm2lmP, llmMag, ulmMag, llm, ulm
   use radial_spectra ! rBrSpec, rBpSpec
#ifndef WITH_SHTNS
   use fft
   use legendre_grid_to_spec, only: legTF2, legTF3
#else
   use shtns, only: spat_to_SH_dist
#endif
   use constants, only: two
   use radial_der, only: get_dr

   implicit none

   private

   !-- Global arrays!!! They are only required for some movie outputs
   !-- but we should definitely try to get rid of them
   complex(cp), public, allocatable :: PstrLM(:,:), TstrLM(:,:), PadvLM(:,:)
   complex(cp), public, allocatable :: TadvLM(:,:), TomeLM(:,:)
   complex(cp), public, allocatable :: PdifLM(:,:), TdifLM(:,:), PadvLMIC(:,:)
   complex(cp), public, allocatable :: PdifLMIC(:,:), TadvLMIC(:,:), TdifLMIC(:,:)

   !-- Container for R to LM MPI transposes
   complex(cp), allocatable, target :: dtB_LMloc_container(:,:,:)
   complex(cp), allocatable, target :: dtB_Rloc_container(:,:,:)

   !-- R-distributed arrays
   complex(cp), pointer :: PstrLM_Rloc(:,:), PadvLM_Rloc(:,:)
   complex(cp), pointer :: TomeRLM_Rloc(:,:), TomeLM_Rloc(:,:)
   complex(cp), pointer :: TstrRLM_Rloc(:,:), TstrLM_Rloc(:,:)
   complex(cp), pointer :: TadvRLM_Rloc(:,:), TadvLM_Rloc(:,:)

   !-- LM-distributed arrays
   complex(cp), public, pointer :: PstrLM_LMloc(:,:), PadvLM_LMloc(:,:)
   complex(cp), public, pointer :: TomeRLM_LMloc(:,:), TomeLM_LMloc(:,:)
   complex(cp), public, pointer :: TstrRLM_LMloc(:,:), TstrLM_LMloc(:,:)
   complex(cp), public, pointer :: TadvRLM_LMloc(:,:), TadvLM_LMloc(:,:)
   complex(cp), public, allocatable :: PdifLM_LMloc(:,:), TdifLM_LMloc(:,:)
   complex(cp), public, allocatable :: PadvLMIC_LMloc(:,:), PdifLMIC_LMloc(:,:)
   complex(cp), public, allocatable :: TadvLMIC_LMloc(:,:), TdifLMIC_LMloc(:,:)

   class(type_mpitransp), pointer :: r2lo_dtB

   public :: initialize_dtB_mod, get_dtBLMfinish, get_dtBLM, get_dH_dtBLM, &
   &         finalize_dtB_mod

contains

   subroutine initialize_dtB_mod
      !
      ! Memory allocation
      !

      !
      ! The remaining global arrays should be suppressed, they are only
      ! needed because of some movie outputs
      !

      if ( l_dtBmovie ) then
         if ( coord_r == 0 ) then
            allocate( PstrLM(lm_max_dtB,n_r_max_dtB) )
            allocate( PadvLM(lm_max_dtB,n_r_max_dtB) )
            allocate( TstrLM(lm_max_dtB,n_r_max_dtB) )
            allocate( TadvLM(lm_max_dtB,n_r_max_dtB) )
            allocate( TomeLM(lm_max_dtB,n_r_max_dtB) )
            allocate( PdifLM(lm_max_dtB,n_r_max_dtB) )
            allocate( TdifLM(lm_max_dtB,n_r_max_dtB) )
            bytes_allocated = bytes_allocated+ &
            &                 7*lm_max_dtB*n_r_max_dtB*SIZEOF_DEF_COMPLEX
         else
            allocate( PstrLM(1,1) )
            allocate( PadvLM(1,1) )
            allocate( PdifLM(1,1) )
            allocate( TdifLM(1,1) )
            allocate( TstrLM(1,1) )
            allocate( TadvLM(1,1) )
            allocate( TomeLM(1,1) )
         end if

         if ( coord_r == 0 ) then
            allocate( PadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
            allocate( PdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
            allocate( TadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
            allocate( TdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
            bytes_allocated = bytes_allocated+ &
            &                 4*lm_max_dtB*n_r_ic_max_dtB*SIZEOF_DEF_COMPLEX
         else
            allocate( PadvLMIC(1,1) )
            allocate( PdifLMIC(1,1) )
            allocate( TadvLMIC(1,1) )
            allocate( TdifLMIC(1,1) )
         end if
      end if

      allocate( PdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
      allocate( TdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
      bytes_allocated = bytes_allocated+ &
                        2*(ulmMag-llmMag+1)*n_r_max_dtB*SIZEOF_DEF_COMPLEX
      allocate( PadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( PdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( TadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( TdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      bytes_allocated = bytes_allocated+ &
      &                 4*(ulmMag-llmMag+1)*n_r_ic_max_dtB*SIZEOF_DEF_COMPLEX

      allocate( dtB_Rloc_container(lm_max_dtB,nRstart:nRstop,8) )
      TomeLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,1)
      TomeRLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,2)
      TstrLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,3)
      TstrRLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,4)
      TadvLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,5)
      TadvRLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,6)
      PstrLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,7)
      PadvLM_Rloc(1:,nRstart:) => dtB_Rloc_container(1:lm_max_dtB,nRstart:nRstop,8)
      bytes_allocated = bytes_allocated+8*(nRstop-nRstart+1)*lm_max_dtB* &
      &                 SIZEOF_DEF_COMPLEX

      allocate( dtB_LMloc_container(llmMag:ulmMag,n_r_max_dtB,8) )
      TomeLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,1)
      TomeRLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,2)
      TstrLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,3)
      TstrRLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,4)
      TadvLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,5)
      TadvRLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,6)
      PstrLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,7)
      PadvLM_LMloc(llmMag:,1:) => dtB_LMloc_container(llmMag:ulmMag,1:n_r_max_dtB,8)
      bytes_allocated = bytes_allocated+8*(ulmMag-llmMag+1)*n_r_max_dtB* &
      &                 SIZEOF_DEF_COMPLEX

      allocate ( type_mpiptop :: r2lo_dtB )

      call r2lo_dtB%create_comm(8)

   end subroutine initialize_dtB_mod
!----------------------------------------------------------------------------
   subroutine finalize_dtB_mod

      if ( l_dtBmovie ) then
         deallocate( PstrLM, PadvLM, TstrLM, TadvLM, TomeLM )
         deallocate( TdifLMIC, TadvLMIC, PdifLMIC, PadvLMIC, TdifLM, PdifLM )
      end if

      deallocate( PdifLM_LMloc, TdifLM_LMloc, PadvLMIC_LMloc, PdifLMIC_LMloc )
      deallocate( TadvLMIC_LMloc, TdifLMIC_LMloc )
      deallocate( dtB_Rloc_container, dtB_LMloc_container )

      call r2lo_dtB%destroy_comm()

   end subroutine finalize_dtB_mod
!----------------------------------------------------------------------------
   subroutine dtb_gather_lo_on_rank0
      !
      ! MPI gather on rank0 for dtBmovie outputs.
      ! This routine should really be suppressed once the movie
      ! outputs have been improved
      !

      call gather_all_from_lo_to_rank0(gt_OC,PstrLM_LMloc,PstrLM)
      call gather_all_from_lo_to_rank0(gt_OC,TstrLM_LMloc,TstrLM)
      call gather_all_from_lo_to_rank0(gt_OC,PadvLM_LMloc,PadvLM)
      call gather_all_from_lo_to_rank0(gt_OC,TadvLM_LMloc,TadvLM)
      call gather_all_from_lo_to_rank0(gt_OC,TomeLM_LMloc,TomeLM)
      call gather_all_from_lo_to_rank0(gt_OC,PdifLM_LMloc,PdifLM)
      call gather_all_from_lo_to_rank0(gt_OC,TdifLM_LMloc,TdifLM)
      call gather_all_from_lo_to_rank0(gt_IC,PadvLMIC_LMloc,PadvLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,TadvLMIC_LMloc,TadvLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,PdifLMIC_LMloc,PdifLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,TdifLMIC_LMloc,TdifLMIC)

   end subroutine dtb_gather_lo_on_rank0
!----------------------------------------------------------------------------
   subroutine  get_dtBLM(nR,vr,vt,vp,br,bt,bp,BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
               &         BtVpLM,BpVtLM,BrVZLM,BtVZLM,BtVpCotLM,BpVtCotLM,  &
               &         BtVZcotLM,BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)

      !
      !  This subroutine calculates non-linear products in grid-space for radial
      !  level n_r and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3
      !
      !  if lvelo >0 velocities are zero only the (vxB)
      !  contributions to bnlr2-3 need to be calculated
      !
      !  vr...sr: (input) velocity, magnetic field comp. and derivs, entropy
      !                   on grid points
      !  i1: (input) range of points in theta for which calculation is done
      !
      !

      !-- Input variables:
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: vp(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: br(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bt(nrp,nThetaStart:nThetaStop)
      real(cp), intent(in) :: bp(nrp,nThetaStart:nThetaStop)

      !-- Output variables:
      complex(cp), intent(out) :: BtVrLM(*),BpVrLM(*)
      complex(cp), intent(out) :: BrVtLM(*),BrVpLM(*)
      complex(cp), intent(out) :: BtVpLM(*),BpVtLM(*)
      complex(cp), intent(out) :: BrVZLM(*),BtVZLM(*)
      complex(cp), intent(out) :: BpVtCotLM(*),BtVpCotLM(*),BtVZcotLM(*)
      complex(cp), intent(out) :: BtVpSn2LM(*),BpVtSn2LM(*)
      complex(cp), intent(out) :: BtVZsn2LM(*)

      !-- Local variables:
      integer :: lm,l, n_theta,n_phi,n_theta_nhs
      real(cp) :: fac,facCot,vpAS
      real(cp) :: BtVr(nrp,nThetaStart:nThetaStop),BpVr(nrp,nThetaStart:nThetaStop)
      real(cp) :: BrVt(nrp,nThetaStart:nThetaStop),BrVp(nrp,nThetaStart:nThetaStop)
      real(cp) :: BtVp(nrp,nThetaStart:nThetaStop),BpVt(nrp,nThetaStart:nThetaStop)
      real(cp) :: BrVZ(nrp,nThetaStart:nThetaStop),BtVZ(nrp,nThetaStart:nThetaStop)
      real(cp) :: BpVtCot(nrp,nThetaStart:nThetaStop)
      real(cp) :: BtVpCot(nrp,nThetaStart:nThetaStop)
      real(cp) :: BpVtSn2(nrp,nThetaStart:nThetaStop)
      real(cp) :: BtVpSn2(nrp,nThetaStart:nThetaStop)
      real(cp) :: BtVZcot(nrp,nThetaStart:nThetaStop)
      real(cp) :: BtVZsn2(nrp,nThetaStart:nThetaStop)

#ifdef WITH_SHTNS
      !$omp parallel do default(shared) &
      !$omp& private(n_theta, n_phi, fac, facCot, n_theta_nhs)
#endif
      do n_theta=nThetaStart,nThetaStop
         n_theta_nhs=(n_theta+1)/2
         fac=osn2(n_theta_nhs)
         facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
         if ( mod(n_theta,2) == 0 ) facCot=-facCot  ! SHS

         do n_phi=1,n_phi_max
            BtVr(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vr(n_phi,n_theta)
            BpVr(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vr(n_phi,n_theta)
         end do
#ifndef WITH_SHTNS
         BtVr(n_phi_max+1,n_theta)=0.0_cp
         BtVr(n_phi_max+2,n_theta)=0.0_cp
         BpVr(n_phi_max+1,n_theta)=0.0_cp
         BpVr(n_phi_max+2,n_theta)=0.0_cp
#endif

         do n_phi=1,n_phi_max
            BrVt(n_phi,n_theta)= fac*orho1(nR)*vt(n_phi,n_theta)*br(n_phi,n_theta)
            BrVp(n_phi,n_theta)= fac*orho1(nR)*vp(n_phi,n_theta)*br(n_phi,n_theta)
         end do
#ifndef WITH_SHTNS
         BrVt(n_phi_max+1,n_theta)=0.0_cp
         BrVt(n_phi_max+2,n_theta)=0.0_cp
         BrVp(n_phi_max+1,n_theta)=0.0_cp
         BrVp(n_phi_max+2,n_theta)=0.0_cp
#endif

         vpAS=0.0_cp
         do n_phi=1,n_phi_max
            BtVp(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVt(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            vpAS=vpAS+orho1(nR)*vp(n_phi,n_theta)
         end do
#ifndef WITH_SHTNS
         BtVp(n_phi_max+1,n_theta)=0.0_cp
         BtVp(n_phi_max+2,n_theta)=0.0_cp
         BpVt(n_phi_max+1,n_theta)=0.0_cp
         BpVt(n_phi_max+2,n_theta)=0.0_cp
#endif
         vpAS=vpAS/real(n_phi_max,kind=cp)

         !---- For toroidal terms that cancel:
         do n_phi=1,n_phi_max
            BpVtCot(n_phi,n_theta)=facCot*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpCot(n_phi,n_theta)=facCot*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVtSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
         end do
#ifndef WITH_SHTNS
         BpVtCot(n_phi_max+1,n_theta)=0.0_cp
         BpVtCot(n_phi_max+2,n_theta)=0.0_cp
         BtVpCot(n_phi_max+1,n_theta)=0.0_cp
         BtVpCot(n_phi_max+2,n_theta)=0.0_cp
         BpVtSn2(n_phi_max+1,n_theta)=0.0_cp
         BpVtSn2(n_phi_max+2,n_theta)=0.0_cp
         BtVpSn2(n_phi_max+1,n_theta)=0.0_cp
         BtVpSn2(n_phi_max+2,n_theta)=0.0_cp
#endif

         !---- For omega effect:
         do n_phi=1,n_phi_max
            BrVZ(n_phi,n_theta)=fac*br(n_phi,n_theta)*vpAS
            BtVZ(n_phi,n_theta)=fac*bt(n_phi,n_theta)*vpAS
            BtVZcot(n_phi,n_theta)=facCot*bt(n_phi,n_theta)*vpAS
            BtVZsn2(n_phi,n_theta)=fac*fac*bt(n_phi,n_theta)*vpAS
         end do
#ifndef WITH_SHTNS
         BrVZ(n_phi_max+1,n_theta)=0.0_cp
         BrVZ(n_phi_max+2,n_theta)=0.0_cp
         BtVZ(n_phi_max+1,n_theta)=0.0_cp
         BtVZ(n_phi_max+2,n_theta)=0.0_cp
         BtVZCot(n_phi_max+1,n_theta)=0.0_cp
         BtVZCot(n_phi_max+2,n_theta)=0.0_cp
         BtVZsn2(n_phi_max+1,n_theta)=0.0_cp
         BtVZsn2(n_phi_max+2,n_theta)=0.0_cp
#endif

      end do
#ifdef WITH_SHTNS
      !$omp end parallel do
#endif

#ifdef WITH_SHTNS
      call shtns_load_cfg(1)

      call spat_to_SH_dist(BtVr, BtVrLM, l_R(nR))
      call spat_to_SH_dist(BpVr, BpVrLM, l_R(nR))
      call spat_to_SH_dist(BrVt, BrVtLM, l_R(nR))

      call spat_to_SH_dist(BrVp, BrVpLM, l_R(nR))
      call spat_to_SH_dist(BtVp, BtVpLM, l_R(nR))
      call spat_to_SH_dist(BpVt, BpVtLM, l_R(nR))

      call spat_to_SH_dist(BtVpCot, BtVpCotLM, l_R(nR))
      call spat_to_SH_dist(BpVtCot, BpVtCotLM, l_R(nR))
      call spat_to_SH_dist(BtVZCot, BtVZCotLM, l_R(nR))

      call spat_to_SH_dist(BrVZ, BrVZLM, l_R(nR))
      call spat_to_SH_dist(BtVZ, BtVZLM, l_R(nR))
      call spat_to_SH_dist(BtVZsn2, BtVZsn2LM, l_R(nR))

      call spat_to_SH_dist(BtVpSn2, BtVpSn2LM, l_R(nR))
      call spat_to_SH_dist(BpVtsn2, BpVtsn2LM, l_R(nR))

      call shtns_load_cfg(0)
#else
      !-- Fourier transform phi2m
      if ( .not. l_axi ) then
         call fft_thetab(BtVpSn2,-1)
         call fft_thetab(BpVtSn2,-1)
         call fft_thetab(BtVpCot,-1)
         call fft_thetab(BpVtCot,-1)
         call fft_thetab(BtVr,-1)
         call fft_thetab(BpVr,-1)
         call fft_thetab(BrVt,-1)
         call fft_thetab(BrVp,-1)
         call fft_thetab(BtVp,-1)
         call fft_thetab(BpVt,-1)
         call fft_thetab(BrVZ,-1)
         call fft_thetab(BtVZ,-1)
         call fft_thetab(BtVZcot,-1)
         call fft_thetab(BtVZsn2,-1)
      end if

      !-- Legendre transform: theta2l
      call legTF3(n_theta_start,BtVrLM,BpVrLM,BrVtLM,BtVr,BpVr,BrVt)
      call legTF3(n_theta_start,BrVpLM,BtVpLM,BpVtLM,BrVp,BtVp,BpVt)
      call legTF3(n_theta_start,BtVpCotLM,BpVtCotLM,BtVZcotLM,BtVpCot,BpVtCot,BtVZcot)
      call legTF3(n_theta_start,BrVZLM,BtVZLM,BtVZsn2LM,BrVZ,BtVZ,BtVZsn2)
      call legTF2(n_theta_start,BtVpSn2LM,BpVtSn2LM,BtVpSn2,BpVtSn2)

      do l=0,l_max
         lm=l2lmAS(l)
         BtVrLM(lm)   =cmplx(real(BtVrLM(lm)),   0.0_cp, kind=cp)
         BpVrLM(lm)   =cmplx(real(BpVrLM(lm)),   0.0_cp, kind=cp)
         BrVtLM(lm)   =cmplx(real(BrVtLM(lm)),   0.0_cp, kind=cp)
         BrVpLM(lm)   =cmplx(real(BrVpLM(lm)),   0.0_cp, kind=cp)
         BtVpLM(lm)   =cmplx(real(BtVpLM(lm)),   0.0_cp, kind=cp)
         BpVtLM(lm)   =cmplx(real(BpVtLM(lm)),   0.0_cp, kind=cp)
         BtVpCotLM(lm)=cmplx(real(BtVpCotLM(lm)),0.0_cp, kind=cp)
         BpVtCotLM(lm)=cmplx(real(BpVtCotLM(lm)),0.0_cp, kind=cp)
         BtVZCotLM(lm)=cmplx(real(BtVZCotLM(lm)),0.0_cp, kind=cp)
         BrVZLM(lm)   =cmplx(real(BrVZLM(lm))   ,0.0_cp, kind=cp)
         BtVZLM(lm)   =cmplx(real(BtVZLM(lm))   ,0.0_cp, kind=cp)
         BtVZSn2LM(lm)=cmplx(real(BtVZSn2LM(lm)),0.0_cp, kind=cp)
         BtVpSn2LM(lm)=cmplx(real(BtVpSn2LM(lm)),0.0_cp, kind=cp)
         BpVtSn2LM(lm)=cmplx(real(BpVtSn2LM(lm)),0.0_cp, kind=cp)
      end do
#endif

   end subroutine get_dtBLM
!-----------------------------------------------------------------------
   subroutine get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
              &            BtVpLM,BpVtLM,BrVZLM,BtVZLM,    &
              &            BtVpCotLM,BpVtCotLM,            &
              &            BtVpSn2LM,BpVtSn2LM)
      !
      !  Purpose of this routine is to calculate theta and phi
      !  derivative related terms of the magnetic production and
      !  advection terms and store them.
      !

      !-- Input variables:
      integer,     intent(in) :: nR
      complex(cp), intent(in) :: BtVrLM(*),BpVrLM(*)
      complex(cp), intent(in) :: BrVtLM(*),BrVpLM(*)
      complex(cp), intent(in) :: BtVpLM(*),BpVtLM(*)
      complex(cp), intent(in) :: BtVpCotLM(*),BpVtCotLM(*)
      complex(cp), intent(in) :: BtVpSn2LM(*),BpVtSn2LM(*)
      complex(cp), intent(in) :: BrVZLM(*),BtVZLM(*)

      !-- Local variables:
      integer :: l,m,lm,lmP,lmPS,lmPA
      real(cp) :: fac

      PstrLM_Rloc(1,nR)=0.0_cp
      PadvLM_Rloc(1,nR)=0.0_cp
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         if ( l > m ) then
            PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   (                     &
            &    dTheta1S(lm)*BtVrLM(lmPS) - dTheta1A(lm)*BtVrLM(lmPA) + &
            &    dPhi(lm)*BpVrLM(lmP)  )
            PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   (                     &
            &    dTheta1S(lm)*BrVtLM(lmPS) - dTheta1A(lm)*BrVtLM(lmPA) + &
            &    dPhi(lm)*BrVpLM(lmP)  )
         else if ( l == m ) then
            PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
            &    - dTheta1A(lm)*BtVrLM(lmPA) + dPhi(lm)*BpVrLM(lmP)  )
            PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
            &    - dTheta1A(lm)*BrVtLM(lmPA) + dPhi(lm)*BrVpLM(lmP) )
         end if
      end do

      !--- Poloidal advection and stretching term finished for radial level nR !

      TstrLM_Rloc(1,nR) =0.0_cp
      TstrRLM_Rloc(1,nR)=0.0_cp
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
            &            fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
            &                                BpVtSn2LM(lmP) )   + &
            &                                             fac * ( &
            &             dTheta1S(lm) * ( or1(nR)*BpVrLM(lmPS) + &
            &                                   BpVtCotLM(lmPS) + &
            &                                 BtVpCotLM(lmPS) ) - &
            &             dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
            &                                   BpVtCotLM(lmPA) + &
            &                               BtVpCotLM(lmPA) ) ) - &
            &                  fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
            TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
            &                        dTheta1S(lm)*BrVpLM(lmPS) - &
            &                        dTheta1A(lm)*BrVpLM(lmPA) - &
            &                            dPhi(lm)*BrVtLM(lmP)  )
         else if ( l == m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
            &            fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
            &                                BpVtSn2LM(lmP) )   + &
            &                                             fac * ( &
            &           - dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
            &                                   BpVtCotLM(lmPA) + &
            &                               BtVpCotLM(lmPA) ) ) - &
            &                  fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
            TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
            &                      - dTheta1A(lm)*BrVpLM(lmPA) - &
            &                            dPhi(lm)*BrVtLM(lmP)  )
         end if
      end do

      TadvLM_Rloc(1,nR) =0.0_cp
      TadvRLM_Rloc(1,nR)=0.0_cp
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
            &           fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
            &                               BtVpSn2LM(lmP) )   + &
            &                                            fac * ( &
            &            dTheta1S(lm) * ( or1(nR)*BrVpLM(lmPS) + &
            &                                  BtVpCotLM(lmPS) + &
            &                                BpVtCotLM(lmPS) ) - &
            &            dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
            &                                  BtVpCotLM(lmPA) + &
            &                              BpVtCotLM(lmPA) ) ) - &
            &    fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
            TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
            &           dTheta1S(lm)*BpVrLM(lmPS) - &
            &           dTheta1A(lm)*BpVrLM(lmPA) - &
            &               dPhi(lm)*BtVrLM(lmP)   )
         else if ( l == m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
            &           fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
            &                               BtVpSn2LM(lmP) )   + &
            &                                            fac * ( &
            &          - dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
            &                                  BtVpCotLM(lmPA) + &
            &                              BpVtCotLM(lmPA) ) ) - &
            &                fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
            TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
            &         - dTheta1A(lm)*BpVrLM(lmPA) - &
            &               dPhi(lm)*BtVrLM(lmP)   )
         end if
      end do

      !--- TomeLM same as TstrLM but where ever Vp appeared
      !    it is replaced by its axisymmetric contribution VZ:
      TomeLM_Rloc(1,nR) =0.0_cp
      TomeRLM_Rloc(1,nR)=0.0_cp
      do lm=2,lm_max
         l  =lm2l(lm)
         m  =lm2m(lm)
         lmP=lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmP)       - &
            &                                 fac*or1(nR)*(     &
            &                       dTheta1S(lm)*BrVZLM(lmPS) - &
            &                       dTheta1A(lm)*BrVZLM(lmPA) )
            TomeRLM_Rloc(lm,nR)=                      fac * ( &
            &                     dTheta1S(lm)*BrVZLM(lmPS) - &
            &                     dTheta1A(lm)*BrVZLM(lmPA) )
         else if ( l == m ) then
            TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmp)       + &
            &         fac*or1(nR)*dTheta1A(lm)*BrVZLM(lmPA)
            TomeRLM_Rloc(lm,nR)=-fac*dTheta1A(lm)*BrVZLM(lmPA)
         end if
      end do

   end subroutine get_dH_dtBLM
!------------------------------------------------------------------------------
   subroutine get_dtBLMfinish(time,n_time_step,omega_ic,         &
     &                        b,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic, &
     &                        aj_ic,dj_ic,ddj_ic,l_frame)

      !-- Input of variables:
      real(cp),    intent(in) :: time
      integer,     intent(in) :: n_time_step
      real(cp),    intent(in) :: omega_ic
      complex(cp), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: ddj(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      logical,     intent(in) :: l_frame

      !-- Local variables:
      integer :: nR
      real(cp) :: dLh
      complex(cp) :: work_LMloc(llmMag:ulmMag,n_r_max)
      integer :: l,m,lm


      !-- Bring some array from rLoc to LMloc
      call r2lo_dtB%transp_r2lm(dtB_Rloc_container, dtB_LMloc_container)

      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               PadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*b_ic(lm,nR)
               TadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*aj_ic(lm,nR)
               PdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
               &    two*real(l+1,cp)*O_r_ic(nR)*db_ic(lm,nR) )
               TdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
               &    two*real(l+1,cp)*O_r_ic(nR)*dj_ic(lm,nR) )
            end do
         end do
      end if

      do nR=1,n_r_max
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            dLh = real(l*(l+1),cp)
            PdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &                   (ddb(lm,nR)-dLh*or2(nR)*b(lm,nR))
            TdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &    ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - dLh*or2(nR)*aj(lm,nR) )
         end do
      end do

      call get_dr(TomeRLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,       &
           &      nocopy=.true.)

      do nR=1,n_r_max
         do lm=llm,ulm
            TomeLM_LMloc(lm,nR)=TomeLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do

      call get_dr(TstrRLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,       &
           &      nocopy=.true.)

      do nR=1,n_r_max
         do lm=llm,ulm
            TstrLM_LMloc(lm,nR)=TstrLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do

      call get_dr(TadvRLM_LMloc(llmMag:ulmMag,:),work_LMloc(llmMag:ulmMag,:), &
           &      ulmMag-llmMag+1,1,ulmMag-llmMag+1,n_r_max,rscheme_oc,       &
           &      nocopy=.true.)

      do nR=1,n_r_max
         do lm=llm,ulm
            TadvLM_LMloc(lm,nR)=TadvLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do

      if ( l_DTrMagSpec .and. n_time_step > 1 ) then
         call rBrSpec(time,PstrLM_LMLoc,PadvLMIC_LMloc,'rBrProSpec',.false.,lo_map)
         call rBrSpec(time,PadvLM_LMLoc,PadvLMIC_LMLoc,'rBrAdvSpec',.true.,lo_map)
         call rBrSpec(time,PdifLM_LMLoc,PdifLMIC_LMLoc,'rBrDifSpec',.true.,lo_map)
         do nR=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,nR)=PstrLM_LMloc(lm,nR)-PadvLM_LMloc(lm,nR)
            end do
         end do
         call rBrSpec(time,work_LMloc,PadvLMIC_LMloc,'rBrDynSpec',.false.,lo_map)

         call rBpSpec(time,TstrLM_LMloc,TadvLMIC_LMloc,'rBpProSpec',.false.,lo_map)
         call rBpSpec(time,TadvLM_LMloc,TadvLMIC_LMloc,'rBpAdvSpec',.true.,lo_map)
         call rBpSpec(time,TdifLM_LMloc,TdifLMIC_LMloc,'rBpDifSpec',.true.,lo_map)
         do nR=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,nR)=TstrLM_LMloc(lm,nR)-TadvLM_LMloc(lm,nR)
            end do
         end do
         call rBpSpec(time,work_LMloc,TadvLMIC_LMloc,'rBpDynSpec',.false.,lo_map)

      end if

      if ( l_dtBmovie .and. l_frame ) then
         !-- If movie is required, let's gather everything on coord_r 0
         call dtb_gather_lo_on_rank0()
      end if

   end subroutine get_dtBLMfinish
!------------------------------------------------------------------------------
end module dtB_mod
