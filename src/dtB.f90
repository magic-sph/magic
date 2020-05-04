module dtB_mod
   !
   !  This module contains magnetic field stretching and advection terms
   !  plus a separate omega-effect.
   !  It is used for movie output.
   !
   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use LMmapping, only: map_dist_st, map_mlo
   use fields, only: work_LMdist
   use truncation, only: nrp, n_r_maxMag, n_r_ic_maxMag, n_r_max, lm_max_dtB, &
       &                 n_r_max_dtB, n_r_ic_max_dtB, lm_max, n_cheb_max,     &
       &                 n_r_ic_max, n_phi_max, l_axi, n_mloMag_loc,          &
       &                 nRstart, nRstop, nThetaStart, nThetaStop, n_lm_loc,  &
       &                 n_lmP_loc
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC
   use mpi_transp, only: type_mpitransp
   use mpi_thetap_mod, only: type_mpisendrecv
   use physical_parameters, only: opm, O_sr
   use radial_functions, only: O_r_ic, lambda, or2, dLlambda, rscheme_oc, &
       &                       or1, orho1, l_R
   use horizontal_data, only: dLh_loc, hdif_B, osn2, cosn2, osn1, dPhi_loc, &
       &                      dTheta1S_loc, dTheta1A_loc
   use logic, only: l_cond_ic, l_DTrMagSpec, l_dtBmovie
   use constants, only: ci, two
   use radial_spectra ! rBrSpec, rBpSpec
   use shtns, only: spat_to_SH_dist
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
   complex(cp), allocatable, target :: dtB_LMdist_container(:,:,:)
   complex(cp), allocatable, target :: dtB_Rdist_container(:,:,:)

   !-- R-distributed arrays
   complex(cp), pointer :: PstrLM_Rdist(:,:), PadvLM_Rdist(:,:)
   complex(cp), pointer :: TomeRLM_Rdist(:,:), TomeLM_Rdist(:,:)
   complex(cp), pointer :: TstrRLM_Rdist(:,:), TstrLM_Rdist(:,:)
   complex(cp), pointer :: TadvRLM_Rdist(:,:), TadvLM_Rdist(:,:)

   !-- LM-distributed arrays
   complex(cp), public, pointer :: PstrLM_LMdist(:,:), PadvLM_LMdist(:,:)
   complex(cp), public, pointer :: TomeRLM_LMdist(:,:), TomeLM_LMdist(:,:)
   complex(cp), public, pointer :: TstrRLM_LMdist(:,:), TstrLM_LMdist(:,:)
   complex(cp), public, pointer :: TadvRLM_LMdist(:,:), TadvLM_LMdist(:,:)
   complex(cp), public, allocatable :: PdifLM_LMdist(:,:), TdifLM_LMdist(:,:)
   complex(cp), public, allocatable :: PadvLMIC_LMdist(:,:), PdifLMIC_LMdist(:,:)
   complex(cp), public, allocatable :: TadvLMIC_LMdist(:,:), TdifLMIC_LMdist(:,:)

   class(type_mpitransp), pointer :: r2lo_dtB_dist

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
         if ( l_master_rank ) then
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

         if ( l_master_rank ) then
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

      allocate( PdifLM_LMdist(n_mloMag_loc,n_r_max_dtB) )
      allocate( TdifLM_LMdist(n_mloMag_loc,n_r_max_dtB) )
      bytes_allocated = bytes_allocated+ &
                        2*n_mloMag_loc*n_r_max_dtB*SIZEOF_DEF_COMPLEX
      allocate( PadvLMIC_LMdist(n_mloMag_loc,n_r_ic_max_dtB) )
      allocate( PdifLMIC_LMdist(n_mloMag_loc,n_r_ic_max_dtB) )
      allocate( TadvLMIC_LMdist(n_mloMag_loc,n_r_ic_max_dtB) )
      allocate( TdifLMIC_LMdist(n_mloMag_loc,n_r_ic_max_dtB) )
      bytes_allocated = bytes_allocated+ &
      &                 4*n_mloMag_loc*n_r_ic_max_dtB*SIZEOF_DEF_COMPLEX

      allocate( dtB_Rdist_container(n_lm_loc,nRstart:nRstop,8) )
      TomeLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,1)
      TomeRLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,2)
      TstrLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,3)
      TstrRLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,4)
      TadvLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,5)
      TadvRLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,6)
      PstrLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,7)
      PadvLM_Rdist(1:,nRstart:) => dtB_Rdist_container(1:n_lm_loc,nRstart:nRstop,8)
      bytes_allocated = bytes_allocated+8*(nRstop-nRstart+1)*n_lm_loc* &
      &                 SIZEOF_DEF_COMPLEX

      allocate( dtB_LMdist_container(n_mloMag_loc,n_r_max_dtB,8) )
      TomeLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,1)
      TomeRLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,2)
      TstrLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,3)
      TstrRLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,4)
      TadvLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,5)
      TadvRLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,6)
      PstrLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,7)
      PadvLM_LMdist(1:,1:) => dtB_LMdist_container(1:n_mloMag_loc,1:n_r_max_dtB,8)
      bytes_allocated = bytes_allocated+8*n_mloMag_loc*n_r_max_dtB* &
      &                 SIZEOF_DEF_COMPLEX

      allocate ( type_mpisendrecv :: r2lo_dtB_dist )

      call r2lo_dtB_dist%create_comm(8)

   end subroutine initialize_dtB_mod
!----------------------------------------------------------------------------
   subroutine finalize_dtB_mod

      if ( l_dtBmovie ) then
         deallocate( PstrLM, PadvLM, TstrLM, TadvLM, TomeLM )
         deallocate( TdifLMIC, TadvLMIC, PdifLMIC, PadvLMIC, TdifLM, PdifLM )
      end if

      deallocate( PdifLM_LMdist, TdifLM_LMdist, PadvLMIC_LMdist, PdifLMIC_LMdist )
      deallocate( TadvLMIC_LMdist, TdifLMIC_LMdist )
      deallocate( dtB_Rdist_container, dtB_LMdist_container )

      call r2lo_dtB_dist%destroy_comm()

   end subroutine finalize_dtB_mod
!----------------------------------------------------------------------------
   subroutine dtb_gather_lo_on_rank0
      !
      ! MPI gather on rank0 for dtBmovie outputs.
      ! This routine should really be suppressed once the movie
      ! outputs have been improved
      !

      use useful, only: abortRun

      call abortRun('Not ported yet! Stop heere...!')

      call gather_all_from_lo_to_rank0(gt_OC,PstrLM_LMdist,PstrLM)
      call gather_all_from_lo_to_rank0(gt_OC,TstrLM_LMdist,TstrLM)
      call gather_all_from_lo_to_rank0(gt_OC,PadvLM_LMdist,PadvLM)
      call gather_all_from_lo_to_rank0(gt_OC,TadvLM_LMdist,TadvLM)
      call gather_all_from_lo_to_rank0(gt_OC,TomeLM_LMdist,TomeLM)
      call gather_all_from_lo_to_rank0(gt_OC,PdifLM_LMdist,PdifLM)
      call gather_all_from_lo_to_rank0(gt_OC,TdifLM_LMdist,TdifLM)
      call gather_all_from_lo_to_rank0(gt_IC,PadvLMIC_LMdist,PadvLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,TadvLMIC_LMdist,TadvLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,PdifLMIC_LMdist,PdifLMIC)
      call gather_all_from_lo_to_rank0(gt_IC,TdifLMIC_LMdist,TdifLMIC)

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
      integer :: n_theta,n_phi,n_theta_nhs
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

      !$omp parallel do default(shared) &
      !$omp& private(n_theta, n_phi, fac, facCot, n_theta_nhs)
      do n_theta=nThetaStart,nThetaStop
         n_theta_nhs=(n_theta+1)/2
         fac=osn2(n_theta_nhs)
         facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
         if ( mod(n_theta,2) == 0 ) facCot=-facCot  ! SHS

         do n_phi=1,n_phi_max
            BtVr(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vr(n_phi,n_theta)
            BpVr(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vr(n_phi,n_theta)
         end do

         do n_phi=1,n_phi_max
            BrVt(n_phi,n_theta)= fac*orho1(nR)*vt(n_phi,n_theta)*br(n_phi,n_theta)
            BrVp(n_phi,n_theta)= fac*orho1(nR)*vp(n_phi,n_theta)*br(n_phi,n_theta)
         end do

         vpAS=0.0_cp
         do n_phi=1,n_phi_max
            BtVp(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVt(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            vpAS=vpAS+orho1(nR)*vp(n_phi,n_theta)
         end do
         vpAS=vpAS/real(n_phi_max,kind=cp)

         !---- For toroidal terms that cancel:
         do n_phi=1,n_phi_max
            BpVtCot(n_phi,n_theta)=facCot*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpCot(n_phi,n_theta)=facCot*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVtSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
         end do

         !---- For omega effect:
         do n_phi=1,n_phi_max
            BrVZ(n_phi,n_theta)=fac*br(n_phi,n_theta)*vpAS
            BtVZ(n_phi,n_theta)=fac*bt(n_phi,n_theta)*vpAS
            BtVZcot(n_phi,n_theta)=facCot*bt(n_phi,n_theta)*vpAS
            BtVZsn2(n_phi,n_theta)=fac*fac*bt(n_phi,n_theta)*vpAS
         end do

      end do
      !$omp end parallel do

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

   end subroutine get_dtBLM
!-----------------------------------------------------------------------
   subroutine get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM,BtVpLM,BpVtLM, &
              &            BrVZLM,BtVZLM,BtVpCotLM,BpVtCotLM,BtVpSn2LM,  &
              &            BpVtSn2LM)
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
      integer :: lm_maybe_skip_first
      integer, pointer :: lm2l(:), lm2m(:), lm2lmP(:), lm2(:,:)
      integer, pointer :: lmP2lmPS(:), lmP2lmPA(:)

      lm2l(1:n_lm_loc) => map_dist_st%lm2l
      lm2m(1:n_lm_loc) => map_dist_st%lm2m
      lmP2lmPS(1:n_lmP_loc) => map_dist_st%lmP2lmPS
      lmP2lmPA(1:n_lmP_loc) => map_dist_st%lmP2lmPA
      lm2lmP(1:n_lm_loc) => map_dist_st%lm2lmP
      lm2(0:,0:) => map_dist_st%lm2

      lm_maybe_skip_first = 1
      if (map_dist_st%lm2(0,0) > 0) lm_maybe_skip_first = 2

      lm =lm2(0,0)   ! This is l=0,m=0
      if ( lm > 0 ) then
         PstrLM_Rdist(lm,nR)=0.0_cp
         PadvLM_Rdist(lm,nR)=0.0_cp
      end if
      do lm=lm_maybe_skip_first,n_lm_loc
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         if ( l > m ) then
            PstrLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) *   (                        &
            &    dTheta1S_loc(lm)*BtVrLM(lmPS) - dTheta1A_loc(lm)*BtVrLM(lmPA) + &
            &    dPhi_loc(lm)*BpVrLM(lmP)  )
            PadvLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) *   (                        &
            &    dTheta1S_loc(lm)*BrVtLM(lmPS) - dTheta1A_loc(lm)*BrVtLM(lmPA) + &
            &    dPhi_loc(lm)*BrVpLM(lmP)  )
         else if ( l == m ) then
            PstrLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) *   ( &
            &    - dTheta1A_loc(lm)*BtVrLM(lmPA) + dPhi_loc(lm)*BpVrLM(lmP)  )
            PadvLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) *   ( &
            &    - dTheta1A_loc(lm)*BrVtLM(lmPA) + dPhi_loc(lm)*BrVpLM(lmP) )
         end if
      end do

      !--- Poloidal advection and stretching term finished for radial level nR !

      lm =lm2(0,0)   ! This is l=0,m=0
      if ( lm > 0 ) then
         TstrLM_Rdist(lm,nR) =0.0_cp
         TstrRLM_Rdist(lm,nR)=0.0_cp
      end if
      do lm=lm_maybe_skip_first,n_lm_loc
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh_loc(lm)
         if ( l > m ) then
            TstrLM_Rdist(lm,nR)=        -or2(nR)*BtVpLM(lmP)     -      &
            &          fac*dPhi_loc(lm)*dPhi_loc(lm)*( BtVpSn2LM(lmP) + &
            &                                BpVtSn2LM(lmP) ) + fac * ( &
            &             dTheta1S_loc(lm) * ( or1(nR)*BpVrLM(lmPS) +   &
            &                                   BpVtCotLM(lmPS) +       &
            &                                 BtVpCotLM(lmPS) ) -       &
            &             dTheta1A_loc(lm) * ( or1(nR)*BpVrLM(lmPA) +   &
            &                                   BpVtCotLM(lmPA) +       &
            &                               BtVpCotLM(lmPA) ) ) -       &
            &                  fac*or1(nR)*dPhi_loc(lm)*BtVrLM(lmP)
            TstrRLM_Rdist(lm,nR)=            or1(nR)/dLh_loc(lm) * ( &
            &                        dTheta1S_loc(lm)*BrVpLM(lmPS) - &
            &                        dTheta1A_loc(lm)*BrVpLM(lmPA) - &
            &                            dPhi_loc(lm)*BrVtLM(lmP)  )
         else if ( l == m ) then
            TstrLM_Rdist(lm,nR)=        -or2(nR)*BtVpLM(lmP)     -      &
            &          fac*dPhi_loc(lm)*dPhi_loc(lm)*( BtVpSn2LM(lmP) + &
            &                                BpVtSn2LM(lmP) ) + fac * ( &
            &           - dTheta1A_loc(lm) * ( or1(nR)*BpVrLM(lmPA) +   &
            &                                   BpVtCotLM(lmPA) +       &
            &                               BtVpCotLM(lmPA) ) ) -       &
            &                  fac*or1(nR)*dPhi_loc(lm)*BtVrLM(lmP)
            TstrRLM_Rdist(lm,nR)=             or1(nR)/dLh_loc(lm) * ( &
            &                      - dTheta1A_loc(lm)*BrVpLM(lmPA) -  &
            &                            dPhi_loc(lm)*BrVtLM(lmP)  )
         end if
      end do

      lm =lm2(0,0)   ! This is l=0,m=0
      if ( lm > 0 ) then
         TadvLM_Rdist(lm,nR) =0.0_cp
         TadvRLM_Rdist(lm,nR)=0.0_cp
      end if
      do lm=lm_maybe_skip_first,n_lm_loc
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh_loc(lm)
         if ( l > m ) then
            TadvLM_Rdist(lm,nR)=       -or2(nR)*BpVtLM(lmP)     -      &
            &       fac*dPhi_loc(lm)*dPhi_loc(lm)*( BpVtSn2LM(lmP) +   &
            &                               BtVpSn2LM(lmP) ) + fac * ( &
            &            dTheta1S_loc(lm) * ( or1(nR)*BrVpLM(lmPS) +   &
            &                                  BtVpCotLM(lmPS) +       &
            &                                BpVtCotLM(lmPS) ) -       &
            &            dTheta1A_loc(lm) * ( or1(nR)*BrVpLM(lmPA) +   &
            &                                  BtVpCotLM(lmPA) +       &
            &                              BpVtCotLM(lmPA) ) ) -       &
            &    fac*or1(nR)*dPhi_loc(lm)*BrVtLM(lmP)
            TadvRLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) * ( &
            &           dTheta1S_loc(lm)*BpVrLM(lmPS) -  &
            &           dTheta1A_loc(lm)*BpVrLM(lmPA) -  &
            &               dPhi_loc(lm)*BtVrLM(lmP)   )
         else if ( l == m ) then
            TadvLM_Rdist(lm,nR)=       -or2(nR)*BpVtLM(lmP)     -        &
            &           fac*dPhi_loc(lm)*dPhi_loc(lm)*( BpVtSn2LM(lmP) + &
            &                               BtVpSn2LM(lmP) )  +  fac * ( &
            &          - dTheta1A_loc(lm) * ( or1(nR)*BrVpLM(lmPA) +     &
            &                                  BtVpCotLM(lmPA) +         &
            &                              BpVtCotLM(lmPA) ) ) -         &
            &                fac*or1(nR)*dPhi_loc(lm)*BrVtLM(lmP)
            TadvRLM_Rdist(lm,nR)=or2(nR)/dLh_loc(lm) * ( &
            &         - dTheta1A_loc(lm)*BpVrLM(lmPA) -  &
            &               dPhi_loc(lm)*BtVrLM(lmP)   )
         end if
      end do

      !--- TomeLM same as TstrLM but where ever Vp appeared
      !    it is replaced by its axisymmetric contribution VZ:
      lm =lm2(0,0)   ! This is l=0,m=0
      if ( lm > 0 ) then
         TomeLM_Rdist(lm,nR) =0.0_cp
         TomeRLM_Rdist(lm,nR)=0.0_cp
      end if
      do lm=lm_maybe_skip_first,n_lm_loc
         l  =lm2l(lm)
         m  =lm2m(lm)
         lmP=lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh_loc(lm)
         if ( l > m ) then
            TomeLM_Rdist(lm,nR)=  -or2(nR)*BtVZLM(lmP)-fac*or1(nR)*( &
            &                     dTheta1S_loc(lm)*BrVZLM(lmPS) -    &
            &                     dTheta1A_loc(lm)*BrVZLM(lmPA) )
            TomeRLM_Rdist(lm,nR)=                      fac * (    &
            &                     dTheta1S_loc(lm)*BrVZLM(lmPS) - &
            &                     dTheta1A_loc(lm)*BrVZLM(lmPA) )
         else if ( l == m ) then
            TomeLM_Rdist(lm,nR)=    -or2(nR)*BtVZLM(lmP)       + &
            &         fac*or1(nR)*dTheta1A_loc(lm)*BrVZLM(lmPA)
            TomeRLM_Rdist(lm,nR)=-fac*dTheta1A_loc(lm)*BrVZLM(lmPA)
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
      complex(cp), intent(in) :: b(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: ddb(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: aj(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: dj(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: ddj(n_mloMag_loc,n_r_maxMag)
      complex(cp), intent(in) :: b_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(in) :: db_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddb_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(in) :: dj_ic(n_mloMag_loc,n_r_ic_maxMag)
      complex(cp), intent(in) :: ddj_ic(n_mloMag_loc,n_r_ic_maxMag)
      logical,     intent(in) :: l_frame

      !-- Local variables:
      integer :: nR,l,m,lm
      real(cp) :: dLh

      !-- Bring some array from rLoc to LMdist
      call r2lo_dtB_dist%transp_r2lm_dist(dtB_Rdist_container, dtB_LMdist_container)

      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            do lm=1,n_mloMag_loc
               l=map_mlo%i2l(lm)
               m=map_mlo%i2m(lm)
               PadvLMIC_LMdist(lm,nR)=-omega_ic*ci*m*b_ic(lm,nR)
               TadvLMIC_LMdist(lm,nR)=-omega_ic*ci*m*aj_ic(lm,nR)
               PdifLMIC_LMdist(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
               &    two*real(l+1,cp)*O_r_ic(nR)*db_ic(lm,nR) )
               TdifLMIC_LMdist(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
               &    two*real(l+1,cp)*O_r_ic(nR)*dj_ic(lm,nR) )
            end do
         end do
      end if

      do nR=1,n_r_max
         do lm=1,n_mloMag_loc
            l=map_mlo%i2l(lm)
            m=map_mlo%i2m(lm)
            dLh = real(l*(l+1),cp)
            PdifLM_LMdist(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &                   (ddb(lm,nR)-dLh*or2(nR)*b(lm,nR))
            TdifLM_LMdist(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &    ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - dLh*or2(nR)*aj(lm,nR) )
         end do
      end do

      call get_dr(TomeRLM_LMdist, work_LMdist, n_mloMag_loc, 1, n_mloMag_loc, &
           &      n_r_max, rscheme_oc,nocopy=.true.)

      do nR=1,n_r_max
         do lm=1,n_mloMag_loc
            TomeLM_LMdist(lm,nR)=TomeLM_LMdist(lm,nR)+or1(nR)*work_LMdist(lm,nR)
         end do
      end do

      call get_dr(TstrRLM_LMdist, work_LMdist, n_mloMag_loc, 1, n_mloMag_loc, &
           &      n_r_max, rscheme_oc, nocopy=.true.)

      do nR=1,n_r_max
         do lm=1,n_mloMag_loc
            TstrLM_LMdist(lm,nR)=TstrLM_LMdist(lm,nR)+or1(nR)*work_LMdist(lm,nR)
         end do
      end do

      call get_dr(TadvRLM_LMdist, work_LMdist, n_mloMag_loc, 1, n_mloMag_loc, &
           &      n_r_max, rscheme_oc, nocopy=.true.)

      do nR=1,n_r_max
         do lm=1,n_mloMag_loc
            TadvLM_LMdist(lm,nR)=TadvLM_LMdist(lm,nR)+or1(nR)*work_LMdist(lm,nR)
         end do
      end do

      if ( l_DTrMagSpec .and. n_time_step > 1 ) then
         call rBrSpec(time,PstrLM_LMdist,PadvLMIC_LMdist,'rBrProSpec',.false.)
         call rBrSpec(time,PadvLM_LMdist,PadvLMIC_LMdist,'rBrAdvSpec',.true.)
         call rBrSpec(time,PdifLM_LMdist,PdifLMIC_LMdist,'rBrDifSpec',.true.)
         do nR=1,n_r_max
            do lm=1,n_mloMag_loc
               work_LMdist(lm,nR)=PstrLM_LMdist(lm,nR)-PadvLM_LMdist(lm,nR)
            end do
         end do
         call rBrSpec(time,work_LMdist,PadvLMIC_LMdist,'rBrDynSpec',.false.)

         call rBpSpec(time,TstrLM_LMdist,TadvLMIC_LMdist,'rBpProSpec',.false.)
         call rBpSpec(time,TadvLM_LMdist,TadvLMIC_LMdist,'rBpAdvSpec',.true.)
         call rBpSpec(time,TdifLM_LMdist,TdifLMIC_LMdist,'rBpDifSpec',.true.)
         do nR=1,n_r_max
            do lm=1,n_mloMag_loc
               work_LMdist(lm,nR)=TstrLM_LMdist(lm,nR)-TadvLM_LMdist(lm,nR)
            end do
         end do
         call rBpSpec(time,work_LMdist,TadvLMIC_LMdist,'rBpDynSpec',.false.)

      end if

      if ( l_dtBmovie .and. l_frame ) then
         !-- If movie is required, let's gather everything on coord_r 0
         call dtb_gather_lo_on_rank0()
      end if

   end subroutine get_dtBLMfinish
!------------------------------------------------------------------------------
end module dtB_mod
