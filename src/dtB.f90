module dtB_mod
   !
   !  This module contains magnetic field stretching and advection terms
   !  plus a separate omega-effect.
   !  It is used for movie output.
   !
   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_maxMag, n_r_ic_maxMag, n_r_max, lm_max_dtB,   &
       &                 n_r_max_dtB, n_r_ic_max_dtB, lm_max, n_cheb_max,  &
       &                 n_r_ic_max, l_max, n_phi_max, ldtBmem, l_axi,     &
       &                 n_theta_max, nlat_padded
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC
   use mpi_transp_mod, only: type_mpitransp
   use mpi_ptop_mod, only: type_mpiptop
   use physical_parameters, only: opm,O_sr
   use radial_functions, only: O_r_ic, lambda, or2, dLlambda, rscheme_oc, &
       &                       or1, orho1, or3
   use radial_data,only: nRstart,nRstop
   use horizontal_data, only: dPhi, dLh, hdif_B, osn2, cosn2, osn1, &
       &                      dTheta1S, dTheta1A
   use logic, only: l_cond_ic, l_DTrMagSpec, l_dtBmovie
   use blocking, only: lo_map, st_map, lm2l, lm2m, lmP2lmPS, lmP2lmPA, &
                       lm2lmP, llmMag, ulmMag, llm, ulm
   use radial_spectra ! rBrSpec, rBpSpec
   use sht, only: scal_to_SH, spat_to_sphertor
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
         if ( rank == 0 ) then
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

         if ( rank == 0 ) then
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
      &                 2*(ulmMag-llmMag+1)*n_r_max_dtB*SIZEOF_DEF_COMPLEX
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
      !
      ! Memory deallocation
      !

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
               &         BtVpLM,BpVtLM,BrVZLM,BtVZLM,BpVtBtVpCotLM,        &
               &         BpVtBtVpSn2LM,BtVZsn2LM)
      !
      !  This subroutine calculates non-linear products in grid-space for radial
      !  level nR.
      !

      !-- Input variables:
      integer,  intent(in) :: nR
      real(cp), intent(in) :: vr(:,:),vt(:,:),vp(:,:)
      real(cp), intent(in) :: br(:,:),bt(:,:),bp(:,:)

      !-- Output variables:
      complex(cp), intent(out) :: BtVrLM(:),BpVrLM(:)
      complex(cp), intent(out) :: BrVtLM(:),BrVpLM(:)
      complex(cp), intent(out) :: BtVpLM(:),BpVtLM(:)
      complex(cp), intent(out) :: BrVZLM(:),BtVZLM(:)
      complex(cp), intent(out) :: BpVtBtVpCotLM(:)
      complex(cp), intent(out) :: BpVtBtVpSn2LM(:)
      complex(cp), intent(out) :: BtVZsn2LM(:)

      !-- Local variables:
      integer :: n_theta,n_phi,n_theta_nhs
      real(cp) :: fac,facCot
      real(cp) :: BtVr(nlat_padded,n_phi_max),BpVr(nlat_padded,n_phi_max)
      real(cp) :: BrVt(nlat_padded,n_phi_max),BrVp(nlat_padded,n_phi_max)
      real(cp) :: BtVp(nlat_padded,n_phi_max),BpVt(nlat_padded,n_phi_max)
      real(cp) :: BrVZ(nlat_padded,n_phi_max),BtVZ(nlat_padded,n_phi_max)
      real(cp) :: BpVtBtVpCot(nlat_padded,n_phi_max)
      real(cp) :: BpVtBtVpSn2(nlat_padded,n_phi_max)
      real(cp) :: BtVZsn2(nlat_padded,n_phi_max)
      real(cp) :: vpAS(n_theta_max)

      vpAS(:)=0.0_cp
      !$omp parallel do default(shared) &
      !$omp& private(n_theta, n_phi, fac, facCot, n_theta_nhs) &
      !$omp& reduction(+:vpAS)
      do n_phi=1,n_phi_max
         do n_theta=1,n_theta_max ! loop over ic-points, alternating north/south
            n_theta_nhs=(n_theta+1)/2
            fac=osn2(n_theta_nhs)
            facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
            if ( mod(n_theta,2) == 0 ) facCot=-facCot  ! SHS

            BtVr(n_theta,n_phi)= orho1(nR)*bt(n_theta,n_phi)*vr(n_theta,n_phi)
            BpVr(n_theta,n_phi)= orho1(nR)*bp(n_theta,n_phi)*vr(n_theta,n_phi)

            BrVt(n_theta,n_phi)= orho1(nR)*vt(n_theta,n_phi)*br(n_theta,n_phi)
            BrVp(n_theta,n_phi)= orho1(nR)*vp(n_theta,n_phi)*br(n_theta,n_phi)

            BtVp(n_theta,n_phi)= fac*orho1(nR)*bt(n_theta,n_phi)*vp(n_theta,n_phi)
            BpVt(n_theta,n_phi)= fac*orho1(nR)*bp(n_theta,n_phi)*vt(n_theta,n_phi)

            BpVtBtVpCot(n_theta,n_phi)=facCot*orho1(nR)*(                        &
            &                             bp(n_theta,n_phi)*vt(n_theta,n_phi) +  &
            &                             bt(n_theta,n_phi)*vp(n_theta,n_phi) )
            BpVtBtVpSn2(n_theta,n_phi)=fac*fac*orho1(nR)*(                       &
            &                             bp(n_theta,n_phi)*vt(n_theta,n_phi) +  &
            &                             bt(n_theta,n_phi)*vp(n_theta,n_phi) )
            vpAS(n_theta)=vpAS(n_theta)+orho1(nR)*vp(n_theta,n_phi)
         end do
      end do
      !$omp end parallel do
      vpAS(:)=vpAS(:)/real(n_phi_max,kind=cp)

      !---- For omega effect:
      !$omp parallel do default(shared) &
      !$omp private(n_phi,n_theta,n_theta_nhs,fac,facCot)
      do n_phi=1,n_phi_max
         do n_theta=1,n_theta_max ! loop over ic-points, alternating north/south
            n_theta_nhs=(n_theta+1)/2
            fac=osn2(n_theta_nhs)
            facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
            if ( mod(n_theta,2) == 0 ) facCot=-facCot  ! SHS
            BrVZ(n_theta,n_phi)=fac*br(n_theta,n_phi)*vpAS(n_theta)
            BtVZ(n_theta,n_phi)=fac*bt(n_theta,n_phi)*vpAS(n_theta)
            BtVZsn2(n_theta,n_phi)=fac*fac*bt(n_theta,n_phi)*vpAS(n_theta)
         end do
      end do
      !$omp end parallel do

      call spat_to_sphertor(BtVr, BpVr, BtVrLM, BpVrLM, l_max)
      call spat_to_sphertor(BrVt, BrVp, BrVtLM, BrVpLM, l_max)

      call scal_to_SH(BtVp, BtVpLM, l_max)
      call scal_to_SH(BpVt, BpVtLM, l_max)

      call scal_to_SH(BpVtBtVpCot, BpVtBtVpCotLM, l_max)
      call scal_to_SH(BpVtBtVpsn2, BpVtBtVpsn2LM, l_max)

      call scal_to_SH(BrVZ, BrVZLM, l_max)
      call scal_to_SH(BtVZ, BtVZLM, l_max)
      call scal_to_SH(BtVZsn2, BtVZsn2LM, l_max)

   end subroutine get_dtBLM
!-----------------------------------------------------------------------
   subroutine get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM,BtVpLM,BpVtLM, &
              &            BrVZLM,BtVZLM,BpVtBtVpCotLM,BpVtBtVpSn2LM)
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
      complex(cp), intent(in) :: BpVtBtVpCotLM(*)
      complex(cp), intent(in) :: BpVtBtVpSn2LM(*)
      complex(cp), intent(in) :: BrVZLM(*),BtVZLM(*)

      !-- Local variables:
      integer :: l,m,lm,lmP,lmPS,lmPA
      real(cp) :: fac

      !$omp parallel default(shared) private(lm,l,m,lmP,lmPS,lmPA,fac)

      PstrLM_Rloc(1,nR)=0.0_cp
      PadvLM_Rloc(1,nR)=0.0_cp
      !$omp do
      do lm=2,lm_max
         lmP =lm2lmP(lm)
         PstrLM_Rloc(lm,nR)=or2(nR) * BtVrLM(lmP)
         PadvLM_Rloc(lm,nR)=or2(nR) * BrVtLM(lmP)
      end do
      !$omp end do

      !--- Poloidal advection and stretching term finished for radial level nR !

      TstrLM_Rloc(1,nR) =0.0_cp
      TstrRLM_Rloc(1,nR)=0.0_cp
      !$omp do
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         TstrRLM_Rloc(lm,nR)=or1(nR) * BrVpLM(lmP)
         if ( l > m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
            &          fac*dPhi(lm)*dPhi(lm)*BpVtBtVpSn2LM(lmP) + &
            &                             or3(nR)* BpVrLM(lmP)  + &
            &                                             fac * ( &
            &             dTheta1S(lm) * BpVtBtVpCotLM(lmPS)  - &
            &             dTheta1A(lm) * BpVtBtVpCotLM(lmPA) )
         else if ( l == m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
            &          fac*dPhi(lm)*dPhi(lm)*BpVtBtVpSn2LM(lmP) + &
            &                             or3(nR)* BpVrLM(lmP)  + &
            &                                             fac * ( &
            &              - dTheta1A(lm) * BpVtBtVpCotLM(lmPA) )
         end if
      end do
      !$omp end do

      TadvLM_Rloc(1,nR) =0.0_cp
      TadvRLM_Rloc(1,nR)=0.0_cp
      !$omp do
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         TadvRLM_Rloc(lm,nR)=or2(nR) * BpVrLM(lmP)
         if ( l > m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
            &         fac*dPhi(lm)*dPhi(lm)*BpVtBtVpSn2LM(lmP) + &
            &                            or3(nR) * BrVpLM(lmP) + &
            &                                            fac * ( &
            &            dTheta1S(lm) * BpVtBtVpCotLM(lmPS)   -  &
            &            dTheta1A(lm) * BpVtBtVpCotLM(lmPA) )
         else if ( l == m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
            &         fac*dPhi(lm)*dPhi(lm)*BpVtBtVpSn2LM(lmP) + &
            &                            or3(nR) * BrVpLM(lmP) + &
            &                                            fac * ( &
            &            - dTheta1A(lm) *  BpVtBtVpCotLM(lmPA) ) 
         end if
      end do
      !$omp end do

      !--- TomeLM same as TstrLM but where ever Vp appeared
      !    it is replaced by its axisymmetric contribution VZ:
      TomeLM_Rloc(1,nR) =0.0_cp
      TomeRLM_Rloc(1,nR)=0.0_cp
      !$omp do
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
      !$omp end do
      !$omp end parallel

   end subroutine get_dH_dtBLM
!------------------------------------------------------------------------------
   subroutine get_dtBLMfinish(time,n_time_step,omega_ic,b,ddb,aj,dj,ddj,b_ic, &
              &               db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic,l_frame)

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
      integer :: nR, start_lm, stop_lm, l, m, lm
      real(cp) :: dL
      complex(cp) :: work_LMloc(llmMag:ulmMag,n_r_max)

      !-- Bring some array from rLoc to LMloc
      call r2lo_dtB%transp_r2lm(dtB_Rloc_container, dtB_LMloc_container)

      !$omp parallel default(shared) private(nR, lm, start_lm, stop_lm, l, m, dL)
      start_lm=llmMag; stop_lm=ulmMag
      call get_openmp_blocks(start_lm,stop_lm)

      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            do lm=start_lm,stop_lm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               PadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*b_ic(lm,nR)
               TadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*aj_ic(lm,nR)
               PdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
               &                     two*real(l+1,cp)*O_r_ic(nR)*db_ic(lm,nR) )
               TdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
               &                     two*real(l+1,cp)*O_r_ic(nR)*dj_ic(lm,nR) )
            end do
         end do
      end if

      do nR=1,n_r_max
         do lm=start_lm, stop_lm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            dL = real(l*(l+1),cp)
            PdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &                   (ddb(lm,nR)-dL*or2(nR)*b(lm,nR))
            TdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(l) * &
            &    ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - dL*or2(nR)*aj(lm,nR) )
         end do
      end do

      call get_dr(TomeRLM_LMloc, work_LMloc, ulmMag-llmMag+1, start_lm-llmMag+1, &
           &      stop_lm-llmMag+1, n_r_max, rscheme_oc, nocopy=.true.)
      !$omp barrier

      do nR=1,n_r_max
         do lm=start_lm, stop_lm
            TomeLM_LMloc(lm,nR)=TomeLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do

      call get_dr(TstrRLM_LMloc, work_LMloc, ulmMag-llmMag+1, start_lm-llmMag+1, &
           &      stop_lm-llmMag+1, n_r_max, rscheme_oc, nocopy=.true.)
      !$omp barrier

      do nR=1,n_r_max
         do lm=start_lm, stop_lm
            TstrLM_LMloc(lm,nR)=TstrLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do

      call get_dr(TadvRLM_LMloc, work_LMloc, ulmMag-llmMag+1, start_lm-llmMag+1, &
           &      stop_lm-llmMag+1, n_r_max, rscheme_oc, nocopy=.true.)
      !$omp barrier

      do nR=1,n_r_max
         do lm=start_lm,stop_lm
            TadvLM_LMloc(lm,nR)=TadvLM_LMloc(lm,nR)+or1(nR)*work_LMloc(lm,nR)
         end do
      end do
      !$omp end parallel

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
         !-- If movie is required, let's gather everything on rank 0
         call dtb_gather_lo_on_rank0()
      end if

   end subroutine get_dtBLMfinish
!------------------------------------------------------------------------------
end module dtB_mod
