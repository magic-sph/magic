!$Id$
module dtB_mod
   !---------------------------------------------------------------------------------
   !  This module contains magnetic field stretching and advection terms 
   !  plus a separate omega-effect.
   !  It is used for movie output.
   !--------------------------------------------------------------------------------
   use truncation, only: nrp, n_r_maxMag, n_r_ic_maxMag, n_r_max, lm_max_dtB, &
                         n_r_max_dtB, n_r_ic_max_dtB, lm_max, n_cheb_max,     &
                         n_r_ic_max, l_max, n_phi_max, ldtBmem, lm_max_real
   use communications, only: gather_all_from_lo_to_rank0, gt_OC, gt_IC
   use physical_parameters, only: opm,O_sr
   use radial_functions, only: O_r_ic, lambda, or2, dLlambda, i_costf_init, &
                               d_costf_init, drx, or1, orho1
   use radial_data,only: nRstart,nRstop
   use parallel_mod,only: nr_per_rank, nr_on_last_rank, MPI_COMM_WORLD, n_procs, &
                          MPI_DOUBLE_COMPLEX, rank
   use horizontal_data, only: dPhi, D_lP1, dLh, hdif_B, osn2, cosn2, osn1, &
                              dTheta1S, dTheta1A
   use logic, only: l_cond_ic, l_DTrMagSpec
   use LMLoop_data, only: llmMag, ulmMag, llm, ulm, llm_real, ulm_real
   use blocking, only: lo_map, st_map, l2lmAS, lm2l, lm2m, lmP2lmPS, lmP2lmPA, &
                       lm2lmP
   use radial_spectra ! rBrSpec, rBpSpec
#if (FFTLIB==JW)
   use fft_JW
#elif (FFTLIB==MKL)
   use fft_MKL
#endif
 
   implicit none
 
   private 
 
   complex(kind=8), public, allocatable :: PstrLM(:,:), TstrLM(:,:), PadvLM(:,:)
   complex(kind=8), public, allocatable :: TadvLM(:,:), TomeLM(:,:)
   complex(kind=8), public, allocatable :: PstrLM_Rloc(:,:),TstrLM_Rloc(:,:)
   complex(kind=8), public, allocatable :: PadvLM_Rloc(:,:), TadvLM_Rloc(:,:)
   complex(kind=8), public, allocatable :: TomeLM_Rloc(:,:)
 
   complex(kind=8), public, allocatable :: PdifLM(:,:), TdifLM(:,:), PadvLMIC(:,:)
   complex(kind=8), public, allocatable :: PdifLMIC(:,:), TadvLMIC(:,:), TdifLMIC(:,:)
   complex(kind=8), public, allocatable :: PdifLM_LMloc(:,:), TdifLM_LMloc(:,:)
   complex(kind=8), public, allocatable :: PadvLMIC_LMloc(:,:), PdifLMIC_LMloc(:,:)
   complex(kind=8), public, allocatable :: TadvLMIC_LMloc(:,:), TdifLMIC_LMloc(:,:)
 
   complex(kind=8), public, allocatable :: TstrRLM(:,:), TadvRLM(:,:), TomeRLM(:,:)
   complex(kind=8), public, allocatable :: TstrRLM_Rloc(:,:), TadvRLM_Rloc(:,:)
   complex(kind=8), public, allocatable :: TomeRLM_Rloc(:,:)
 
   real(kind=8), public :: PstrRms, PstrAsRms
   real(kind=8), public :: PadvRms, PadvAsRms
   real(kind=8), public :: PdifRms, PdifAsRms
   real(kind=8), public :: TstrRms, TstrAsRms
   real(kind=8), public :: TadvRms, TadvAsRms
   real(kind=8), public :: TdifRms, TdifAsRms
   real(kind=8), public :: TomeRms, TomeAsRms

   public :: initialize_dtB_mod, dtb_gather_Rloc_on_rank0, get_dtBLMfinish, &
             get_dtBLM, get_dH_dtBLM

contains

   subroutine initialize_dtB_mod

      if ( rank == 0 ) then
         allocate( PstrLM(lm_max_dtB,n_r_max_dtB) )
         allocate( PadvLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TstrLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TadvLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TomeLM(lm_max_dtB,n_r_max_dtB) )
      else
         allocate( PstrLM(1,1) )
         allocate( PadvLM(1,1) )
         allocate( TstrLM(1,1) )
         allocate( TadvLM(1,1) )
         allocate( TomeLM(1,1) )
      end if
      allocate( PstrLM_Rloc(lm_max_dtB,nRstart:nRstop) )
      allocate( PadvLM_Rloc(lm_max_dtB,nRstart:nRstop) )
      allocate( TstrLM_Rloc(lm_max_dtB,nRstart:nRstop) )
      allocate( TadvLM_Rloc(lm_max_dtB,nRstart:nRstop) )
      allocate( TomeLM_Rloc(lm_max_dtB,nRstart:nRstop) )

      if ( rank == 0 ) then
         allocate( PdifLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TdifLM(lm_max_dtB,n_r_max_dtB) )
         allocate( PadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
         allocate( PdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
         allocate( TadvLMIC(lm_max_dtB,n_r_ic_max_dtB) )
         allocate( TdifLMIC(lm_max_dtB,n_r_ic_max_dtB) )
      else
         allocate( PdifLM(1,1) )
         allocate( TdifLM(1,1) )
         allocate( PadvLMIC(1,1) )
         allocate( PdifLMIC(1,1) )
         allocate( TadvLMIC(1,1) )
         allocate( TdifLMIC(1,1) )
      end if
      allocate( PdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
      allocate( TdifLM_LMloc(llmMag:ulmMag,n_r_max_dtB) )
      allocate( PadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( PdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( TadvLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )
      allocate( TdifLMIC_LMloc(llmMag:ulmMag,n_r_ic_max_dtB) )

      if ( rank == 0 ) then
         allocate( TstrRLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TadvRLM(lm_max_dtB,n_r_max_dtB) )
         allocate( TomeRLM(lm_max_dtB,n_r_max_dtB) )
      else
         allocate( TstrRLM(1,1) )
         allocate( TadvRLM(1,1) )
         allocate( TomeRLM(1,1) )
      end if
      allocate(TstrRLM_Rloc(lm_max_dtB,nRstart:nRstop))
      allocate(TadvRLM_Rloc(lm_max_dtB,nRstart:nRstop))
      allocate(TomeRLM_Rloc(lm_max_dtB,nRstart:nRstop))

   end subroutine initialize_dtB_mod
!----------------------------------------------------------------------------
   subroutine dtb_gather_Rloc_on_rank0

      integer :: sendcount,recvcounts(0:n_procs-1),displs(0:n_procs-1)
      integer :: i,ierr
    
      if ( ldtBmem == 1 ) then
         sendcount  = (nRstop-nRstart+1)*lm_max_dtB
         recvcounts = nr_per_rank*lm_max_dtB
         recvcounts(n_procs-1) = nr_on_last_rank*lm_max_dtB
         do i=0,n_procs-1
            displs(i) = i*nr_per_rank*lm_max_dtB
         end do
         call MPI_GatherV(TstrRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TstrRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(TadvRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TadvRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(TomeRLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TomeRLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    
         call MPI_GatherV(TstrLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TstrLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(TadvLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TadvLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(PstrLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & PstrLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         call MPI_GatherV(PadvLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & PadvLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    
         call MPI_GatherV(TomeLM_Rloc,sendcount,MPI_DOUBLE_COMPLEX,&
              & TomeLM,recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      end if

   end subroutine dtb_gather_Rloc_on_rank0
!----------------------------------------------------------------------------
   subroutine  get_dtBLM(nR,vr,vt,vp,br,bt,bp,n_theta_start,n_theta_block, &
     &                   BtVrLM,BpVrLM,BrVtLM,BrVpLM,BtVpLM,BpVtLM,BrVZLM, &
     &                   BtVZLM,BtVpCotLM,BpVtCotLM,BtVZcotLM,BtVpSn2LM,   &
     &                   BpVtSn2LM,BtVZsn2LM)
      !-----------------------------------------------------------------------

      !  calculates non-linear products in grid-space for radial
      !  level n_r and returns them in arrays wnlr1-3, snlr1-3, bnlr1-3

      !  if lvelo >0 velocities are zero only the (vxB)
      !  contributions to bnlr2-3 need to be calculated

      !  vr...sr: (input) velocity, magnetic field comp. and derivs, entropy
      !                   on grid points
      !  i1: (input) range of points in theta for which calculation is done

      !-----------------------------------------------------------------------

      !-- Input variables:
      integer,      intent(in) :: n_theta_start,n_theta_block,nR
      real(kind=8), intent(in) :: vr(nrp,nfs),vt(nrp,nfs),vp(nrp,nfs)
      real(kind=8), intent(in) :: br(nrp,nfs),bt(nrp,nfs),bp(nrp,nfs)
    
      !-- Output variables:
      complex(kind=8), intent(out) :: BtVrLM(*),BpVrLM(*)
      complex(kind=8), intent(out) :: BrVtLM(*),BrVpLM(*)
      complex(kind=8), intent(out) :: BtVpLM(*),BpVtLM(*)
      complex(kind=8), intent(out) :: BrVZLM(*),BtVZLM(*)
      complex(kind=8), intent(out) :: BpVtCotLM(*),BtVpCotLM(*),BtVZcotLM(*)
      complex(kind=8), intent(out) :: BtVpSn2LM(*),BpVtSn2LM(*)
      complex(kind=8), intent(out) :: BtVZsn2LM(*)
    
      !-- Local variables:
      integer :: n_theta,n_phi,n_theta_nhs
      real(kind=8) :: fac,facCot
      real(kind=8) :: BtVr(nrp,nfs),BpVr(nrp,nfs)
      real(kind=8) :: BrVt(nrp,nfs),BrVp(nrp,nfs)
      real(kind=8) :: BtVp(nrp,nfs),BpVt(nrp,nfs)
      real(kind=8) :: BrVZ(nrp,nfs),BtVZ(nrp,nfs)
      real(kind=8) :: BpVtCot(nrp,nfs),BtVpCot(nrp,nfs)
      real(kind=8) :: BpVtSn2(nrp,nfs),BtVpSn2(nrp,nfs)
      real(kind=8) :: BtVZcot(nrp,nfs),BtVZsn2(nrp,nfs)
      real(kind=8) :: vpAS
    
      integer :: lm,l
    
      do n_theta=1,n_theta_block ! loop over ic-points, alternating north/south
    
         n_theta_nhs=(n_theta_start+n_theta)/2
         fac=osn2(n_theta_nhs)
         facCot=cosn2(n_theta_nhs)*osn1(n_theta_nhs)
         if ( mod(n_theta,2) == 0 ) facCot=-facCot  ! SHS
    
         do n_phi=1,n_phi_max
            BtVr(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vr(n_phi,n_theta)
            BpVr(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vr(n_phi,n_theta)
         end do
         BtVr(n_phi_max+1,n_theta)=0.d0
         BtVr(n_phi_max+2,n_theta)=0.d0
         BpVr(n_phi_max+1,n_theta)=0.d0
         BpVr(n_phi_max+2,n_theta)=0.d0
    
         do n_phi=1,n_phi_max
            BrVt(n_phi,n_theta)= fac*orho1(nR)*vt(n_phi,n_theta)*br(n_phi,n_theta)
            BrVp(n_phi,n_theta)= fac*orho1(nR)*vp(n_phi,n_theta)*br(n_phi,n_theta)
         end do
         BrVt(n_phi_max+1,n_theta)=0.d0
         BrVt(n_phi_max+2,n_theta)=0.d0
         BrVp(n_phi_max+1,n_theta)=0.d0
         BrVp(n_phi_max+2,n_theta)=0.d0
    
         vpAS=0.d0
         do n_phi=1,n_phi_max
            BtVp(n_phi,n_theta)= fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVt(n_phi,n_theta)= fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            vpAS=vpAS+orho1(nR)*vp(n_phi,n_theta)
         end do
         BtVp(n_phi_max+1,n_theta)=0.d0
         BtVp(n_phi_max+2,n_theta)=0.d0
         BpVt(n_phi_max+1,n_theta)=0.d0
         BpVt(n_phi_max+2,n_theta)=0.d0
         vpAS=vpAS/dble(n_phi_max)
    
         !---- For toroidal terms that cancel:
         do n_phi=1,n_phi_max
            BpVtCot(n_phi,n_theta)=facCot*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpCot(n_phi,n_theta)=facCot*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
            BpVtSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bp(n_phi,n_theta)*vt(n_phi,n_theta)
            BtVpSn2(n_phi,n_theta)=fac*fac*orho1(nR)*bt(n_phi,n_theta)*vp(n_phi,n_theta)
         end do
         BpVtCot(n_phi_max+1,n_theta)=0.d0
         BpVtCot(n_phi_max+2,n_theta)=0.d0
         BtVpCot(n_phi_max+1,n_theta)=0.d0
         BtVpCot(n_phi_max+2,n_theta)=0.d0
         BpVtSn2(n_phi_max+1,n_theta)=0.d0
         BpVtSn2(n_phi_max+2,n_theta)=0.d0
         BtVpSn2(n_phi_max+1,n_theta)=0.d0
         BtVpSn2(n_phi_max+2,n_theta)=0.d0
    
         !---- For omega effect:
         do n_phi=1,n_phi_max
            BrVZ(n_phi,n_theta)=fac*br(n_phi,n_theta)*vpAS
            BtVZ(n_phi,n_theta)=fac*bt(n_phi,n_theta)*vpAS
            BtVZcot(n_phi,n_theta)=facCot*bt(n_phi,n_theta)*vpAS
            BtVZsn2(n_phi,n_theta)=fac*fac*bt(n_phi,n_theta)*vpAS
         end do
         BrVZ(n_phi_max+1,n_theta)=0.d0
         BrVZ(n_phi_max+2,n_theta)=0.d0
         BtVZ(n_phi_max+1,n_theta)=0.d0
         BtVZ(n_phi_max+2,n_theta)=0.d0
         BtVZCot(n_phi_max+1,n_theta)=0.d0
         BtVZCot(n_phi_max+2,n_theta)=0.d0
         BtVZsn2(n_phi_max+1,n_theta)=0.d0
         BtVZsn2(n_phi_max+2,n_theta)=0.d0
    
      end do
    
      !-- Fourier transform phi2m
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
    
      !-- Legendre transform: theta2l
      call legTF3(n_theta_start,BtVrLM,BpVrLM,BrVtLM,BtVr,BpVr,BrVt)
      call legTF3(n_theta_start,BrVpLM,BtVpLM,BpVtLM,BrVp,BtVp,BpVt)
      call legTF3(n_theta_start,BtVpCotLM,BpVtCotLM,BtVZcotLM,BtVpCot,BpVtCot,BtVZcot)
      call legTF3(n_theta_start,BrVZLM,BtVZLM,BtVZsn2LM,BrVZ,BtVZ,BtVZsn2)
      call legTF2(n_theta_start,BtVpSn2LM,BpVtSn2LM,BtVpSn2,BpVtSn2)
    
      do l=0,l_max
         lm=l2lmAS(l)
         BtVrLM(lm)   =cmplx(real(BtVrLM(lm)),   0.D0,kind=kind(0.d0))
         BpVrLM(lm)   =cmplx(real(BpVrLM(lm)),   0.D0,kind=kind(0.d0))
         BrVtLM(lm)   =cmplx(real(BrVtLM(lm)),   0.D0,kind=kind(0.d0))
         BrVpLM(lm)   =cmplx(real(BrVpLM(lm)),   0.D0,kind=kind(0.d0))
         BtVpLM(lm)   =cmplx(real(BtVpLM(lm)),   0.D0,kind=kind(0.d0))
         BpVtLM(lm)   =cmplx(real(BpVtLM(lm)),   0.D0,kind=kind(0.d0))
         BtVpCotLM(lm)=cmplx(real(BtVpCotLM(lm)),0.D0,kind=kind(0.d0))
         BpVtCotLM(lm)=cmplx(real(BpVtCotLM(lm)),0.D0,kind=kind(0.d0))
         BtVZCotLM(lm)=cmplx(real(BtVZCotLM(lm)),0.D0,kind=kind(0.d0))
         BrVZLM(lm)   =cmplx(real(BrVZLM(lm))   ,0.D0,kind=kind(0.d0))
         BtVZLM(lm)   =cmplx(real(BtVZLM(lm))   ,0.D0,kind=kind(0.d0))
         BtVZSn2LM(lm)=cmplx(real(BtVZSn2LM(lm)),0.D0,kind=kind(0.d0))
         BtVpSn2LM(lm)=cmplx(real(BtVpSn2LM(lm)),0.D0,kind=kind(0.d0))
         BpVtSn2LM(lm)=cmplx(real(BpVtSn2LM(lm)),0.D0,kind=kind(0.d0))
      end do
    
   end subroutine get_dtBLM
!------------------------------------------------------------------------------------
   subroutine get_dtBLMfinish(time,n_time_step,omega_ic, &
     &                     b,ddb,aj,dj,ddj,b_ic,db_ic,ddb_ic, &
     &                     aj_ic,dj_ic,ddj_ic)

      !-- Input of variables:
      real(kind=8),    intent(in) :: time
      integer,         intent(in) :: n_time_step
      real(kind=8),    intent(in) :: omega_ic
      complex(kind=8), intent(in) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: ddj(llmMag:ulmMag,n_r_maxMag)
      complex(kind=8), intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: db_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: ddb_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: dj_ic(llmMag:ulmMag,n_r_ic_maxMag)
      complex(kind=8), intent(in) :: ddj_ic(llmMag:ulmMag,n_r_ic_maxMag)
    
      !-- Local variables:
      integer :: nR
    
      complex(kind=8) :: workA(lm_max,n_r_max), workB(lm_max,n_r_max)
    
      integer :: l,m,lm
      
      ! gathering TstrRLM,TadvRLM and TomeRLM on rank0,
      ! they are then in st_map order
      call dtB_gather_Rloc_on_rank0
    
      if ( l_cond_ic ) then
         do nR=1,n_r_ic_max
            do lm=llm,ulm
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               PadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*b_ic(lm,nR)
               TadvLMIC_LMloc(lm,nR)=-omega_ic*dPhi(st_map%lm2(l,m))*aj_ic(lm,nR)
               PdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddb_ic(lm,nR) + &
                    2.D0*D_lP1(st_map%lm2(l,m))*O_r_ic(nR)*db_ic(lm,nR) )
               TdifLMIC_LMloc(lm,nR)=opm*O_sr * ( ddj_ic(lm,nR) + &
                    2.D0*D_lP1(st_map%lm2(l,m))*O_r_ic(nR)*dj_ic(lm,nR) )
            end do
         end do
      end if
    
      do nR=1,n_r_max
         do lm=llm,ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            PdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(st_map%lm2(l,m)) * &
                 (ddb(lm,nR)-dLh(st_map%lm2(l,m))*or2(nR)*b(lm,nR))
            TdifLM_LMloc(lm,nR)= opm*lambda(nR)*hdif_B(st_map%lm2(l,m)) * &
                 ( ddj(lm,nR) + dLlambda(nR)*dj(lm,nR) - &
                 dLh(st_map%lm2(l,m))*or2(nR)*aj(lm,nR) )
         end do
      end do
    
      if ( rank == 0 ) then
         call get_drNS(TstrRLM,workA,lm_max_real,1,lm_max_real, &
                       n_r_max,n_cheb_max,workB,i_costf_init,d_costf_init,drx)
    
         do nR=1,n_r_max
            do lm=1,lm_max
               TstrLM(lm,nR)=TstrLM(lm,nR)+or1(nR)*workA(lm,nR)
            end do
         end do
         
         call get_drNS(TomeRLM,workA,lm_max_real,1,lm_max_real, &
                       n_r_max,n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         
         do nR=1,n_r_max
            do lm=1,lm_max
               TomeLM(lm,nR)=TomeLM(lm,nR)+or1(nR)*workA(lm,nR)
            end do
         end do
         
         call get_drNS(TadvRLM,workA,lm_max_real,1,lm_max_real, &
                       n_r_max,n_cheb_max,workB,i_costf_init,d_costf_init,drx)
         
         do nR=1,n_r_max
            do lm=1,lm_max
               TadvLM(lm,nR)=TadvLM(lm,nR)+or1(nR)*workA(lm,nR)
            end do
         end do
      end if
      !end do
    
      ! PdifLM and TdifLM need to be gathered over lm
      call gather_all_from_lo_to_rank0(gt_OC,PdifLM_LMloc,PdifLM)
      call gather_all_from_lo_to_rank0(gt_OC,TdifLM_LMloc,TdifLM)
         
      if ( l_DTrMagSpec .and. n_time_step > 1 ) then
    
         ! also gather PadvLMIC,TadvLMIC,PdifLMIC and TdifLMIC
         call gather_all_from_lo_to_rank0(gt_IC,PadvLMIC_LMloc,PadvLMIC)
         call gather_all_from_lo_to_rank0(gt_IC,TadvLMIC_LMloc,TadvLMIC)
         call gather_all_from_lo_to_rank0(gt_IC,PdifLMIC_LMloc,PdifLMIC)
         call gather_all_from_lo_to_rank0(gt_IC,TdifLMIC_LMloc,TdifLMIC)
    
         if ( rank == 0 ) then
            call rBrSpec(time,PstrLM,PadvLMIC,'rBrProSpec',.false.,st_map)
            call rBrSpec(time,PadvLM,PadvLMIC,'rBrAdvSpec',.true.,st_map)
            call rBrSpec(time,PdifLM,PdifLMIC,'rBrDifSpec',.true.,st_map)
            do nR=1,n_r_max
               do lm=1,lm_max
                  PstrLM(lm,nR)=PstrLM(lm,nR)-PadvLM(lm,nR)
               end do
            end do
            call rBrSpec(time,PstrLM,PadvLMIC,'rBrDynSpec',.false.,st_map)
    
            call rBpSpec(time,TstrLM,TadvLMIC,'rBpProSpec',.false.,st_map)
            call rBpSpec(time,TadvLM,TadvLMIC,'rBpAdvSpec',.true.,st_map)
            call rBpSpec(time,TdifLM,TdifLMIC,'rBpDifSpec',.true.,st_map)
            do nR=1,n_r_max
               do lm=1,lm_max
                  TstrLM(lm,nR)=TstrLM(lm,nR)-TadvLM(lm,nR)
               end do
            end do
            call rBpSpec(time,TstrLM,TadvLMIC,'rBpDynSpec',.false.,st_map)
         end if
      end if

   end subroutine get_dtBLMfinish
!-----------------------------------------------------------------------
   subroutine get_dH_dtBLM(nR,BtVrLM,BpVrLM,BrVtLM,BrVpLM, &
                           BtVpLM,BpVtLM,BrVZLM,BtVZLM,    &
                           BtVpCotLM,BpVtCotLM,BtVZcotLM,  &
                           BtVpSn2LM,BpVtSn2LM,BtVZsn2LM)
      !  +-------------------------------------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this routine is to calculate theta and phi            |
      !  |  derivative related terms of the magnetic production and          |
      !  |  advection terms and store them.                                  |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      integer,         intent(in) :: nR
      complex(kind=8), intent(in) :: BtVrLM(*),BpVrLM(*)
      complex(kind=8), intent(in) :: BrVtLM(*),BrVpLM(*)
      complex(kind=8), intent(in) :: BtVpLM(*),BpVtLM(*)
      complex(kind=8), intent(in) :: BtVpCotLM(*),BpVtCotLM(*)
      complex(kind=8), intent(in) :: BtVpSn2LM(*),BpVtSn2LM(*)
      complex(kind=8), intent(in) :: BtVZcotLM(*),BtVZsn2LM(*)
      complex(kind=8), intent(in) :: BrVZLM(*),BtVZLM(*)
    
      !-- Local variables:
      integer :: l,m,lm,lmP,lmPS,lmPA
      real(kind=8) :: fac
    
      PstrLM_Rloc(1,nR)=0.d0
      PadvLM_Rloc(1,nR)=0.D0
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         if ( l > m ) then
            PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   (                     &
                 dTheta1S(lm)*BtVrLM(lmPS) - dTheta1A(lm)*BtVrLM(lmPA) + &
                 dPhi(lm)*BpVrLM(lmP)  )
            PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   (                     &
                 dTheta1S(lm)*BrVtLM(lmPS) - dTheta1A(lm)*BrVtLM(lmPA) + &
                 dPhi(lm)*BrVpLM(lmP)  )
         else if ( l == m ) then
            PstrLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
                 - dTheta1A(lm)*BtVrLM(lmPA) + dPhi(lm)*BpVrLM(lmP)  )
            PadvLM_Rloc(lm,nR)=or2(nR)/dLh(lm) *   ( &
                 - dTheta1A(lm)*BrVtLM(lmPA) + dPhi(lm)*BrVpLM(lmP) )
         end if
      end do
    
      !--- Poloidal advection and stretching term finished for radial level nR !
    
      TstrLM_Rloc(1,nR) =0.D0
      TstrRLM_Rloc(1,nR)=0.D0
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
                         fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
                                             BpVtSn2LM(lmP) )   + &
                                                          fac * ( &
                          dTheta1S(lm) * ( or1(nR)*BpVrLM(lmPS) + &
                                                BpVtCotLM(lmPS) + &
                                              BtVpCotLM(lmPS) ) - &
                          dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
                                                BpVtCotLM(lmPA) + &
                                            BtVpCotLM(lmPA) ) ) - &
                               fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
            TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
                                     dTheta1S(lm)*BrVpLM(lmPS) - &
                                     dTheta1A(lm)*BrVpLM(lmPA) - &
                                         dPhi(lm)*BrVtLM(lmP)  )
         else if ( l == m ) then
            TstrLM_Rloc(lm,nR)=        -or2(nR)*BtVpLM(lmP)     - &
                         fac*dPhi(lm)*dPhi(lm)*( BtVpSn2LM(lmP) + &
                                             BpVtSn2LM(lmP) )   + &
                                                          fac * ( &
                        - dTheta1A(lm) * ( or1(nR)*BpVrLM(lmPA) + &
                                                BpVtCotLM(lmPA) + &
                                            BtVpCotLM(lmPA) ) ) - &
                               fac*or1(nR)*dPhi(lm)*BtVrLM(lmP)
            TstrRLM_Rloc(lm,nR)=             or1(nR)/dLh(lm) * ( &
                                   - dTheta1A(lm)*BrVpLM(lmPA) - &
                                         dPhi(lm)*BrVtLM(lmP)  )
         end if
      end do
    
      TadvLM_Rloc(1,nR)=0.D0
      TadvRLM_Rloc(1,nR)  =0.D0
      do lm=2,lm_max
         l   =lm2l(lm)
         m   =lm2m(lm)
         lmP =lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
                        fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
                                            BtVpSn2LM(lmP) )   + &
                                                         fac * ( &
                         dTheta1S(lm) * ( or1(nR)*BrVpLM(lmPS) + &
                                               BtVpCotLM(lmPS) + &
                                             BpVtCotLM(lmPS) ) - &
                         dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
                                               BtVpCotLM(lmPA) + &
                                           BpVtCotLM(lmPA) ) ) - &
                 fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
            TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
                        dTheta1S(lm)*BpVrLM(lmPS) - &
                        dTheta1A(lm)*BpVrLM(lmPA) - &
                            dPhi(lm)*BtVrLM(lmP)   )
         else if ( l == m ) then
            TadvLM_Rloc(lm,nR)=       -or2(nR)*BpVtLM(lmP)     - &
                        fac*dPhi(lm)*dPhi(lm)*( BpVtSn2LM(lmP) + &
                                            BtVpSn2LM(lmP) )   + &
                                                         fac * ( &
                       - dTheta1A(lm) * ( or1(nR)*BrVpLM(lmPA) + &
                                               BtVpCotLM(lmPA) + &
                                           BpVtCotLM(lmPA) ) ) - &
                             fac*or1(nR)*dPhi(lm)*BrVtLM(lmP)
            TadvRLM_Rloc(lm,nR)=or2(nR)/dLh(lm) * ( &
                      - dTheta1A(lm)*BpVrLM(lmPA) - &
                            dPhi(lm)*BtVrLM(lmP)   )
         end if
      end do
    
      !--- TomeLM same as TstrLM but where ever Vp appeared
      !    it is replaced by its axisymmetric contribution VZ:
      TomeLM_Rloc(1,nR) =0.D0
      TomeRLM_Rloc(1,nR)=0.D0
      do lm=2,lm_max
         l  =lm2l(lm)
         m  =lm2m(lm)
         lmP=lm2lmP(lm)
         lmPS=lmP2lmPS(lmP)
         lmPA=lmP2lmPA(lmP)
         fac=or2(nR)/dLh(lm)
         if ( l > m ) then
            TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
                      fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    + &
                           fac*( dTheta1S(lm)*BtVZCotLM(lmPS) - &
                                 dTheta1A(lm)*BtVZCotLM(lmPA) )
            TomeRLM_Rloc(lm,nR)=          or2(nR)/dLh(lm) * ( &
                                  dTheta1S(lm)*BrVZLM(lmPS) - &
                                  dTheta1A(lm)*BrVZLM(lmPA) )
         else if ( l == m ) then
            TomeLM_Rloc(lm,nR)=    -or2(nR)*BtVZLM(lmp)       - &
                      fac*dPhi(lm)*dPhi(lm)*BtVZsn2LM(lmP)    - &
                           fac*dTheta1A(lm)*BtVZCotLM(lmPA)
            TomeRLM_Rloc(lm,nR)=- or2(nR)/dLh(lm) * dTheta1A(lm)*BrVZLM(lmPA)
         end if
      end do
    
   end subroutine get_dH_dtBLM
!-----------------------------------------------------------------------
end module dtB_mod
