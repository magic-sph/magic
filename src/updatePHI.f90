module updatePhi_mod
   !
   ! This module handles the time advance of the phase field phi.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: or2, rscheme_oc, r, or1
   use num_param, only: dct_counter, solve_counter
   use physical_parameters, only: pr, phaseDiffFac, stef, ktopphi, kbotphi
   use init_fields, only: phi_top, phi_bot
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use logic, only: l_finite_diff, l_full_sphere, l_parallel_solve
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
   use radial_der, only: get_ddr, get_ddr_ghost, exch_ghosts, bulk_to_ghost
   use constants, only: zero, one, two
   use fields, only: work_LMloc
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use real_matrices
   use band_matrices
   use parallel_solvers, only: type_tri_par

   implicit none

   private

   !-- Local variables
   real(cp), allocatable :: rhs1(:,:,:)
   integer :: maxThreads
   class(type_realmat), pointer :: phiMat(:), phi0Mat
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: phiMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: phi0Mat_fac(:)
#endif
   logical, public, allocatable :: lPhimat(:)
   type(type_tri_par), public :: phiMat_FD
   complex(cp), public, allocatable :: phi_ghost(:,:)

   public :: initialize_updatePhi, finalize_updatePhi, updatePhi, assemble_phase,  &
   &         get_phase_rhs_imp, get_phase_rhs_imp_ghost, updatePhase_FD,           &
   &         preparePhase_FD, fill_ghosts_Phi, assemble_phase_Rloc

contains

   subroutine initialize_updatePhi

      integer :: ll, n_bands
      integer, pointer :: nLMBs2(:)
#ifdef WITH_OMP_GPU
      logical :: use_gpu, use_pivot
      use_gpu = .false.; use_pivot = .true.
#endif

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( l_finite_diff ) then
            allocate( type_bandmat :: phiMat(nLMBs2(1+rank)) )
            allocate( type_bandmat :: phi0Mat )

            if ( rscheme_oc%order == 2 .and. rscheme_oc%order_boundary <= 2 ) then ! Dirichelt BCs
               n_bands = rscheme_oc%order+1
            else
               n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
            end if

#ifdef WITH_OMP_GPU
            call phi0Mat%initialize(n_bands,n_r_max,use_pivot,use_gpu)
            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_bands,n_r_max,use_pivot,use_gpu)
            end do
#else
            call phi0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
#endif
         else
            allocate( type_densemat :: phiMat(nLMBs2(1+rank)) )
            allocate( type_densemat :: phi0Mat )

#ifdef WITH_OMP_GPU
            use_gpu = .true.
            call phi0Mat%initialize(n_r_max,n_r_max,use_pivot,use_gpu)
            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_r_max,n_r_max,use_pivot,use_gpu)
            end do
#else
            call phi0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
#endif
         end if

#ifdef WITH_PRECOND_S
         allocate(phiMat_fac(n_r_max,nLMBs2(1+rank)))
         phiMat_fac(:,:) = 0.0_cp
         bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: phiMat_fac)
         !$omp target update to(phiMat_fac)
         gpu_bytes_allocated = gpu_bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif
#endif
#ifdef WITH_PRECOND_S0
         allocate(phi0Mat_fac(n_r_max))
         phi0Mat_fac(:) = 0.0_cp
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: phi0Mat_fac)
         !$omp target update to(phi0Mat_fac)
         gpu_bytes_allocated = gpu_bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
#endif

#ifdef WITHOMP
         maxThreads=omp_get_max_threads()
#else
         maxThreads=1
#endif
         allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         rhs1 = 0.0_cp
         bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
         &                 maxThreads*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: rhs1)
         !$omp target update to(rhs1)
         gpu_bytes_allocated = gpu_bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
         &                 maxThreads*SIZEOF_DEF_COMPLEX
#endif
      else ! Parallel solvers are requested

         !-- Create matrix
         call phiMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Allocate an array with ghost zones
         allocate( phi_ghost(lm_max, nRstart-1:nRstop+1) )
         bytes_allocated=bytes_allocated + lm_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         phi_ghost(:,:)=zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: phi_ghost)
         !$omp target update to(phi_ghost)
         gpu_bytes_allocated=gpu_bytes_allocated + lm_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
#endif

      end if

      allocate( lPhimat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

   end subroutine initialize_updatePhi
!------------------------------------------------------------------------------
   subroutine finalize_updatePhi
      !
      ! This subroutine deallocates the matrices involved in the time-advance of
      ! phi.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lPhimat )
      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call phiMat(ll)%finalize()
         end do
         call phi0Mat%finalize()

#ifdef WITH_PRECOND_S
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: phiMat_fac)
#endif
         deallocate(phiMat_fac)
#endif
#ifdef WITH_PRECOND_S0
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: phi0Mat_fac)
#endif
         deallocate(phi0Mat_fac)
#endif
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: rhs1)
#endif
         deallocate( rhs1 )
      else
         call phiMat_FD%finalize()
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: phi_ghost)
#endif
         deallocate(phi_ghost)
      end if

   end subroutine finalize_updatePhi
!------------------------------------------------------------------------------
   subroutine updatePhi(phi, dphidt, tscheme)
      !
      !  Updates the phase field
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max) ! Chemical composition
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables:
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out            ! counts cheb modes
      real(cp), allocatable :: rhs(:) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid,iChunk,nChunks,size_of_last_chunk,lmB0

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank

      allocate(rhs(n_r_max))
      rhs = 0.0_cp

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dphidt)

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: rhs)
      !$omp target update to(rhs)
      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! one subblock is linked to one l value and needs therefore once the matrix
      do nLMB2=1,nLMBs2(nLMB)
         lmB=0

         !-- LU factorisation (big loop but hardly any work because of lPhimat)
         do lm=1,sizeLMB2(nLMB2,nLMB)
            l1=lm22l(lm,nLMB2,nLMB)
            if ( .not. lPhimat(l1) ) then

               if ( l1 == 0 ) then
#ifdef WITH_PRECOND_S0
                  call get_phi0Mat(tscheme,phi0Mat,phi0Mat_fac)
#else
                  call get_phi0Mat(tscheme,phi0Mat)
#endif
               else ! l /= 0
#ifdef WITH_PRECOND_S
                  call get_phiMat(tscheme,l1,phiMat(nLMB2),phiMat_fac(:,nLMB2))
#else
                  call get_phiMat(tscheme,l1,phiMat(nLMB2))
#endif
               end if
               lPhimat(l1)=.true.
            end if
         end do

         !-- Assemble RHS
         !$omp target map(tofrom: lmB) &
         !$omp& private(lm1, l1, nR)
         do lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1=lm22l(lm,nLMB2,nLMB)
            if ( l1 == 0 ) then
               rhs(1)      =phi_top
               rhs(n_r_max)=phi_bot
               do nR=2,n_r_max-1
                  rhs(nR)=real(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S0
               rhs(:) = phi0Mat_fac(:)*rhs(:)
#endif
            else ! l1  /=  0
               lmB=lmB+1

               rhs1(1,2*lmB-1,0)      =0.0_cp
               rhs1(1,2*lmB,0)        =0.0_cp
               rhs1(n_r_max,2*lmB-1,0)=0.0_cp
               rhs1(n_r_max,2*lmB,0)  =0.0_cp
               do nR=2,n_r_max-1
                  rhs1(nR,2*lmB-1,0)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lmB,0)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lmB-1,0)=phiMat_fac(:,nLMB2)*rhs1(:,2*lmB-1,0)
               rhs1(:,2*lmB,0)  =phiMat_fac(:,nLMB2)*rhs1(:,2*lmB,0)
#endif
            end if
         end do
         !$omp end target

         !-- Solve matrices with batched RHS (hipsolver)
         if ( lmB == 0 ) then
            if(.not. phi0Mat%gpu_is_used) then
               !$omp target update from(rhs)
            end if
            call phi0Mat%solve(rhs)
            if(.not. phi0Mat%gpu_is_used) then
               !$omp target update to(rhs)
            end if
         else
            if(.not. phiMat(nLMB2)%gpu_is_used) then
               !$omp target update from(rhs1)
            end if
            call phiMat(nLMB2)%solve(rhs1(:,:,0),2*lmB)
            if(.not. phiMat(nLMB2)%gpu_is_used) then
               !$omp target update to(rhs1)
            end if
         end if

         lmB=0
         !-- Loop to reassemble fields
         !$omp target map(tofrom: lmB) &
         !$omp& private(lm1, l1, m1, n_r_out)
         do lm=1,sizeLMB2(nLMB2,nLMB)
            lm1=lm22lm(lm,nLMB2,nLMB)
            l1=lm22l(lm,nLMB2,nLMB)
            m1=lm22m(lm,nLMB2,nLMB)

            if ( l1 == 0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  phi(lm1,n_r_out)=rhs(n_r_out)
               end do
            else
               lmB=lmB+1
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lmB-1,0), &
                     &                      rhs1(n_r_out,2*lmB,0),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lmB-1,0), &
                     &                       0.0_cp,kind=cp)
                  end do
               end if
            end if
         end do
         !$omp end target

      end do     ! loop over lm blocks
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single
      !$omp target exit data map(delete: rhs)
#else
      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! one subblock is linked to one l value and needs therefore once the matrix
      !$omp single
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lm,lm1,l1,m1,lmB,threadid) &
         !$omp private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 == 0 ) then
            if ( .not. lPhimat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_phi0Mat(tscheme,phi0Mat,phi0Mat_fac)
#else
               call get_phi0Mat(tscheme,phi0Mat)
#endif
               lPhimat(l1)=.true.
            end if
         else
            if ( .not. lPhimat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_phiMat(tscheme,l1,phiMat(nLMB2),phiMat_fac(:,nLMB2))
#else
               call get_phiMat(tscheme,l1,phiMat(nLMB2))
#endif
                lPhimat(l1)=.true.
            end if
         end if

         do iChunk=1,nChunks
            !$omp task default(shared) &
            !$omp firstprivate(iChunk) &
            !$omp private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$omp private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  rhs(1)      =phi_top
                  rhs(n_r_max)=phi_bot
                  do nR=2,n_r_max-1
                     rhs(nR)=real(work_LMloc(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs(:) = phi0Mat_fac(:)*rhs(:)
#endif

                  call phi0Mat%solve(rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,2*lmB-1,threadid)      =0.0_cp
                  rhs1(1,2*lmB,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lmB-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lmB,threadid)  =0.0_cp
                  do nR=2,n_r_max-1
                     rhs1(nR,2*lmB-1,threadid)= real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lmB,threadid)  =aimag(work_LMloc(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S
                  rhs1(:,2*lmB-1,threadid)=phiMat_fac(:,nLMB2)*rhs1(:,2*lmB-1,threadid)
                  rhs1(:,2*lmB,threadid)  =phiMat_fac(:,nLMB2)*rhs1(:,2*lmB,threadid)
#endif

               end if
            end do

            if ( lmB  >  lmB0 ) then
               call phiMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*lmB,threadid), &
                    &                   2*(lmB-lmB0))
            end if

            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        phi(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                      rhs1(n_r_out,2*lmB,threadid),kind=cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        phi(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                       0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !$omp end task
         end do
         !$omp taskwait
         !$omp end task
      end do     ! loop over lm blocks
      !$omp end single
      !$omp taskwait
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single
#endif

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
#ifdef WITH_OMP_GPU
      !$omp target
#else
      !$omp do private(n_r_out,lm1) collapse(2)
#endif
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            phi(lm1,n_r_out)=zero
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target
#else
      !$omp end do

      !$omp end parallel
#endif

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dphidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_phase_rhs_imp(phi, dphidt, 1, tscheme%l_imp_calc_rhs(1),  &
              &                l_in_cheb_space=.true.)
      else
         call get_phase_rhs_imp(phi, dphidt, tscheme%istage+1,            &
              &                tscheme%l_imp_calc_rhs(tscheme%istage+1),  &
              &                l_in_cheb_space=.true.)
      end if

      deallocate(rhs)

   end subroutine updatePhi
!------------------------------------------------------------------------------
   subroutine preparePhase_FD(tscheme, dphidt)
      !
      ! This subroutine is used to assemble the r.h.s. of the phase field equation
      ! when parallel F.D solvers are used. Boundary values are set here.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm, l

      !-- LU factorisation of the matrix if needed
      if ( .not. lPhimat(0) ) then
         call get_phiMat_Rdist(tscheme,phiMat_FD)
         lPhimat(:)=.true.
      end if

#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
#else
      !$omp parallel default(shared) private(lm_start,lm_stop, nR, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier
#endif

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(phi_ghost, dphidt, lm_start, lm_stop, 1)

      !-- Set boundary conditions
      if ( nRstart == n_r_cmb ) then
         nR=n_r_cmb
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( ktopphi == 1 ) then ! Dirichlet
               if ( l == 0 ) then
                  phi_ghost(lm,nR)=phi_top
               else
                  phi_ghost(lm,nR)=zero
               end if
            !else ! Neuman
            end if
            phi_ghost(lm,nR-1)=zero ! Set ghost zone to zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)

            if ( l_full_sphere ) then
               if ( l == 1 ) then
                  phi_ghost(lm,nR)=0.0_cp
               end if
            else
               if ( kbotphi == 1 ) then
                  if ( l == 0 ) then
                     phi_ghost(lm,nR)=phi_bot
                  else
                     phi_ghost(lm,nR)=zero
                  end if
               end if
            end if
            phi_ghost(lm,nR+1)=zero ! Set ghost zone to zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

   end subroutine preparePhase_FD
!------------------------------------------------------------------------------
   subroutine fill_ghosts_Phi(phig)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Sherwood
      ! numbers).
      !
      complex(cp), intent(inout) :: phig(lm_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: lm, l, lm_start, lm_stop
      real(cp) :: dr

#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
#else
      !$omp parallel default(shared) private(lm_start, lm_stop, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier
#endif
      !-- Handle upper boundary
      dr = r(2)-r(1)
      if ( nRstart == n_r_cmb ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=lm_start,lm_stop
            if ( ktopphi == 1 ) then
               phig(lm,nRstart-1)=two*phig(lm,nRstart)-phig(lm,nRstart+1)
            else
               phig(lm,nRstart-1)=phig(lm,nRstart+1)
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

      !-- Handle Lower boundary
      dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_icb ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( l_full_sphere ) then
               if ( l == 1 ) then
                  phig(lm,nRstop+1)=two*phig(lm,nRstop)-phig(lm,nRstop-1)
               else
                  if ( l == 0 ) then
                     phig(lm,nRstop+1)=phig(lm,nRstop-1)+two*dr*phi_bot
                  else
                     phig(lm,nRstop+1)=phig(lm,nRstop-1)
                  end if
               end if
            else ! Not a full sphere
               if ( kbotphi == 1 ) then
                  phig(lm,nRstop+1)=two*phig(lm,nRstop)-phig(lm,nRstop-1)
               else
                  phig(lm,nRstop+1)=phig(lm,nRstop-1)
               end if
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

   end subroutine fill_ghosts_Phi
!------------------------------------------------------------------------------
   subroutine updatePhase_FD(phi, dphidt, tscheme)
      !
      ! This subroutine is called after the linear solves have been completed.
      ! This is then assembling the linear terms that will be used in the r.h.s.
      ! for the next iteration.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dphidt
      complex(cp),       intent(inout) :: phi(lm_max,nRstart:nRstop) ! Phase field

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dphidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_phase_rhs_imp_ghost(phi_ghost, dphidt, 1, tscheme%l_imp_calc_rhs(1))
      else
         call get_phase_rhs_imp_ghost(phi_ghost, dphidt, tscheme%istage+1,       &
              &                       tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !-- Array copy from phi_ghost to phi
#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared) private(lm_start,lm_stop,nR,lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
#endif
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            phi(lm,nR)=phi_ghost(lm,nR)
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end parallel
#endif

   end subroutine updatePhase_FD
!------------------------------------------------------------------------------
   subroutine get_phase_rhs_imp(phi, dphidt, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for phase field. This is the LM-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      complex(cp), allocatable :: dphi(:,:)
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l
      real(cp) :: dL
      integer, pointer :: lm2l(:)

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      allocate(dphi(llm:ulm,n_r_max))
      dphi = zero

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: dphi)
      !$omp target update to(dphi)
      start_lm=llm; stop_lm=ulm
      call dct_counter%start_count()
#else
      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
#endif

#ifdef WITH_OMP_GPU
      call get_ddr(phi, dphi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) then
         call rscheme_oc%costf1(phi,ulm-llm+1,start_lm-llm+1, &
                               &                 stop_lm-llm+1,.true.)
      end if
#else
      call get_ddr(phi, dphi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(phi,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
#endif

#ifdef WITH_OMP_GPU
      call dct_counter%stop_count(l_increment=.false.)
#else
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
#endif

      if ( istage == 1 ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
         do n_r=1,n_r_max
            do lm=llm,ulm
               dphidt%old(lm,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phi(lm,n_r)
            end do
         end do
         !$omp end target teams distribute parallel do
#else
         !$omp do
         do n_r=1,n_r_max
            dphidt%old(:,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phi(:,n_r)
         end do
         !$omp end do
#endif
      end if

      if ( l_calc_lin ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#else
         !$omp do private(n_r,lm,l,dL)
#endif
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dphidt%impl(lm,n_r,istage)=phaseDiffFac*( work_LMloc(lm,n_r) + &
               &                                two*or1(n_r) * dphi(lm,n_r) - &
               &                                 dL*or2(n_r) *  phi(lm,n_r) )
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end do
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dphi)
#endif
      deallocate(dphi)

   end subroutine get_phase_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_phase_rhs_imp_ghost(phig, dphidt, istage, l_calc_lin)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for phase field. This is the R-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin

      !-- Output variable
      complex(cp),       intent(inout) :: phig(lm_max,nRstart-1:nRstop+1)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      complex(cp), allocatable :: dphi(:,:) ! Radial derivative of phase field
      complex(cp), allocatable :: work_Rloc(:,:)
      integer :: n_r, lm, start_lm, stop_lm, l
      real(cp) :: dL
      integer, pointer :: lm2l(:)

      allocate(dphi(lm_max,nRstart:nRstop), work_Rloc(lm_max,nRstart:nRstop))
      dphi = zero; work_Rloc = zero

      lm2l(1:lm_max) => st_map%lm2l
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: dphi, work_Rloc)
      !$omp target update to(dphi, work_Rloc)
#endif

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
      call dct_counter%start_count()
      call get_ddr_ghost(phig, dphi, work_Rloc, lm_max, start_lm, stop_lm,  nRstart, &
           &             nRstop, rscheme_oc)
      call dct_counter%stop_count(l_increment=.false.)
#else
      !$omp parallel default(shared)  private(start_lm, stop_lm, n_r, lm, l, dL)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr_ghost(phig, dphi, work_Rloc, lm_max, start_lm, stop_lm,  nRstart, &
           &             nRstop, rscheme_oc)
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
      !$omp barrier
#endif

      if ( istage == 1 ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               dphidt%old(lm,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phig(lm,n_r)
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

      if ( l_calc_lin ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dphidt%impl(lm,n_r,istage)= phaseDiffFac * ( work_Rloc(lm,n_r) + &
               &                                 two*or1(n_r) *  dphi(lm,n_r) - &
               &                                  dL*or2(n_r) *  phig(lm,n_r) )
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dphi, work_Rloc)
#endif
      deallocate(dphi, work_Rloc)

   end subroutine get_phase_rhs_imp_ghost
!------------------------------------------------------------------------------
   subroutine assemble_phase(phi, dphidt, tscheme)
      !
      ! This subroutine is used to assemble the phase field when an
      ! IMEX-RK with an assembly stage is employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: lm, l, m, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dphidt)

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m)
#endif
      do n_r=2,n_r_max
         do lm=llm,ulm
            m = lm2m(lm)
            if ( m == 0 ) then
               phi(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp) * &
               &           6.0_cp/5.0_cp/stef/pr
            else
               phi(lm,n_r)=work_LMloc(lm,n_r)*6.0_cp/5.0_cp/stef/pr
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
      !$omp target update from(phi) !-- TODO: Mandatory as robin_bc is on CPU currently
#else
      !$omp end do
#endif

      !-- Boundary conditions
      if ( l_full_sphere) then
         if ( ktopphi == 1 ) then ! Dirichlet
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm,l)
#endif
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 1 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               else
                  if ( l == 0 ) then
                     call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                          &                   one, 0.0_cp, cmplx(phi_bot,0.0_cp,cp), &
                          &                   phi(lm,:))
                  else
                     call rscheme_oc%robin_bc(0.0_cp, one, zero, one, 0.0_cp, &
                          &                   zero, phi(lm,:))
                  end if
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         else ! Neummann
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm,l)
#endif
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 1 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, &
                       &                   zero, phi(lm,:))
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         end if

      else ! Spherical shell

         if ( ktopphi==1 .and. kbotphi==1 ) then
            !-- Boundary conditions: Dirichlet on both sides
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm,l)
#endif
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   0.0_cp, one, cmplx(phi_bot,0.0_cp,cp), &
                       &                   phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         else if ( ktopphi==1 .and. kbotphi /= 1 ) then
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm,l)
#endif
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   one, 0.0_cp, zero, phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, one, 0.0_cp, &
                       &                   zero, phi(lm,:))
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         else if ( ktopphi/=1 .and. kbotphi == 1 ) then
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm,l)
#endif
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   cmplx(phi_bot,0.0_cp,cp), phi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         else if ( ktopphi/=1 .and. kbotphi /= 1 ) then
            !-- Boundary conditions: Neuman on both sides
#ifdef WITH_OMP_GPU
            !$omp parallel do default(shared) private(lm,l)
#else
            !$omp do private(lm)
#endif
            do lm=llm,ulm
               call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, zero, phi(lm,:))
            end do
#ifdef WITH_OMP_GPU
            !$omp end parallel do
#else
            !$omp end do
#endif
         end if

      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

#ifdef WITH_OMP_GPU
      !$omp target update to(phi)  !-- TODO: Update on GPU after robin_bc
#endif

      call get_phase_rhs_imp(phi, dphidt, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_phase
!------------------------------------------------------------------------------
   subroutine assemble_phase_Rloc(phi, dphidt, tscheme)
      !
      ! This subroutine is used when an IMEX Runge-Kutta time scheme with an assembly
      ! stage is used. This is used when R is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: phi(lm_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: lm, l, m, n_r, start_lm, stop_lm
      complex(cp), allocatable :: work_Rloc(:,:)

      allocate(work_Rloc(lm_max,nRstart:nRstop))
      work_Rloc = zero
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: work_Rloc)
      !$omp target update to(work_Rloc)
#endif

      call tscheme%assemble_imex(work_Rloc, dphidt)

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
#else
      !$omp parallel default(shared) private(start_lm, stop_lm, l, m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier
#endif

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            m = st_map%lm2m(lm)
            if ( m == 0 ) then
               phi(lm,n_r)=cmplx(real(work_Rloc(lm,n_r)),0.0_cp,cp)* &
               &           6.0_cp/5.0_cp/stef/pr
            else
               phi(lm,n_r)=work_Rloc(lm,n_r)*6.0_cp/5.0_cp/stef/pr
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      if ( ktopphi==1 .and. nRstart==n_r_cmb ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l == 0 ) then
               phi(lm,nRstart)=phi_top
            else
               phi(lm,nRstart)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      end if

      if ( kbotphi==1 .and. nRstop==n_r_icb ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do
#endif
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l == 0 ) then
               phi(lm,nRstop)=phi_bot
            else
               phi(lm,nRstop)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

#ifdef WITH_OMP_GPU
      call bulk_to_ghost(phi, phi_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm, .true.)
      !$omp target update from(phi_ghost)
#else
      call bulk_to_ghost(phi, phi_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm)
#endif

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

      call exch_ghosts(phi_ghost, lm_max, nRstart, nRstop, 1) !-- Run on CPU (MPI comm)

#ifdef WITH_OMP_GPU
      !$omp target update to(phi_ghost) !-- Update on GPU after MPI exchange
#endif
      call fill_ghosts_Phi(phi_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_phase_rhs_imp_ghost(phi_ghost, dphidt, 1, tscheme%l_imp_calc_rhs(1))

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: work_Rloc)
#endif
      deallocate(work_Rloc)

   end subroutine assemble_phase_Rloc
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_phi0Mat(tscheme,phiMat,phiMat_fac)
#else
   subroutine get_phi0Mat(tscheme,phiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  phiMat0
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      class(type_realmat), intent(inout) :: phiMat
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: phiMat_fac(n_r_max)
#endif

      !-- Local variables:
      real(cp), allocatable :: dat(:,:)
      integer :: info, nR_out, nR

      allocate(dat(n_r_max,n_r_max))
      dat(:,:) = 0.0_cp

      !----- Boundary condition:
      if ( ktopphi == 1 ) then
         !--------- Constant phase field at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         !--------- dphi/dr=0 at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      else
         if ( kbotphi == 1 ) then
            !--------- Constant phase field at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            !--------- dphi/dr=0 at ICB
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: dat)
      !$omp target update to(dat)
#endif

      !-- Fill bulk points
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                                   &
            &                 5.0_cp/6.0_cp*stef*pr*  rscheme_oc%rMat(nR,nR_out) - &
            &  tscheme%wimp_lin(1)*phaseDiffFac*(   rscheme_oc%d2rMat(nR,nR_out) + &
            &                         two*or1(nR)*   rscheme_oc%drMat(nR,nR_out) ) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !----- Factors for highest and lowest cheb mode:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

#ifdef WITH_PRECOND_S0
      ! compute the linesum of each line
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         phiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each line by the linesum to regularize the matrix
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*phiMat_fac(nR)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
#endif

#ifdef WITH_OMP_GPU
      if(.not. phiMat%gpu_is_used) then
         !$omp target update from(dat)
      end if
#endif

      !-- Array copy
      call phiMat%set_data(dat)

      !---- LU decomposition:
      call phiMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix phiMat0!')

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dat)
#endif
      deallocate(dat)

   end subroutine get_phi0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_phiMat(tscheme,l,phiMat,phiMat_fac)
#else
   subroutine get_phiMat(tscheme,l,phiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  phiMat(i,j) for the equation for phase field.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: phiMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: phiMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: dLh
      real(cp), allocatable :: dat(:,:)

      allocate(dat(n_r_max,n_r_max))
      dat(:,:) = 0.0_cp

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktopphi == 1 ) then ! Dirichlet
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else ! Neumann
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         if ( l == 1 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      else
         if ( kbotphi == 1 ) then ! Dirichlet
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: dat)
      !$omp target update to(dat)
#endif

      !----- Bulk points
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                                &
            &               5.0_cp/6.0_cp*stef*pr* rscheme_oc%rMat(nR,nR_out) - &
            &  tscheme%wimp_lin(1)*phaseDiffFac*(rscheme_oc%d2rMat(nR,nR_out) + &
            &                     two*or1(nR)*    rscheme_oc%drMat(nR,nR_out) - &
            &                     dLh*or2(nR)*     rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !----- Factor for highest and lowest cheb:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         phiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each line by the linesum to regularize the matrix
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*phiMat_fac(nR)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
#endif

#ifdef WITH_OMP_GPU
      if(.not. phiMat%gpu_is_used) then
         !$omp target update from(dat)
      end if
#endif

      !-- Array copy
      call phiMat%set_data(dat)

      !----- LU decomposition:
      call phiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix phiMat!')

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dat)
#endif
      deallocate(dat)

   end subroutine get_phiMat
!-----------------------------------------------------------------------------
   subroutine get_phiMat_Rdist(tscheme,phiMat)
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  phiMat(i,j) for the equation for the phase field. This is
      !  used when parallel F.D. solvers are employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      type(type_tri_par), intent(inout) :: phiMat

      !-- Local variables:
      integer :: nR, l
      real(cp) :: dLh
      real(cp) :: wimp_lin

      !-- Copie into local variable
      wimp_lin = tscheme%wimp_lin(1)

      !----- Bulk points
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared) private(nR,l,dLh)
      !$omp do
#endif
      do nR=1,n_r_max
         do l=0,l_max
            dLh=real(l*(l+1),kind=cp)
            phiMat%diag(l,nR)=    5.0_cp/6.0_cp*stef*pr-            &
            &                  wimp_lin*phaseDiffFac*(   &
            &                                rscheme_oc%ddr(nR,1) + &
            &                     two*or1(nR)*rscheme_oc%dr(nR,1) - &
            &                                        dLh*or2(nR) )
            phiMat%low(l,nR)=-wimp_lin*phaseDiffFac*(    &
            &                                rscheme_oc%ddr(nR,0) + &
            &                 two*or1(nR)*    rscheme_oc%dr(nR,0) )
            phiMat%up(l,nR) =-wimp_lin*phaseDiffFac*(    &
            &                                rscheme_oc%ddr(nR,2) + &
            &                 two*or1(nR)*    rscheme_oc%dr(nR,2) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end do
#endif

      !----- Boundary conditions:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp do
#endif
      do l=0,l_max
         if ( ktopphi == 1 ) then ! Dirichlet
            phiMat%diag(l,1)=one
            phiMat%up(l,1)  =0.0_cp
            phiMat%low(l,1) =0.0_cp
         else
            phiMat%up(l,1)=phiMat%up(l,1)+phiMat%low(l,1)
         end if

         if ( l_full_sphere ) then
            !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
            if ( l == 1 ) then
               phiMat%diag(l,n_r_max)=one
               phiMat%up(l,n_r_max)  =0.0_cp
               phiMat%low(l,n_r_max) =0.0_cp
            else
               phiMat%low(l,n_r_max)=phiMat%up(l,n_r_max)+phiMat%low(l,n_r_max)
               !fd_fac_bot(l)=two*(r(n_r_max-1)-r(n_r_max))*phiMat%up(l,n_r_max)
            end if
         else
            if ( kbotphi == 1 ) then ! Dirichlet
               phiMat%diag(l,n_r_max)=one
               phiMat%up(l,n_r_max)  =0.0_cp
               phiMat%low(l,n_r_max) =0.0_cp
            else
               phiMat%low(l,n_r_max)=phiMat%up(l,n_r_max)+phiMat%low(l,n_r_max)
            end if
         end if
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end do
      !$omp end parallel
#endif

      !----- LU decomposition:
      call phiMat%prepare_mat()

   end subroutine get_phiMat_Rdist
!-----------------------------------------------------------------------------
end module updatePhi_mod
