module updateWP_mod
   !
   ! This module handles the time advance of the poloidal potential w and the pressure p.
   ! It contains the computation of the implicit terms and the linear solves.
   !

   use omp_lib
   use precision_mod
#ifdef WITH_OMP_GPU
   use mem_alloc, only: bytes_allocated, gpu_bytes_allocated
#else
   use mem_alloc, only: bytes_allocated
#endif
   use truncation, only: lm_max, n_r_max, l_max, m_min
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: or1, or2, rho0, rgrav, visc, dLvisc, r, &
       &                       alpha0, temp0, beta, dbeta, ogrun,      &
       &                       rscheme_oc, ddLvisc, ddbeta, orho1
   use physical_parameters, only: kbotv, ktopv, ra, BuoFac, ChemFac,   &
       &                          ViscHeatFac, ThExpNb, ktopp
   use num_param, only: dct_counter, solve_counter
   use blocking, only: lo_sub_map, lo_map, st_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_V
   use logic, only: l_update_v, l_chemical_conv, l_RMS, l_double_curl, &
       &            l_fluxProfs, l_finite_diff, l_full_sphere, l_heat, &
       &            l_parallel_solve
   use RMS, only: DifPol2hInt, DifPolLMr
   use communications, only: get_global_sum
   use parallel_mod
   use radial_der, only: get_dddr, get_ddr, get_dr, get_dr_Rloc, get_ddddr_ghost, &
       &                 bulk_to_ghost, exch_ghosts
   use integration, only: rInt_R
   use fields, only: work_LMloc, s_Rloc, xi_Rloc !TODO> pass directly
   use constants, only: zero, one, two, three, four, third, half
   use useful, only: abortRun, cc2real
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use parallel_solvers
   use dense_matrices
   use real_matrices
   use band_matrices
#ifdef WITH_OMP_GPU
   use iso_c_binding
   use hipfort_check
   use hipfort_hipblas
#endif

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: ddddw(:,:)
   complex(cp), allocatable :: dwold(:,:)
   real(cp), allocatable :: work(:)
#ifdef WITH_OMP_GPU
   real(cp), allocatable, target :: rhs1(:,:,:), rhs(:)
#else
   real(cp), allocatable :: rhs1(:,:,:), rhs(:)
#endif
   real(cp), allocatable :: rhs0(:,:,:)
   real(cp), allocatable :: wpMat_fac(:,:,:)
   class(type_realmat), allocatable :: wpMat(:), p0Mat, ellMat(:)
   logical, public, allocatable :: lWPmat(:)
   logical, allocatable :: l_ellMat(:)
   type(type_penta_par), public :: wMat_FD
   type(type_tri_par), public :: ellMat_FD, p0Mat_FD
   complex(cp), public, allocatable :: w_ghost(:,:), p0_ghost(:)
   integer :: maxThreads, size_rhs1

#ifdef WITH_OMP_GPU
   real(cp), allocatable :: dat(:,:)
   type(c_ptr) :: handle = c_null_ptr
   integer, allocatable, target :: devInfo(:)
#endif

   public :: initialize_updateWP, finalize_updateWP, updateWP, assemble_pol,       &
   &         finish_exp_pol, get_pol_rhs_imp, finish_exp_pol_Rdist, fill_ghosts_W, &
   &         prepareW_FD, updateW_FD, get_pol_rhs_imp_ghost, assemble_pol_Rloc

contains

   subroutine initialize_updateWP(tscheme)
      !
      ! Purpose of this subroutine is to allocate the matrices needed
      ! to time advance the poloidal/pressure equations. Depending on the
      ! radial scheme, it can be either full or band matrices.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer, pointer :: nLMBs2(:)
      integer :: ll, n_bands
#ifdef WITH_OMP_GPU
      logical :: use_gpu, use_pivot
      use_gpu = .false.; use_pivot = .true.
#endif

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

#ifdef WITHOMP
         maxThreads=omp_get_max_threads()
#else
         maxThreads=1
#endif

         if ( l_finite_diff ) then
            allocate( type_bandmat :: wpMat(nLMBs2(1+rank)) )

            if ( rscheme_oc%order <= 2 .and. rscheme_oc%order_boundary <= 2 ) then
               n_bands =rscheme_oc%order+3
            else
               n_bands = max(rscheme_oc%order+3,2*rscheme_oc%order_boundary+3)
            end if
            do ll=1,nLMBs2(1+rank)
               call wpMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
            allocate( wpMat_fac(n_r_max,2,nLMBs2(1+rank)) )
            wpMat_fac(:,:,:) = 0.0_cp
            bytes_allocated=bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
            &               SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: wpMat_fac)
            !$omp target update to(wpMat_fac)
            gpu_bytes_allocated=gpu_bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
            &                   SIZEOF_DEF_REAL
#endif
            allocate( type_bandmat :: p0Mat )
            n_bands = rscheme_oc%order+1
            call p0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         else
            allocate( type_densemat :: wpMat(nLMBs2(1+rank)) )
            if ( l_double_curl ) then
#ifdef WITH_OMP_GPU
               use_gpu = .true.
               do ll=1,nLMBs2(1+rank)
                  call wpMat(ll)%initialize(n_r_max,n_r_max,use_pivot,use_gpu)
               end do
#else
               do ll=1,nLMBs2(1+rank)
                  call wpMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
               end do
#endif
               allocate( wpMat_fac(n_r_max,2,nLMBs2(1+rank)) )
               bytes_allocated=bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
               &               SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: wpMat_fac)
            !$omp target update to(wpMat_fac)
            gpu_bytes_allocated=gpu_bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
            &                   SIZEOF_DEF_REAL
#endif
            else
#ifdef WITH_OMP_GPU
               use_gpu = .true.
               do ll=1,nLMBs2(1+rank)
                  call wpMat(ll)%initialize(2*n_r_max,2*n_r_max,use_pivot,use_gpu)
               end do
#else
               do ll=1,nLMBs2(1+rank)
                  call wpMat(ll)%initialize(2*n_r_max,2*n_r_max,l_pivot=.true.)
               end do
#endif
               allocate( wpMat_fac(2*n_r_max,2,nLMBs2(1+rank)) )
               bytes_allocated=bytes_allocated+4*n_r_max*nLMBs2(1+rank)*    &
               &               SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
               !$omp target enter data map(alloc: wpMat_fac)
               !$omp target update to(wpMat_fac)
               gpu_bytes_allocated=bytes_allocated+4*n_r_max*nLMBs2(1+rank)*    &
               &                   SIZEOF_DEF_REAL
#endif
            end if

            allocate( type_densemat :: p0Mat )
#ifdef WITH_OMP_GPU
            use_gpu = .true.
            call p0Mat%initialize(n_r_max,n_r_max,use_pivot,use_gpu)
#else
            call p0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
#endif
         end if

         if ( l_double_curl ) then
            allocate( ddddw(llm:ulm,n_r_max) )
            ddddw = zero
            bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: ddddw)
            !$omp target update to(ddddw)
            gpu_bytes_allocated = gpu_bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
#endif
            if ( l_RMS .or. l_FluxProfs ) then
               allocate( dwold(llm:ulm,n_r_max) )
               bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
               dwold(:,:)=zero
            end if
         end if

         allocate( work(n_r_max) )
         bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL

         allocate( rhs(n_r_max) )
         bytes_allocated=bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX
         rhs(:)=0.0_cp
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: rhs)
         !$omp target update to(rhs)
         gpu_bytes_allocated=gpu_bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
         if ( l_double_curl ) then
            size_rhs1 = n_r_max
            allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
            rhs1 = zero
            bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
            &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: rhs1)
            !$omp target update to(rhs1)
            gpu_bytes_allocated=gpu_bytes_allocated+n_r_max*maxThreads* &
            &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
#endif
         else
            size_rhs1 = 2*n_r_max
            allocate( rhs1(2*n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
            rhs1 = zero
            bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
            &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: rhs1)
            !$omp target update to(rhs1)
            gpu_bytes_allocated=gpu_bytes_allocated+2*n_r_max*maxThreads* &
            &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
#endif
         end if

         if ( tscheme%l_assembly .and. l_double_curl ) then
            allocate( type_bandmat :: ellMat(nLMBs2(1+rank)) )
            if ( rscheme_oc%order <= 2 .and. rscheme_oc%order_boundary <= 2 .and. &
            &    ktopv /=1 .and. kbotv /=1 ) then
               !n_bands =rscheme_oc%order+1 # should be that but yield matrix singularity?
               n_bands = max(rscheme_oc%order+1,2*rscheme_oc%order_boundary+1)
            else
               n_bands = max(rscheme_oc%order+1,2*rscheme_oc%order_boundary+1)
            end if
            do ll=1,nLMBs2(1+rank)
               call ellMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
            allocate( rhs0(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
            rhs0(:,:,:)=zero
            bytes_allocated = bytes_allocated+n_r_max*maxThreads*2* &
            &                 lo_sub_map%sizeLMB2max*SIZEOF_DEF_REAL
#ifdef WITH_OMP_GPU
            !$omp target enter data map(alloc: rhs0)
            !$omp target update to(rhs0)
            gpu_bytes_allocated = gpu_bytes_allocated+n_r_max*maxThreads*2* &
            &                 lo_sub_map%sizeLMB2max*SIZEOF_DEF_REAL
#endif
         end if

      else ! Parallel solver

         call p0Mat_FD%initialize(1,n_r_max,1,1)
         call wMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Allocate an array with ghost zones
         allocate( w_ghost(lm_max,nRstart-2:nRstop+2), p0_ghost(nRstart-1:nRstop+1) )
         bytes_allocated=bytes_allocated+lm_max*(nRstop-nRstart+5)*SIZEOF_DEF_COMPLEX &
         &               +(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         w_ghost(:,:)=zero
         p0_ghost(:) =zero
#ifdef WITH_OMP_GPU
         !$omp target enter data map(alloc: w_ghost, p0_ghost)
         !$omp target update to(w_ghost, p0_ghost)
         gpu_bytes_allocated=gpu_bytes_allocated+lm_max*(nRstop-nRstart+5)*SIZEOF_DEF_COMPLEX &
         &               +(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
#endif

         if ( l_RMS .or. l_FluxProfs ) then
            allocate( dwold(lm_max,nRstart:nRstop) )
            bytes_allocated = bytes_allocated+lm_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
            dwold(:,:)=zero
         end if
         if ( tscheme%l_assembly .and. l_double_curl ) then
            call ellMat_FD%initialize(1,n_r_max,0,l_max)
         end if
      end if

      if ( tscheme%l_assembly .and. l_double_curl ) then
         allocate( l_ellMat(0:l_max) )
         l_ellMat(:) = .false.
         bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_LOGICAL
      end if

      allocate( lWPmat(0:l_max) )
      bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITH_OMP_GPU
      allocate(dat(n_r_max,n_r_max))
      dat(:,:) = 0.0_cp
      !$omp target enter data map(alloc: dat)
      !$omp target update to(dat)
      if ( ( .not. l_parallel_solve ) .and. ( .not. l_finite_diff) ) then
         call hipblasCheck(hipblasCreate(handle))
         allocate(devInfo(1))
         devInfo(1) = 0
         !$omp target enter data map(alloc: devInfo)
         !$omp target update to(devInfo)
      end if
#endif

   end subroutine initialize_updateWP
!-----------------------------------------------------------------------------
   subroutine finalize_updateWP(tscheme)
      !
      ! Deallocation of the matrices used to time-advance the poloidal/pressure
      ! equations.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer, pointer :: nLMBs2(:)
      integer :: ll

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( tscheme%l_assembly .and. l_double_curl ) then
            do ll=1,nLMBs2(1+rank)
               call ellMat(ll)%finalize()
            end do
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: rhs0)
#endif
            deallocate( rhs0 )
         end if

         do ll=1,nLMBs2(1+rank)
            call wpMat(ll)%finalize()
         end do
         call p0Mat%finalize()

#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: rhs1, rhs, wpMat_fac)
#endif
         deallocate( wpMat_fac, rhs1, work, rhs )
         if ( l_double_curl ) then
#ifdef WITH_OMP_GPU
            !$omp target exit data map(delete: ddddw)
#endif
            deallocate( ddddw )
            if ( l_RMS .or. l_FluxProfs ) deallocate( dwold )
         end if
      else ! Parallel solver
         call p0Mat_FD%finalize()
         call wMat_FD%finalize()
#ifdef WITH_OMP_GPU
         !$omp target exit data map(delete: w_ghost, p0_ghost)
#endif
         deallocate( w_ghost, p0_ghost )
         if ( l_RMS .or. l_FluxProfs ) deallocate( dwold )
         if ( tscheme%l_assembly .and. l_double_curl ) call ellMat_FD%finalize()
      end if
      if ( tscheme%l_assembly .and. l_double_curl ) deallocate(l_ellMat)
      deallocate( lWPmat )

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: dat)
      deallocate(dat)
      if ( ( .not. l_parallel_solve ) .and. ( .not. l_finite_diff) ) then
         call hipblasCheck(hipblasDestroy(handle))
         !$omp target exit data map(delete: devInfo)
         deallocate(devInfo)
      end if
#endif

   end subroutine finalize_updateWP
!-----------------------------------------------------------------------------
   subroutine updateWP(s, xi, w, dw, ddw, dwdt, p, dp, dpdt, tscheme, &
              &        lRmsNext, lPressNext)
      !
      !  updates the poloidal velocity potential w, the pressure p, and
      !  their radial derivatives.
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext

      type(type_tarray), intent(inout) :: dpdt
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: ddw(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)

      complex(cp),       intent(inout) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1      ! degree and order
      integer :: lm1,lm     ! position of (l,m) in array
      integer :: nLMB2
      integer :: nR         ! counts radial grid points
      integer :: n_r_out    ! counts cheb modes
      integer :: nLMB

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid
#ifdef WITH_OMP_GPU
      complex(cp), pointer :: ptr_expl(:,:)
      real(cp) :: wimp_lin
      integer :: nRows
      integer :: n_max_rSchemeOc
      ptr_expl => dwdt%expl(:,:,tscheme%istage)
      wimp_lin = tscheme%wimp_lin(1)
      if( l_double_curl ) then
         nRows = n_r_max
      else
         nRows = 2*n_r_max
      end if
      n_max_rSchemeOc = rscheme_oc%n_max
#endif

      if ( .not. l_update_v ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dwdt)
      if ( .not. l_double_curl ) then
         call tscheme%set_imex_rhs(ddw, dpdt)
      end if

#ifndef WITH_OMP_GPU
      !$omp parallel default(shared)
#endif

      if ( l_double_curl .and. lPressNext .and. tscheme%istage == 1) then
         ! Store old dw
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#else
         !$omp do private(lm)
#endif
         do nR=1,n_r_max
            do lm=llm,ulm
               dwold(lm,nR)=dw(lm,nR)
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#else
         !$omp end do
#endif
      end if

#ifdef WITH_OMP_GPU
      !$omp single
      call solve_counter%start_count()
      !$omp end single

      ! each of the nLMBs2(nLMB) subblocks have one l value
      !-- MPI Level
      do nLMB2=1,nLMBs2(nLMB)

         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 == 0 ) then ! Spherically-symmetric modes
            lm1=lm2(0,0)
            if ( .not. lWPmat(l1) ) then
               call get_p0Mat(p0Mat)
               lWPmat(l1)=.true.
            end if

            !-- The integral of rho' r^2 dr vanishes
            if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then
               if ( rscheme_oc%version == 'cheb' ) then
                  !$omp target teams distribute parallel do
                  do nR=1,n_r_max
                     work(nR)=ThExpNb*alpha0(nR)*temp0(nR)*rho0(nR)*r(nR)*&
                     &        r(nR)*real(s(lm1,nR))
                  end do
                  !$omp end target teams distribute parallel do
                  rhs(1)=rInt_R(work,r,rscheme_oc)
               else
                  rhs(1)=0.0_cp
               end if
               !$omp target update to(rhs(1))
            else
               rhs(1)=0.0_cp
               !$omp target update to(rhs(1))
            end if

            if ( l_chemical_conv ) then
               !$omp target teams distribute parallel do
               do nR=2,n_r_max
                  rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))+   &
                  &       rho0(nR)*ChemFac*rgrav(nR)*real(xi(lm1,nR))+ &
                  &       real(ptr_expl(lm1,nR))
               end do
               !$omp end target teams distribute parallel do
            else
               !$omp target teams distribute parallel do
               do nR=2,n_r_max
                  rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))+  &
                  &       real(ptr_expl(lm1,nR))
               end do
               !$omp end target teams distribute parallel do
            end if

            !-- Matrix solve for l=m=0
            if (.not. p0Mat%gpu_is_used) then
               !$omp target update from(rhs)
               call p0Mat%solve(rhs)
               !$omp target update to(rhs)
            else
               call p0Mat%solve(rhs, handle, devInfo)
            end if

            !$omp target teams distribute parallel do
            do n_r_out=1,n_max_rSchemeOc
               p(lm1,n_r_out)=rhs(n_r_out)
               w(lm1,n_r_out)=zero
            end do
            !$omp end target teams distribute parallel do

         else ! l /= 0 

            if ( .not. lWPmat(l1) ) then
               if ( l_double_curl ) then
                  call get_wMat(tscheme,l1,hdif_V(l1),wpMat(nLMB2), &
                       &        wpMat_fac(:,:,nLMB2))
               else
                  call get_wpMat(tscheme,l1,hdif_V(l1),wpMat(nLMB2), &
                       &         wpMat_fac(:,:,nLMB2))
               end if
               lWPmat(l1)=.true.
            end if

            !-- Assemble RHS
            !$omp target teams distribute parallel do private(lm1,m1,nR)
            do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1=lm22m(lm,nLMB2,nLMB)

               rhs1(1,2*lm-1,0)      =0.0_cp
               rhs1(1,2*lm,0)        =0.0_cp
               rhs1(n_r_max,2*lm-1,0)=0.0_cp
               rhs1(n_r_max,2*lm,0)  =0.0_cp
               if ( l_double_curl ) then
                  rhs1(2,2*lm-1,0)        =0.0_cp
                  rhs1(2,2*lm,0)          =0.0_cp
                  rhs1(n_r_max-1,2*lm-1,0)=0.0_cp
                  rhs1(n_r_max-1,2*lm,0)  =0.0_cp
                  do nR=3,n_r_max-2
                     rhs1(nR,2*lm-1,0)= real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lm,0)  =aimag(work_LMloc(lm1,nR))
                  end do

                  if ( l_heat .and. (.not. l_parallel_solve) ) then
                     do nR=3,n_r_max-2
                        rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)+             &
                        &      wimp_lin*real(l1*(l1+1),cp) *             &
                        &      or2(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))
                        rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)+               &
                        &      wimp_lin*real(l1*(l1+1),cp) *             &
                        &      or2(nR)*BuoFac*rgrav(nR)*aimag(s(lm1,nR))
                     end do
                  end if

                  if ( l_chemical_conv .and. ( .not. l_parallel_solve ) ) then
                     do nR=3,n_r_max-2
                        rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)+            &
                        &      wimp_lin*real(l1*(l1+1),cp) *            &
                        &      or2(nR)*ChemFac*rgrav(nR)*real(xi(lm1,nR))
                        rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)+            &
                        &      wimp_lin*real(l1*(l1+1),cp) *          &
                        &      or2(nR)*ChemFac*rgrav(nR)*aimag(xi(lm1,nR))
                     end do
                  end if
               else
                  rhs1(n_r_max+1,2*lm-1,0)=0.0_cp
                  rhs1(n_r_max+1,2*lm,0)  =0.0_cp
                  rhs1(2*n_r_max,2*lm-1,0)=0.0_cp
                  rhs1(2*n_r_max,2*lm,0)  =0.0_cp
                  do nR=2,n_r_max-1
                     rhs1(nR,2*lm-1,0)        = real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lm,0)          =aimag(work_LMloc(lm1,nR))
                     rhs1(nR+n_r_max,2*lm-1,0)= real(ddw(lm1,nR)) ! ddw is a work array
                     rhs1(nR+n_r_max,2*lm,0)  =aimag(ddw(lm1,nR))
                  end do

                  if ( l_heat ) then
                     do nR=2,n_r_max-1
                        rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)+              &
                        &         wimp_lin*rho0(nR)*BuoFac*               &
                        &                      rgrav(nR)*real(s(lm1,nR))
                        rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)+                &
                        &         wimp_lin*rho0(nR)*BuoFac*               &
                        &                      rgrav(nR)*aimag(s(lm1,nR))
                     end do
                  end if

                  if ( l_chemical_conv ) then
                     do nR=2,n_r_max-1
                        rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)+               &
                        &         wimp_lin*rho0(nR)*ChemFac*               &
                        &                      rgrav(nR)*real(xi(lm1,nR))
                        rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)+                 &
                        &         wimp_lin*rho0(nR)*ChemFac*               &
                        &                      rgrav(nR)*aimag(xi(lm1,nR))
                     end do
                  end if
               end if

               ! use the mat_fac(:,1) to scale the rhs
               do nR=1,size_rhs1
                  rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)*wpMat_fac(nR,1,nLMB2)
                  rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)*wpMat_fac(nR,1,nLMB2)
               end do
            end do
            !$omp end target teams distribute parallel do

            !-- Solve matrices with batched RHS (hipsolver)
            lm=sizeLMB2(nLMB2,nLMB)
            if (.not. wpMat(nLMB2)%gpu_is_used) then
               !$omp target update from(rhs1)
               call wpMat(nLMB2)%solve(rhs1(:,:,0),2*lm)
               !$omp target update to(rhs1)
            else
               call wpMat(nLMB2)%solve(rhs1(:,:,0),2*lm,handle,devInfo)
            end if

            !-- Loop to reassemble fields
            !$omp target teams distribute parallel do private(lm1,m1,nR,n_r_out)
            do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1=lm22m(lm,nLMB2,nLMB)

               do nR=1,size_rhs1 ! rescale solution
                  rhs1(nR,2*lm-1,0)=rhs1(nR,2*lm-1,0)*wpMat_fac(nR,2,nLMB2)
                  rhs1(nR,2*lm,0)  =rhs1(nR,2*lm,0)*wpMat_fac(nR,2,nLMB2)
               end do

               if ( l_double_curl ) then
                  if ( m1 > 0 ) then
                     do n_r_out=1,n_max_rSchemeOc
                        w(lm1,n_r_out)  =cmplx(rhs1(n_r_out,2*lm-1,0), &
                        &                      rhs1(n_r_out,2*lm,0),cp)
                     end do
                  else
                     do n_r_out=1,n_max_rSchemeOc
                        w(lm1,n_r_out)  = cmplx(rhs1(n_r_out,2*lm-1,0),&
                        &                       0.0_cp,kind=cp)
                     end do
                  end if
               else
                  if ( m1 > 0 ) then
                     do n_r_out=1,n_max_rSchemeOc
                        w(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,0), &
                        &                    rhs1(n_r_out,2*lm,0),cp)
                        p(lm1,n_r_out)=cmplx(rhs1(n_r_max+n_r_out,2*lm-1,0), &
                        &                    rhs1(n_r_max+n_r_out,2*lm,0),cp)
                     end do
                  else
                     do n_r_out=1,n_max_rSchemeOc
                        w(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,0), &
                        &                     0.0_cp,kind=cp)
                        p(lm1,n_r_out)= cmplx(rhs1(n_r_max+n_r_out,2*lm-1,0), &
                        &                     0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !$omp end target teams distribute parallel do

         end if ! Test of l /= 0
      end do

      !$omp single
      call solve_counter%stop_count()
      !$omp end single
#else
      !$omp single
      call solve_counter%start_count()
      !$omp end single

      !$omp single
      ! each of the nLMBs2(nLMB) subblocks have one l value
      !-- MPI Level
      do nLMB2=1,nLMBs2(nLMB)
         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lm,lm1,l1,m1,iChunk,nChunks,size_of_last_chunk,threadid) &
         !$omp shared(dwold,nLMB,nLMBs2,rhs1)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 == 0 ) then
            if ( .not. lWPmat(l1) ) then
               call get_p0Mat(p0Mat)
               lWPmat(l1)=.true.
            end if
         else
            if ( .not. lWPmat(l1) ) then
               if ( l_double_curl ) then
                  call get_wMat(tscheme,l1,hdif_V(l1),wpMat(nLMB2),wpMat_fac(:,:,nLMB2))
               else
                  call get_wpMat(tscheme,l1,hdif_V(l1),wpMat(nLMB2),wpMat_fac(:,:,nLMB2))
               end if
               lWPmat(l1)=.true.
            end if
         end if

         if ( l1 == 0 ) then
            lm1 = lm2(0,0)

            !-- The integral of rho' r^2 dr vanishes
            if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then
               if ( rscheme_oc%version == 'cheb' ) then
                  do nR=1,n_r_max
                     work(nR)=ThExpNb*alpha0(nR)*temp0(nR)*rho0(nR)*r(nR)*&
                     &        r(nR)*real(s(lm1,nR))
                  end do
                  rhs(1)=rInt_R(work,r,rscheme_oc)
               else
                  rhs(1)=0.0_cp
               end if
            else
               rhs(1)=0.0_cp
            end if

            if ( l_chemical_conv ) then
               do nR=2,n_r_max
                  rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))+   &
                  &       rho0(nR)*ChemFac*rgrav(nR)*real(xi(lm1,nR))+ &
                  &       real(dwdt%expl(lm1,nR,tscheme%istage))
               end do
            else
               do nR=2,n_r_max
                  rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))+  &
                  &       real(dwdt%expl(lm1,nR,tscheme%istage))
               end do
            end if

            call p0Mat%solve(rhs)

            do n_r_out=1,rscheme_oc%n_max
               p(lm1,n_r_out)=rhs(n_r_out)
               w(lm1,n_r_out)=zero
            end do

         else ! l1 > 0

            do iChunk=1,nChunks
               !$omp task if (nChunks>1) default(shared) &
               !$omp firstprivate(iChunk) &
               !$omp private(lmB0,lm,lm1,m1,nR,n_r_out) &
               !$omp private(threadid)

#ifdef WITHOMP
               threadid = omp_get_thread_num()
#else
               threadid = 0
#endif
               lmB0=(iChunk-1)*chunksize
               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)

                  rhs1(1,2*lm-1,threadid)      =0.0_cp
                  rhs1(1,2*lm,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lm-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lm,threadid)  =0.0_cp
                  if ( l_double_curl ) then
                     rhs1(2,2*lm-1,threadid)        =0.0_cp
                     rhs1(2,2*lm,threadid)          =0.0_cp
                     rhs1(n_r_max-1,2*lm-1,threadid)=0.0_cp
                     rhs1(n_r_max-1,2*lm,threadid)  =0.0_cp
                     do nR=3,n_r_max-2
                        rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                        rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
                     end do

                     if ( l_heat .and. (.not. l_parallel_solve) ) then
                        do nR=3,n_r_max-2
                           rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)+ &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *    &
                           &      or2(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))
                           rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)+   &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *    &
                           &      or2(nR)*BuoFac*rgrav(nR)*aimag(s(lm1,nR))
                        end do
                     end if

                     if ( l_chemical_conv .and. ( .not. l_parallel_solve ) ) then
                        do nR=3,n_r_max-2
                           rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)+ &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *    &
                           &      or2(nR)*ChemFac*rgrav(nR)*real(xi(lm1,nR))
                           rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)+   &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *    &
                           &      or2(nR)*ChemFac*rgrav(nR)*aimag(xi(lm1,nR))
                        end do
                     end if
                  else
                     rhs1(n_r_max+1,2*lm-1,threadid)=0.0_cp
                     rhs1(n_r_max+1,2*lm,threadid)  =0.0_cp
                     rhs1(2*n_r_max,2*lm-1,threadid)=0.0_cp
                     rhs1(2*n_r_max,2*lm,threadid)  =0.0_cp
                     do nR=2,n_r_max-1
                        rhs1(nR,2*lm-1,threadid)        = real(work_LMloc(lm1,nR))
                        rhs1(nR,2*lm,threadid)          =aimag(work_LMloc(lm1,nR))
                        rhs1(nR+n_r_max,2*lm-1,threadid)= real(ddw(lm1,nR)) ! ddw is a work array
                        rhs1(nR+n_r_max,2*lm,threadid)  =aimag(ddw(lm1,nR))
                     end do

                     if ( l_heat ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)+ &
                           &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*     &
                           &                      rgrav(nR)*real(s(lm1,nR))
                           rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)+   &
                           &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*     &
                           &                      rgrav(nR)*aimag(s(lm1,nR))
                        end do
                     end if

                     if ( l_chemical_conv ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)+ &
                           &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*    &
                           &                      rgrav(nR)*real(xi(lm1,nR))
                           rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)+   &
                           &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*    &
                           &                      rgrav(nR)*aimag(xi(lm1,nR))
                        end do
                     end if

                  end if

                  do nR=1,size_rhs1
                     rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)* &
                     &                        wpMat_fac(nR,1,nLMB2)
                     rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)* &
                     &                        wpMat_fac(nR,1,nLMB2)
                  end do
               end do

               call wpMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                    &                  2*(lm-1-lmB0))

               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)

                  do nR=1,size_rhs1
                     rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)* &
                     &                        wpMat_fac(nR,2,nLMB2)
                     rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)* &
                     &                        wpMat_fac(nR,2,nLMB2)
                  end do

                  if ( l_double_curl ) then
                     if ( m1 > 0 ) then
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  =cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                           &                      rhs1(n_r_out,2*lm,threadid),cp)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  = cmplx(rhs1(n_r_out,2*lm-1,threadid),&
                           &                       0.0_cp,kind=cp)
                        end do
                     end if
                  else
                     if ( m1 > 0 ) then
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                           &                    rhs1(n_r_out,2*lm,threadid),cp)
                           p(lm1,n_r_out)=cmplx(rhs1(n_r_max+n_r_out,2*lm-1,    &
                           &                    threadid),rhs1(n_r_max+n_r_out, &
                           &                    2*lm,threadid),cp)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                           &                     0.0_cp,kind=cp)
                           p(lm1,n_r_out)= cmplx(rhs1(n_r_max+n_r_out,2*lm-1, &
                           &                     threadid),0.0_cp,kind=cp)
                        end do
                     end if
                  end if
               end do
               !$omp end task
            end do
         end if
         !$omp taskwait
         !$omp end task
      end do   ! end of loop over l1 subblocks
      !$omp end single
      !$omp taskwait

      !$omp single
      call solve_counter%stop_count()
      !$omp end single
#endif

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
#ifdef WITH_OMP_GPU
      !$omp target teams
      do n_r_out=n_max_rSchemeOc+1,n_r_max
         !$omp distribute parallel do
#else
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
#endif
         do lm1=llm,ulm
            w(lm1,n_r_out)=zero
            p(lm1,n_r_out)=zero
         end do
#ifdef WITH_OMP_GPU
         !$omp end distribute parallel do
#endif
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams
#else
      !$omp end do
      !$omp end parallel
#endif

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dwdt)
      if ( .not. l_double_curl ) call tscheme%rotate_imex(dpdt)

      if ( tscheme%istage == tscheme%nstages ) then
         call get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt,       &
              &               tscheme, 1, tscheme%l_imp_calc_rhs(1),      &
              &               lPressNext, lRmsNext,                       &
              &               dpdt%expl(:,:,1), l_in_cheb_space=.true.)
         ! dpdt%expl(:,:,1) needed for RMS calc: first stage explicit term
      else
         call get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt,       &
              &               tscheme, tscheme%istage+1,                  &
              &               tscheme%l_imp_calc_rhs(tscheme%istage+1),   &
              &               lPressNext, lRmsNext,                       &
              &               dpdt%expl(:,:,1), l_in_cheb_space=.true.)
      end if

   end subroutine updateWP
!------------------------------------------------------------------------------
   subroutine prepareW_FD(tscheme, dwdt, lPressNext)

      !-- Input of variable
      logical,             intent(in) :: lPressNext
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dwdt

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm, l, lm00

      if ( .not. l_update_v ) return

      !-- LU factorisation of the matrix if needed
      if ( .not. lWPmat(1) ) then
         call get_wMat_Rdist(tscheme, hdif_V, wMat_FD)
         call get_p0Mat_Rdist(p0Mat_FD)
         lWPmat(:)=.true.
      end if

      if ( lPressNext .and. (m_min==0) ) then
         lm00=st_map%lm2(0,0)
         do nR=nRstart,nRstop
            p0_ghost(nR)=dwdt%expl(lm00,nR,tscheme%istage)
            if ( l_heat ) then
               p0_ghost(nR)=p0_ghost(nR)+rho0(nR)*BuoFac*rgrav(nR)*s_Rloc(lm00,nR)
            end if
            if ( l_chemical_conv ) then
               p0_ghost(nR)=p0_ghost(nR)+rho0(nR)*ChemFac*rgrav(nR)*xi_Rloc(lm00,nR)
            end if
         end do
         if ( nRstart == n_r_cmb ) p0_ghost(nRstart)=zero
      end if

#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
#else
      !$omp parallel default(shared) private(lm_start,lm_stop, nR, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier
#endif

      !-- Now assemble the right hand side and store it in work_LMloc
#ifdef WITH_OMP_GPU
      !$omp target update to(dwdt)
#endif
      call tscheme%set_imex_rhs_ghost(w_ghost, dwdt, lm_start, lm_stop, 2)
#ifdef WITH_OMP_GPU
      !$omp target update from(w_ghost)
#endif

      !-- Ensure that l=m=0 is zero
#ifdef WITH_OMP_GPU
#else
#endif
      if ( m_min == 0 ) then
         do nR=nRstart,nRstop
            w_ghost(1,nR)=zero
         end do
      end if
#ifdef WITH_OMP_GPU
#else
#endif

      !-- Set boundary conditions
      if ( nRstart == n_r_cmb ) then
         nR=n_r_cmb
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=lm_start,lm_stop
            l=st_map%lm2l(lm)
            if ( l /= 0 ) then
               w_ghost(lm,nR)  =zero ! Non-penetration condition
               w_ghost(lm,nR-1)=zero ! Ghost zones set to zero
               w_ghost(lm,nR-2)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=lm_start,lm_stop
            l=st_map%lm2l(lm)
            if ( l /= 0 ) then
               w_ghost(lm,nR)=zero ! Non-penetration condition
               w_ghost(lm,nR+1)=zero ! Ghost zones set to zero
               w_ghost(lm,nR+2)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if
#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

   end subroutine prepareW_FD
!------------------------------------------------------------------------------
   subroutine fill_ghosts_W(wg,p0g,lPressNext)
      !
      ! This subroutine is used to fill the ghost zones.
      !

      logical,     intent(in)    :: lPressNext
      complex(cp), intent(inout) :: p0g(nRstart-1:nRstop+1)
      complex(cp), intent(inout) :: wg(lm_max, nRstart-2:nRstop+2)

      !-- Local variables
      integer :: lm, l, lm_start, lm_stop
      real(cp) :: dr

      if ( .not. l_update_v ) return

      if ( lPressNext ) then
         if ( nRstart == n_r_cmb ) then
            p0g(nRstart-1)=two*p0g(nRstart)-p0g(nRstart+1)
         end if
         if ( nRstop == n_r_icb ) then
            p0g(nRstop+1)=two*p0g(nRstop)-p0g(nRstop-1)
         end if
      end if

#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
#else
      !$omp parallel default(shared) private(lm_start, lm_stop, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier
#endif

      !-- Upper boundary
      dr = r(2)-r(1)
      if ( nRstart == n_r_cmb ) then ! Rank with n_r_mcb
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=lm_start,lm_stop
            l=st_map%lm2l(lm)
            if ( l /= 0 ) then
               if ( ktopv == 1 ) then  ! Stress-free
                  wg(lm,nRstart-1)=-(one-half*(two*or1(1)+beta(1))*dr)/ &
                  &                 (one+half*(two*or1(1)+beta(1))*dr) * wg(lm,nRstart+1)
               else ! Rigid boundary condition
                  wg(lm,nRstart-1)=wg(lm,nRstart+1) ! dw=0
               end if
               wg(lm,nRstart-2)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if

      !-- Lower boundary
      dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_icb ) then
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=lm_start,lm_stop
            l=st_map%lm2l(lm)
            if ( l /= 0 ) then
               if ( l_full_sphere ) then
                  if ( l == 1 ) then
                     wg(lm,nRstop+1)=wg(lm,nRstop-1) ! dw=0
                  else
                     wg(lm,nRstop+1)=-wg(lm,nRstop-1) ! ddw=0
                  end if
               else
                  if ( kbotv == 1 ) then ! Stress-free
                     wg(lm,nRstop+1)=-(one+half*(two*or1(n_r_max)+beta(n_r_max))*dr)/ &
                     &                (one-half*(two*or1(n_r_max)+beta(n_r_max))*dr) *&
                     &                wg(lm,nRstop-1)
                  else
                     wg(lm,nRstop+1)=wg(lm,nRstop-1) ! dw=0
                  end if
               end if
               wg(lm,nRstop+2)=zero
            end if
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

   end subroutine fill_ghosts_W
!------------------------------------------------------------------------------
   subroutine updateW_FD(w, dw, ddw, dwdt, p, dp, dpdt, tscheme, lRmsNext, &
              &          lPressNext, lP00Next)

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lP00Next
      type(type_tarray),   intent(in) :: dpdt

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(inout) :: w(lm_max,nRstart:nRstop) ! Poloidal potential
      !-- Output: ds
      complex(cp),       intent(inout) :: dw(lm_max,nRstart:nRstop) ! Radial derivative of w
      complex(cp),       intent(out) :: ddw(lm_max,nRstart:nRstop) ! Radial derivative of dw
      complex(cp),       intent(inout) :: p(lm_max,nRstart:nRstop) ! Pressure
      complex(cp),       intent(out) :: dp(lm_max,nRstart:nRstop) ! Radial derivative of p

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm, l

      if ( .not. l_update_v ) return

      if ( lPressNext .and. tscheme%istage == 1) then
         ! Store old dw
#ifdef WITH_OMP_GPU
#else
         !$omp parallel do collapse(2)
#endif
         do nR=nRstart,nRstop
            do lm=1,lm_max
               dwold(lm,nR)=dw(lm,nR)
            end do
         end do
#ifdef WITH_OMP_GPU
#else
         !$omp end parallel do
#endif
      end if

      !-- Roll the arrays before filling again the first block
#ifdef WITH_OMP_GPU
      !$omp target update to(dwdt)
#endif
      call tscheme%rotate_imex(dwdt)
#ifdef WITH_OMP_GPU
      !$omp target update from(dwdt)
#endif

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_pol_rhs_imp_ghost(w_ghost, dw, ddw, p, dp, dwdt, tscheme, 1, &
              &                     tscheme%l_imp_calc_rhs(1), lPressNext,     &
              &                     lRmsNext, dpdt%expl(:,:,1))
      else
         call get_pol_rhs_imp_ghost(w_ghost, dw, ddw, p, dp, dwdt, tscheme,   &
              &                     tscheme%istage+1,                         &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                     lPressNext, lRmsNext, dpdt%expl(:,:,1))
      end if

#ifdef WITH_OMP_GPU
      lm_start=1; lm_stop=lm_max
#else
      !$omp parallel default(shared) private(lm_start,lm_stop,nR,lm,l)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier
#endif
      !-- Array copy from w_ghost to w
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( l /= 0 ) then
               w(lm,nR)=w_ghost(lm,nR)
            else
               if ( lPressNext .or. lP00Next ) p(lm,nR)=p0_ghost(nR)
               cycle
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
#else
      !$omp end parallel
#endif

   end subroutine updateW_FD
!------------------------------------------------------------------------------
   subroutine get_pol(w, work)
      !
      !  Get the poloidal potential from the solve of an elliptic equation.
      !  Careful: the output is in Chebyshev space!
      !

      !-- Input field
      complex(cp), intent(in) :: work(llm:ulm,n_r_max)

      !-- Output field
      complex(cp), intent(out) :: w(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1      ! degree and order
      integer :: lm1,lm     ! position of (l,m) in array
      integer :: nLMB2
      integer :: nR         ! counts radial grid points
      integer :: n_r_out    ! counts cheb modes
      integer :: nLMB

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      if ( .not. l_update_v ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank

#ifndef WITH_OMP_GPU
      !$omp parallel default(shared)
#endif

      !-- Compute the right hand side
#ifdef WITH_OMP_GPU
      !$omp single
      call solve_counter%start_count()
      !$omp end single

      ! each of the nLMBs2(nLMB) subblocks have one l value
      !-- MPI Level
      do nLMB2=1,nLMBs2(nLMB)

         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 > 0 ) then
            !-- LU factorisation (big loop but hardly any work because of l_ellMat)
            if ( .not. l_ellMat(l1) ) then
               call get_elliptic_mat(l1, ellMat(nLMB2))
               l_ellMat(l1) = .true.
            end if

            !-- Assemble RHS
            do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)

               rhs0(1,2*lm-1,0)        =0.0_cp
               rhs0(1,2*lm,0)          =0.0_cp
               rhs0(2,2*lm-1,0)        =0.0_cp
               rhs0(2,2*lm,0)          =0.0_cp
               rhs0(n_r_max-1,2*lm-1,0)=0.0_cp
               rhs0(n_r_max-1,2*lm,0)  =0.0_cp
               rhs0(n_r_max,2*lm-1,0)  =0.0_cp
               rhs0(n_r_max,2*lm,0)    =0.0_cp
               do nR=3,n_r_max-2
                  rhs0(nR,2*lm-1,0)= real(work(lm1,nR))
                  rhs0(nR,2*lm,0)  =aimag(work(lm1,nR))
               end do
            end do

            !-- Solve matrices with batched RHS (hipsolver)
            lm=sizeLMB2(nLMB2,nLMB)
            call ellMat(nLMB2)%solve(rhs0(:,:,0),2*lm)

            !-- Loop to reassemble fields
            do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     w(lm1,n_r_out)  =cmplx(rhs0(n_r_out,2*lm-1,0), &
                     &                      rhs0(n_r_out,2*lm,0),cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     w(lm1,n_r_out)  = cmplx(rhs0(n_r_out,2*lm-1,0),&
                     &                       0.0_cp,kind=cp)
                  end do
               end if
            end do

         else ! l=0: ensure spherically-symmetric part is zero

            lm1=lm2(0,0)
            do n_r_out=1,rscheme_oc%n_max
               w(lm1,n_r_out) =zero
            end do

         end if
      end do   ! end of loop over l1 subblocks

      !$omp single
      call solve_counter%stop_count()
      !$omp end single
#else
      !$omp single
      call solve_counter%start_count()
      !$omp end single

      !$omp single
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)
         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lm,lm1,l1,m1,iChunk,nChunks,size_of_last_chunk,threadid) &
         !$omp shared(dwold,nLMB,nLMBs2,rhs0)

         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 > 0 ) then

            if ( .not. l_ellMat(l1) ) then
               call get_elliptic_mat(l1, ellMat(nLMB2))
               l_ellMat(l1) = .true.
            end if

            do iChunk=1,nChunks
               !$omp task if (nChunks>1) default(shared) &
               !$omp firstprivate(iChunk) &
               !$omp private(lmB0,lm,lm1,m1,nR,n_r_out) &
               !$omp private(threadid)

#ifdef WITHOMP
               threadid = omp_get_thread_num()
#else
               threadid = 0
#endif
               lmB0=(iChunk-1)*chunksize
               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  rhs0(1,2*lm-1,threadid)        =0.0_cp
                  rhs0(1,2*lm,threadid)          =0.0_cp
                  rhs0(2,2*lm-1,threadid)        =0.0_cp
                  rhs0(2,2*lm,threadid)          =0.0_cp
                  rhs0(n_r_max-1,2*lm-1,threadid)=0.0_cp
                  rhs0(n_r_max-1,2*lm,threadid)  =0.0_cp
                  rhs0(n_r_max,2*lm-1,threadid)  =0.0_cp
                  rhs0(n_r_max,2*lm,threadid)    =0.0_cp
                  do nR=3,n_r_max-2
                     rhs0(nR,2*lm-1,threadid)= real(work(lm1,nR))
                     rhs0(nR,2*lm,threadid)  =aimag(work(lm1,nR))
                  end do
               end do

               call ellMat(nLMB2)%solve(rhs0(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                    &                   2*(lm-1-lmB0))

               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        w(lm1,n_r_out)  =cmplx(rhs0(n_r_out,2*lm-1,threadid), &
                        &                      rhs0(n_r_out,2*lm,threadid),cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        w(lm1,n_r_out)  = cmplx(rhs0(n_r_out,2*lm-1,threadid),&
                        &                       0.0_cp,kind=cp)
                     end do
                  end if
               end do
               !$omp end task
            end do

         else ! l1 == 0

            lm1=lm2(0,0)
            do n_r_out=1,rscheme_oc%n_max
               w(lm1,n_r_out) =zero
            end do

         end if
         !$omp taskwait
         !$omp end task

      end do   ! end of loop over l1 subblocks
      !$omp end single
      !$omp taskwait
      !$omp single
      call solve_counter%stop_count()
      !$omp end single
#endif

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
#ifdef WITH_OMP_GPU
      !$omp parallel do private(n_r_out,lm1) collapse(2)
#else
      !$omp do private(n_r_out,lm1) collapse(2)
#endif
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            w(lm1,n_r_out)=zero
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end parallel do
#else
      !$omp end do
      !$omp end parallel
#endif

   end subroutine get_pol
!------------------------------------------------------------------------------
   subroutine finish_exp_pol(dVxVhLM, dw_exp_last)

#ifdef WITH_OMP_GPU
      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(:,:)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(:,:)
#else
      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(llm:ulm,n_r_max)
#endif
      !-- Local variables
      integer :: n_r, start_lm, stop_lm, lm

#ifdef WITH_OMP_GPU
      start_lm=llm; stop_lm=ulm
      call get_dr( dVxVhLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
           &       stop_lm-llm+1,n_r_max,rscheme_oc, .true., .true., .true. )
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVxVhLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
           &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      !$omp barrier
      !$omp do
#endif
      do n_r=1,n_r_max
         do lm=llm,ulm
            dw_exp_last(lm,n_r)= dw_exp_last(lm,n_r)+or2(n_r)*work_LMloc(lm,n_r)
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end do
      !$omp end parallel
#endif

   end subroutine finish_exp_pol
!------------------------------------------------------------------------------
   subroutine finish_exp_pol_Rdist(dVxVhLM, dw_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(lm_max,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(lm_max,nRstart:nRstop)

      !-- Local variables
      complex(cp), allocatable :: work_Rloc(:,:)
      integer :: n_r, start_lm, stop_lm, l, lm
      real(cp) :: dLh

      allocate(work_Rloc(lm_max,nRstart:nRstop))
      work_Rloc = zero
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: work_Rloc)
      !$omp target update to(work_Rloc)
#endif

      call get_dr_Rloc(dVxVhLM, work_Rloc, lm_max, nRstart, nRstop, n_r_max, &
           &           rscheme_oc)

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared) private(n_r, lm, l, dLh, start_lm, stop_lm)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm, stop_lm)
      !$omp barrier
#endif
      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l /= 0 ) then
               dw_exp_last(lm,n_r)=dw_exp_last(lm,n_r)+or2(n_r)*work_Rloc(lm,n_r)
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      if ( l_heat .and. l_parallel_solve ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l /= 0 ) then
                  dLh = real(l*(l+1),cp)
                  dw_exp_last(lm,n_r)=dw_exp_last(lm,n_r)+dLh*or2(n_r)*BuoFac* &
                  &                   rgrav(n_r)*s_Rloc(lm,n_r)
                end if
            end do
         end do
#ifdef WITH_OMP_GPU
         !$omp end target teams distribute parallel do
#endif
      end if

      if ( l_chemical_conv .and. l_parallel_solve ) then
#ifdef WITH_OMP_GPU
         !$omp target teams distribute parallel do collapse(2)
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l /= 0 ) then
                  dLh = real(l*(l+1),cp)
                  dw_exp_last(lm,n_r)=dw_exp_last(lm,n_r)+dLh*or2(n_r)*ChemFac* &
                  &                   rgrav(n_r)*xi_Rloc(lm,n_r)
               end if
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
      !$omp target exit data map(delete: work_Rloc)
#endif
      deallocate(work_Rloc)

   end subroutine finish_exp_pol_Rdist
!------------------------------------------------------------------------------
   subroutine get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, tscheme,     &
              &               istage, l_calc_lin, lPressNext, lRmsNext, dp_expl, &
              &               l_in_cheb_space)
      !
      ! This subroutine computes the derivatives of w and p and assemble the
      ! implicit stage if needed.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: xi(llm:ulm,n_r_max)
      logical,             intent(in) :: l_calc_lin
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      logical, optional,   intent(in) :: l_in_cheb_space
      complex(cp),         intent(in) :: dp_expl(llm:ulm,n_r_max)

      !-- Output variables
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dp(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp) :: Dif, Buo, Pre
      logical :: l_in_cheb
      integer :: n_r_top, n_r_bot, l1, m1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:),lm2m(:)
      real(cp) :: dL
      integer :: istage_loc, nstages
      real(cp) :: dt_loc

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      if ( m_min == 0 ) then
         lmStart_00=max(2,llm)
      else
         lmStart_00=llm
      end if

      istage_loc  = tscheme%istage
      nstages     = tscheme%nstages
      dt_loc      = tscheme%dt(1)

#ifdef WITH_OMP_GPU
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
      if ( l_double_curl ) then
         call get_ddr( w, dw, ddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc,      &
              &       l_dct_in=.not. l_in_cheb )
         call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc )
      else
         call get_dddr( w, dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &         stop_lm-llm+1, n_r_max, rscheme_oc,                &
              &         l_dct_in=.not. l_in_cheb)
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb, use_gpu=.true.)
         if ( l_in_cheb ) then
            call rscheme_oc%costf1(p,ulm-llm+1,start_lm-llm+1, &
                                  &                 stop_lm-llm+1,.true.)
         end if
      end if
      if ( l_in_cheb ) then
         call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1, &
                               &                 stop_lm-llm+1,.true.)
      end if
#else
      if ( l_double_curl ) then
         call get_ddr( w, dw, ddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc,      &
              &       l_dct_in=.not. l_in_cheb )
         call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc )
      else
         call get_dddr( w, dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &         stop_lm-llm+1, n_r_max, rscheme_oc,                &
              &         l_dct_in=.not. l_in_cheb)
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
         if ( l_in_cheb ) call rscheme_oc%costf1(p,ulm-llm+1,start_lm-llm+1, &
                               &                 stop_lm-llm+1)
      end if
      if ( l_in_cheb ) call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
#endif

#ifdef WITH_OMP_GPU
      call dct_counter%stop_count()
#else
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single
#endif

      if ( istage == 1 ) then
         if ( l_double_curl ) then
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2)
#else
            !$omp do private(n_r,lm,l1,dL)
#endif
            do n_r=2,n_r_max-1
               do lm=lmStart_00,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dwdt%old(lm,n_r,istage)=dL*or2(n_r)* ( -orho1(n_r)*(  &
                  &                   ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
                  &                             dL*or2(n_r)* w(lm,n_r) ) )
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end do
#endif
         else
#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2)
#else
            !$omp do private(n_r,lm,l1,dL)
#endif
            do n_r=2,n_r_max-1
               do lm=lmStart_00,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dwdt%old(lm,n_r,istage)= dL*or2(n_r)*w(lm,n_r)
                  dpdt%old(lm,n_r,istage)=-dL*or2(n_r)*dw(lm,n_r)
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end do
#endif
         end if
      end if

      if ( l_calc_lin .or. (tscheme%istage==tscheme%nstages .and. lRmsNext)) then

         if ( lRmsNext .and. tscheme%istage == tscheme%nstages ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !-- Calculate explicit time step part:
         if ( l_double_curl ) then

            if ( lPressNext ) then
               n_r_top=n_r_cmb
               n_r_bot=n_r_icb
            end if

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2) private(l1,m1,Dif,Buo,dL)
#else
            !$omp do private(n_r,lm,l1,m1,dL,Dif,Buo)
#endif
            do n_r=n_r_top,n_r_bot
               do lm=lmStart_00,ulm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif=-hdif_V(l1)*dL*or2(n_r)*visc(n_r)*orho1(n_r)*      (      &
                  &                                              ddddw(lm,n_r)  &
                  &        +two*( dLvisc(n_r)-beta(n_r) ) * work_LMloc(lm,n_r)  &
                  &    +( ddLvisc(n_r)-two*dbeta(n_r)+dLvisc(n_r)*dLvisc(n_r)+  &
                  &       beta(n_r)*beta(n_r)-three*dLvisc(n_r)*beta(n_r)-two*  &
                  &       or1(n_r)*(dLvisc(n_r)+beta(n_r))-two*or2(n_r)*dL ) *  &
                  &                                                ddw(lm,n_r)  &
                  &    +( -ddbeta(n_r)-dbeta(n_r)*(two*dLvisc(n_r)-beta(n_r)+   &
                  &       two*or1(n_r))-ddLvisc(n_r)*(beta(n_r)+two*or1(n_r))+  &
                  &       beta(n_r)*beta(n_r)*(dLvisc(n_r)+two*or1(n_r))-       &
                  &       beta(n_r)*(dLvisc(n_r)*dLvisc(n_r)-two*or2(n_r))-     &
                  &       two*dLvisc(n_r)*or1(n_r)*(dLvisc(n_r)-or1(n_r))+      &
                  &       two*(two*or1(n_r)+beta(n_r)-dLvisc(n_r))*or2(n_r)*dL) &
                  &                                *                dw(lm,n_r)  &
                  &    + dL*or2(n_r)* ( two*dbeta(n_r)+ddLvisc(n_r)+            &
                  &      dLvisc(n_r)*dLvisc(n_r)-two*third*beta(n_r)*beta(n_r)+ &
                  &      dLvisc(n_r)*beta(n_r)+two*or1(n_r)*(two*dLvisc(n_r)-   &
                  &      beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   w(lm,n_r) )

                  Buo=zero
                  if ( l_heat ) Buo=BuoFac*dL*or2(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo=Buo+ChemFac*dL*or2(n_r)*&
                  &                                rgrav(n_r)*xi(lm,n_r)

                  if ( l_parallel_solve ) then
                     dwdt%impl(lm,n_r,istage)=Dif
                  else
                     dwdt%impl(lm,n_r,istage)=Dif+Buo
                  end if

                  if ( l1 /= 0 .and. lPressNext .and. &
                  &    istage_loc==nstages) then
                     ! In the double curl formulation, we can estimate the pressure
                     ! if required.
                     p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                     &            -one/dt_loc*(dw(lm,n_r)-dwold(lm,n_r))+   &
                     &              hdif_V(l1)*visc(n_r)* ( work_LMloc(lm,n_r)     &
                     &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                     &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                     &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                     &                                              ) * dw(lm,n_r) &
                     &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                     &                     +dLvisc(n_r) )   *            w(lm,n_r) )
                  end if

                  if ( lRmsNext .and. istage_loc==nstages ) then
                     !-- In case RMS force balance is required, one needs to also
                     !-- compute the classical diffusivity that is used in the non
                     !-- double-curl version
                     Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                     &    +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                     &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                     &       beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                     &                                     *       w(lm,n_r) )
                     DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                     DifPolLMr(lm,n_r)  =r(n_r)**2/dL * Dif
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end do
#endif

         else

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do collapse(2) private(l1,m1,Dif,Buo,Pre,dL)
#else
            !$omp do private(n_r,lm,l1,m1,dL,Dif,Buo,Pre)
#endif
            do n_r=n_r_top,n_r_bot
               do lm=lmStart_00,ulm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r)*(      ddw(lm,n_r)   &
                  &    +(two*dLvisc(n_r)-third*beta(n_r))*       dw(lm,n_r)   &
                  &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &      beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                               w(lm,n_r)  )
                  Pre = -dp(lm,n_r)+beta(n_r)*p(lm,n_r)
                  Buo = zero
                  if ( l_heat )  Buo=BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo=Buo+ChemFac*rho0(n_r)* &
                  &                              rgrav(n_r)*xi(lm,n_r)
                  if ( l_parallel_solve ) then
                     dwdt%impl(lm,n_r,istage)=Pre+Dif
                  else
                     dwdt%impl(lm,n_r,istage)=Pre+Dif+Buo
                  end if
                  dpdt%impl(lm,n_r,istage)=               dL*or2(n_r)*p(lm,n_r) &
                  &            + hdif_V(l1)*visc(n_r)*dL*or2(n_r)               &
                  &                                     * ( -work_LMloc(lm,n_r) &
                  &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                           ) *    dw(lm,n_r) &
                  &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
                  &                     +dLvisc(n_r) )   *           w(lm,n_r)  )
                  if ( lRmsNext .and. istage_loc==nstages ) then
                     DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                     DifPolLMr(lm,n_r)  =r(n_r)**2/dL * Dif
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#else
            !$omp end do
#endif

         end if

      end if

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext .and. l_double_curl ) then
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc, use_gpu=.true.)
#ifndef WITH_OMP_GPU
         !$omp barrier
#endif
      end if

#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

   end subroutine get_pol_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_pol_rhs_imp_ghost(wg, dw, ddw, p, dp, dwdt, tscheme, istage, &
              &                    l_calc_lin, lPressNext, lRmsNext, dp_expl)
      !
      ! This subroutine computes the derivatives of w and p and assemble the
      ! implicit stage if needed.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: l_calc_lin
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      complex(cp),         intent(in) :: dp_expl(lm_max,nRstart:nRstop)

      !-- Output variables
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(inout) :: wg(lm_max,nRstart-2:nRstop+2)
      complex(cp),       intent(inout) :: p(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: dp(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: dw(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: ddw(lm_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: Dif
      complex(cp), allocatable :: work_Rloc(:,:), dddw_Rloc(:,:)
      integer :: n_r, l, m, lm, start_lm, stop_lm
      real(cp) :: dL

      allocate(work_Rloc(lm_max,nRstart:nRstop), dddw_Rloc(lm_max,nRstart:nRstop))
      work_Rloc = zero; dddw_Rloc = zero
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: work_Rloc, dddw_Rloc)
      !$omp target update to(work_Rloc, dddw_Rloc)
#endif

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
      call dct_counter%start_count()
      call get_ddddr_ghost(wg, dw, ddw, dddw_Rloc, work_Rloc, lm_max, start_lm, &
           &               stop_lm, nRstart, nRstop, rscheme_oc)
      call dct_counter%stop_count()
#else
      !$omp parallel default(shared)  private(start_lm, stop_lm, n_r, lm, l, dL)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddddr_ghost(wg, dw, ddw, dddw_Rloc, work_Rloc, lm_max, start_lm, &
           &               stop_lm, nRstart, nRstop, rscheme_oc)
      !$omp single
      call dct_counter%stop_count()
      !$omp end single
      !$omp barrier
#endif

      if ( istage == 1 ) then
#ifdef WITH_OMP_GPU
#else
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               dL = real(l*(l+1),cp)
               dwdt%old(lm,n_r,istage)=dL*or2(n_r)* ( -orho1(n_r)*(  &
               &                   ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
               &                            dL*or2(n_r)* wg(lm,n_r) ) )
            end do
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if

      if ( l_calc_lin ) then

#ifdef WITH_OMP_GPU
#else
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l=st_map%lm2l(lm)
               if ( l /= 0 ) then
                  dL=real(l*(l+1),cp)

                  dwdt%impl(lm,n_r,istage)=-hdif_V(l)*dL*or2(n_r)*visc(n_r)         &
                  &                                                   *orho1(n_r)*( &
                  &                                              work_Rloc(lm,n_r)  &
                  &             +two*( dLvisc(n_r)-beta(n_r) ) * dddw_Rloc(lm,n_r)  &
                  &        +( ddLvisc(n_r)-two*dbeta(n_r)+dLvisc(n_r)*dLvisc(n_r)+  &
                  &           beta(n_r)*beta(n_r)-three*dLvisc(n_r)*beta(n_r)-two*  &
                  &           or1(n_r)*(dLvisc(n_r)+beta(n_r))-two*or2(n_r)*dL ) *  &
                  &                                                    ddw(lm,n_r)  &
                  &        +( -ddbeta(n_r)-dbeta(n_r)*(two*dLvisc(n_r)-beta(n_r)+   &
                  &           two*or1(n_r))-ddLvisc(n_r)*(beta(n_r)+two*or1(n_r))+  &
                  &           beta(n_r)*beta(n_r)*(dLvisc(n_r)+two*or1(n_r))-       &
                  &           beta(n_r)*(dLvisc(n_r)*dLvisc(n_r)-two*or2(n_r))-     &
                  &           two*dLvisc(n_r)*or1(n_r)*(dLvisc(n_r)-or1(n_r))+      &
                  &           two*(two*or1(n_r)+beta(n_r)-dLvisc(n_r))*or2(n_r)*dL) &
                  &                                    *                dw(lm,n_r)  &
                  &        + dL*or2(n_r)* ( two*dbeta(n_r)+ddLvisc(n_r)+            &
                  &          dLvisc(n_r)*dLvisc(n_r)-two*third*beta(n_r)*beta(n_r)+ &
                  &          dLvisc(n_r)*beta(n_r)+two*or1(n_r)*(two*dLvisc(n_r)-   &
                  &          beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   wg(lm,n_r) )
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if
#ifndef WITH_OMP_GPU
      !$omp end parallel
#endif

      if ( (tscheme%istage==tscheme%nstages .and. lRmsNext) .or. lPressNext ) then
         !-- Recompute third derivative to have the boundary point right
#ifdef WITH_OMP_GPU
         !$omp target update to(ddw, dddw_Rloc)
#endif
         call get_dr_Rloc(ddw, dddw_Rloc, lm_max, nRstart, nRstop, n_r_max, rscheme_oc )
#ifdef WITH_OMP_GPU
         !$omp target update from(dddw_Rloc)
#endif

#ifdef WITH_OMP_GPU
         start_lm=1; stop_lm=lm_max
#else
         !$omp parallel default(shared) private(start_lm, stop_lm, n_r, lm, l, m, dL, Dif) &
         !$omp reduction(+:DifPol2hInt)
         start_lm=1; stop_lm=lm_max
         call get_openmp_blocks(start_lm,stop_lm)
#endif
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l=st_map%lm2l(lm)
               m=st_map%lm2m(lm)
               dL=real(l*(l+1),cp)
               if ( l == 0 ) cycle

               if ( lPressNext ) then
                  ! In the double curl formulation, we can estimate the pressure
                  ! if required.
                  p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                  &            -one/tscheme%dt(1)*(dw(lm,n_r)-dwold(lm,n_r))+   &
                  &              hdif_V(l)*visc(n_r)* (   dddw_Rloc(lm,n_r)     &
                  &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                              ) * dw(lm,n_r) &
                  &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                  &                     +dLvisc(n_r) )  *           wg(lm,n_r) )
               end if

               if ( lRmsNext ) then
                  !-- In case RMS force balance is required, one needs to also
                  !-- compute the classical diffusion that is used in the non
                  !-- double-curl version
                  Dif =  hdif_V(l)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                  &    +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                  &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                  &       beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                  &                                     *       wg(lm,n_r) )
               end if

               if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                  DifPol2hInt(l,n_r)=DifPol2hInt(l,n_r)+r(n_r)**2*cc2real(Dif,m)
                  DifPolLMr(lm,n_r) =r(n_r)**2/dL * Dif
               end if
            end do
         end do
#ifndef WITH_OMP_GPU
         !$omp end parallel
#endif
      end if

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext ) then
#ifdef WITH_OMP_GPU
         !$omp target update to(p, dp)
#endif
         call get_dr_Rloc(p, dp, lm_max, nRstart, nRstop, n_r_max, rscheme_oc )
#ifdef WITH_OMP_GPU
         !$omp target update from(dp)
#endif
      end if

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: work_Rloc, dddw_Rloc)
#endif
      deallocate(work_Rloc, dddw_Rloc)

   end subroutine get_pol_rhs_imp_ghost
!------------------------------------------------------------------------------
   subroutine assemble_pol(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, dp_expl, &
              &            tscheme, lPressNext, lRmsNext)
      !
      ! This subroutine is used to assemble w and dw/dr when IMEX RK time schemes
      ! which necessitate an assembly stage are employed. Robin-type boundary
      ! conditions are enforced using Canuto (1986) approach.
      !

      !-- Input variables
      complex(cp),         intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: dp_expl(llm:ulm,n_r_max)
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext

      !-- Output variable
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: dp(llm:ulm,n_r_max)

      !-- Local variables
      complex(cp) :: Dif, Buo
      real(cp) :: fac_top, fac_bot
      integer :: n_r_top, n_r_bot, l1, m1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:), lm2m(:)
      real(cp) :: dL

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      if ( m_min == 0 ) then
         lmStart_00=max(2,llm)
      else
         lmStart_00=llm
      end if

#ifdef WITH_OMP_GPU
      !$omp target update to(dwdt)
#endif
      call tscheme%assemble_imex(work_LMloc, dwdt)
#ifdef WITH_OMP_GPU
      !$omp target update from(work_LMloc)
#endif
      if ( l_double_curl) then
         call get_pol(w, work_LMloc)
      else
#ifdef WITH_OMP_GPU
         !$omp target update to(dpdt)
#endif
         call tscheme%assemble_imex(ddw, dpdt) ! Use ddw as a work array
#ifdef WITH_OMP_GPU
      !$omp target update from(ddw)
#endif
      end if

      if ( l_double_curl) then
#ifdef WITH_OMP_GPU
         start_lm=lmStart_00; stop_lm=ulm
#else
         !$omp parallel default(shared)  private(start_lm,stop_lm,n_r,lm,l1,dL,m1,Dif,Buo)
         start_lm=lmStart_00; stop_lm=ulm
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier
#endif

#ifdef WITH_OMP_GPU
         !$omp target update to(dw, ddw, work_LMloc, ddddw)
         !$omp target update to(w)
#endif

#ifdef WITH_OMP_GPU
         call dct_counter%start_count()
#else
         !$omp single
         call dct_counter%start_count()
         !$omp end single
#endif

#ifdef WITH_OMP_GPU
         call get_ddr( w, dw, ddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc,      &
              &       l_dct_in=.false. )
         call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc )
         call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,.true.)
#else
         call get_ddr( w, dw, ddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc,      &
              &       l_dct_in=.false. )
         call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc )
         call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
#endif

#ifdef WITH_OMP_GPU
         call dct_counter%stop_count()
#else
         !$omp barrier
         !$omp single
         call dct_counter%stop_count()
         !$omp end single
#endif

#ifdef WITH_OMP_GPU
         !$omp target update from(dw, ddw, work_LMloc, ddddw)
         !$omp target update from(w)
#endif

#ifdef WITH_OMP_GPU
#else
#endif
         do n_r=2,n_r_max-1
            do lm=start_lm,stop_lm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dwdt%old(lm,n_r,1)=-dL*or2(n_r)*orho1(n_r)*(         &
               &                  ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
               &                         dL*or2(n_r)* w(lm,n_r) )
            end do
         end do
#ifdef WITH_OMP_GPU
#else
#endif

         if ( tscheme%l_imp_calc_rhs(1) .or. lPressNext ) then
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1

#ifdef WITH_OMP_GPU
#else
#endif
            do n_r=n_r_top,n_r_bot
               do lm=start_lm,stop_lm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif=-hdif_V(l1)*dL*or2(n_r)*visc(n_r)*orho1(n_r)*      (      &
                  &                                              ddddw(lm,n_r)  &
                  &        +two*( dLvisc(n_r)-beta(n_r) ) * work_LMloc(lm,n_r)  &
                  &    +( ddLvisc(n_r)-two*dbeta(n_r)+dLvisc(n_r)*dLvisc(n_r)+  &
                  &       beta(n_r)*beta(n_r)-three*dLvisc(n_r)*beta(n_r)-two*  &
                  &       or1(n_r)*(dLvisc(n_r)+beta(n_r))-two*or2(n_r)*dL ) *  &
                  &                                                ddw(lm,n_r)  &
                  &    +( -ddbeta(n_r)-dbeta(n_r)*(two*dLvisc(n_r)-beta(n_r)+   &
                  &       two*or1(n_r))-ddLvisc(n_r)*(beta(n_r)+two*or1(n_r))+  &
                  &       beta(n_r)*beta(n_r)*(dLvisc(n_r)+two*or1(n_r))-       &
                  &       beta(n_r)*(dLvisc(n_r)*dLvisc(n_r)-two*or2(n_r))-     &
                  &       two*dLvisc(n_r)*or1(n_r)*(dLvisc(n_r)-or1(n_r))+      &
                  &       two*(two*or1(n_r)+beta(n_r)-dLvisc(n_r))*or2(n_r)*dL) &
                  &                                *                dw(lm,n_r)  &
                  &    + dL*or2(n_r)* ( two*dbeta(n_r)+ddLvisc(n_r)+            &
                  &      dLvisc(n_r)*dLvisc(n_r)-two*third*beta(n_r)*beta(n_r)+ &
                  &      dLvisc(n_r)*beta(n_r)+two*or1(n_r)*(two*dLvisc(n_r)-   &
                  &      beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   w(lm,n_r) )

                  Buo=zero
                  if ( l_heat ) Buo=BuoFac*dL*or2(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo=Buo+ChemFac*dL*or2(n_r)*&
                  &                              rgrav(n_r)*xi(lm,n_r)

                  if ( l_parallel_solve ) then
                     dwdt%impl(lm,n_r,1)=Dif
                  else
                     dwdt%impl(lm,n_r,1)=Dif+Buo
                  end if

                  if ( lPressNext ) then
                     ! In the double curl formulation, we can estimate the pressure
                     ! if required.
                     p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                     &            -one/tscheme%dt(1)*(dw(lm,n_r)-dwold(lm,n_r))+   &
                     &              hdif_V(l1)*visc(n_r)* ( work_LMloc(lm,n_r)     &
                     &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                     &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                     &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                     &                                              ) * dw(lm,n_r) &
                     &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                     &                     +dLvisc(n_r) )   *            w(lm,n_r) )
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
#else
#endif
         end if

         ! In case pressure is needed in the double curl formulation
         ! we also have to compute the radial derivative of p
         if ( lPressNext .and. l_double_curl ) then
#ifdef WITH_OMP_GPU
            !$omp target update to(p)
#endif
            call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
                 &       n_r_max, rscheme_oc, use_gpu=.true.)
#ifdef WITH_OMP_GPU
            !$omp target update from(p, dp)
#endif
#ifndef WITH_OMP_GPU
            !$omp barrier
#endif
         end if
#ifndef WITH_OMP_GPU
         !$omp end parallel
#endif

         if ( lRmsNext ) then
#ifdef WITH_OMP_GPU
            start_lm=lmStart_00; stop_lm=ulm
#else
            !$omp parallel default(shared)  private(start_lm,stop_lm,n_r,lm,l1,dL,m1,Dif) &
            !$omp reduction(+:DifPol2hInt)
            start_lm=lmStart_00; stop_lm=ulm
            call get_openmp_blocks(start_lm,stop_lm)
            !$omp barrier
#endif
            do n_r=n_r_cmb,n_r_icb
               do lm=start_lm,stop_lm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)
                  Dif = hdif_V(l1)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                  &    +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                  &    -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                  &       beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                  &                                     *       w(lm,n_r) )
                  DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                  DifPolLMr(lm,n_r) =r(n_r)**2/dL * Dif
               end do
            end do
#ifndef WITH_OMP_GPU
            !$omp end parallel
#endif
         end if

      else

#ifdef WITH_OMP_GPU
         start_lm=lmStart_00; stop_lm=ulm
#else
         !$omp parallel default(shared)  private(start_lm,stop_lm,n_r,lm,l1,dL,m1,Dif,Buo)
         start_lm=lmStart_00; stop_lm=ulm
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier
#endif

         !-- Now get the poloidal from the assembly
#ifdef WITH_OMP_GPU
#else
#endif
         do n_r=2,n_r_max-1
            do lm=start_lm,stop_lm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               dL = real(l1*(l1+1),cp)
               if ( m1 == 0 ) then
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*cmplx(real(ddw(lm,n_r)),0.0_cp,cp)
               else
                  w(lm,n_r) = r(n_r)*r(n_r)/dL*work_LMloc(lm,n_r)
                  dw(lm,n_r)=-r(n_r)*r(n_r)/dL*ddw(lm,n_r)
               end if
            end do
         end do
#ifdef WITH_OMP_GPU
#else
#endif

         !-- Non-penetration: u_r=0 -> w_lm=0 on both boundaries
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=start_lm,stop_lm
            w(lm,1)      =zero
            w(lm,n_r_max)=zero
         end do
#ifdef WITH_OMP_GPU
#else
#endif

         !-- Other boundary condition: stress-free or rigid
         if ( l_full_sphere ) then
            if ( ktopv == 1 ) then ! Stress-free
               fac_top=-two*or1(1)-beta(1)
#ifdef WITH_OMP_GPU
#else
#endif
               do lm=start_lm,stop_lm
                  l1 = lm2l(lm)
                  if ( l1 == 1 ) then
                     call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
                  else
                     call rscheme_oc%robin_bc(one, fac_top, zero, one, 0.0_cp, zero, dw(lm,:))
                  end if
               end do
#ifdef WITH_OMP_GPU
#else
#endif
            else
#ifdef WITH_OMP_GPU
#else
#endif
               do lm=start_lm,stop_lm
                  l1 = lm2l(lm)
                  if ( l1 == 1 ) then
                     dw(lm,1)      =zero
                     dw(lm,n_r_max)=zero
                  else
                     call rscheme_oc%robin_bc(0.0_cp, one, zero, one, 0.0_cp, zero, dw(lm,:))
                  end if
               end do
            end if
         else ! Spherical shell
            if ( ktopv /= 1 .and. kbotv /= 1 ) then ! Rigid at both boundaries
               do lm=start_lm,stop_lm
                  dw(lm,1)      =zero
                  dw(lm,n_r_max)=zero
               end do
            else if ( ktopv /= 1 .and. kbotv == 1 ) then ! Rigid top/Stress-free bottom
               fac_bot=-two*or1(n_r_max)-beta(n_r_max)
               do lm=start_lm,stop_lm
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, one, fac_bot, zero, dw(lm,:))
               end do
            else if ( ktopv == 1 .and. kbotv /= 1 ) then ! Rigid bottom/Stress-free top
               fac_top=-two*or1(1)-beta(1)
               do lm=start_lm,stop_lm
                  call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
               end do
            else if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free at both boundaries
               fac_bot=-two*or1(n_r_max)-beta(n_r_max)
               fac_top=-two*or1(1)-beta(1)
               do lm=start_lm,stop_lm
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, zero, dw(lm,:))
               end do
            end if

         end if

#ifdef WITH_OMP_GPU
         !$omp target update to(ddw, work_LMloc)
         !$omp target update to(dw)
         call dct_counter%start_count()
         call get_ddr( dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &         stop_lm-llm+1, n_r_max, rscheme_oc)
         call dct_counter%stop_count()
         !$omp target update from(ddw, work_LMloc)
#else
         !$omp single
         call dct_counter%start_count()
         !$omp end single
         call get_ddr( dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &         stop_lm-llm+1, n_r_max, rscheme_oc)
         !$omp barrier
         !$omp single
         call dct_counter%stop_count()
         !$omp end single
#endif

#ifdef WITH_OMP_GPU
#else
#endif
         do n_r=2,n_r_max-1
            do lm=start_lm,stop_lm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dwdt%old(lm,n_r,1)= dL*or2(n_r)*w(lm,n_r)
               dpdt%old(lm,n_r,1)=-dL*or2(n_r)*dw(lm,n_r)
            end do
         end do
#ifdef WITH_OMP_GPU
#else
#endif

         if ( tscheme%l_imp_calc_rhs(1) ) then
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1

#ifdef WITH_OMP_GPU
#else
#endif
            do n_r=n_r_top,n_r_bot
               do lm=start_lm,stop_lm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif=hdif_V(l1)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)   &
                  &   +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
                  &   -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &     beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                              w(lm,n_r)  )
                  Buo=zero
                  if ( l_heat )  Buo=BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo=Buo+ChemFac*rho0(n_r)* &
                  &                              rgrav(n_r)*xi(lm,n_r)
                  if ( l_parallel_solve ) then
                     dwdt%impl(lm,n_r,1)=Dif+Buo
                  else
                     dwdt%impl(lm,n_r,1)=Dif+Buo
                  end if
                  dpdt%impl(lm,n_r,1)=hdif_V(l1)*visc(n_r)*dL*or2(n_r)*         &
                  &                                       ( -work_LMloc(lm,n_r) &
                  &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                           ) *    dw(lm,n_r) &
                  &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
                  &                     +dLvisc(n_r) )   *            w(lm,n_r) )
               end do
            end do
#ifdef WITH_OMP_GPU
#else
#endif
         end if
#ifndef WITH_OMP_GPU
         !$omp end parallel
#endif

         if ( lRmsNext ) then
#ifdef WITH_OMP_GPU
            start_lm=lmStart_00; stop_lm=ulm
#else
            !$omp parallel default(shared) private(start_lm,stop_lm,n_r,lm,l1,dL,m1,Dif)
            start_lm=lmStart_00; stop_lm=ulm
            call get_openmp_blocks(start_lm,stop_lm)
            !$omp barrier
#endif
            do n_r=n_r_cmb,n_r_icb
               do lm=start_lm,stop_lm
                  l1=lm2l(lm)
                  m1=lm2m(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif=hdif_V(l1)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)   &
                  &   +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
                  &   -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &     beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                              w(lm,n_r)  )
                  DifPol2hInt(l1,n_r)=DifPol2hInt(l1,n_r)+r(n_r)**2*cc2real(Dif,m1)
                  DifPolLMr(lm,n_r) =r(n_r)**2/dL * Dif
               end do
            end do
#ifndef WITH_OMP_GPU
            !$omp end parallel
#endif
         end if

      end if

   end subroutine assemble_pol
!------------------------------------------------------------------------------
   subroutine assemble_pol_Rloc(block_sze, nblocks, w, dw, ddw, p, dp, dwdt, dp_expl, &
              &                 tscheme, lPressNext, lRmsNext)
      !
      ! This subroutine is used to assemble w and dw/dr when IMEX RK time schemes
      ! which necessitate an assembly stage are employed. Robin-type boundary
      ! conditions are enforced using Canuto (1986) approach.
      !

      !-- Input variables
      integer,             intent(in) :: block_sze, nblocks
      complex(cp),         intent(in) :: dp_expl(lm_max,nRstart:nRstop)
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext

      !-- Output variables
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(inout) :: w(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: dw(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: ddw(lm_max,nRstart:nRstop)
      complex(cp),       intent(inout) :: p(lm_max,nRstart:nRstop)
      complex(cp),       intent(inout) :: dp(lm_max,nRstart:nRstop)

      !-- Local variables
      integer :: nlm_block, start_lm, stop_lm, req, tag, lms_block
      integer :: n_r, lm, l
      complex(cp), allocatable :: work_Rloc(:,:)
      complex(cp) :: work_ghost(lm_max, nRstart-1:nRstop+1)
      integer :: array_of_requests(4*nblocks)

      !-- LU factorisation of the matrix if needed
      if ( .not. l_ellMat(1) ) then
         call get_elliptic_mat_Rdist(ellMat_FD)
         l_ellMat(:)=.true.
      end if

      allocate(work_Rloc(lm_max,nRstart:nRstop))
      work_Rloc(:,:) = zero
#ifdef WITH_OMP_GPU
      !$omp target enter data map(alloc: work_Rloc)
      !$omp target update to(work_Rloc)
#endif

      !-- First assemble IMEX to get an r.h.s. stored in work_Rloc
#ifdef WITH_OMP_GPU
      !$omp target update to(dwdt)
#endif
      call tscheme%assemble_imex(work_Rloc, dwdt)
#ifdef WITH_OMP_GPU
      !$omp target update from(work_Rloc)
#endif

#ifdef WITH_MPI
      array_of_requests(:)=MPI_REQUEST_NULL
#endif

      !-- Now solve to finally get w
#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
#else
      !$omp parallel default(shared) private(tag, req, start_lm, stop_lm, lm, n_r, nlm_block, lms_block)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier
#endif

      !-- Non-penetration boundary condition
      if ( nRstart==n_r_cmb ) then
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=start_lm,stop_lm
            work_Rloc(lm,n_r_cmb)=zero
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if
      if ( nRstop==n_r_icb ) then
#ifdef WITH_OMP_GPU
#else
#endif
         do lm=start_lm,stop_lm
            work_Rloc(lm,n_r_icb)=zero
         end do
#ifdef WITH_OMP_GPU
#else
#endif
      end if

      !-- Now copy into an array with proper ghost zones
      call bulk_to_ghost(work_Rloc, work_ghost, 1, nRstart, nRstop, lm_max, start_lm, &
           &             stop_lm)

      tag = 0
      req=1

      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
#ifdef WITH_OMP_GPU
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
#else
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier
#endif
         call ellMat_FD%solver_up(work_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
              &                   array_of_requests, req, lms_block, nlm_block)
         tag = tag+1
      end do

      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
#ifdef WITH_OMP_GPU
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
#else
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier
#endif
         call ellMat_FD%solver_dn(work_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
              &                   array_of_requests, req, lms_block, nlm_block)
         tag = tag+1
      end do

#ifdef WITH_OMP_GPU
#else
      !$omp master
#endif
      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze

         call ellMat_FD%solver_finish(work_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                    tag, array_of_requests, req)
         tag = tag+1
      end do

#ifdef WITH_MPI
      call MPI_Waitall(req-1, array_of_requests(1:req-1), MPI_STATUSES_IGNORE, ierr)
      if ( ierr /= MPI_SUCCESS ) call abortRun('MPI_Waitall failed in assemble_pol_Rloc')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
#ifdef WITH_OMP_GPU
#else
      !$omp end master
      !$omp barrier
#endif

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
#else
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier
#endif
      do n_r=nRstart-1,nRstop+1
         do lm=start_lm,stop_lm
            w_ghost(lm,n_r)=work_ghost(lm,n_r)
         end do
      end do
#ifdef WITH_OMP_GPU
#else
      !$omp end parallel
#endif

      ! nRstart-1 and nRstop+1 are already known, only the next one is not known
      ! call exch_ghosts(w_ghost, lm_max, nRstart-1, nRstop+1, 1)
      ! Apparently it yields some problems, not sure why yet
      call exch_ghosts(w_ghost, lm_max, nRstart, nRstop, 2)
      call fill_ghosts_W(w_ghost, p0_ghost, .false.)
      call get_pol_rhs_imp_ghost(w_ghost, dw, ddw, p, dp, dwdt, tscheme, 1, &
           &                     tscheme%l_imp_calc_rhs(1), lPressNext,     &
           &                     lRmsNext, dp_expl)

#ifdef WITH_OMP_GPU
      start_lm=1; stop_lm=lm_max
#else
      !$omp parallel default(shared) private(start_lm,stop_lm,n_r,lm,l)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier
#endif
      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l /= 0 ) then
               w(lm,n_r)=w_ghost(lm,n_r)
            end if
         end do
      end do
#ifdef WITH_OMP_GPU
#else
      !$omp end parallel
#endif

#ifdef WITH_OMP_GPU
      !$omp target exit data map(delete: work_Rloc)
#endif
      deallocate(work_Rloc)

   end subroutine assemble_pol_Rloc
!------------------------------------------------------------------------------
   subroutine get_wpMat(tscheme,l,hdif,wpMat,wpMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  ``wpMat`` for the Navier-Stokes equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme
      real(cp),            intent(in) :: hdif    ! hyperdiffusion
      integer,             intent(in) :: l       ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: wpMat
      real(cp), intent(out) :: wpMat_fac(2*n_r_max,2)

      !-- local variables:
      integer :: nR,nR_out,nR_p,nR_out_p,n_r_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: wimp_lin

      dLh =real(l*(l+1),kind=cp)

      wimp_lin = tscheme%wimp_lin(1)

      !-- Now mode l>0

      !----- Boundary conditions, see above:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR_out=1,rscheme_oc%n_max
         nR_out_p=nR_out+n_r_max

         wpMat%dat(1,nR_out)        =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         wpMat%dat(1,nR_out_p)      =0.0_cp
         wpMat%dat(n_r_max,nR_out)  =rscheme_oc%rnorm* &
         &                       rscheme_oc%rMat(n_r_max,nR_out)
         wpMat%dat(n_r_max,nR_out_p)=0.0_cp

         if ( ktopv == 1 ) then  ! free slip !
            wpMat%dat(n_r_max+1,nR_out)= rscheme_oc%rnorm * (        &
            &                        rscheme_oc%d2rMat(1,nR_out) -   &
            &    (two*or1(1)+beta(1))*rscheme_oc%drMat(1,nR_out) )
         else                    ! no slip, note exception for l=1,m=0
            wpMat%dat(n_r_max+1,nR_out)=rscheme_oc%rnorm*    &
            &                       rscheme_oc%drMat(1,nR_out)
         end if
         wpMat%dat(n_r_max+1,nR_out_p)=0.0_cp

         if ( l_full_sphere ) then
            if ( l == 1 ) then
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%drMat(n_r_max,nR_out)
            else
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%d2rMat(n_r_max,nR_out)
            end if
         else
            if ( kbotv == 1 ) then  ! free slip !
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * (       &
               &                  rscheme_oc%d2rMat(n_r_max,nR_out) - &
               &      ( two*or1(n_r_max)+beta(n_r_max))*              &
               &                  rscheme_oc%drMat(n_r_max,nR_out))
            else                 ! no slip, note exception for l=1,m=0
               wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
               &                       rscheme_oc%drMat(n_r_max,nR_out)
            end if
         end if
         wpMat%dat(2*n_r_max,nR_out_p)=0.0_cp

      end do   !  loop over nR_out
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
#ifdef WITH_OMP_GPU
         !$omp target
#endif
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_p=nR_out+n_r_max
            wpMat%dat(1,nR_out)          =0.0_cp
            wpMat%dat(n_r_max,nR_out)    =0.0_cp
            wpMat%dat(n_r_max+1,nR_out)  =0.0_cp
            wpMat%dat(2*n_r_max,nR_out)  =0.0_cp
            wpMat%dat(1,nR_out_p)        =0.0_cp
            wpMat%dat(n_r_max,nR_out_p)  =0.0_cp
            wpMat%dat(n_r_max+1,nR_out_p)=0.0_cp
            wpMat%dat(2*n_r_max,nR_out_p)=0.0_cp
         end do
#ifdef WITH_OMP_GPU
         !$omp end target
#endif
      end if

      !----- Other points:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR_out=1,n_r_max
         nR_out_p=nR_out+n_r_max
         do nR=2,n_r_max-1
            ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
            nR_p=nR+n_r_max
            wpMat%dat(nR,nR_out)= rscheme_oc%rnorm *  (                     &
            &               dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)          &
            &         - wimp_lin*hdif*visc(nR)*dLh*or2(nR) * (              &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &        +(two*dLvisc(nR)-third*beta(nR))*                      &
            &                               rscheme_oc%drMat(nR,nR_out)     &
            &        -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)        &
            &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
            &          )                    *rscheme_oc%rMat(nR,nR_out)  )  )

            wpMat%dat(nR,nR_out_p)= rscheme_oc%rnorm*wimp_lin*(             &
            &                            rscheme_oc%drMat(nR,nR_out)        &
            &                  -beta(nR)* rscheme_oc%rMat(nR,nR_out))
            ! the following part gives sometimes very large
            ! matrix entries
            wpMat%dat(nR_p,nR_out)=rscheme_oc%rnorm * (                       &
            &                  -dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)       &
            &         -wimp_lin*hdif*visc(nR)*dLh*or2(nR)      *(             &
            &                                 - rscheme_oc%d3rMat(nR,nR_out)  &
            &          +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out)  &
            &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
            &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
            &                                    rscheme_oc%drMat(nR,nR_out)  &
            &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
            &          +two*third*beta(nR)   )*   rscheme_oc%rMat(nR,nR_out) ) )

            wpMat%dat(nR_p,nR_out_p)= -rscheme_oc%rnorm*wimp_lin*           &
            &                          dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)
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
         nR_p=nR+n_r_max
         wpMat%dat(nR,1)          =rscheme_oc%boundary_fac*wpMat%dat(nR,1)
         wpMat%dat(nR,n_r_max)    =rscheme_oc%boundary_fac*wpMat%dat(nR,n_r_max)
         wpMat%dat(nR,n_r_max+1)  =rscheme_oc%boundary_fac*wpMat%dat(nR,n_r_max+1)
         wpMat%dat(nR,2*n_r_max)  =rscheme_oc%boundary_fac*wpMat%dat(nR,2*n_r_max)
         wpMat%dat(nR_p,1)        =rscheme_oc%boundary_fac*wpMat%dat(nR_p,1)
         wpMat%dat(nR_p,n_r_max)  =rscheme_oc%boundary_fac*wpMat%dat(nR_p,n_r_max)
         wpMat%dat(nR_p,n_r_max+1)=rscheme_oc%boundary_fac*wpMat%dat(nR_p,n_r_max+1)
         wpMat%dat(nR_p,2*n_r_max)=rscheme_oc%boundary_fac*wpMat%dat(nR_p,2*n_r_max)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      ! compute the linesum of each line
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,2*n_r_max
         wpMat_fac(nR,1)=one/maxval(abs(wpMat%dat(nR,:)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each line by the linesum to regularize the matrix
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
      do n_r_out=1,2*n_r_max
         do nr=1,2*n_r_max
            wpMat%dat(nR,n_r_out) = wpMat%dat(nR,n_r_out)*wpMat_fac(nR,1)
         end do
      end do
      !$omp end target teams distribute parallel do
#else
      do nr=1,2*n_r_max
         wpMat%dat(nR,:) = wpMat%dat(nR,:)*wpMat_fac(nR,1)
      end do
#endif

      ! also compute the rowsum of each column
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,2*n_r_max
         wpMat_fac(nR,2)=one/maxval(abs(wpMat%dat(:,nR)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each row by the rowsum
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
      do nR=1,2*n_r_max
         do n_r_out=1,2*n_r_max
            wpMat%dat(n_r_out,nR) = wpMat%dat(n_r_out,nR)*wpMat_fac(nR,2)
         end do
      end do
      !$omp end target teams distribute parallel do
#else
      do nR=1,2*n_r_max
         wpMat%dat(:,nR) = wpMat%dat(:,nR)*wpMat_fac(nR,2)
      end do
#endif

#ifdef MATRIX_CHECK
      block

      integer ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
      real(cp) :: work(8*n_r_max),anorm,linesum,rcond
      real(cp) :: temp_wpMat(2*n_r_max,2*n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
      logical :: first_run=.true.

      write(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

#ifdef WITH_OMP_GPU
      !$omp target update from(wpMat%dat)
#endif

      do i=1,2*n_r_max
         do j=1,2*n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") wpMat%dat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_wpMat=wpMat%dat
      anorm = 0.0_cp
      do i=1,2*n_r_max
         linesum = 0.0_cp
         do j=1,2*n_r_max
            linesum = linesum + abs(temp_wpMat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(2*n_r_max,2*n_r_max,temp_wpMat,2*n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',2*n_r_max,temp_wpMat,2*n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of wpMat for l=",l," is ",rcond

      end block
#endif

#ifdef WITH_OMP_GPU
      if(.not. wpMat%gpu_is_used) then
         !$omp target update from(wpMat%dat)
      end if
#endif

#ifdef WITH_OMP_GPU
      if(.not. wpMat%gpu_is_used) then
         call wpMat%prepare(info)
      else
         call wpMat%prepare(info, handle, devInfo)
      end if
#else
      call wpMat%prepare(info)
#endif

      if ( info /= 0 ) call abortRun('Singular matrix wpMat!')

   end subroutine get_wpMat
!------------------------------------------------------------------------------
   subroutine get_elliptic_mat(l, ellMat)
      !
      !  Purpose of this subroutine is to contruct the matrix needed
      !  for the derivation of w for the time advance of the poloidal equation
      !  if the double curl form is used.
      !

      !-- Input variables:
      integer, intent(in) :: l       ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: ellMat

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: dLh
#ifndef WITH_OMP_GPU
      real(cp) :: dat(n_r_max,n_r_max)
#endif

      dLh =real(l*(l+1),kind=cp)

      dat(:,:) = 0.0_cp

      !----- Boundary conditions:
      dat(1,:)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
      if ( ktopv == 1 ) then
         dat(2,:)=rscheme_oc%rnorm*(rscheme_oc%d2rMat(1,:)-(beta(1)+two*or1(1))* &
         &                          rscheme_oc%drMat(1,:))
      else
         dat(2,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if
      if ( kbotv == 1 ) then
         dat(n_r_max-1,:)=rscheme_oc%rnorm*(rscheme_oc%d2rMat(n_r_max,:)-&
         &                (beta(n_r_max)+two*or1(n_r_max))*rscheme_oc%drMat(n_r_max,:))
      else
         dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      end if

#ifdef WITH_OMP_GPU
      !$omp target update to(dat)
#endif

      !----- Bulk points:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR_out=1,n_r_max
         do nR=3,n_r_max-2
            dat(nR,nR_out)= -rscheme_oc%rnorm*dLh*orho1(nR)*or2(nR)* (   &
            &                              rscheme_oc%d2rMat(nR,nR_out)  &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out)  &
            &                   -dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )
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

#ifdef WITH_OMP_GPU
      if( .not. ellMat%gpu_is_used) then
         !$omp target update from(dat)
      end if
#endif

      !-- Array copy
      call ellMat%set_data(dat)

      !-- Lu factorisation
      call ellMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix ellMat!')

   end subroutine get_elliptic_mat
!-----------------------------------------------------------------------------
   subroutine get_wMat(tscheme,l,hdif,wMat,wMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation. This matrix corresponds here to the
      !  radial component of the double-curl of the Navier-Stokes equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme  ! time scheme
      real(cp),            intent(in) :: hdif     ! hyperdiffusion
      integer,             intent(in) :: l        ! degree :math:`\ell`

      !-- Output variables:
      class(type_realmat), intent(inout) :: wMat
      real(cp), intent(out) :: wMat_fac(n_r_max,2)

      !-- local variables:
      integer :: nR, nR_out,n_r_out
      integer :: info
      real(cp) :: dLh
#ifndef WITH_OMP_GPU
      real(cp) :: dat(n_r_max,n_r_max)
#endif
      real(cp) :: wimp_lin

      dLh =real(l*(l+1),kind=cp)

      wimp_lin = tscheme%wimp_lin(1)

#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
      do nR=1,n_r_max
         !----- Boundary conditions:
         !-- Non-penetration condition at both boundaries
         dat(1,nR)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR)
         dat(n_r_max,nR)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR)

         !-- Second line and n_r_max-1 lines are used for the second BCs
         if ( ktopv == 1 ) then  ! free slip
            dat(2,nR)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(1,nR)-   &
            &               (two*or1(1)+beta(1))*rscheme_oc%drMat(1,nR) )
         else                    ! no slip
            dat(2,nR)=rscheme_oc%rnorm*rscheme_oc%drMat(1,nR)
         end if

         if ( l_full_sphere ) then
            if ( l == 1 ) then
               dat(n_r_max-1,nR)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,nR)
            else
               dat(n_r_max-1,nR)=rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,nR)
            end if
         else
            if ( kbotv == 1 ) then  ! free slip
               dat(n_r_max-1,nR)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(n_r_max,nR)-  &
               &                      (two*or1(n_r_max)+beta(n_r_max))*           &
               &                                    rscheme_oc%drMat(n_r_max,nR) )
            else                 ! no slip
               dat(n_r_max-1,nR)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,nR)
            end if
         end if
      end do
      !$omp end target teams distribute parallel do
#else
      !----- Boundary conditions:
      !-- Non-penetration condition at both boundaries
      dat(1,:)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)

      !-- Second line and n_r_max-1 lines are used for the second BCs
      if ( ktopv == 1 ) then  ! free slip
         dat(2,:)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(1,:)-   &
         &               (two*or1(1)+beta(1))*rscheme_oc%drMat(1,:) )
      else                    ! no slip
         dat(2,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         if ( l == 1 ) then
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%d2rMat(n_r_max,:)
         end if
      else
         if ( kbotv == 1 ) then  ! free slip
            dat(n_r_max-1,:)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(n_r_max,:)-  &
            &                      (two*or1(n_r_max)+beta(n_r_max))*           &
            &                                    rscheme_oc%drMat(n_r_max,:) )
         else                 ! no slip
            dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if
#endif

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
#ifdef WITH_OMP_GPU
         !$omp target
#endif
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)        =0.0_cp
            dat(2,nR_out)        =0.0_cp
            dat(n_r_max-1,nR_out)=0.0_cp
            dat(n_r_max,nR_out)  =0.0_cp
         end do
#ifdef WITH_OMP_GPU
         !$omp end target
#endif
      end if

      !----- Bulk points:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR_out=1,n_r_max
         do nR=3,n_r_max-2
            dat(nR,nR_out)=rscheme_oc%rnorm* (                              &
            &          -dLh*or2(nR)*orho1(nR)*(                             &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out) -   &
            &                    dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )   &
            &  + wimp_lin*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (           &
            &                               rscheme_oc%d4rMat(nR,nR_out)    &
            &   +two*(dLvisc(nR)-beta(nR))* rscheme_oc%d3rMat(nR,nR_out)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                               rscheme_oc%d2rMat(nR,nR_out)    &
            &    +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
            &       two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+        &
            &       beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-beta(nR)*    &
            &       (dLvisc(nR)*dLvisc(nR)-two*or2(nR))-two*dLvisc(nR)*     &
            &       or1(nR)*(dLvisc(nR)-or1(nR))+two*(two*or1(nR)+          &
            &       beta(nR)-dLvisc(nR))*dLh*or2(nR) ) *                    &
            &                                rscheme_oc%drMat(nR,nR_out)    &
            &    + dLh*or2(nR)*( two*dbeta(nR)+ddLvisc(nR)+dLvisc(nR)*      &
            &      dLvisc(nR)-two*third*beta(nR)*beta(nR)+dLvisc(nR)*       &
            &      beta(nR)+two*or1(nR)*(two*dLvisc(nR)-beta(nR)-three*     &
            &      or1(nR) ) + dLh*or2(nR) ) *                              &
            &                                rscheme_oc%rMat(nR,nR_out) ) )
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

      ! compute the linesum of each line
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         wMat_fac(nR,1)=one/maxval(abs(dat(nR,:)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each line by the linesum to regularize the matrix
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
      do n_r_out=1,n_r_max
         do nr=1,n_r_max
            dat(nR,n_r_out) = dat(nR,n_r_out)*wMat_fac(nR,1)
         end do
      end do
      !$omp end target teams distribute parallel do
#else
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*wMat_fac(nR,1)
      end do
#endif

      ! also compute the rowsum of each column
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=1,n_r_max
         wMat_fac(nR,2)=one/maxval(abs(dat(:,nR)))
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif
      ! now divide each row by the rowsum
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
      do nr=1,n_r_max
         do n_r_out=1,n_r_max
            dat(n_r_out,nR) = dat(n_r_out,nR)*wMat_fac(nR,2)
         end do
      end do
      !$omp end target teams distribute parallel do
#else
      do nR=1,n_r_max
         dat(:,nR) = dat(:,nR)*wMat_fac(nR,2)
      end do
#endif

#ifdef WITH_OMP_GPU
      if( .not. wMat%gpu_is_used) then
         !$omp target update from(dat)
      end if
#endif

      !-- Array copy
      call wMat%set_data(dat)

      !-- LU decomposition
#ifdef WITH_OMP_GPU
      if(.not. wMat%gpu_is_used) then
         call wMat%prepare(info)
      else
         call wMat%prepare(info, handle, devInfo)
      end if
#else
      call wMat%prepare(info)
#endif

      if ( info /= 0 ) call abortRun('Singular matrix wMat!')

   end subroutine get_wMat
!-----------------------------------------------------------------------------
   subroutine get_elliptic_mat_Rdist(ellMat)
      !
      !  Purpose of this subroutine is to contruct the matrix needed
      !  for the derivation of w for the time advance of the poloidal equation
      !  if the double curl form is used. This is the R-dist version.
      !

      !-- Output variables:
      type(type_tri_par), intent(inout) :: ellMat

      !-- local variables:
      integer :: nR, l
      real(cp) :: dLh

      !----- Bulk points:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR=2,n_r_max-1
         do l=1,l_max
            dLh =real(l*(l+1),kind=cp)
            ellMat%diag(l,nR)=-dLh*orho1(nR)*or2(nR)* ( rscheme_oc%ddr(nR,1) - &
            &                                   beta(nR)*rscheme_oc%dr(nR,1) - &
            &                                           dLh*or2(nR) )
            ellMat%up(l,nR)  =-dLh*orho1(nR)*or2(nR)* ( rscheme_oc%ddr(nR,2) - &
            &                                   beta(nR)*rscheme_oc%dr(nR,2) )
            ellMat%low(l,nR) =-dLh*orho1(nR)*or2(nR)* ( rscheme_oc%ddr(nR,0) - &
            &                                   beta(nR)*rscheme_oc%dr(nR,0) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !-- Non penetrative boundary condition
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do l=1,l_max
         ellMat%diag(l,1)      =one
         ellMat%up(l,1)        =0.0_cp
         ellMat%low(l,1)       =0.0_cp
         ellMat%diag(l,n_r_max)=one
         ellMat%up(l,n_r_max)  =0.0_cp
         ellMat%low(l,n_r_max) =0.0_cp
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !-- Lu factorisation
      call ellMat%prepare_mat()

   end subroutine get_elliptic_mat_Rdist
!-----------------------------------------------------------------------------
   subroutine get_wMat_Rdist(tscheme,hdif,wMat)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wMat_FD  for the NS equation. This matrix corresponds here to the
      !  radial component of the double-curl of the Navier-Stokes equation.
      !  This routine is used when parallel F.D. solvers are employed.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme  ! time scheme
      real(cp),            intent(in) :: hdif(0:l_max)     ! hyperdiffusion

      !-- Output variables:
      type(type_penta_par), intent(inout) :: wMat

      !-- local variables:
      integer :: nR, l
      real(cp) :: dLh, dr, fac
      real(cp) :: wimp_lin

      !-- Copie into local variable
      wimp_lin = tscheme%wimp_lin(1)

      !----- Bulk points (first and last lines always set for non-penetration condition)
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#else
      !$omp parallel default(shared) private(nR,l,dLh,dr,fac)
      !$omp do
#endif
      do nR=2,n_r_max-1
         do l=1,l_max
            dLh=real(l*(l+1),cp)
            wMat%diag(l,nR)=-dLh*or2(nR)*orho1(nR)*( rscheme_oc%ddr(nR,1)   &
            &              -beta(nR)*rscheme_oc%dr(nR,1)-   dLh*or2(nR) )   &
            &  +wimp_lin*orho1(nR)*hdif(l)*visc(nR)*dLh*or2(nR)*(&
            &                                     rscheme_oc%ddddr(nR,2)    &
            &          +two*(dLvisc(nR)-beta(nR))* rscheme_oc%dddr(nR,2)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                                       rscheme_oc%ddr(nR,1)    &
            &    +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
            &       two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+        &
            &       beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-beta(nR)*    &
            &       (dLvisc(nR)*dLvisc(nR)-two*or2(nR))-two*dLvisc(nR)*     &
            &       or1(nR)*(dLvisc(nR)-or1(nR))+two*(two*or1(nR)+          &
            &       beta(nR)-dLvisc(nR))*dLh*or2(nR) ) *                    &
            &                                        rscheme_oc%dr(nR,1)    &
            &    + dLh*or2(nR)*( two*dbeta(nR)+ddLvisc(nR)+dLvisc(nR)*      &
            &      dLvisc(nR)-two*third*beta(nR)*beta(nR)+dLvisc(nR)*       &
            &      beta(nR)+two*or1(nR)*(two*dLvisc(nR)-beta(nR)-three*     &
            &      or1(nR) ) + dLh*or2(nR) ) )
            wMat%low1(l,nR)=-dLh*or2(nR)*orho1(nR)*( rscheme_oc%ddr(nR,0)   &
            &                             -beta(nR)*rscheme_oc%dr(nR,0) )   &
            &  +wimp_lin*orho1(nR)*hdif(l)*visc(nR)*dLh*or2(nR)*(&
            &                                     rscheme_oc%ddddr(nR,1)    &
            &          +two*(dLvisc(nR)-beta(nR))* rscheme_oc%dddr(nR,1)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                                       rscheme_oc%ddr(nR,0)    &
            &    +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
            &       two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+        &
            &       beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-beta(nR)*    &
            &       (dLvisc(nR)*dLvisc(nR)-two*or2(nR))-two*dLvisc(nR)*     &
            &       or1(nR)*(dLvisc(nR)-or1(nR))+two*(two*or1(nR)+          &
            &       beta(nR)-dLvisc(nR))*dLh*or2(nR) ) *                    &
            &                                        rscheme_oc%dr(nR,0) )
            wMat%up1(l,nR)= -dLh*or2(nR)*orho1(nR)*( rscheme_oc%ddr(nR,2)   &
            &                             -beta(nR)*rscheme_oc%dr(nR,2) )   &
            &  +wimp_lin*orho1(nR)*hdif(l)*visc(nR)*dLh*or2(nR)*(&
            &                                     rscheme_oc%ddddr(nR,3)    &
            &          +two*(dLvisc(nR)-beta(nR))* rscheme_oc%dddr(nR,3)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                                       rscheme_oc%ddr(nR,2)    &
            &    +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
            &       two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+        &
            &       beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-beta(nR)*    &
            &       (dLvisc(nR)*dLvisc(nR)-two*or2(nR))-two*dLvisc(nR)*     &
            &       or1(nR)*(dLvisc(nR)-or1(nR))+two*(two*or1(nR)+          &
            &       beta(nR)-dLvisc(nR))*dLh*or2(nR) ) *                    &
            &                                        rscheme_oc%dr(nR,2) )
            wMat%low2(l,nR)=wimp_lin*orho1(nR)*hdif(l)*visc(nR)* &
            &               dLh*or2(nR) * (       rscheme_oc%ddddr(nR,0)    &
            &          +two*(dLvisc(nR)-beta(nR))* rscheme_oc%dddr(nR,0) )
            wMat%up2(l,nR)=wimp_lin*orho1(nR)*hdif(l)*visc(nR)* &
            &               dLh*or2(nR) * (       rscheme_oc%ddddr(nR,4)   &
            &          +two*(dLvisc(nR)-beta(nR))* rscheme_oc%dddr(nR,4) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !----- Boundary conditions:
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#else
      !$omp do
#endif
      do l=1,l_max
         !-- Non-penetration condition at both boundaries
         wMat%diag(l,1)=one
         wMat%low1(l,1)=0.0_cp
         wMat%low2(l,1)=0.0_cp
         wMat%up1(l,1) =0.0_cp
         wMat%up2(l,1) =0.0_cp

         wMat%diag(l,n_r_max)=one
         wMat%low1(l,n_r_max)=0.0_cp
         wMat%low2(l,n_r_max)=0.0_cp
         wMat%up1(l,n_r_max) =0.0_cp
         wMat%up2(l,n_r_max) =0.0_cp

         !-- Second part of B.C.
         if ( ktopv == 1 ) then ! free slip
            dr=r(2)-r(1)
            fac=(one-half*(two*or1(1)+beta(1))*dr)/(one+half*(two*or1(1)+beta(1))*dr)
            wMat%diag(l,2)=wMat%diag(l,2)-fac*wMat%low2(l,2)
            wMat%low1(l,2)=wMat%low1(l,2)+two/(one+half*(two*or1(1)+beta(1))*dr)*&
            &              wMat%low2(l,2)
         else ! No slip
            wMat%diag(l,2)=wMat%diag(l,2)+wMat%low2(l,2)
         end if

         if ( l_full_sphere ) then
            if ( l== 1 ) then ! dw=0
               wMat%diag(l,n_r_max-1)=wMat%diag(l,n_r_max-1)+wMat%up2(l,n_r_max-1)
            else ! ddw=0
               wMat%diag(l,n_r_max-1)=wMat%diag(l,n_r_max-1)-wMat%up2(l,n_r_max-1)
            end if
         else
            if ( kbotv == 1 ) then ! free-slip
               dr=r(n_r_max)-r(n_r_max-1)
               fac=(one+half*(two*or1(n_r_max)+beta(n_r_max))*dr)/ &
               &   (one-half*(two*or1(n_r_max)+beta(n_r_max))*dr)
               wMat%diag(l,n_r_max-1)=wMat%diag(l,n_r_max-1)-fac*wMat%up2(l,n_r_max-1)
               wMat%up1(l,n_r_max-1)=wMat%up1(l,n_r_max-1)+two*wMat%up2(l,n_r_max-1) &
               &                     /(one-half*(two*or1(n_r_max)+beta(n_r_max))*dr)
            else ! no slip
               wMat%diag(l,n_r_max-1)=wMat%diag(l,n_r_max-1)+wMat%up2(l,n_r_max-1)
            end if
         end if
      end do ! Loop over \ell
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#else
      !$omp end do
      !$omp end parallel
#endif

      call wMat%prepare_mat()

   end subroutine get_wMat_Rdist
!-----------------------------------------------------------------------------
   subroutine get_p0Mat(pMat)
      !
      ! This subroutine solves the linear problem of the spherically-symmetric
      ! pressure
      !

      !-- Output variables:
      class(type_realmat), intent(inout) :: pMat ! matrix

      !-- Local variables:
#ifndef WITH_OMP_GPU
      real(cp) :: dat(n_r_max,n_r_max), delr, work(n_r_max)
#else
      real(cp) :: delr, work(n_r_max)
#endif
      integer :: info, nCheb, nR_out, nR, nCheb_in

      !-- Bulk points
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do collapse(2)
#endif
      do nR_out=1,n_r_max
         do nR=2,n_r_max
            dat(nR,nR_out)= rscheme_oc%rnorm * (              &
            &                    rscheme_oc%drMat(nR,nR_out)- &
            &            beta(nR)*rscheme_oc%rMat(nR,nR_out) )
         end do
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
#endif

      !-- Boundary condition for spherically-symmetric pressure
      !-- The integral of rho' r^2 dr vanishes
      if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then

         work(:) = ThExpNb*ViscHeatFac*ogrun(:)*alpha0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work)
         work(:)      =work(:)*rscheme_oc%rnorm
         work(1)      =rscheme_oc%boundary_fac*work(1)
         work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)

         if ( rscheme_oc%version == 'cheb' ) then

#ifdef WITH_OMP_GPU
            !$omp target teams distribute parallel do
#endif
            do nCheb=1,rscheme_oc%n_max
               dat(1,nCheb)=0.0_cp
               do nCheb_in=1,rscheme_oc%n_max
                  if ( mod(nCheb+nCheb_in-2,2)==0 ) then
                     dat(1,nCheb)=dat(1,nCheb)+ &
                     &             ( one/(one-real(nCheb_in-nCheb,cp)**2)    + &
                     &               one/(one-real(nCheb_in+nCheb-2,cp)**2) )* &
                     &               work(nCheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do
#ifdef WITH_OMP_GPU
            !$omp end target teams distribute parallel do
#endif

         else

            !!-- In the finite differences case, we restrict the integral boundary
            !!-- condition to a trapezoidal rule of integration
            !do nR_out=2,rscheme_oc%n_max-1
            !   dat(1,nR_out)=half*work(nR_out)*( r(nR_out+1)-r(nR_out-1) )
            !end do
            !dat(1,1)=half*work(1)*(r(2)-r(1))
#ifdef WITH_OMP_GPU
            !$omp target
#endif
            dat(1,:)=rscheme_oc%rMat(1,:)
            delr = r(n_r_max)-r(n_r_max-1)

            !-- First-order on the last point
            dat(n_r_max,n_r_max)    =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)  =-one/delr
            dat(n_r_max,1:n_r_max-2)=0.0_cp
#ifdef WITH_OMP_GPU
            !$omp end target
#endif
         end if

      else

#ifdef WITH_OMP_GPU
         !$omp target
#endif
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
         if ( rscheme_oc%version == 'fd' ) then
            delr = r(n_r_max)-r(n_r_max-1)
            !-- First-order on the last point
            dat(n_r_max,n_r_max)  =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)=-one/delr
         end if
#ifdef WITH_OMP_GPU
         !$omp end target
#endif

      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros
#ifdef WITH_OMP_GPU
         !$omp target
#endif
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)=0.0_cp
         end do
#ifdef WITH_OMP_GPU
         !$omp end target
#endif
      end if

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

#ifdef WITH_OMP_GPU
      if( .not. pMat%gpu_is_used) then
         !$omp target update from(dat)
      end if
#endif

      !-- Array copy
      call pMat%set_data(dat)

      !---- LU decomposition:
#ifdef WITH_OMP_GPU
      if(.not. pMat%gpu_is_used) then
         call pMat%prepare(info)
      else
         call pMat%prepare(info, handle, devInfo)
      end if
#else
      call pMat%prepare(info)
#endif
      if ( info /= 0 ) call abortRun('! Singular matrix p0Mat!')

   end subroutine get_p0Mat
!-----------------------------------------------------------------------------
   subroutine get_p0Mat_Rdist(pMat)
      !
      ! This subroutine solves the linear problem of the spherically-symmetric
      ! pressure. This is the R-distributed variant of the function used when
      ! parallel F.D. solvers are employed
      !

      !-- Output variables:
      type(type_tri_par), intent(inout) :: pMat ! matrix

      !-- Local variables:
      real(cp) :: dr
      integer :: l, nR

      l=1
      !-- Bulk points
#ifdef WITH_OMP_GPU
      !$omp target teams distribute parallel do
#endif
      do nR=2,n_r_max-1
         pMat%diag(l,nR)=rscheme_oc%dr(nR,1)-beta(nR)
         pMat%low(l,nR) =rscheme_oc%dr(nR,0)
         pMat%up(l,nR)  =rscheme_oc%dr(nR,2)
      end do
#ifdef WITH_OMP_GPU
      !$omp end target teams distribute parallel do
      !$omp target update from(pMat%diag, pMat%low, pMat%up)
#endif

      !-- Boundary conditions for spherically-symmetric pressure
      pMat%diag(l,1)=one
      pMat%low(l,1) =0.0_cp
      pMat%up(l,1)  =0.0_cp

      !-- First order on the last point (no way around)
      dr = r(n_r_max)-r(n_r_max-1)
      pMat%diag(l,n_r_max)=one/dr-beta(n_r_max)
      pMat%low(l,n_r_max) =-one/dr
      pMat%up(l,n_r_max)  =0.0_cp

#ifdef WITH_OMP_GPU
      !$omp target update to(pMat%diag, pMat%low, pMat%up)
#endif

      !---- LU decomposition:
      call pMat%prepare_mat()

   end subroutine get_p0Mat_Rdist
!-----------------------------------------------------------------------------
end module updateWP_mod
