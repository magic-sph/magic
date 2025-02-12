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
   use constants, only: zero, one, two, three
   use fields, only: work_LMloc
   use mem_alloc, only: bytes_allocated
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
   class(type_realmat), pointer :: phiMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: phiMat_fac(:,:)
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

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( l_finite_diff ) then
            allocate( type_bandmat :: phiMat(nLMBs2(1+rank)) )

            if ( rscheme_oc%order == 2 .and. rscheme_oc%order_boundary <= 2 ) then ! Dirichelt BCs
               n_bands = rscheme_oc%order+1
            else
               n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
            end if

            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
         else
            allocate( type_densemat :: phiMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
         end if

#ifdef WITH_PRECOND_S
         allocate(phiMat_fac(n_r_max,nLMBs2(1+rank)))
         bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif

#ifdef WITHOMP
         maxThreads=omp_get_max_threads()
#else
         maxThreads=1
#endif
         allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
         &                 maxThreads*SIZEOF_DEF_COMPLEX
      else ! Parallel solvers are requested

         !-- Create matrix
         call phiMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Allocate an array with ghost zones
         allocate( phi_ghost(lm_max, nRstart-1:nRstop+1) )
         bytes_allocated=bytes_allocated + lm_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         phi_ghost(:,:)=zero

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

#ifdef WITH_PRECOND_S
         deallocate(phiMat_fac)
#endif
         deallocate( rhs1 )
      else
         call phiMat_FD%finalize()
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
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR             ! counts radial grid points
      integer :: n_r_out        ! counts cheb modes

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

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dphidt)

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
         !$omp private(lm,lm1,l1,m1,threadid) &
         !$omp private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)

         if ( .not. lPhimat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_phiMat(tscheme,l1,phiMat(nLMB2),phiMat_fac(:,nLMB2))
#else
            call get_phiMat(tscheme,l1,phiMat(nLMB2))
#endif
             lPhimat(l1)=.true.
         end if

         do iChunk=1,nChunks
            !$omp task default(shared) &
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

               if ( l1 == 0 ) then
                  rhs1(1,2*lm-1,threadid)      =phi_top
                  rhs1(1,2*lm,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lm-1,threadid)=phi_bot
                  rhs1(n_r_max,2*lm,threadid)  =0.0_cp
               else
                  rhs1(1,2*lm-1,threadid)      =0.0_cp
                  rhs1(1,2*lm,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lm-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lm,threadid)  =0.0_cp
               end if
               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=phiMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =phiMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif

            end do

            call phiMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                   2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                      rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                       0.0_cp,kind=cp)
                  end do
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

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            phi(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

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

      !$omp parallel default(shared) private(lm_start,lm_stop, nR, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(phi_ghost, dphidt, lm_start, lm_stop, 1)

      !-- Set boundary conditions
      if ( nRstart == n_r_cmb ) then
         nR=n_r_cmb
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
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)

            if ( l_full_sphere ) then
               if ( l > 0 ) then
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
      end if
      !$omp end parallel


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

      !$omp parallel default(shared) private(lm_start, lm_stop, l, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Handle upper boundary
      dr = r(2)-r(1)
      if ( nRstart == n_r_cmb ) then
         do lm=lm_start,lm_stop
            if ( ktopphi == 1 ) then
               phig(lm,nRstart-1)=two*phig(lm,nRstart)-phig(lm,nRstart+1)
            else
               phig(lm,nRstart-1)=phig(lm,nRstart+1)
            end if
         end do
      end if

      !-- Handle Lower boundary
      dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_icb ) then
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            if ( l_full_sphere ) then
               if ( l == 0 ) then
                  phig(lm,nRstop+1)=phig(lm,nRstop-1)
               else
                  phig(lm,nRstop+1)=two*phig(lm,nRstop)-phig(lm,nRstop-1)
               end if
            else ! Not a full sphere
               if ( kbotphi == 1 ) then
                  phig(lm,nRstop+1)=two*phig(lm,nRstop)-phig(lm,nRstop-1)
               else
                  phig(lm,nRstop+1)=phig(lm,nRstop-1)
               end if
            end if
         end do
      end if
      !$omp end parallel

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

      !$omp parallel default(shared) private(lm_start,lm_stop,nR,lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)

      !-- Array copy from phi_ghost to phi
      !!$omp parallel do simd collapse(2) schedule(simd:static)
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            phi(lm,nR)=phi_ghost(lm,nR)
         end do
      end do
      !!$omp end parallel do simd
      !$omp end parallel

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
      complex(cp) :: dphi(llm:ulm,n_r_max)
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

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(phi, dphi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(phi,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do
         do n_r=1,n_r_max
            dphidt%old(:,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phi(:,n_r)
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then

         !$omp do private(n_r,lm,l,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dphidt%impl(lm,n_r,istage)=phaseDiffFac*( work_LMloc(lm,n_r) + &
               &                                two*or1(n_r) * dphi(lm,n_r) - &
               &                                 dL*or2(n_r) *  phi(lm,n_r) )
            end do
         end do
         !$omp end do

      end if

      !$omp end parallel

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
      complex(cp) :: dphi(lm_max,nRstart:nRstop) ! Radial derivative of phase field
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      integer :: n_r, lm, start_lm, stop_lm, l
      real(cp) :: dL
      integer, pointer :: lm2l(:)

      lm2l(1:lm_max) => st_map%lm2l

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

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               dphidt%old(lm,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phig(lm,n_r)
            end do
         end do
      end if

      if ( l_calc_lin ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dphidt%impl(lm,n_r,istage)= phaseDiffFac * ( work_Rloc(lm,n_r) + &
               &                                 two*or1(n_r) *  dphi(lm,n_r) - &
               &                                  dL*or2(n_r) *  phig(lm,n_r) )
               if ( l_full_sphere .and. l==0 .and. n_r==n_r_icb ) then
                  dphidt%impl(lm,n_r,istage)=phaseDiffFac*three*work_Rloc(lm,n_r)
               end if
            end do
         end do
      end if
      !$omp end parallel

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

      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m)
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
      !$omp end do

      !-- Boundary conditions
      if ( l_full_sphere) then
         if ( ktopphi == 1 ) then ! Dirichlet
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   one, 0.0_cp, cmplx(phi_bot,0.0_cp,cp), &
                       &                   phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         else ! Neummann
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, &
                       &                   zero, phi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         end if

      else ! Spherical shell

         if ( ktopphi==1 .and. kbotphi==1 ) then
            !-- Boundary conditions: Dirichlet on both sides
            !$omp do private(lm,l)
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
            !$omp end do
         else if ( ktopphi==1 .and. kbotphi /= 1 ) then
            !$omp do private(lm,l)
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
            !$omp end do
         else if ( ktopphi/=1 .and. kbotphi == 1 ) then
            !$omp do private(lm,l)
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
            !$omp end do
         else if ( ktopphi/=1 .and. kbotphi /= 1 ) then
            !-- Boundary conditions: Neuman on both sides
            !$omp do private(lm)
            do lm=llm,ulm
               call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, zero, phi(lm,:))
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

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
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)

      call tscheme%assemble_imex(work_Rloc, dphidt)

      !$omp parallel default(shared) private(start_lm, stop_lm, l, m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier

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

      if ( ktopphi==1 .and. nRstart==n_r_cmb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l == 0 ) then
               phi(lm,nRstart)=phi_top
            else
               phi(lm,nRstart)=zero
            end if
         end do
      end if

      if ( kbotphi==1 .and. nRstop==n_r_icb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            if ( l == 0 ) then
               phi(lm,nRstop)=phi_bot
            else
               phi(lm,nRstop)=zero
            end if
         end do
      end if

      call bulk_to_ghost(phi, phi_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm)
      !$omp end parallel

      call exch_ghosts(phi_ghost, lm_max, nRstart, nRstop, 1)
      call fill_ghosts_Phi(phi_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_phase_rhs_imp_ghost(phi_ghost, dphidt, 1, tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_phase_Rloc
!------------------------------------------------------------------------------
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
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktopphi == 1 ) then ! Dirichlet
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else ! Neumann
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         if ( l == 0 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
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

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                                &
            &               5.0_cp/6.0_cp*stef*pr* rscheme_oc%rMat(nR,nR_out) - &
            &  tscheme%wimp_lin(1)*phaseDiffFac*(rscheme_oc%d2rMat(nR,nR_out) + &
            &                     two*or1(nR)*    rscheme_oc%drMat(nR,nR_out) - &
            &                     dLh*or2(nR)*     rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         phiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*phiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call phiMat%set_data(dat)

      !----- LU decomposition:
      call phiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix phiMat!')

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

      !----- Bulk points
      !$omp parallel default(shared) private(nR,l,dLh)
      !$omp do
      do nR=1,n_r_max
         do l=0,l_max
            dLh=real(l*(l+1),kind=cp)
            phiMat%diag(l,nR)=    5.0_cp/6.0_cp*stef*pr-            &
            &                  tscheme%wimp_lin(1)*phaseDiffFac*(   &
            &                                rscheme_oc%ddr(nR,1) + &
            &                     two*or1(nR)*rscheme_oc%dr(nR,1) - &
            &                                        dLh*or2(nR) )
            phiMat%low(l,nR)=-tscheme%wimp_lin(1)*phaseDiffFac*(    &
            &                                rscheme_oc%ddr(nR,0) + &
            &                 two*or1(nR)*    rscheme_oc%dr(nR,0) )
            phiMat%up(l,nR) =-tscheme%wimp_lin(1)*phaseDiffFac*(    &
            &                                rscheme_oc%ddr(nR,2) + &
            &                 two*or1(nR)*    rscheme_oc%dr(nR,2) )

            if ( l==0 .and. l_full_sphere .and. nR==n_r_icb ) then
               !-- Use L'Hôpital's rule to replace 2/r d/dr at the center
               !-- by 2*d^2/dr^2
               phiMat%diag(l,nR)=    5.0_cp/6.0_cp*stef*pr-           &
               &                  tscheme%wimp_lin(1)*phaseDiffFac*   &
               &                       three*   rscheme_oc%ddr(nR,1)
               phiMat%low(l,nR)=-tscheme%wimp_lin(1)*phaseDiffFac*    &
               &                       three*   rscheme_oc%ddr(nR,0)
               phiMat%up(l,nR) =-tscheme%wimp_lin(1)*phaseDiffFac*    &
               &                       three*   rscheme_oc%ddr(nR,2)
            end if
         end do
      end do
      !$omp end do

      !----- Boundary conditions:
      !$omp do
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
            if ( l == 0 ) then
               phiMat%low(l,n_r_max)=phiMat%up(l,n_r_max)+phiMat%low(l,n_r_max)
            else
               phiMat%diag(l,n_r_max)=one
               phiMat%up(l,n_r_max)  =0.0_cp
               phiMat%low(l,n_r_max) =0.0_cp
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
      !$omp end do
      !$omp end parallel

      !----- LU decomposition:
      call phiMat%prepare_mat()

   end subroutine get_phiMat_Rdist
!-----------------------------------------------------------------------------
end module updatePhi_mod
