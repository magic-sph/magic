module updateS_mod
   !
   ! This module handles the time advance of the entropy s.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: orho1, or1, or2, beta, dentropy0, rscheme_oc,  &
       &                       kappa, dLkappa, dLtemp0, temp0, r, l_R
   use physical_parameters, only: opr, kbots, ktops, stef
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: tops, bots
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_S
   use logic, only: l_anelastic_liquid, l_finite_diff, l_phase_field, &
       &            l_full_sphere, l_parallel_solve
   use parallel_mod
   use radial_der, only: get_ddr, get_dr, get_dr_Rloc, get_ddr_ghost, &
       &                 exch_ghosts, bulk_to_ghost
   use fields, only:  work_LMloc
   use constants, only: zero, one, two, three
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
   complex(cp), allocatable, public :: s_ghost(:,:)
   class(type_realmat), pointer :: sMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac(:,:)
#endif
   logical, public, allocatable :: lSmat(:)

   type(type_tri_par), public :: sMat_FD
   real(cp), allocatable :: fd_fac_top(:), fd_fac_bot(:)

   integer :: maxThreads

   public :: initialize_updateS, updateS, finalize_updateS, assemble_entropy,  &
   &         finish_exp_entropy, get_entropy_rhs_imp, finish_exp_entropy_Rdist,&
   &         prepareS_FD, updateS_FD, get_entropy_rhs_imp_ghost, fill_ghosts_S,&
   &         assemble_entropy_Rloc

contains

   subroutine initialize_updateS
      !
      ! This subroutine allocates the arrays involved in the time-advance of the
      ! entropy/temperature equation.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll,n_bands

      if ( .not. l_parallel_solve ) then

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( l_finite_diff ) then
            allocate( type_bandmat :: sMat(nLMBs2(1+rank)) )

            if ( ktops == 1 .and. kbots == 1 .and. rscheme_oc%order == 2 &
             &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed entropy at both boundaries
               n_bands = rscheme_oc%order+1
            else
               n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
            end if

            do ll=1,nLMBs2(1+rank)
               call sMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
         else
            allocate( type_densemat :: sMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call sMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
         end if

#ifdef WITH_PRECOND_S
         allocate(sMat_fac(n_r_max,nLMBs2(1+rank)))
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
      else

         !-- Create matrix
         call sMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Allocate an array with ghost zones
         allocate( s_ghost(lm_max, nRstart-1:nRstop+1) )
         bytes_allocated=bytes_allocated + lm_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         s_ghost(:,:)=zero

         allocate( fd_fac_top(0:l_max), fd_fac_bot(0:l_max) )
         bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_DEF_REAL
         fd_fac_top(:)=0.0_cp
         fd_fac_bot(:)=0.0_cp
      end if

      allocate( lSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

   end subroutine initialize_updateS
!------------------------------------------------------------------------------
   subroutine finalize_updateS
      !
      ! Memory deallocation of updateS module
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lSmat)
      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call sMat(ll)%finalize()
         end do
         deallocate(rhs1)

#ifdef WITH_PRECOND_S
         deallocate( sMat_fac )
#endif
      else
         deallocate( s_ghost, fd_fac_top, fd_fac_bot )
         call sMat_FD%finalize()
      end if

   end subroutine finalize_updateS
!------------------------------------------------------------------------------
   subroutine updateS(s, ds, dsdt, phi, tscheme)
      !
      !  Updates the entropy field s and its radial derivative.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(llm:ulm,n_r_max) ! Phase field

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max) ! Entropy
      type(type_tarray), intent(inout) :: dsdt
      !-- Output: ds
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max) ! Radial derivative of entropy

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
      call tscheme%set_imex_rhs(work_LMloc, dsdt)

      !$omp parallel default(shared)

      if ( l_phase_field ) then
         !-- Add the last remaining term to assemble St*\partial \phi/\partial t
         !$omp do private(nR,lm)
         do nR=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,nR)=work_LMloc(lm,nR)+stef*phi(lm,nR)
            end do
         end do
         !$omp end do
      end if

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

         if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_sMat(tscheme,l1,hdif_S(l1),sMat(nLMB2),sMat_fac(:,nLMB2))
#else
            call get_sMat(tscheme,l1,hdif_S(l1),sMat(nLMB2))
#endif
            lSmat(l1)=.true.
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

               rhs1(1,2*lm-1,threadid)      = real(tops(l1,m1))
               rhs1(1,2*lm,threadid)        =aimag(tops(l1,m1))
               rhs1(n_r_max,2*lm-1,threadid)= real(bots(l1,m1))
               rhs1(n_r_max,2*lm,threadid)  =aimag(bots(l1,m1))

               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=sMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =sMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif
            end do

            call sMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                 2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     0.0_cp,kind=cp)
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
            s(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do
      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dsdt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_entropy_rhs_imp(s, ds, dsdt, phi, 1, tscheme%l_imp_calc_rhs(1), &
              &                   l_in_cheb_space=.true.)
      else
         call get_entropy_rhs_imp(s, ds, dsdt, phi, tscheme%istage+1,       &
              &                   tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                   l_in_cheb_space=.true.)
      end if

   end subroutine updateS
!------------------------------------------------------------------------------
   subroutine prepareS_FD(tscheme, dsdt, phi)
      !
      ! This subroutine is used to assemble the r.h.s. of the entropy equation
      ! when parallel F.D solvers are used. Boundary values are set here.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(lm_max,nRstart:nRstop)

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm, l, m

      !-- LU factorisation of the matrix if needed
      if ( .not. lSmat(0) ) then
         call get_sMat_Rdist(tscheme,hdif_S,sMat_FD)
         lSmat(:)=.true.
      end if

      !$omp parallel default(shared) private(lm_start,lm_stop, nR, l, m, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(s_ghost, dsdt, lm_start, lm_stop, 1)

      if ( l_phase_field ) then
         !-- Add the last remaining term to assemble St*\partial \phi/\partial t
         do nR=nRstart,nRstop
            do lm=lm_start,lm_stop
               s_ghost(lm,nR)=s_ghost(lm,nR)+stef*phi(lm,nR)
            end do
         end do
      end if

      !-- Set boundary conditions
      if ( nRstart == n_r_cmb ) then
         nR=n_r_cmb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            if ( ktops == 1 ) then ! Fixed temperature
               s_ghost(lm,nR)=tops(l,m)
            else ! Fixed flux
               !TBD
               s_ghost(lm,nR)=s_ghost(lm,nR)+fd_fac_top(l)*tops(l,m)
               !s_ghost(lm,nR)=tops(l,m)
            end if
            s_ghost(lm,nR-1)=zero ! Set ghost zone to zero
         end do
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)

            if ( l_full_sphere ) then
               if ( l == 0 ) then
                  s_ghost(lm,nR)=s_ghost(lm,nR)+fd_fac_bot(l)*bots(l,m)
               else
                  ! TBD
                  s_ghost(lm,nR)=bots(l,m)
               end if
            else
               if ( kbots == 1 ) then ! Fixed temperature
                  s_ghost(lm,nR)=bots(l,m)
               else
                  ! TBD
                  s_ghost(lm,nR)=s_ghost(lm,nR)+fd_fac_bot(l)*bots(l,m)
               end if
            end if
            s_ghost(lm,nR+1)=zero ! Set ghost zone to zero
         end do
      end if
      !$omp end parallel

   end subroutine prepareS_FD
!------------------------------------------------------------------------------
   subroutine fill_ghosts_S(sg)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Nusselt
      ! numbers).
      !
      complex(cp), intent(inout) :: sg(lm_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: lm, l, m, lm_start, lm_stop
      real(cp) :: dr

      !$omp parallel default(shared) private(lm_start, lm_stop, l, m, lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Handle upper boundary
      dr = r(2)-r(1)
      if ( nRstart == n_r_cmb ) then
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            if ( ktops == 1 ) then
               sg(lm,nRstart-1)=two*sg(lm,nRstart)-sg(lm,nRstart+1)
            else
               sg(lm,nRstart-1)=sg(lm,nRstart+1)-two*dr*tops(l,m)
            end if
         end do
      end if

      !-- Handle Lower boundary
      dr = r(n_r_max)-r(n_r_max-1)
      if ( nRstop == n_r_icb ) then
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            if ( l_full_sphere ) then
               if (l == 0 ) then
                  sg(lm,nRstop+1)=sg(lm,nRstop-1)+two*dr*bots(l,m)
               else
                  sg(lm,nRstop+1)=two*sg(lm,nRstop)-sg(lm,nRstop-1)
               end if
            else ! Not a full sphere
               if (kbots == 1) then ! Fixed temperature at bottom
                  sg(lm,nRstop+1)=two*sg(lm,nRstop)-sg(lm,nRstop-1)
               else
                  sg(lm,nRstop+1)=sg(lm,nRstop-1)+two*dr*bots(l,m)
               end if
            end if
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_S
!------------------------------------------------------------------------------
   subroutine updateS_FD(s, ds, dsdt, phi, tscheme)
      !
      ! This subroutine is called after the linear solves have been completed.
      ! This is then assembling the linear terms that will be used in the r.h.s.
      ! for the next iteration.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(lm_max,nRstart:nRstop) ! Phase field

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dsdt
      complex(cp),       intent(inout) :: s(lm_max,nRstart:nRstop) ! Entropy
      !-- Output: ds
      complex(cp),       intent(out) :: ds(lm_max,nRstart:nRstop) ! Radial derivative of entropy

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dsdt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_entropy_rhs_imp_ghost(s_ghost, ds, dsdt, phi, 1, &
              &                         tscheme%l_imp_calc_rhs(1))
      else
         call get_entropy_rhs_imp_ghost(s_ghost, ds, dsdt, phi, tscheme%istage+1, &
              &                         tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !$omp parallel default(shared) private(lm_start,lm_stop,nR,lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Array copy from s_ghost to s
      !!$omp parallel do simd collapse(2) schedule(simd:static)
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            s(lm,nR)=s_ghost(lm,nR)
         end do
      end do
      !!$omp end parallel do simd
      !$omp end parallel

   end subroutine updateS_FD
!------------------------------------------------------------------------------
   subroutine finish_exp_entropy(w, dVSrLM, ds_exp_last)
      !
      ! This subroutine completes the computation of the advection term by
      ! computing the radial derivative (LM-distributed variant).
      !

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      real(cp) :: dL
      integer :: n_r, lm, start_lm, stop_lm, l
      integer, pointer :: lm2l(:)

      lm2l(1:lm_max) => lo_map%lm2l

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      if ( l_anelastic_liquid ) then
         !$omp do private(n_r,l,lm,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*     ds_exp_last(lm,n_r) - &
               &        or2(n_r)*orho1(n_r)*        work_LMloc(lm,n_r) + &
               &       or2(n_r)*orho1(n_r)*dLtemp0(n_r)*dVSrLM(lm,n_r) - &
               &        dL*or2(n_r)*orho1(n_r)*temp0(n_r)*dentropy0(n_r)*&
               &                                             w(lm,n_r)
            end do
         end do
         !$omp end do
      else
         !$omp do private(n_r,l,dL,lm)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*(      ds_exp_last(lm,n_r)- &
               &                             or2(n_r)*work_LMloc(lm,n_r)- &
               &                    dL*or2(n_r)*dentropy0(n_r)*w(lm,n_r))
            end do
         end do
         !$omp end do
      end if
      !$omp end parallel

   end subroutine finish_exp_entropy
!-----------------------------------------------------------------------------
   subroutine finish_exp_entropy_Rdist(w, dVSrLM, ds_exp_last)
      !
      ! This subroutine completes the computation of the advection term by
      ! computing the radial derivative (R-distributed variant).
      !

      !-- Input variables
      complex(cp), intent(in) :: w(lm_max,nRstart:nRstop)
      complex(cp), intent(inout) :: dVSrLM(lm_max,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(lm_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      real(cp) :: dL
      integer :: n_r, lm, l, start_lm, stop_lm

      call get_dr_Rloc(dVSrLM, work_Rloc, lm_max, nRstart, nRstop, n_r_max, &
           &           rscheme_oc)

      !$omp parallel default(shared) private(n_r, lm, l, dL, start_lm, stop_lm)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm, stop_lm)
      !$omp barrier

      if ( l_anelastic_liquid ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*     ds_exp_last(lm,n_r) - &
               &         or2(n_r)*orho1(n_r)*        work_Rloc(lm,n_r) + &
               &       or2(n_r)*orho1(n_r)*dLtemp0(n_r)*dVSrLM(lm,n_r) - &
               &        dL*or2(n_r)*orho1(n_r)*temp0(n_r)*dentropy0(n_r)*&
               &                                             w(lm,n_r)
            end do
         end do
      else
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*(      ds_exp_last(lm,n_r)- &
               &                              or2(n_r)*work_Rloc(lm,n_r)- &
               &                    dL*or2(n_r)*dentropy0(n_r)*w(lm,n_r))
            end do
         end do
      end if
      !$omp end parallel

   end subroutine finish_exp_entropy_Rdist
!-----------------------------------------------------------------------------
   subroutine get_entropy_rhs_imp(s, ds, dsdt, phi, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms that enters the r.h.s.. This is
      ! used with LM-distributed
      !

      !-- Input variables
      integer,           intent(in) :: istage
      logical,           intent(in) :: l_calc_lin
      logical, optional, intent(in) :: l_in_cheb_space
      complex(cp),       intent(in) :: phi(llm:ulm,n_r_max)

      !-- Output variable
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm, l1
      logical :: l_in_cheb
      integer, pointer :: lm2l(:),lm2m(:)
      real(cp) :: dL

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(s, ds, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do private(n_r)
         do n_r=1,n_r_max
            dsdt%old(:,n_r,istage)=s(:,n_r)
         end do
         !$omp end do
         if ( l_phase_field ) then
            !$omp do private(n_r)
            do n_r=1,n_r_max
               dsdt%old(:,n_r,istage)=dsdt%old(:,n_r,istage)-stef*phi(:,n_r)
            end do
            !$omp end do
         end if
      end if

      if ( l_calc_lin ) then

         !-- Calculate explicit time step part:
         if ( l_anelastic_liquid ) then
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=llm,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l1)* kappa(n_r) *  (    &
                  &                                          work_LMloc(lm,n_r) &
                  &     + ( beta(n_r)+two*or1(n_r)+dLkappa(n_r) ) *  ds(lm,n_r) &
                  &                                     - dL*or2(n_r)*s(lm,n_r) )
               end do
            end do
            !$omp end do
         else
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=llm,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l1)*kappa(n_r) *   (       &
                  &                                        work_LMloc(lm,n_r)      &
                  &        + ( beta(n_r)+dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )  &
                  &                                              * ds(lm,n_r)      &
                  &        - dL*or2(n_r)                         *  s(lm,n_r) )
               end do
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

   end subroutine get_entropy_rhs_imp
!-----------------------------------------------------------------------------
   subroutine get_entropy_rhs_imp_ghost(sg, ds, dsdt, phi, istage, l_calc_lin)
      !
      ! This subroutine computes the linear terms that enters the r.h.s.. This is
      ! used with R-distributed
      !

      !-- Input variables
      integer,     intent(in) :: istage
      logical,     intent(in) :: l_calc_lin
      complex(cp), intent(in) :: sg(lm_max,nRstart-1:nRstop+1)
      complex(cp), intent(in) :: phi(lm_max,nRstart:nRstop)

      !-- Output variable
      complex(cp),       intent(out) :: ds(lm_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      integer :: n_r, lm, start_lm, stop_lm, l
      real(cp) :: dL

      !$omp parallel default(shared)  private(start_lm, stop_lm, n_r, lm, l, dL)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr_ghost(sg, ds, work_Rloc, lm_max,start_lm, stop_lm,  nRstart, nRstop, &
           &             rscheme_oc)
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
      !$omp barrier

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               dsdt%old(lm,n_r,istage)=sg(lm,n_r)
            end do
            if ( l_phase_field ) then
               do lm=start_lm,stop_lm
                  dsdt%old(lm,n_r,istage)=dsdt%old(lm,n_r,istage)-stef*phi(lm,n_r)
               end do
            end if
         end do
      end if

      !-- Calculate explicit time step part:
      if ( l_calc_lin ) then
         if ( l_anelastic_liquid ) then
            do n_r=nRstart,nRstop
               do lm=start_lm,stop_lm
                  l = st_map%lm2l(lm)
                  dL = real(l*(l+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l)* kappa(n_r) *  (     &
                  &                                           work_Rloc(lm,n_r) &
                  &     + ( beta(n_r)+two*or1(n_r)+dLkappa(n_r) ) *  ds(lm,n_r) &
                  &                                    - dL*or2(n_r)*sg(lm,n_r) )
                  if ( l_full_sphere .and. l==0 .and. n_r==n_r_icb ) then
                     dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l)* kappa(n_r) *  (  &
                     &                              three *   work_Rloc(lm,n_r) &
                     &               + ( beta(n_r)+dLkappa(n_r) ) *  ds(lm,n_r) )
                  end if
               end do
            end do
         else
            do n_r=nRstart,nRstop
               do lm=start_lm,stop_lm
                  l = st_map%lm2l(lm)
                  dL = real(l*(l+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l)*kappa(n_r) *   (        &
                  &                                         work_Rloc(lm,n_r)      &
                  &        + ( beta(n_r)+dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )  &
                  &                                              * ds(lm,n_r)      &
                  &        - dL*or2(n_r)                        *  sg(lm,n_r) )
                  if ( l_full_sphere .and. l==0 .and. n_r==n_r_icb ) then
                     dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l)*kappa(n_r) *   (     &
                     &                                  three* work_Rloc(lm,n_r)   &
                     &   + ( beta(n_r)+dLtemp0(n_r)+dLkappa(n_r) )  * ds(lm,n_r) )
                  end if
               end do
            end do
         end if
      end if
      !$omp end parallel

   end subroutine get_entropy_rhs_imp_ghost
!-----------------------------------------------------------------------------
   subroutine assemble_entropy_Rloc(s, ds, dsdt, phi, tscheme)
      !
      ! This subroutine is used when an IMEX Runge-Kutta time scheme with an assembly
      ! stage is used. This is used when R is distributed.
      !

      !-- Input variable
      complex(cp),         intent(in) :: phi(lm_max,nRstart:nRstop)
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: s(lm_max,nRstart:nRstop)
      complex(cp),       intent(out) :: ds(lm_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: lm, l, m, n_r, start_lm, stop_lm
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)

      call tscheme%assemble_imex(work_Rloc, dsdt)

      !$omp parallel default(shared) private(start_lm, stop_lm, l, m)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier

      !-- In case phase field is used it needs to be substracted from work_LMloc
      !-- since time advance handles \partial/\partial t (T-St*Phi)
      if ( l_phase_field ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               work_Rloc(lm,n_r)=work_Rloc(lm,n_r)+stef*phi(lm,n_r)
            end do
         end do
      end if

      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            m = st_map%lm2m(lm)
            if ( m == 0 ) then
               s(lm,n_r)=cmplx(real(work_Rloc(lm,n_r)),0.0_cp,cp)
            else
               s(lm,n_r)=work_Rloc(lm,n_r)
            end if
         end do
      end do

      if ( ktops == 1 .and. nRstart==n_r_cmb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            s(lm,nRstart)=tops(l,m)
         end do
      end if

      if ( kbots == 1 .and. nRstop==n_r_icb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            s(lm,nRstop)=bots(l,m)
         end do
      end if

      call bulk_to_ghost(s, s_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm)
      !$omp end parallel

      call exch_ghosts(s_ghost, lm_max, nRstart, nRstop, 1)
      call fill_ghosts_S(s_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_entropy_rhs_imp_ghost(s_ghost, ds, dsdt, phi, 1, &
           &                         tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_entropy_Rloc
!-----------------------------------------------------------------------------
   subroutine assemble_entropy(s, ds, dsdt, phi, tscheme)
      !
      ! This subroutine is used to assemble the entropy/temperature at assembly
      ! stages of IMEX-RK time schemes. This is used when LM is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dsdt)

      !$omp parallel default(shared)

      !-- In case phase field is used it needs to be substracted from work_LMloc
      !-- since time advance handles \partial/\partial t (T-St*Phi)
      if ( l_phase_field ) then
         !$omp do private(n_r,lm)
         do n_r=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,n_r)=work_LMloc(lm,n_r)+stef*phi(lm,n_r)
            end do
         end do
         !$omp end do
      end if

      !$omp do private(n_r,lm,m1)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m1 = lm2m(lm)
            if ( m1 == 0 ) then
               s(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               s(lm,n_r)=work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !-- Get the boundary points using Canuto (1986) approach
      if ( l_full_sphere) then
         if ( ktops == 1 ) then ! Fixed entropy at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), one, 0.0_cp, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               end if
            end do
            !$omp end do
         else ! Fixed flux at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), one, 0.0_cp, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               end if
            end do
            !$omp end do
         end if
      else ! Spherical shell
         !-- Boundary conditions
         if ( ktops==1 .and. kbots==1 ) then ! Dirichlet on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( ktops==1 .and. kbots /= 1 ) then ! Dirichlet: top and Neumann: bot
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), one, 0.0_cp, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( kbots==1 .and. ktops /= 1 ) then ! Dirichlet: bot and Neumann: top
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), 0.0_cp, one, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( kbots /=1 .and. kbots /= 1 ) then ! Neumann on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), one, 0.0_cp, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         end if
      end if
      !$omp end parallel

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_entropy_rhs_imp(s, ds, dsdt, phi, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_entropy
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_sMat(tscheme,l,hdif,sMat,sMat_fac)
#else
   subroutine get_sMat(tscheme,l,hdif,sMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  sMat(i,j) and s0mat for the entropy equation.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: sMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktops == 1 ) then
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         if ( l == 0 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         end if
      else
         if ( kbots == 1 ) then
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

      !----- Bulk points:
      if ( l_anelastic_liquid ) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &   rscheme_oc%rMat(nR,nR_out)-tscheme%wimp_lin(1)*opr*hdif* &
               &                kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*                       &
               &                              rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &   rscheme_oc%rMat(nR,nR_out)-tscheme%wimp_lin(1)*opr*hdif* &
               &                kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               & ( beta(nR)+dLtemp0(nR)+                                    &
               &   two*or1(nR)+dLkappa(nR) )* rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      end if

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         sMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*sMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      block

      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename

      ! copy the sMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "sMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_max
         do j=1,n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") dat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=dat
      anorm = 0.0_cp
      do i=1,n_r_max
         linesum = 0.0_cp
         do j=1,n_r_max
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_max,n_r_max,temp_Mat,n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_max,temp_Mat,n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of sMat for l=",l," is ",rcond

      end block
#endif

      !-- Array copy
      call sMat%set_data(dat)

      !-- LU decomposition:
      call sMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix sMat!')

   end subroutine get_sMat
!-----------------------------------------------------------------------------
   subroutine get_sMat_Rdist(tscheme,hdif,sMat)
      !
      !  This subroutine is used to construct the matrices when the parallel
      !  solver for F.D. is employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif(0:l_max)

      !-- Output variables
      type(type_tri_par), intent(inout) :: sMat

      !-- Local variables:
      real(cp) :: dLh
      integer :: nR, l

      !-- Bulk points: we fill all the points: this is then easier to handle
      !-- Neumann boundary conditions
      !$omp parallel default(shared) private(nR,l,dLh)
      !$omp do
      do nR=1,n_r_max
         do l=0,l_max
            dLh = real(l*(l+1),cp)
            if ( l_anelastic_liquid ) then
               sMat%diag(l,nR)=  one -     tscheme%wimp_lin(1)*opr*hdif(l)*&
               &                kappa(nR)*(         rscheme_oc%ddr(nR,1) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,1) - &
               &                dLh*or2(nR)    )
               sMat%low(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*         &
               &                kappa(nR)*(         rscheme_oc%ddr(nR,0) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,0) )
               sMat%up(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*          &
               &                kappa(nR)*(         rscheme_oc%ddr(nR,2) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,2) )
            else
               sMat%diag(l,nR)=one-tscheme%wimp_lin(1)*opr*hdif(l)* &
               &                kappa(nR)*(  rscheme_oc%ddr(nR,1) + &
               & ( beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR) )*  &
               &                              rscheme_oc%dr(nR,1) - &
               &                                        dLh*or2(nR) )
               sMat%low(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*  &
               &                kappa(nR)*(  rscheme_oc%ddr(nR,0) + &
               & ( beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR) )*  &
               &                              rscheme_oc%dr(nR,0) )
               sMat%up(l,nR) =   -tscheme%wimp_lin(1)*opr*hdif(l)*  &
               &                kappa(nR)*(  rscheme_oc%ddr(nR,2) + &
               & ( beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR) )*  &
               &                              rscheme_oc%dr(nR,2) )
            end if

            if ( l==0 .and. l_full_sphere .and. nR==n_r_icb ) then
               !-- Use L'Hôpital's rule to replace 2/r d/dr at the center
               !-- by 2*d^2/dr^2
               if ( l_anelastic_liquid ) then
                  sMat%diag(l,nR)=  one - tscheme%wimp_lin(1)*opr*hdif(l)*&
                  &            kappa(nR)*(  three* rscheme_oc%ddr(nR,1) + &
                  &        ( beta(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,1))
                  sMat%low(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*     &
                  &            kappa(nR)*( three*  rscheme_oc%ddr(nR,0) + &
                  &        ( beta(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,0) )
                  sMat%up(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*      &
                  &            kappa(nR)*( three*  rscheme_oc%ddr(nR,2) + &
                  &        ( beta(nR)+dLkappa(nR) )*rscheme_oc%dr(nR,2) )
               else
                  sMat%diag(l,nR)=one-tscheme%wimp_lin(1)*opr*hdif(l)* &
                  &           kappa(nR)*( three*rscheme_oc%ddr(nR,1) + &
                  &             ( beta(nR)+dLtemp0(nR)+dLkappa(nR) )*  &
                  &                              rscheme_oc%dr(nR,1) )
                  sMat%low(l,nR)=   -tscheme%wimp_lin(1)*opr*hdif(l)*  &
                  &          kappa(nR)*(  three*rscheme_oc%ddr(nR,0) + &
                  &             ( beta(nR)+dLtemp0(nR)+dLkappa(nR) )*  &
                  &                              rscheme_oc%dr(nR,0) )
                  sMat%up(l,nR) =   -tscheme%wimp_lin(1)*opr*hdif(l)*  &
                  &           kappa(nR)*( three*rscheme_oc%ddr(nR,2) + &
                  & ( beta(nR)+dLtemp0(nR)+dLkappa(nR) )*  &
                  &                              rscheme_oc%dr(nR,2) )
               end if
            end if
         end do
      end do
      !$omp end do

      !----- Boundary conditions:
      !$omp do
      do l=0,l_max
         if ( ktops == 1 ) then
            sMat%diag(l,1)=one
            sMat%up(l,1)  =0.0_cp
            sMat%low(l,1) =0.0_cp
         else
            sMat%up(l,1)=sMat%up(l,1)+sMat%low(l,1)
            fd_fac_top(l)=two*(r(2)-r(1))*sMat%low(l,1)
            !dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
         end if

         if ( l_full_sphere ) then
            if ( l == 0 ) then
               !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
               sMat%low(l,n_r_max)=sMat%up(l,n_r_max)+sMat%low(l,n_r_max)
               fd_fac_bot(l)=two*(r(n_r_max-1)-r(n_r_max))*sMat%up(l,n_r_max)
            else
               sMat%diag(l,n_r_max)=one
               sMat%up(l,n_r_max)  =0.0_cp
               sMat%low(l,n_r_max) =0.0_cp
            end if
         else
            if ( kbots == 1 ) then
               sMat%diag(l,n_r_max)=one
               sMat%up(l,n_r_max)  =0.0_cp
               sMat%low(l,n_r_max) =0.0_cp
            else
               !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
               sMat%low(l,n_r_max)=sMat%up(l,n_r_max)+sMat%low(l,n_r_max)
               fd_fac_bot(l)=two*(r(n_r_max-1)-r(n_r_max))*sMat%up(l,n_r_max)
            end if
         end if
      end do
      !$omp end do
      !$omp end parallel

      !-- LU decomposition:
      call sMat%prepare_mat()

   end subroutine get_sMat_Rdist
!-----------------------------------------------------------------------------
end module updateS_mod
