module updateXi_mod
   !
   ! This module handles the time advance of the chemical composition xi.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max, l_axi
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: orho1, or1, or2, beta, rscheme_oc, r, dxicond, l_R, r_cmb
   use physical_parameters, only: osc, kbotxi, ktopxi
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: topxi, botxi
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_Xi
   use logic, only: l_finite_diff, l_full_sphere, l_parallel_solve, l_onset
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
#ifdef WITH_MPI
   use mpimod
#endif
   use special, only: amp_tide, omega_tide, tide_fac20, tide_fac22p, tide_fac22n
   use radial_der, only: get_robin_ghost_coeffs, get_ddr, get_dr, get_dr_Rloc, get_ddr_ghost, exch_ghosts,&
       &                 bulk_to_ghost
   use constants, only: zero, one, two, three, osq4pi
   use fields, only: work_LMloc, dxi_LMloc, dxi_Rloc
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
   class(type_realmat), pointer :: xiMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: xiMat_fac(:,:)
#endif
   logical, public, allocatable :: lXimat(:)
   type(type_tri_par), public :: xiMat_FD
   complex(cp), allocatable :: fd_fac_top(:), fd_fac_bot(:)
   complex(cp), public, allocatable :: xi_ghost(:,:)

   !-- Perturbation of the outer boundary value/flux due to the tidal
   !-- displacement of the boundary (standard lm order, see update_xibc)
   complex(cp), allocatable :: xibc_cmb(:)
   logical :: l_xibc_tide ! Tidal perturbation of a fixed composition at r_cmb (see update_xibc)

   public :: initialize_updateXi, finalize_updateXi, updateXi, assemble_comp,  &
   &         finish_exp_comp, get_comp_rhs_imp, finish_exp_comp_Rdist,         &
   &         get_comp_rhs_imp_ghost, updateXi_FD, prepareXi_FD, fill_ghosts_Xi,&
   &         assemble_comp_Rloc, &
   &         update_xibc

contains

   subroutine initialize_updateXi

      integer :: ll, n_bands
      integer, pointer :: nLMBs2(:)
      real(cp) :: a, b

      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         if ( l_finite_diff ) then
            allocate( type_bandmat :: xiMat(nLMBs2(1+rank)) )

            if ( ktopxi == 1 .and. kbotxi == 1 .and. rscheme_oc%order == 2 &
             &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed composition at both boundaries
               n_bands = rscheme_oc%order+1
            else
               n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
            end if

            do ll=1,nLMBs2(1+rank)
               call xiMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
            end do
         else
            allocate( type_densemat :: xiMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call xiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
         end if

#ifdef WITH_PRECOND_S
         allocate(xiMat_fac(n_r_max,nLMBs2(1+rank)))
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
         call xiMat_FD%initialize(1,n_r_max,0,l_max)

         !-- Allocate an array with ghost zones
         allocate( xi_ghost(lm_max, nRstart-1:nRstop+1) )
         bytes_allocated=bytes_allocated + lm_max*(nRstop-nRstart+3)*SIZEOF_DEF_COMPLEX
         xi_ghost(:,:)=zero

         allocate( fd_fac_top(0:l_max), fd_fac_bot(0:l_max) )
         bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_DEF_REAL
         fd_fac_top(:)=0.0_cp
         fd_fac_bot(:)=0.0_cp

      end if

      allocate( lXimat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL


      allocate( xibc_cmb(lm_max) )
      bytes_allocated=bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
      xibc_cmb(:)=zero

      !-- The tidal correction only applies to a fixed value: for a fixed flux it
      !-- is O(delta/H) and omitted (see update_xibc)
      call get_xibc_coeffs(n_r_cmb, 2, a, b)
      l_xibc_tide=amp_tide /= 0.0_cp .and. a == 0.0_cp

   end subroutine initialize_updateXi
!------------------------------------------------------------------------------
   subroutine finalize_updateXI
      !
      ! This subroutine deallocates the matrices involved in the time-advance of
      ! xi.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lXimat )
      if ( .not. l_parallel_solve ) then
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call xiMat(ll)%finalize()
         end do

#ifdef WITH_PRECOND_S
         deallocate(xiMat_fac)
#endif
         deallocate( rhs1 )
      else
         call xiMat_FD%finalize()
         deallocate( fd_fac_top, fd_fac_bot, xi_ghost )
      end if

      deallocate( xibc_cmb )

   end subroutine finalize_updateXI
!------------------------------------------------------------------------------
   subroutine updateXi(xi, dxi, dxidt, tscheme)
      !
      !  Updates the chemical composition field s and its radial derivative.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max) ! Chemical composition
      type(type_tarray), intent(inout) :: dxidt
      !-- Output: dxi
      complex(cp),       intent(out) :: dxi(llm:ulm,n_r_max) ! Radial derivative of xi

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
      call tscheme%set_imex_rhs(work_LMloc, dxidt)

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

         if ( .not. lXimat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_xiMat(tscheme,l1,hdif_Xi(l1),xiMat(nLMB2),xiMat_fac(:,nLMB2))
#else
            call get_xiMat(tscheme,l1,hdif_Xi(l1),xiMat(nLMB2))
#endif
             lXimat(l1)=.true.
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

               rhs1(1,2*lm-1,threadid)      = real(topxi(l1,m1)+xibc_cmb(st_map%lm2(l1,m1)))
               rhs1(1,2*lm,threadid)        =aimag(topxi(l1,m1)+xibc_cmb(st_map%lm2(l1,m1)))
               rhs1(n_r_max,2*lm-1,threadid)= real(botxi(l1,m1))
               rhs1(n_r_max,2*lm,threadid)  =aimag(botxi(l1,m1))
               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=xiMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =xiMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif

            end do

            call xiMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                  2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     xi(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     xi(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                      0.0_cp,kind=cp)
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
            xi(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dxidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_comp_rhs_imp(xi, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1),  &
              &                l_in_cheb_space=.true.)
      else
         call get_comp_rhs_imp(xi, dxi, dxidt, tscheme%istage+1,          &
              &                tscheme%l_imp_calc_rhs(tscheme%istage+1),  &
              &                l_in_cheb_space=.true.)
      end if

   end subroutine updateXi
!------------------------------------------------------------------------------
   subroutine prepareXi_FD(tscheme, dxidt)
      !
      ! This subroutine is used to assemble the r.h.s. of the composition equation
      ! when parallel F.D solvers are used. Boundary values are set here.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm, l, m
      real(cp) :: a, b

      !-- LU factorisation of the matrix if needed
      if ( .not. lXimat(0) ) then
         call get_xiMat_Rdist(tscheme,hdif_Xi,xiMat_FD)
         lXimat(:)=.true.
      end if

      !$omp parallel default(shared) private(lm_start,lm_stop, nR, l, m, lm, a, b)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Now assemble the right hand side
      call tscheme%set_imex_rhs_ghost(xi_ghost, dxidt, lm_start, lm_stop, 1)

      !-- Set boundary conditions: a fixed value replaces the equation, a flux
      !-- enters through the ghost point folded into the matrix (get_xiMat_Rdist)
      if ( nRstart == n_r_cmb ) then
         nR=n_r_cmb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_coeffs(nR, l, a, b)
            if ( a == 0.0_cp ) then
               xi_ghost(lm,nR)=(topxi(l,m)+xibc_cmb(lm))
            else
               xi_ghost(lm,nR)=xi_ghost(lm,nR)+fd_fac_top(l)*(topxi(l,m)+xibc_cmb(lm))
            end if
            xi_ghost(lm,nR-1)=zero ! Set ghost zone to zero
         end do
      end if

      if ( nRstop == n_r_icb ) then
         nR=n_r_icb
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_coeffs(nR, l, a, b)
            if ( a == 0.0_cp ) then
               xi_ghost(lm,nR)=botxi(l,m)
            else
               xi_ghost(lm,nR)=xi_ghost(lm,nR)+fd_fac_bot(l)*botxi(l,m)
            end if
            xi_ghost(lm,nR+1)=zero ! Set ghost zone to zero
         end do
      end if
      !$omp end parallel

   end subroutine prepareXi_FD
!------------------------------------------------------------------------------
   subroutine fill_ghosts_Xi(xig)
      !
      ! This subroutine is used to fill the ghosts zones that are located at
      ! nR=n_r_cmb-1 and nR=n_r_icb+1. This is used to properly set the Neuman
      ! boundary conditions. In case Dirichlet BCs are used, a simple first order
      ! extrapolation is employed. This is anyway only used for outputs (like Sherwood
      ! numbers).
      !
      complex(cp), intent(inout) :: xig(lm_max,nRstart-1:nRstop+1)

      !-- Local variables
      integer :: lm, l, m, lm_start, lm_stop
      real(cp) :: c_in, c_bc, c_rhs

      !$omp parallel default(shared) private(lm_start, lm_stop, l, m, lm, c_in, c_bc, c_rhs)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)
      !$omp barrier

      !-- Handle upper boundary
      if ( nRstart == n_r_cmb ) then
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_ghost_coeffs(n_r_cmb, l, c_in, c_bc, c_rhs)
            xig(lm,nRstart-1)=c_in*xig(lm,nRstart+1)+c_bc*xig(lm,nRstart)+ &
            &                 c_rhs*(topxi(l,m)+xibc_cmb(lm))
         end do
      end if

      !-- Handle Lower boundary
      if ( nRstop == n_r_icb ) then
         do lm=lm_start,lm_stop
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_ghost_coeffs(n_r_icb, l, c_in, c_bc, c_rhs)
            xig(lm,nRstop+1)=c_in*xig(lm,nRstop-1)+c_bc*xig(lm,nRstop)+ &
            &                c_rhs*botxi(l,m)
         end do
      end if
      !$omp end parallel

   end subroutine fill_ghosts_Xi
!------------------------------------------------------------------------------
   subroutine updateXi_FD(xi, dxidt, tscheme)
      !
      ! This subroutine is called after the linear solves have been completed.
      ! This is then assembling the linear terms that will be used in the r.h.s.
      ! for the next iteration.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      type(type_tarray), intent(inout) :: dxidt
      complex(cp),       intent(inout) :: xi(lm_max,nRstart:nRstop) ! Composition

      !-- Local variables
      integer :: nR, lm_start, lm_stop, lm

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dxidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_comp_rhs_imp_ghost(xi_ghost, dxidt, 1, tscheme%l_imp_calc_rhs(1))
      else
         call get_comp_rhs_imp_ghost(xi_ghost, dxidt, tscheme%istage+1,       &
              &                      tscheme%l_imp_calc_rhs(tscheme%istage+1))
      end if

      !$omp parallel default(shared) private(lm_start,lm_stop,nR,lm)
      lm_start=1; lm_stop=lm_max
      call get_openmp_blocks(lm_start,lm_stop)

      !-- Array copy from xi_ghost to xi
      !!$omp parallel do simd collapse(2) schedule(simd:static)
      do nR=nRstart,nRstop
         do lm=lm_start,lm_stop
            xi(lm,nR)=xi_ghost(lm,nR)
         end do
      end do
      !!$omp end parallel do simd
      !$omp end parallel

   end subroutine updateXi_FD
!------------------------------------------------------------------------------
   subroutine finish_exp_comp(w, dVXirLM, dxi_exp_last)
      !
      ! This subroutine completes the computation of the advection term which
      ! enters the composition equation by taking the radial derivative. This is
      ! the LM-distributed version.
      !

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVXirLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dxi_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      real(cp) :: dLh
      integer :: n_r, start_lm, stop_lm, l, lm

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVXirLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(lm,l)
      do n_r=1,n_r_max
         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            if ( l > l_R(n_r) ) cycle
            dxi_exp_last(lm,n_r)=orho1(n_r)*( dxi_exp_last(lm,n_r)-   &
            &                         or2(n_r)*work_LMloc(lm,n_r) )
         end do
      end do
      !$omp end do

      if ( l_onset ) then
         !$omp do private(lm,l,dLh)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lo_map%lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dLh = real(l*(l+1),cp)
               dxi_exp_last(lm,n_r)=-dLh*or2(n_r)*orho1(n_r)*w(lm,n_r)*dxicond(n_r)
            end do
         end do
         !$omp end do
      end if
      !$omp end parallel

   end subroutine finish_exp_comp
!------------------------------------------------------------------------------
   subroutine finish_exp_comp_Rdist(w, dVXirLM, dxi_exp_last)
      !
      ! This subroutine completes the computation of the advection term which
      ! enters the composition equation by taking the radial derivative. This is
      ! the R-distributed version.
      !

      !-- Input variables
      complex(cp), intent(in) :: w(lm_max,nRstart:nRstop)
      complex(cp), intent(inout) :: dVXirLM(lm_max,nRstart:nRstop)

      !-- Output variables
      complex(cp), intent(inout) :: dxi_exp_last(lm_max,nRstart:nRstop)

      !-- Local variables
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)
      real(cp) :: dLh
      integer :: n_r, lm, start_lm, stop_lm, l

      call get_dr_Rloc(dVXirLM, work_Rloc, lm_max, nRstart, nRstop, n_r_max, &
           &           rscheme_oc)

      !$omp parallel default(shared) private(n_r, lm, start_lm, stop_lm, l, dLh)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm, stop_lm)
      !$omp barrier
      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            l=st_map%lm2l(lm)
            if ( l > l_R(n_r) ) cycle
            dxi_exp_last(lm,n_r)=orho1(n_r)*( dxi_exp_last(lm,n_r)-   &
            &                           or2(n_r)*work_Rloc(lm,n_r) )
         end do
      end do

      if ( l_onset ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = st_map%lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dLh = real(l*(l+1),cp)
               dxi_exp_last(lm,n_r)=-dLh*or2(n_r)*orho1(n_r)*w(lm,n_r)*dxicond(n_r)
            end do
         end do
      end if
      !$omp end parallel

   end subroutine finish_exp_comp_Rdist
!------------------------------------------------------------------------------
   subroutine get_comp_rhs_imp(xi, dxi, dxidt, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for composition. This is the LM-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dxi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l1
      real(cp) :: dL
      integer, pointer :: lm2l(:),lm2m(:)

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
      call get_ddr(xi, dxi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(xi,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do
         do n_r=1,n_r_max
            dxidt%old(:,n_r,istage) = xi(:,n_r)
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then

         !$omp do private(n_r,lm,l1,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dxidt%impl(lm,n_r,istage)=                    osc*hdif_Xi(l1) *   &
               &     ( work_LMloc(lm,n_r)+(beta(n_r)+two*or1(n_r)) * dxi(lm,n_r) &
               &                                       - dL*or2(n_r)* xi(lm,n_r) )
            end do
         end do
         !$omp end do

      end if

      !$omp end parallel

   end subroutine get_comp_rhs_imp
!------------------------------------------------------------------------------
   subroutine get_comp_rhs_imp_ghost(xig, dxidt, istage, l_calc_lin)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for composition. This is the R-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin

      !-- Output variable
      complex(cp),       intent(inout) :: xig(lm_max,nRstart-1:nRstop+1)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
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
      !-- The radial derivative is kept in dxi_Rloc (used by update_xibc)
      call get_ddr_ghost(xig, dxi_Rloc, work_Rloc, lm_max, start_lm, stop_lm,  nRstart, &
           &             nRstop, rscheme_oc)
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single
      !$omp barrier

      if ( istage == 1 ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               dxidt%old(lm,n_r,istage) = xig(lm,n_r)
            end do
         end do
      end if

      if ( l_calc_lin ) then
         do n_r=nRstart,nRstop
            do lm=start_lm,stop_lm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dxidt%impl(lm,n_r,istage)=                     osc*hdif_Xi(l) *   &
               &     ( work_Rloc(lm,n_r)+(beta(n_r)+two*or1(n_r)) *  dxi_Rloc(lm,n_r) &
               &                                       - dL*or2(n_r)* xig(lm,n_r) )

               if ( l==0 .and. l_full_sphere .and. n_r==n_r_icb ) then
                  dxidt%impl(lm,n_r,istage)= osc*hdif_Xi(l) *   &
                  &               ( three*work_Rloc(lm,n_r)+beta(n_r)*dxi_Rloc(lm,n_r) )
               end if
            end do
         end do
      end if
      !$omp end parallel

   end subroutine get_comp_rhs_imp_ghost
!------------------------------------------------------------------------------
   subroutine assemble_comp(xi, dxi, dxidt, tscheme)
      !
      ! This subroutine is used to assemble the chemical composition when an
      ! IMEX-RK with an assembly stage is employed. Non-Dirichlet boundary
      ! conditions are handled using Canuto (1986) approach. This is the LM
      ! distributed version.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dxi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      real(cp) :: a_top, b_top, a_bot, b_bot
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dxidt)

      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m1)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m1 = lm2m(lm)
            if ( m1 == 0 ) then
               xi(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               xi(lm,n_r)=work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !-- Boundary conditions
      !$omp do private(lm,l1,m1,a_top,b_top,a_bot,b_bot)
      do lm=llm,ulm
         l1 = lm2l(lm)
         m1 = lm2m(lm)
         call get_xibc_coeffs(n_r_cmb, l1, a_top, b_top)
         call get_xibc_coeffs(n_r_icb, l1, a_bot, b_bot)
         call rscheme_oc%robin_bc(a_top, b_top, topxi(l1,m1)+xibc_cmb(st_map%lm2(l1,m1)), &
              &                   a_bot, b_bot, &
              &                   botxi(l1,m1), xi(lm,:))
      end do
      !$omp end do
      !$omp end parallel

      call get_comp_rhs_imp(xi, dxi, dxidt, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_comp
!------------------------------------------------------------------------------
   subroutine assemble_comp_Rloc(xi, dxidt, tscheme)
      !
      ! This subroutine is used when an IMEX Runge-Kutta time scheme with an assembly
      ! stage is used. This is used when R is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: xi(lm_max,nRstart:nRstop)
      type(type_tarray), intent(inout) :: dxidt

      !-- Local variables
      integer :: lm, l, m, n_r, start_lm, stop_lm
      real(cp) :: a, b
      complex(cp) :: work_Rloc(lm_max,nRstart:nRstop)

      call tscheme%assemble_imex(work_Rloc, dxidt)

      !$omp parallel default(shared) private(start_lm, stop_lm, l, m, a, b)
      start_lm=1; stop_lm=lm_max
      call get_openmp_blocks(start_lm,stop_lm)
      !$omp barrier

      do n_r=nRstart,nRstop
         do lm=start_lm,stop_lm
            m = st_map%lm2m(lm)
            if ( m == 0 ) then
               xi(lm,n_r)=cmplx(real(work_Rloc(lm,n_r)),0.0_cp,cp)
            else
               xi(lm,n_r)=work_Rloc(lm,n_r)
            end if
         end do
      end do

      !-- Fixed values at the boundaries (fluxes are handled by the ghost points)
      if ( nRstart==n_r_cmb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_coeffs(n_r_cmb, l, a, b)
            if ( a == 0.0_cp ) xi(lm,nRstart)=(topxi(l,m)+xibc_cmb(lm))
         end do
      end if

      if ( nRstop==n_r_icb ) then
         do lm=start_lm,stop_lm
            l = st_map%lm2l(lm)
            m = st_map%lm2m(lm)
            call get_xibc_coeffs(n_r_icb, l, a, b)
            if ( a == 0.0_cp ) xi(lm,nRstop)=botxi(l,m)
         end do
      end if

      call bulk_to_ghost(xi, xi_ghost, 1, nRstart, nRstop, lm_max, start_lm, stop_lm)
      !$omp end parallel

      call exch_ghosts(xi_ghost, lm_max, nRstart, nRstop, 1)
      call fill_ghosts_Xi(xi_ghost)

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_comp_rhs_imp_ghost(xi_ghost, dxidt, 1, tscheme%l_imp_calc_rhs(1))

   end subroutine assemble_comp_Rloc
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_xiMat(tscheme,l,hdif,xiMat,xiMat_fac)
#else
   subroutine get_xiMat(tscheme,l,hdif,xiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  xiMat(i,j) for the equation for the chemical composition.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: xiMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: xiMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: dLh, a_top, b_top, a_bot, b_bot
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      call get_xibc_coeffs(n_r_cmb, l, a_top, b_top)
      call get_xibc_coeffs(n_r_icb, l, a_bot, b_bot)
      dat(1,:)      =rscheme_oc%rnorm*(a_top*rscheme_oc%drMat(1,:)+ &
      &                                b_top*rscheme_oc%rMat(1,:))
      dat(n_r_max,:)=rscheme_oc%rnorm*(a_bot*rscheme_oc%drMat(n_r_max,:)+ &
      &                                b_bot*rscheme_oc%rMat(n_r_max,:))

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                           &
            &                                 rscheme_oc%rMat(nR,nR_out) - &
            & tscheme%wimp_lin(1)*osc*hdif*(rscheme_oc%d2rMat(nR,nR_out) + &
            &   ( beta(nR)+two*or1(nR) )*    rscheme_oc%drMat(nR,nR_out) - &
            &        dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
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
         xiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*xiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call xiMat%set_data(dat)

      !----- LU decomposition:
      call xiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix xiMat!')

   end subroutine get_xiMat
!-----------------------------------------------------------------------------
   subroutine get_xiMat_Rdist(tscheme,hdif,xiMat)
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  xiMat(i,j) for the equation for the chemical composition. This is
      !  used when parallel F.D. solvers are employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif(0:l_max)

      !-- Output variables
      type(type_tri_par), intent(inout) :: xiMat

      !-- Local variables:
      integer :: nR, l
      real(cp) :: dLh, a, b, c_in, c_bc, c_rhs

      !----- Bulk points
      !$omp parallel default(shared) private(nR,l,dLh,a,b,c_in,c_bc,c_rhs)
      !$omp do
      do nR=1,n_r_max
         do l=0,l_max
            dLh=real(l*(l+1),kind=cp)
            xiMat%diag(l,nR)=one-tscheme%wimp_lin(1)*osc*hdif(l)*(         &
            &                                       rscheme_oc%ddr(nR,1) + &
            &           ( beta(nR)+two*or1(nR) )*    rscheme_oc%dr(nR,1) - &
            &                                        dLh*or2(nR) )
            xiMat%low(l,nR)=-tscheme%wimp_lin(1)*osc*hdif(l)*(             &
            &                                       rscheme_oc%ddr(nR,0) + &
            &           ( beta(nR)+two*or1(nR) )*    rscheme_oc%dr(nR,0) )
            xiMat%up(l,nR) =-tscheme%wimp_lin(1)*osc*hdif(l)*(             &
            &                                       rscheme_oc%ddr(nR,2) + &
            &           ( beta(nR)+two*or1(nR) )*    rscheme_oc%dr(nR,2) )

            if ( l==0 .and. l_full_sphere .and. nR==n_r_icb ) then
               !-- Use L'Hôpital's rule to replace 2/r d/dr at the center
               !-- by 2*d^2/dr^2
               xiMat%diag(l,nR)=one-tscheme%wimp_lin(1)*osc*hdif(l)*(         &
               &                               three*  rscheme_oc%ddr(nR,1) + &
               &                            beta(nR)*   rscheme_oc%dr(nR,1) )
               xiMat%low(l,nR)=-tscheme%wimp_lin(1)*osc*hdif(l)*(             &
               &                               three*  rscheme_oc%ddr(nR,0) + &
               &                            beta(nR)*   rscheme_oc%dr(nR,0) )
               xiMat%up(l,nR) =-tscheme%wimp_lin(1)*osc*hdif(l)*(             &
               &                               three*  rscheme_oc%ddr(nR,2) + &
               &                            beta(nR)*   rscheme_oc%dr(nR,2) )
            end if
         end do
      end do
      !$omp end do

      !----- Boundary conditions: a fixed value replaces the equation, a flux
      !----- condition gives the ghost point that is folded into the matrix
      !$omp do
      do l=0,l_max
         call get_xibc_coeffs(n_r_cmb, l, a, b)
         if ( a == 0.0_cp ) then
            xiMat%diag(l,1)=one
            xiMat%up(l,1)  =0.0_cp
            xiMat%low(l,1) =0.0_cp
         else
            call get_xibc_ghost_coeffs(n_r_cmb, l, c_in, c_bc, c_rhs)
            xiMat%diag(l,1)=xiMat%diag(l,1)+c_bc*xiMat%low(l,1)
            xiMat%up(l,1)  =xiMat%up(l,1)+c_in*xiMat%low(l,1)
            fd_fac_top(l)  =-c_rhs*xiMat%low(l,1)
         end if

         call get_xibc_coeffs(n_r_icb, l, a, b)
         if ( a == 0.0_cp ) then
            xiMat%diag(l,n_r_max)=one
            xiMat%up(l,n_r_max)  =0.0_cp
            xiMat%low(l,n_r_max) =0.0_cp
         else
            call get_xibc_ghost_coeffs(n_r_icb, l, c_in, c_bc, c_rhs)
            xiMat%diag(l,n_r_max)=xiMat%diag(l,n_r_max)+c_bc*xiMat%up(l,n_r_max)
            xiMat%low(l,n_r_max) =xiMat%low(l,n_r_max)+c_in*xiMat%up(l,n_r_max)
            fd_fac_bot(l)        =-c_rhs*xiMat%up(l,n_r_max)
         end if
      end do
      !$omp end do
      !$omp end parallel

      !----- LU decomposition:
      call xiMat%prepare_mat()

   end subroutine get_xiMat_Rdist
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   subroutine get_xibc_coeffs(nR, l, a, b)
      !
      ! Boundary condition on the composition at the boundary nR:
      ! a*dxi/dr+b*xi=rhs. ktopxi/kbotxi=1 means a fixed value, any other
      ! value a fixed flux. At the centre of a full sphere, regularity requires
      ! dxi/dr=0 for l=0 and xi=0 otherwise.
      ! A fixed value always has (a,b)=(0,1): the finite-difference paths
      ! replace the boundary row by an identity row and set the value directly.
      !

      !-- Input variables:
      integer, intent(in) :: nR ! n_r_cmb or n_r_icb
      integer, intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: a, b

      !-- Local variable:
      integer :: k

      if ( nR == n_r_icb .and. l_full_sphere ) then
         if ( l == 0 ) then
            k=2
         else
            k=1
         end if
      else if ( nR == n_r_cmb ) then
         k=ktopxi
      else
         k=kbotxi
      end if

      if ( k == 1 ) then ! Fixed value
         a=0.0_cp
         b=one
      else ! Fixed flux
         a=one
         b=0.0_cp
      end if

   end subroutine get_xibc_coeffs
!-----------------------------------------------------------------------------
   subroutine get_xibc_ghost_coeffs(nR, l, c_in, c_bc, c_rhs)
      !
      ! Ghost point beyond the boundary nR for the composition (see
      ! get_robin_ghost_coeffs). For a flux condition this enforces the
      ! condition, for a fixed value the ghost point is only used for outputs.
      !

      !-- Input variables:
      integer, intent(in) :: nR ! n_r_cmb or n_r_icb
      integer, intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: c_in, c_bc, c_rhs

      !-- Local variables:
      real(cp) :: a, b, h

      call get_xibc_coeffs(nR, l, a, b)
      if ( nR == n_r_cmb ) then
         h=r(1)-r(2)
      else
         h=r(n_r_max)-r(n_r_max-1)
      end if
      call get_robin_ghost_coeffs(a, b, h, c_in, c_bc, c_rhs)

   end subroutine get_xibc_ghost_coeffs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   subroutine update_xibc(time)
      !
      ! Perturbation of the composition b.c. at the outer boundary due to its tidal
      ! displacement delta (Taylor expansion around r_cmb, see the notes on
      ! tidal forcing). For a fixed value, the condition f=rhs holds on the
      ! displaced boundary r_cmb+delta, which to first order adds
      !
      ! .. code-block:: fortran
      !
      !    -delta*df0/dr
      !
      ! to its r.h.s. at r_cmb. f0 is the current spherically-symmetric
      ! composition (recorded when the implicit terms were last computed), so the
      ! correction follows the boundary layer (df0/dr ~ Nu times the conducting
      ! value). For a fixed flux, the corresponding correction -delta*d2f0/dr2
      ! is omitted: at a no-slip, impermeable wall the l=0 balance gives
      ! d2f0/dr2 = -2/r*df0/dr + (df0/dt - sources)/kappa, so it is only
      ! O(delta/H) relative to the flux. It is called once per stage, before
      ! the boundary values are used.
      !

      !-- Input variable:
      real(cp), intent(in) :: time

      !-- Local variables:
      complex(cp) :: delta20, delta22
      real(cp) :: fac, f0_dr
      integer :: lm20, lm22
#ifdef WITH_MPI
      integer :: ierr
#endif

      if ( .not. l_xibc_tide ) return

      !-- d(f00)/dr at r_cmb, from the derivatives computed with the implicit
      !-- terms (only known on the rank that holds l=m=0 at r_cmb)
      f0_dr=0.0_cp
      if ( l_parallel_solve ) then
         if ( nRstart == n_r_cmb ) f0_dr=real(dxi_Rloc(st_map%lm2(0,0),n_r_cmb))
      else
         if ( llm <= lo_map%lm2(0,0) .and. ulm >= lo_map%lm2(0,0) ) &
         &    f0_dr=real(dxi_LMloc(lo_map%lm2(0,0),n_r_cmb))
      end if
#ifdef WITH_MPI
      call MPI_Allreduce(MPI_IN_PLACE, f0_dr, 1, MPI_DEF_REAL, MPI_SUM, &
           &             MPI_COMM_WORLD, ierr)
#endif

      !-- Radial displacement of the outer boundary by the tide: time integral
      !-- of u_r=l(l+1)/r^2*W with W from the tidal flow (see update_wbc)
      lm20=st_map%lm2(2,0)
      lm22=st_map%lm2(2,2)
      fac=6.0_cp/r_cmb**2/omega_tide
      delta20=fac*tide_fac20*cmplx(sin(omega_tide*time), 0.0_cp, cp)
      if ( .not. l_axi ) delta22=fac*cmplx( (tide_fac22p+tide_fac22n)*sin(omega_tide*time), &
      &                                    (-tide_fac22p+tide_fac22n)*cos(omega_tide*time), cp)

      fac=-osq4pi*f0_dr ! Physical f0=osq4pi*f00
      xibc_cmb(lm20)=fac*delta20 ! Other entries stay zero
      if ( .not. l_axi ) xibc_cmb(lm22)=fac*delta22 ! m=0 tide only if axisymmetric

   end subroutine update_xibc
!-----------------------------------------------------------------------------
end module updateXi_mod
