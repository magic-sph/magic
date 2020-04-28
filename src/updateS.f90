#include "perflib_preproc.cpp"
module updateS_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max, n_r_cmb, n_r_icb, &
       &                 get_openmp_blocks, n_lo_loc, n_mlo_loc
   use radial_functions, only: orho1, or1, or2, beta, dentropy0, rscheme_oc,  &
       &                       kappa, dLkappa, dLtemp0, temp0, r
   use physical_parameters, only: opr, kbots, ktops
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: tops, bots
   use blocking, only: lo_map, lo_sub_map, llm, ulm
   use horizontal_data, only: hdif_S
   use logic, only: l_update_s, l_anelastic_liquid, l_finite_diff, &
       &            l_full_sphere
   use parallel_mod, only: coord_r, chunksize, n_ranks_r
   use radial_der, only: get_ddr, get_dr
   use fields, only:  work_LMloc, work_LMdist
   use constants, only: zero, one, two
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use real_matrices
   use band_matrices
   use LMmapping, only: map_mlo

   implicit none

   private

   !-- Local variables
   real(cp), allocatable :: rhs1(:,:,:)
   class(type_realmat), pointer :: sMat(:), s0Mat
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: s0Mat_fac(:)
#endif
   logical, public, allocatable :: lSmat(:)
   
   
   real(cp), allocatable :: rhs1_dist(:,:,:)
   class(type_realmat), pointer :: sMat_dist(:), s0Mat_dist
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac_dist(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: s0Mat_fac_dist(:)
#endif
   logical, public, allocatable :: lSmat_dist(:)

   integer :: maxThreads

   public :: initialize_updateS, updateS, finalize_updateS, assemble_entropy, &
   &         finish_exp_entropy, get_entropy_rhs_imp, updateS_dist,           &
   &         initialize_updateS_dist, finish_exp_entropy_dist

contains

   subroutine initialize_updateS

      integer, pointer :: nLMBs2(:)
      integer :: ll,n_bands

      nLMBs2(1:n_ranks_r) => lo_sub_map%nLMBs2

      if ( l_finite_diff ) then
         allocate( type_bandmat :: sMat(nLMBs2(1+coord_r)) )
         allocate( type_bandmat :: s0Mat )

         if ( ktops == 1 .and. kbots == 1 .and. rscheme_oc%order == 2 &
          &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed entropy at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         !print*, 'S', n_bands
         call s0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+coord_r)
            call sMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: sMat(nLMBs2(1+coord_r)) )
         allocate( type_densemat :: s0Mat )

         call s0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+coord_r)
            call sMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_S
      allocate(sMat_fac(n_r_max,nLMBs2(1+coord_r)))
      bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+coord_r)*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(s0Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      allocate( lSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateS
!------------------------------------------------------------------------------
   subroutine initialize_updateS_dist

      integer :: ll, n_bands

      if ( l_finite_diff ) then
         allocate( type_bandmat :: sMat_dist(n_lo_loc) )
         allocate( type_bandmat :: s0Mat_dist )

         if ( ktops == 1 .and. kbots == 1 .and. rscheme_oc%order == 2 &
          &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed entropy at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         !print*, 'S', n_bands
         call s0Mat_dist%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,n_lo_loc
            call sMat_dist(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: sMat_dist(n_lo_loc) )
         allocate( type_densemat :: s0Mat_dist )

         call s0Mat_dist%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,n_lo_loc
            call sMat_dist(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_S
      allocate(sMat_fac_dist(n_r_max,n_lo_loc))
      bytes_allocated = bytes_allocated+n_r_max*n_lo_loc*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(s0Mat_fac_dist(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      allocate( lSmat_dist(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate( rhs1_dist(n_r_max,2*maxval(map_mlo%n_mi(:)), 1))
      bytes_allocated = bytes_allocated + n_r_max*maxval(map_mlo%n_mi(:))*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateS_dist
!------------------------------------------------------------------------------
   subroutine finalize_updateS

      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_ranks_r) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+coord_r)
         call sMat(ll)%finalize()
      end do
      call s0Mat%finalize()

      deallocate( lSmat )
#ifdef WITH_PRECOND_S
      deallocate( sMat_fac )
#endif
#ifdef WITH_PRECOND_S0
      deallocate( s0Mat_fac )
#endif
      deallocate( rhs1 )

   end subroutine finalize_updateS
!------------------------------------------------------------------------------
   subroutine updateS_dist(s, ds, dsdt, tscheme)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: s(n_mlo_loc,n_r_max)
      type(type_tarray), intent(inout) :: dsdt
      !-- Output: ds
      complex(cp),       intent(out) :: ds(n_mlo_loc,n_r_max)

      !-- Local variables:
      integer :: l, m               ! degree and order corresponding
      integer :: lj, mi, i          ! l, m and ml counter
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out            ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0
      
      if ( .not. l_update_s ) return

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMdist, dsdt, 1, n_mlo_loc, n_r_max)

      call solve_counter%start_count()
      
      ! Loop over local l
      do lj=1, n_lo_loc
         l = map_mlo%lj2l(lj)

         ! Builds matrices if needed
         if ( .not. lSmat_dist(lj) ) then
            if ( l == 0 ) then
#ifdef WITH_PRECOND_S0
               call get_s0Mat(tscheme,s0Mat_dist,s0Mat_fac_dist)
#else
               call get_s0Mat(tscheme,s0Mat_dist)
#endif
            else
#ifdef WITH_PRECOND_S
               call get_sMat(tscheme,l,hdif_S(l), &
                    &        sMat_dist(lj),sMat_fac_dist(:,lj))
#else
               call get_sMat(tscheme,l,hdif_S(l),sMat_dist(lj))
#endif
            end if
            lSmat_dist(lj)=.true.
         end if
         
         ! Loop over m corresponding to current l
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)
            
            if ( l==0 ) then
            
               rhs(1)      =real(tops(0,0))
               do nR=2,n_r_max-1
                  rhs(nR)=real(work_LMdist(i,nR))
               end do
               rhs(n_r_max)=real(bots(0,0))
#ifdef WITH_PRECOND_S0
               rhs(:) = s0Mat_fac_dist(:)*rhs(:)
#endif
               call s0Mat_dist%solve(rhs)

            else ! l  /=  0

               ! Builds real RHS
               rhs1_dist(1,2*mi-1,1)      = real(tops(l,m))
               do nR=2,n_r_max-1
                  rhs1_dist(nR,2*mi-1,1)= real(work_LMdist(i,nR))
               end do
               rhs1_dist(n_r_max,2*mi-1,1)= real(bots(l,m))
               
               ! Builds imag RHS
               rhs1_dist(1,2*mi,1)        =aimag(tops(l,m))
               do nR=2,n_r_max-1
                  rhs1_dist(nR,2*mi,1)  =aimag(work_LMdist(i,nR))
               end do
               rhs1_dist(n_r_max,2*mi,1)  =aimag(bots(l,m))

#ifdef WITH_PRECOND_S
               rhs1_dist(:,2*mi-1,1)=sMat_fac_dist(:,lj)*rhs1_dist(:,2*mi-1,1)
               rhs1_dist(:,2*mi,1)  =sMat_fac_dist(:,lj)*rhs1_dist(:,2*mi,1)
#endif
            end if
         end do ! loop over m
         
         ! Solve for all RHS at once
         if ( l>0 ) then
            call sMat_dist(lj)%solve(rhs1_dist(:,1:2*map_mlo%n_mi(lj),1), &
                  &                 2*map_mlo%n_mi(lj))
         end if
         
      
         ! Loop over m corresponding to current l (again)
         do mi=1,map_mlo%n_mi(lj)
            m = map_mlo%milj2m(mi,lj)
            i = map_mlo%milj2i(mi,lj)
            
            if ( l == 0 ) then
               do n_r_out=1,rscheme_oc%n_max
                  s(i,n_r_out)=rhs(n_r_out)
               end do
            else
               if ( m > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(i,n_r_out)= cmplx(rhs1_dist(n_r_out,2*mi-1,1), &
                     &                     rhs1_dist(n_r_out,2*mi,1),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     s(i,n_r_out)= cmplx(rhs1_dist(n_r_out,2*mi-1,1), &
                     &                     0.0_cp,kind=cp)
                  end do
               end if
            end if
         end do  ! loop over m (again)
      end do     ! loop over l
      
      call solve_counter%stop_count(l_increment=.false.)

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do i=1,n_mlo_loc
            s(i,n_r_out)=zero
         end do
      end do
 
      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dsdt, 1, n_mlo_loc, n_r_max)
 
      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_entropy_rhs_imp_dist(s, ds, dsdt, 1, tscheme%l_imp_calc_rhs(1), &
              &                   l_in_cheb_space=.true.)
      else
         call get_entropy_rhs_imp_dist(s, ds, dsdt, tscheme%istage+1,            &
              &                   tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                   l_in_cheb_space=.true.)
      end if

   end subroutine updateS_dist
!------------------------------------------------------------------------------
   subroutine updateS(s, ds, dsdt, tscheme)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt
      !-- Output: ds
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid,iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_s ) return

      nLMBs2(1:n_ranks_r) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+coord_r

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dsdt, llm, ulm, n_r_max)

      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! one subblock is linked to one l value and needs therefore once the matrix
      !$omp single
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,threadid) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 == 0 ) then
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_s0Mat(tscheme,s0Mat,s0Mat_fac)
#else
               call get_s0Mat(tscheme,s0Mat)
#endif
               lSmat(l1)=.true.
            end if
         else
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_sMat(tscheme,l1,hdif_S(l1), &
                    &        sMat(nLMB2),sMat_fac(:,nLMB2))
#else
               call get_sMat(tscheme,l1,hdif_S(l1),sMat(nLMB2))
#endif
               lSmat(l1)=.true.
            end if
          end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  rhs(1)      =real(tops(0,0))
                  rhs(n_r_max)=real(bots(0,0))
                  do nR=2,n_r_max-1
                     rhs(nR)=real(work_LMloc(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs(:) = s0Mat_fac(:)*rhs(:)
#endif

                  call s0Mat%solve(rhs)

               else ! l1  /=  0

                  lmB=lmB+1

                  rhs1(1,2*lmB-1,threadid)      = real(tops(l1,m1))
                  rhs1(1,2*lmB,threadid)        =aimag(tops(l1,m1))
                  rhs1(n_r_max,2*lmB-1,threadid)= real(bots(l1,m1))
                  rhs1(n_r_max,2*lmB,threadid)  =aimag(bots(l1,m1))

                  do nR=2,n_r_max-1
                     rhs1(nR,2*lmB-1,threadid)= real(work_LMloc(lm1,nR))
                     rhs1(nR,2*lmB,threadid)  =aimag(work_LMloc(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S
                  rhs1(:,2*lmB-1,threadid)=sMat_fac(:,nLMB2)*rhs1(:,2*lmB-1,threadid)
                  rhs1(:,2*lmB,threadid)  =sMat_fac(:,nLMB2)*rhs1(:,2*lmB,threadid)
#endif
               end if
            end do
            !PERFOFF

            !PERFON('upS_sol')
            if ( lmB  >  lmB0 ) then
               call sMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*lmB,threadid), &
                    &                 2*(lmB-lmB0))
            end if
            !PERFOFF

            lmB=lmB0
            !PERFON('upS_af')
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                     rhs1(n_r_out,2*lmB,threadid),kind=cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                        &                     0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do     ! loop over lm blocks
      !$OMP END SINGLE
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
      call tscheme%rotate_imex(dsdt, llm, ulm, n_r_max)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_entropy_rhs_imp(s, ds, dsdt, 1, tscheme%l_imp_calc_rhs(1), &
              &                   l_in_cheb_space=.true.)
      else
         call get_entropy_rhs_imp(s, ds, dsdt, tscheme%istage+1,            &
              &                   tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                   l_in_cheb_space=.true.)
      end if

   end subroutine updateS
!------------------------------------------------------------------------------
   subroutine finish_exp_entropy(w, dVSrLM, ds_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      real(cp) :: dL
      integer :: n_r, lm, start_lm, stop_lm, l1
      integer, pointer :: lm2l(:),lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      if ( l_anelastic_liquid ) then
         !$omp do private(n_r,l1,lm,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*     ds_exp_last(lm,n_r) - &
               &        or2(n_r)*orho1(n_r)*        work_LMloc(lm,n_r) + &
               &       or2(n_r)*orho1(n_r)*dLtemp0(n_r)*dVSrLM(lm,n_r) - &
               &        dL*or2(n_r)*orho1(n_r)*temp0(n_r)*dentropy0(n_r)*&
               &                                             w(lm,n_r)
            end do
         end do
         !$omp end do
      else
         !$omp do private(n_r,l1,dL,lm)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
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
   subroutine finish_exp_entropy_dist(w, dVSrLM, ds_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: w(n_mlo_loc,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(n_mlo_loc,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(n_mlo_loc,n_r_max)

      !-- Local variables
      real(cp) :: dL
      integer :: n_r, lm, start_lm, stop_lm, l

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM, work_LMdist, n_mlo_loc, start_lm,  &
           &       stop_lm, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      if ( l_anelastic_liquid ) then
         !$omp do private(n_r,l,lm,dL)
         do n_r=1,n_r_max
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*     ds_exp_last(lm,n_r) - &
               &        or2(n_r)*orho1(n_r)*       work_LMdist(lm,n_r) + &
               &       or2(n_r)*orho1(n_r)*dLtemp0(n_r)*dVSrLM(lm,n_r) - &
               &        dL*or2(n_r)*orho1(n_r)*temp0(n_r)*dentropy0(n_r)*&
               &                                             w(lm,n_r)
            end do
         end do
         !$omp end do
      else
         !$omp do private(n_r,l,dL,lm)
         do n_r=1,n_r_max
            do lm=1,n_mlo_loc
               l = map_mlo%i2l(lm)
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*(      ds_exp_last(lm,n_r)- &
               &                            or2(n_r)*work_LMdist(lm,n_r)- &
               &                    dL*or2(n_r)*dentropy0(n_r)*w(lm,n_r))
            end do
         end do
         !$omp end do
      end if
      !$omp end parallel

   end subroutine finish_exp_entropy_dist
!-----------------------------------------------------------------------------
   subroutine get_entropy_rhs_imp_dist(s, ds, dsdt, istage, l_calc_lin, l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: s(n_mlo_loc,n_r_max)
      complex(cp),       intent(out) :: ds(n_mlo_loc,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm, l1
      logical :: l_in_cheb
      real(cp) :: dL

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=1; stop_lm=n_mlo_loc
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(s, ds, work_LMdist, n_mlo_loc, start_lm, stop_lm, n_r_max, &
           &       rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(s, n_mlo_loc, start_lm, stop_lm)
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
      end if

      if ( l_calc_lin ) then

         !-- Calculate explicit time step part:
         if ( l_anelastic_liquid ) then
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=1,n_mlo_loc
                  l1 = map_mlo%i2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=opr*hdif_S(l1)* kappa(n_r) *  ( &
                  &                                         work_LMdist(lm,n_r) &
                  &     + ( beta(n_r)+two*or1(n_r)+dLkappa(n_r) ) *  ds(lm,n_r) &
                  &                                     - dL*or2(n_r)*s(lm,n_r) )
               end do
            end do
            !$omp end do
         else
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=1,n_mlo_loc
                  l1 = map_mlo%i2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=opr*hdif_S(l1)*kappa(n_r) *   (    &
                  &                                       work_LMdist(lm,n_r)      &
                  &        + ( beta(n_r)+dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )  &
                  &                                              * ds(lm,n_r)      &
                  &        - dL*or2(n_r)                         *  s(lm,n_r) )
               end do
            end do
            !$omp end do
         end if

      end if

      !$omp end parallel

   end subroutine get_entropy_rhs_imp_dist
!-----------------------------------------------------------------------------
   subroutine get_entropy_rhs_imp(s, ds, dsdt, istage, l_calc_lin, l_in_cheb_space)

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

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
   subroutine assemble_entropy(s, ds, dsdt, tscheme)
      !
      ! This subroutine is used to assemble the entropy/temperature at assembly
      ! stages of IMEX-RK time schemes
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dsdt, llm, ulm, n_r_max)

      !$omp parallel default(shared)
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
               if ( l1 == 1 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), one, 0.0_cp, &
                       &                   bots(l1,m1), s(lm,:))
               end if
            end do
            !$omp end do
         else ! Fixed flux at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 1 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), one, 0.0_cp, &
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
               call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         end if
      end if
      !$omp end parallel

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_entropy_rhs_imp(s, ds, dsdt, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_entropy
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_s0Mat(tscheme,sMat,sMat_fac)
#else
   subroutine get_s0Mat(tscheme,sMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  sMat0
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step

      !-- Output variables
      class(type_realmat), intent(inout) :: sMat
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max)
      integer :: info,nR_out,nR

      !----- Boundary conditions:
      if ( ktops == 1 ) then
         !--------- Constant entropy at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         !--------- Constant flux at CMB:
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      else
         if ( kbots == 1 ) then
            !--------- Constant entropy at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            !--------- Constant flux at ICB:
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !-- Fill bulk points:
      if ( l_anelastic_liquid ) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                      &
               &      rscheme_oc%rMat(nR,nR_out) - tscheme%wimp_lin(1)*  &
               &       opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
               & (beta(nR)+two*or1(nR)+dLkappa(nR))*                     &
               &                           rscheme_oc%drMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                     &
               &      rscheme_oc%rMat(nR,nR_out) - tscheme%wimp_lin(1)* &
               &      opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
               & (beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR))*        &
               &                          rscheme_oc%drMat(nR,nR_out) ) )
            end do
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S0
      ! compute the linesum of each line
      do nR=1,n_r_max
         sMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*sMat_fac(nR)
      end do
#endif

      !-- Array copy
      call sMat%set_data(dat)

      !---- LU decomposition:
      call sMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix sMat0!')

   end subroutine get_s0Mat
!-----------------------------------------------------------------------------
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
         if ( l == 1 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
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

   end subroutine get_Smat
!-----------------------------------------------------------------------------
end module updateS_mod
