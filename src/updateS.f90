#include "perflib_preproc.cpp"
module updateS_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: orho1, or1, or2, beta, dentropy0, rscheme_oc,  &
       &                       kappa, dLkappa, dLtemp0, temp0, r
   use physical_parameters, only: opr, kbots, ktops
   use num_param, only: alpha, dct_counter, solve_counter
   use init_fields, only: tops,bots
   use blocking, only: st_map, lo_map, lo_sub_map, llm, ulm
   use horizontal_data, only: dLh,hdif_S
   use logic, only: l_update_s, l_anelastic_liquid, l_finite_diff, &
       &            l_full_sphere
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
   use radial_der, only: get_ddr, get_dr
   use fields, only:  work_LMloc
   use constants, only: zero, one, two
   use useful, only: abortRun
   use dense_matrices
   use real_matrices
   use band_matrices

   implicit none

   private

   !-- Local variables
   complex(cp), allocatable :: rhs1(:,:,:)
   class(type_realmat), pointer :: sMat(:), s0Mat
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: s0Mat_fac(:)
#endif
   logical, public, allocatable :: lSmat(:)

   integer :: maxThreads

   public :: initialize_updateS, updateS, updateS_ala, finalize_updateS

contains

   subroutine initialize_updateS

      integer, pointer :: nLMBs2(:)
      integer :: ll,n_bands

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      if ( l_finite_diff ) then
         allocate( type_bandmat :: sMat(nLMBs2(1+rank)) )
         allocate( type_bandmat :: s0Mat )

         if ( ktops == 1 .and. kbots == 1 .and. rscheme_oc%order == 2 &
          &   .and. rscheme_oc%order_boundary <= 2 ) then ! Fixed entropy at both boundaries
            n_bands = rscheme_oc%order+1
         else
            n_bands = max(2*rscheme_oc%order_boundary+1,rscheme_oc%order+1)
         end if

         !print*, 'S', n_bands
         call s0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+rank)
            call sMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
      else
         allocate( type_densemat :: sMat(nLMBs2(1+rank)) )
         allocate( type_densemat :: s0Mat )

         call s0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
         do ll=1,nLMBs2(1+rank)
            call sMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
         end do
      end if

#ifdef WITH_PRECOND_S
      allocate(sMat_fac(n_r_max,nLMBs2(1+rank)))
      bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
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
      allocate( rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateS
!------------------------------------------------------------------------------
   subroutine finalize_updateS

      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+rank)
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
   subroutine updateS(s,ds,w,dVSrLM,dsdt,dsdtLast,w1,coex,dt)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1        ! weight for time step !
      real(cp),    intent(in) :: coex      ! factor depending on alpha
      real(cp),    intent(in) :: dt        ! time step
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdtLast(llm:ulm,n_r_max)
      !-- Output: udpated s,ds,dsdtLast
      complex(cp), intent(out) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: start_lm,stop_lm,threadid
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_s ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank
      w2  =one-w1
      O_dt=one/dt

      !PERFON('upS_fin')
      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM,work_LMloc,ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(nR,lm)
      do nR=1,n_r_max
         do lm=llm,ulm
            dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*work_LMloc(lm,nR)- &
            &           dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(nR)*       &
            &           dentropy0(nR)*w(lm,nR))
         end do
      end do
      !$omp end do
      !PERFOFF

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
               call get_s0Mat(dt,s0Mat,s0Mat_fac)
#else
               call get_s0Mat(dt,s0Mat)
#endif
               lSmat(l1)=.true.
            end if
         else
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(nLMB2),sMat_fac(:,nLMB2))
#else
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)),sMat(nLMB2))
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
                     rhs(nR)=real(s(lm1,nR))*O_dt+w1*real(dsdt(lm1,nR))  + &
                     &       w2*real(dsdtLast(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs(:) = s0Mat_fac(:)*rhs(:)
#endif

                  call s0Mat%solve(rhs)

#ifdef TOTO
                  block
                     use radial_functions, only: r
                     integer :: ff, ii

                     open(newunit=ff, file='s0.txt')
                     do nR=1,n_r_max
                        write(ff, '(2ES20.12)') r(nR), rhs(nR)
                     end do
                     close(ff)
                     stop

                  end block
#endif

               else ! l1  /=  0

                  lmB=lmB+1

                  rhs1(1,lmB,threadid)      =tops(l1,m1)
                  rhs1(n_r_max,lmB,threadid)=bots(l1,m1)

                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt+w1*dsdt(lm1,nR)  &
                     &                     +w2*dsdtLast(lm1,nR)
                  end do

#ifdef WITH_PRECOND_S
                  rhs1(:,lmB,threadid) = sMat_fac(:,nLMB2)*rhs1(:,lmB,threadid)
#endif
               end if
            end do
            !PERFOFF

            !PERFON('upS_sol')
            if ( lmB  >  lmB0 ) then
               call sMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
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
                        s(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                    0.0_cp,kind=cp)
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
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            s(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      !$omp single
      call dct_counter%start_count()
      !$omp end single

      call get_ddr(s, ds, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
      call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      !$omp barrier

      !-- Calculate explicit time step part:
      !$omp do private(nR,lm1)
      do nR=n_r_cmb+1,n_r_icb-1
         do lm1=llm,ulm
            dsdtLast(lm1,nR)=                              dsdt(lm1,nR) &
            &      - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1))) * &
            &        kappa(nR) *                   ( work_LMloc(lm1,nR) &
            &        + ( beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR) ) &
            &                                              * ds(lm1,nR) &
            &        - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)     &
            &                                              *  s(lm1,nR) )
         end do
      end do
      !$omp end do

      !PERFOFF
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      !$omp end parallel

   end subroutine updateS
!------------------------------------------------------------------------------
   subroutine updateS_ala(s,ds,w,dVSrLM,dsdt,dsdtLast,w1,coex,dt)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1        ! weight for time step !
      real(cp),    intent(in) :: coex      ! factor depending on alpha
      real(cp),    intent(in) :: dt        ! time step
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdtLast(llm:ulm,n_r_max)
      !-- Output: udpated s,ds,dsdtLast
      complex(cp), intent(out) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid, start_lm, stop_lm
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_s ) return

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m


      nLMB=1+rank
      w2  =one-w1
      O_dt=one/dt


      !PERFON('upS_fin')
      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm, stop_lm)

      !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
      call get_dr( dVSrLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
           &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(nR,lm)
      do nR=1,n_r_max
         do lm=llm,ulm
            dsdt(lm,nR)=           orho1(nR)*        dsdt(lm,nR) - &
            &        or2(nR)*orho1(nR)*        work_LMloc(lm,nR) + &
            &        or2(nR)*orho1(nR)*dLtemp0(nR)*dVSrLM(lm,nR) - &
            &        dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(nR)*   &
            &        orho1(nR)*temp0(nR)*dentropy0(nR)*w(lm,nR)

         end do
      end do
      !$omp end do
      !PERFOFF

      !$omp single
      call solve_counter%start_count()
      !$omp end single

      ! one subblock is linked to one l value and needs therefore once the matrix
      !$OMP SINGLE
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
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 == 0 ) then
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_s0Mat(dt,s0Mat,s0Mat_fac)
#else
               call get_s0Mat(dt,s0Mat)
#endif
               lSmat(l1)=.true.
            end if
         else
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(nLMB2),sMat_fac(:,nLMB2))
#else
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)),sMat(nLMB2))
#endif
               lSmat(l1)=.true.
             !write(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
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
                  rhs(1)=real(tops(0,0))
                  rhs(n_r_max)=real(bots(0,0))
                  do nR=2,n_r_max-1
                     rhs(nR)=real(s(lm1,nR))*O_dt+w1*real(dsdt(lm1,nR)) + &
                     &       w2*real(dsdtLast(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs(:) = s0Mat_fac(:)*rhs(:)
#endif

                  call s0Mat%solve(rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)=tops(l1,m1)
                  rhs1(n_r_max,lmB,threadid)=bots(l1,m1)

                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt+w1*dsdt(lm1,nR)  &
                     &                     +w2*dsdtLast(lm1,nR)
                  end do

#ifdef WITH_PRECOND_S
                  do nR=1,n_r_max
                     rhs1(nR,lmB,threadid) = sMat_fac(nR,nLMB2)*rhs1(nR,lmB,threadid)
                  end do
#endif
               end if
            end do
            !PERFOFF

            !PERFON('upS_sol')
            if ( lmB  >  lmB0 ) then
               call sMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
            end if
            !PERFOFF

            lmB=lmB0
            !PERFON('upS_af')
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
             !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                    0.0_cp,kind=cp)
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
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            s(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      !$omp single
      call dct_counter%start_count()
      !$omp end single

      call get_ddr(s, ds, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
      call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      !$omp barrier

      !-- Calculate explicit time step part:
      !$omp do private(nR,lm1)
      do nR=n_r_cmb+1,n_r_icb-1
         do lm1=llm,ulm
           dsdtLast(lm1,nR)=dsdt(lm1,nR) &
           &      - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1)))*kappa(nR) * &
           &        (                                         work_LMloc(lm1,nR) &
           &                + ( beta(nR)+two*or1(nR)+dLkappa(nR) ) *  ds(lm1,nR) &
           &          - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*  s(lm1,nR) )
         end do
      end do
      !$omp end do
      !PERFOFF
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      !$omp end parallel

   end subroutine updateS_ala
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_s0Mat(dt,sMat,sMat_fac)
#else
   subroutine get_s0Mat(dt,sMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  sMat0
      !

      !-- Input variables
      real(cp), intent(in) :: dt

      !-- Output variables
      class(type_realmat), intent(inout) :: sMat
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max)
      integer :: info,nR_out,nR
      real(cp) :: O_dt

      O_dt=one/dt

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
               &                       O_dt*rscheme_oc%rMat(nR,nR_out) - &
               & alpha*opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
               & (beta(nR)+two*or1(nR)+dLkappa(nR))*                     &
               &                           rscheme_oc%drMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                     &
               &                      O_dt*rscheme_oc%rMat(nR,nR_out) - &
               &alpha*opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
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

#ifdef TOTO
      block

         integer :: fileHandle

         open(newunit=fileHandle, file='s0Mat', form='unformatted', access='stream')

         write(fileHandle) dat
         write(fileHandle) sMat%dat

         close(fileHandle)

      end block
#endif

      !---- LU decomposition:
      call sMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix sMat0!')

   end subroutine get_s0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_sMat(dt,l,hdif,sMat,sMat_fac)
#else
   subroutine get_sMat(dt,l,hdif,sMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
      !  sMat(i,j) and s0mat for the entropy equation.
      !

      !-- Input variables
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: sMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: O_dt,dLh
      real(cp) :: dat(n_r_max,n_r_max)

#ifdef MATRIX_CHECK
      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

      O_dt=one/dt

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktops == 1 ) then
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
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
               &                          O_dt*rscheme_oc%rMat(nR,nR_out) - &
               & alpha*opr*hdif*kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*                       &
               &                              rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &                          O_dt*rscheme_oc%rMat(nR,nR_out) - &
               & alpha*opr*hdif*kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
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
#endif

      !-- Array copy
      call sMat%set_data(dat)

      !-- LU decomposition:
      call sMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix sMat!')

   end subroutine get_Smat
!-----------------------------------------------------------------------------
end module updateS_mod
