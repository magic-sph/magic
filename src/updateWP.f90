#include "perflib_preproc.cpp"
module updateWP_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: or1, or2, rho0, rgrav, visc, dLvisc, r, &
       &                       alpha0, temp0, beta, dbeta, ogrun,      &
       &                       rscheme_oc, ddLvisc, ddbeta, orho1
   use physical_parameters, only: kbotv, ktopv, ra, BuoFac, ChemFac,    &
       &                          ViscHeatFac, ThExpNb, ktopp
   use num_param, only: alpha, dct_counter, solve_counter
   use blocking, only: lo_sub_map, lo_map, st_map, st_sub_map, llm, ulm
   use horizontal_data, only: hdif_V, dLh
   use logic, only: l_update_v, l_chemical_conv, l_RMS, l_double_curl, &
       &            l_fluxProfs, l_finite_diff
   use RMS, only: DifPol2hInt, DifPolLMr
   use algebra, only: prepare_mat, solve_mat
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize, rank, n_procs, get_openmp_blocks
   use RMS_helpers, only:  hInt2Pol
   use radial_der, only: get_dddr, get_ddr, get_dr
   use integration, only: rInt_R
   use fields, only: work_LMloc
   use constants, only: zero, one, two, three, four, third, half
   use useful, only: abortRun
   use dense_matrices
   use real_matrices
   use band_matrices

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: ddddw(:,:)
   complex(cp), allocatable :: dwold(:,:)
   real(cp), allocatable :: work(:)
   complex(cp), allocatable :: Dif(:),Pre(:),Buo(:)
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: wpMat_fac(:,:,:)
   class(type_realmat), allocatable :: wpMat(:), p0Mat
   logical, public, allocatable :: lWPmat(:)
   integer :: maxThreads, size_rhs1

   public :: initialize_updateWP, finalize_updateWP, updateWP

contains

   subroutine initialize_updateWP

      integer, pointer :: nLMBs2(:)
      integer :: ll, n_bands

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      if ( l_finite_diff ) then
         allocate( type_bandmat :: wpMat(nLMBs2(1+rank)) )

         if ( rscheme_oc%order <= 2 .and. rscheme_oc%order_boundary <= 2 ) then
            n_bands =rscheme_oc%order+3
         else
            n_bands = max(rscheme_oc%order+3,2*rscheme_oc%order_boundary+3)
         end if
         !print*, 'WP', n_bands
         do ll=1,nLMBs2(1+rank)
            call wpMat(ll)%initialize(n_bands,n_r_max,l_pivot=.true.)
         end do
         allocate( wpMat_fac(n_r_max,2,nLMBs2(1+rank)) )
         bytes_allocated=bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
         &               SIZEOF_DEF_REAL

         allocate( type_bandmat :: p0Mat )

         n_bands = rscheme_oc%order+1
         call p0Mat%initialize(n_bands,n_r_max,l_pivot=.true.)
      else
         allocate( type_densemat :: wpMat(nLMBs2(1+rank)) )
         if ( l_double_curl ) then
            do ll=1,nLMBs2(1+rank)
               call wpMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
            allocate( wpMat_fac(n_r_max,2,nLMBs2(1+rank)) )
            bytes_allocated=bytes_allocated+2*n_r_max*nLMBs2(1+rank)*    &
            &               SIZEOF_DEF_REAL
         else
            do ll=1,nLMBs2(1+rank)
               call wpMat(ll)%initialize(2*n_r_max,2*n_r_max,l_pivot=.true.)
            end do
            allocate( wpMat_fac(2*n_r_max,2,nLMBs2(1+rank)) )
            bytes_allocated=bytes_allocated+4*n_r_max*nLMBs2(1+rank)*    &
            &               SIZEOF_DEF_REAL
         end if

         allocate( type_densemat :: p0Mat )
         call p0Mat%initialize(n_r_max,n_r_max,l_pivot=.true.)
      end if

      allocate( lWPmat(0:l_max) )
      bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      if ( l_double_curl ) then
         allocate( ddddw(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         if ( l_RMS .or. l_FluxProfs ) then
            allocate( dwold(llm:ulm,n_r_max) )
            bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         end if
      end if

      allocate( work(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL

      allocate( Dif(llm:ulm) )
      allocate( Pre(llm:ulm) )
      allocate( Buo(llm:ulm) )
      bytes_allocated = bytes_allocated+3*(ulm-llm+1)*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      if ( l_double_curl ) then
         size_rhs1 = n_r_max
         allocate( rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
         &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
      else
         size_rhs1 = 2*n_r_max
         allocate( rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
         &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_updateWP
!-----------------------------------------------------------------------------
   subroutine finalize_updateWP

      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+rank)
         call wpMat(ll)%finalize()
      end do
      call p0Mat%finalize()

      deallocate( wpMat_fac,lWPmat )
      deallocate( rhs1, work )
      deallocate( Dif, Pre, Buo )
      if ( l_double_curl ) then
         deallocate( ddddw )
         if ( l_RMS .or. l_FluxProfs ) then
            deallocate( dwold )
         end if
      end if

   end subroutine finalize_updateWP
!-----------------------------------------------------------------------------
   subroutine updateWP(w,dw,ddw,dVxVhLM,dwdt,dwdtLast,p,dp,dpdt,dpdtLast,s,xi, &
              &        w1,coex,dt,lRmsNext,lPressNext)
      !
      !  updates the poloidal velocity potential w, the pressure p,  and
      !  their derivatives
      !  adds explicit part to time derivatives of w and p
      !

      !-- Input/output of scalar fields:
      real(cp),    intent(in) :: w1       ! weight for time step !
      real(cp),    intent(in) :: coex     ! factor depending on alpha
      real(cp),    intent(in) :: dt       ! time step
      logical,     intent(in) :: lRmsNext
      logical,     intent(in) :: lPressNext
      complex(cp), intent(in) :: dpdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: s(llm:ulm,n_r_max)
      complex(cp), intent(in) :: xi(llm:ulm,n_r_max)

      complex(cp), intent(inout) :: dwdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVxVhLM(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dwdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dpdtLast(llm:ulm,n_r_max)

      complex(cp), intent(out) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm,lmB     ! position of (l,m) in array
      integer :: lmStart_00     ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR             ! counts radial grid points
      integer :: n_r_out         ! counts cheb modes
      real(cp) :: rhs(n_r_max)  ! real RHS for l=m=0
      integer :: n_r_top, n_r_bot, nLMB

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: start_lm, stop_lm
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

      nLMB       =1+rank
      lmStart_00 =max(2,llm)

      w2  =one-w1
      O_dt=one/dt

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      if ( l_double_curl ) then
         call get_dr( dVxVhLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
              &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
         !$omp barrier

         !$omp do private(nR,lm)
         do nR=1,n_r_max
            do lm=llm,ulm
               dwdt(lm,nR)= dwdt(lm,nR)+or2(nR)*work_LMloc(lm,nR)
            end do
         end do
         !$omp end do
      end if

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      !PERFON('upWP_ssol')
      !$OMP SINGLE
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)
         !write(*,"(2(A,I3))") "Constructing next task for ",nLMB2,"/",nLMBs2(nLMB)

         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,iChunk,nChunks,size_of_last_chunk,threadid) &
         !$OMP shared(dwold,nLMB,nLMBs2,rhs1)

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
               !PERFON('upWP_mat')
               if ( l_double_curl ) then
                  call get_wMat(dt,l1,hdif_V(st_map%lm2(l1,0)),    &
                       &        wpMat(nLMB2),wpMat_fac(:,:,nLMB2))
               else
                  call get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,0)),   &
                       &         wpMat(nLMB2),wpMat_fac(:,:,nLMB2))
               end if
               lWPmat(l1)=.true.
               !PERFOFF
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(threadid)

            !PERFON('upWP_set')
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
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  !-- The integral of rho' r^2 dr vanishes
                  if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then
                     if ( rscheme_oc%version == 'cheb' ) then
                        do nR=1,n_r_max
                           work(nR)=ThExpNb*alpha0(nR)*temp0(nR)*rho0(nR)*r(nR)*&
                           &        r(nR)*real(s(st_map%lm2(0,0),nR))
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
                        rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*    &
                        &       real(s(st_map%lm2(0,0),nR))+  &
                        &       rho0(nR)*ChemFac*rgrav(nR)*   &
                        &       real(xi(st_map%lm2(0,0),nR))+ &
                        &       real(dwdt(st_map%lm2(0,0),nR))
                     end do
                  else
                     do nR=2,n_r_max
                        rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*    &
                        &       real(s(st_map%lm2(0,0),nR))+  &
                        &       real(dwdt(st_map%lm2(0,0),nR))
                     end do
                  end if

                  call p0Mat%solve(rhs)

               else ! l1 /= 0
                  lmB=lmB+1
                  rhs1(1,lmB,threadid)      =0.0_cp
                  rhs1(n_r_max,lmB,threadid)=0.0_cp
                  if ( l_double_curl ) then
                     rhs1(2,lmB,threadid)        =0.0_cp
                     rhs1(n_r_max-1,lmB,threadid)=0.0_cp
                     if ( l_chemical_conv ) then
                        do nR=3,n_r_max-2
                           rhs1(nR,lmB,threadid)=dLh(st_map%lm2(l1,m1))*or2(nR)* (   &
                           &                     -orho1(nR)*O_dt*(    ddw(lm1,nR)    &
                           &                     -beta(nR)*dw(lm1,nR)-               &
                           &                     dLh(st_map%lm2(l1,m1))*or2(nR)*     &
                           &                                w(lm1,nR) ) +            &
                           &                     alpha*BuoFac *rgrav(nR)* s(lm1,nR)+ &
                           &                     alpha*ChemFac*rgrav(nR)*xi(lm1,nR) )&
                           &                     +w1*dwdt(lm1,nR) +                  &
                           &                     w2*dwdtLast(lm1,nR)
                        end do
                     else
                        do nR=3,n_r_max-2
                           rhs1(nR,lmB,threadid)=dLh(st_map%lm2(l1,m1))*or2(nR)* (   &
                           &                     -orho1(nR)*O_dt*(    ddw(lm1,nR)    &
                           &                     -beta(nR)*dw(lm1,nR)-               &
                           &                     dLh(st_map%lm2(l1,m1))*or2(nR)*     &
                           &                                w(lm1,nR) ) +            &
                           &                     alpha*BuoFac *rgrav(nR)* s(lm1,nR) )&
                           &                     +w1*dwdt(lm1,nR) +                  &
                           &                     w2*dwdtLast(lm1,nR)
                        end do
                     end if
                  else
                     rhs1(n_r_max+1,lmB,threadid)=0.0_cp
                     rhs1(2*n_r_max,lmB,threadid)=0.0_cp
                     if ( l_chemical_conv ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,lmB,threadid)=O_dt*dLh(st_map%lm2(l1,m1))*    &
                           &                     or2(nR)*w(lm1,nR) +             &
                           &                     rho0(nR)*alpha*BuoFac*rgrav(nR)*&
                           &                     s(lm1,nR) + rho0(nR)*alpha*     &
                           &                     ChemFac*rgrav(nR)*xi(lm1,nR) +  &
                           &                     w1*dwdt(lm1,nR) +               &
                           &                     w2*dwdtLast(lm1,nR)
                           rhs1(nR+n_r_max,lmB,threadid)=-O_dt*dLh(st_map%lm2(l1,&
                           &                           m1))*or2(nR)*dw(lm1,nR) + &
                           &                              w1*dpdt(lm1,nR) +      &
                           &                              w2*dpdtLast(lm1,nR)
                        end do
                     else
                        do nR=2,n_r_max-1
                           rhs1(nR,lmB,threadid)=O_dt*dLh(st_map%lm2(l1,m1))* &
                           &                     or2(nR)*w(lm1,nR) +          &
                           &                     rho0(nR)*alpha*BuoFac*       &
                           &                     rgrav(nR)*s(lm1,nR) +        &
                           &                     w1*dwdt(lm1,nR) +            &
                           &                     w2*dwdtLast(lm1,nR)
                           rhs1(nR+n_r_max,lmB,threadid)=-O_dt*dLh(st_map%lm2(l1,&
                           &                           m1))*or2(nR)*dw(lm1,nR) + &
                           &                             w1*dpdt(lm1,nR) +       &
                           &                             w2*dpdtLast(lm1,nR)
                        end do
                     end if
                  end if
               end if
            end do
            !PERFOFF

            !PERFON('upWP_sol')
            if ( lmB > 0 ) then

               ! use the mat_fac(:,1) to scale the rhs
               do lm=lmB0+1,lmB
                  do nR=1,size_rhs1
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,1,nLMB2)
                  end do
               end do
               call wpMat(nLMB2)%solve(rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,size_rhs1
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,2,nLMB2)
                  end do
               end do
            end if
            !PERFOFF

            if ( l_double_curl .and. lPressNext ) then ! Store old dw
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     dwold(lm1,nR)=dw(lm1,nR)
                  end do
               end do
            end if

            !PERFON('upWP_aft')
            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     p(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( l_double_curl ) then
                     if ( m1 > 0 ) then
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  =rhs1(n_r_out,lmB,threadid)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  = cmplx(real(rhs1(n_r_out,lmB,threadid)),&
                           &                       0.0_cp,kind=cp)
                        end do
                     end if
                  else
                     if ( m1 > 0 ) then
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                           p(lm1,n_r_out)=rhs1(n_r_max+n_r_out,lmB,threadid)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                           &                    0.0_cp,kind=cp)
                           p(lm1,n_r_out)= cmplx(real(rhs1(n_r_max+n_r_out,lmB, &
                           &                    threadid)),0.0_cp,kind=cp)
                        end do
                     end if
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do   ! end of loop over l1 subblocks
      !$OMP END SINGLE
      !PERFOFF
      !$omp single
      call solve_counter%stop_count()
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp single
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            w(lm1,n_r_out)=zero
            p(lm1,n_r_out)=zero
         end do
      end do
      !$omp end single

      !$omp single
      call dct_counter%start_count()
      !$omp end single

      if ( l_double_curl ) then
         call get_ddr( w, dw, ddw, ulm-llm+1, start_lm-llm+1,  &
              &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
         call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1,                  &
              &        start_lm-llm+1, stop_lm-llm+1, n_r_max, rscheme_oc, &
              &        l_dct_in=.false. )
         call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)

      else
         call get_dddr( w, dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
              &         stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc, l_dct_in=.false. )
      end if
      call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      call rscheme_oc%costf1(p,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      !$omp barrier

      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      if ( lRmsNext ) then
         n_r_top=n_r_cmb
         n_r_bot=n_r_icb
      else
         n_r_top=n_r_cmb+1
         n_r_bot=n_r_icb-1
      end if

      !PERFON('upWP_ex')
      !-- Calculate explicit time step part:
      if ( l_double_curl ) then

         if ( lPressNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         end if

         !$omp do private(nR,lm1,l1,m1,Dif,Buo)
         do nR=n_r_top,n_r_bot
            do lm1=lmStart_00,ulm
               l1=lm2l(lm1)
               m1=lm2m(lm1)

               Dif(lm1) = -hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*  &
               &          or2(nR)*visc(nR) * orho1(nR)*       ( ddddw(lm1,nR) &
               &            +two*( dLvisc(nR)-beta(nR) ) * work_LMloc(lm1,nR) &
               &        +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+   &
               &           beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-two*   &
               &           or1(nR)*(dLvisc(nR)+beta(nR))-two*or2(nR)*         &
               &           dLh(st_map%lm2(l1,m1)) ) *             ddw(lm1,nR) &
               &        +( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+    &
               &           two*or1(nR))-ddLvisc(nR)*(beta(nR)+two*or1(nR))+   &
               &           beta(nR)*beta(nR)*(dLvisc(nR)+two*or1(nR))-        &
               &           beta(nR)*(dLvisc(nR)*dLvisc(nR)-two*or2(nR))-      &
               &           two*dLvisc(nR)*or1(nR)*(dLvisc(nR)-or1(nR))+       &
               &           two*(two*or1(nR)+beta(nR)-dLvisc(nR))*or2(nR)*     &
               &           dLh(st_map%lm2(l1,m1)) ) *              dw(lm1,nR) &
               &        + dLh(st_map%lm2(l1,m1))*or2(nR)* ( two*dbeta(nR)+    &
               &           ddLvisc(nR)+dLvisc(nR)*dLvisc(nR)-two*third*       &
               &           beta(nR)*beta(nR)+dLvisc(nR)*beta(nR)+two*or1(nR)* &
               &           (two*dLvisc(nR)-beta(nR)-three*or1(nR))+           &
               &           dLh(st_map%lm2(l1,m1))*or2(nR) ) *       w(lm1,nR) )

               Buo(lm1) = BuoFac*dLh(st_map%lm2(l1,m1))*or2(nR)*rgrav(nR)*s(lm1,nR)
               if ( l_chemical_conv ) then
                  Buo(lm1) = Buo(lm1)+ChemFac*dLh(st_map%lm2(l1,m1))*or2(nR)*&
                  &          rgrav(nR)*xi(lm1,nR)
               end if

               dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Buo(lm1)+Dif(lm1))

               if ( l1 /= 0 .and. lPressNext ) then
                  ! In the double curl formulation, we can estimate the pressure
                  ! if required.
                  p(lm1,nR)=-r(nR)*r(nR)/dLh(st_map%lm2(l1,m1))*dpdt(lm1,nR) &
                  &                -O_dt*(dw(lm1,nR)-dwold(lm1,nR))+         &
                  &                 hdif_V(st_map%lm2(l1,m1))*visc(nR)*      &
                  &                                    ( work_LMloc(lm1,nR)  &
                  &                       - (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                  &               - ( dLh(st_map%lm2(l1,m1))*or2(nR)         &
                  &                  + dLvisc(nR)*beta(nR)+ dbeta(nR)        &
                  &                  + two*(dLvisc(nR)+beta(nR))*or1(nR)     &
                  &                                           ) * dw(lm1,nR) &
                  &               + dLh(st_map%lm2(l1,m1))*or2(nR)           &
                  &                  * ( two*or1(nR)+two*third*beta(nR)      &
                  &                     +dLvisc(nR) )   *         w(lm1,nR)  &
                  &                                         )
               end if

               if ( lRmsNext ) then
                  !-- In case RMS force balance is required, one needs to also
                  !-- compute the classical diffusivity that is used in the non
                  !-- double-curl version
                  Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
                  &          or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                  &        +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
                  &        -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
                  &             dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                  &             +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
                  &                                                 w(lm1,nR)  )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,ulm, &
                    &        DifPolLMr(llm:ulm,nR),         &
                    &        DifPol2hInt(:,nR),lo_map)
            end if
         end do
         !$omp end do

      else

         !$omp do private(nR,lm1,l1,m1,Dif,Buo,Pre)
         do nR=n_r_top,n_r_bot
            do lm1=lmStart_00,ulm
               l1=lm2l(lm1)
               m1=lm2m(lm1)

               Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
               &          or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
               &        +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
               &        -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
               &             dbeta(nR)+dLvisc(nR)*beta(nR)                  &
               &             +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
               &                                                 w(lm1,nR)  )
               Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
               Buo(lm1) = BuoFac*rho0(nR)*rgrav(nR)*s(lm1,nR)
               if ( l_chemical_conv ) then
                  Buo(lm1) = Buo(lm1)+ChemFac*rho0(nR)*rgrav(nR)*xi(lm1,nR)
               end if
               dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
               dpdtLast(lm1,nR)= dpdt(lm1,nR) - coex*(                    &
               &                 dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
               &               + hdif_V(st_map%lm2(l1,m1))*               &
               &                 visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR)  &
               &                                  * ( -work_LMloc(lm1,nR) &
               &                       + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
               &               + ( dLh(st_map%lm2(l1,m1))*or2(nR)         &
               &                  + dLvisc(nR)*beta(nR)+ dbeta(nR)        &
               &                  + two*(dLvisc(nR)+beta(nR))*or1(nR)     &
               &                                           ) * dw(lm1,nR) &
               &               - dLh(st_map%lm2(l1,m1))*or2(nR)           &
               &                  * ( two*or1(nR)+two*third*beta(nR)      &
               &                     +dLvisc(nR) )   *         w(lm1,nR)  &
               &                                         ) )
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,ulm, &
                    &        DifPolLMr(llm:ulm,nR),         &
                    &        DifPol2hInt(:,nR),lo_map)
            end if
         end do
         !$omp end do

      end if
      !PERFOFF

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext .and. l_double_curl ) then
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc)
         !$omp barrier
      end if

      !$omp end parallel

   end subroutine updateWP
!------------------------------------------------------------------------------
   subroutine get_wpMat(dt,l,hdif,wpMat,wpMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation.
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables:
      class(type_realmat), intent(inout) :: wpMat
      real(cp), intent(out) :: wpMat_fac(2*n_r_max,2)

      !-- local variables:
      integer :: nR,nR_out,nR_p,nR_out_p
      integer :: info
      real(cp) :: O_dt,dLh

#ifdef MATRIX_CHECK
      integer ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
      real(cp) :: work(8*n_r_max),anorm,linesum,rcond
      real(cp) :: temp_wpMat(2*n_r_max,2*n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
      logical :: first_run=.true.
#endif

      O_dt=one/dt
      dLh =real(l*(l+1),kind=cp)

      !-- Now mode l>0

      !----- Boundary conditions, see above:
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

         if ( kbotv == 1 ) then  ! free slip !
            wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * (       &
            &                  rscheme_oc%d2rMat(n_r_max,nR_out) - &
            &      ( two*or1(n_r_max)+beta(n_r_max))*              &
            &                  rscheme_oc%drMat(n_r_max,nR_out))
         else                 ! no slip, note exception for l=1,m=0
            wpMat%dat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
            &                       rscheme_oc%drMat(n_r_max,nR_out)
         end if
         wpMat%dat(2*n_r_max,nR_out_p)=0.0_cp

      end do   !  loop over nR_out

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
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
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         nR_out_p=nR_out+n_r_max
         do nR=2,n_r_max-1
            !write(*,"(I3,A,6ES11.3)") nR,", visc,beta,dLvisc,dbeta = ",&
            !     & visc(nR),beta(nR),dLvisc(nR),dbeta(nR),hdif,alpha
            ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
            nR_p=nR+n_r_max
            wpMat%dat(nR,nR_out)= rscheme_oc%rnorm *  (                     &
            &               O_dt*dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)     &
            &            - alpha*hdif*visc(nR)*dLh*or2(nR) * (              &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &        +(two*dLvisc(nR)-third*beta(nR))*                      &
            &                               rscheme_oc%drMat(nR,nR_out)     &
            &        -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)        &
            &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
            &          )                    *rscheme_oc%rMat(nR,nR_out)  )  )

            wpMat%dat(nR,nR_out_p)= rscheme_oc%rnorm*alpha*(             &
            &                            rscheme_oc%drMat(nR,nR_out)     &
            &                  -beta(nR)* rscheme_oc%rMat(nR,nR_out))
            ! the following part gives sometimes very large
            ! matrix entries
            wpMat%dat(nR_p,nR_out)=rscheme_oc%rnorm * (                       &
            &                  -O_dt*dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)  &
            &         -alpha*hdif*visc(nR)*dLh*or2(nR)      *(                &
            &                                 - rscheme_oc%d3rMat(nR,nR_out)  &
            &          +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out)  &
            &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
            &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
            &                                    rscheme_oc%drMat(nR,nR_out)  &
            &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
            &          +two*third*beta(nR)   )*   rscheme_oc%rMat(nR,nR_out) ) )

            wpMat%dat(nR_p,nR_out_p)= -rscheme_oc%rnorm*alpha*dLh*or2(nR)* &
            &                      rscheme_oc%rMat(nR,nR_out)
         end do
      end do

      !----- Factor for highest and lowest cheb:
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

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         wpMat_fac(nR,1)=one/maxval(abs(wpMat%dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         wpMat%dat(nR,:) = wpMat%dat(nR,:)*wpMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         wpMat_fac(nR,2)=one/maxval(abs(wpMat%dat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         wpMat%dat(:,nR) = wpMat%dat(:,nR)*wpMat_fac(nR,2)
      end do

#ifdef MATRIX_CHECK
      ! copy the wpMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

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
#endif

      call wpMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix wpMat!')

   end subroutine get_wpMat
!-----------------------------------------------------------------------------
   subroutine get_wMat(dt,l,hdif,wMat,wMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation. This matrix corresponds here to the
      !  radial component of the double-curl of the Navier-Stokes equation.
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables:
      class(type_realmat), intent(inout) :: wMat
      real(cp), intent(out) :: wMat_fac(n_r_max,2)

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: O_dt, dLh
      real(cp) :: dat(n_r_max,n_r_max)

      O_dt=one/dt
      dLh =real(l*(l+1),kind=cp)

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

      if ( kbotv == 1 ) then  ! free slip
         dat(n_r_max-1,:)=rscheme_oc%rnorm *(rscheme_oc%d2rMat(n_r_max,:)-  &
         &                      (two*or1(n_r_max)+beta(n_r_max))*           &
         &                                    rscheme_oc%drMat(n_r_max,:) )
      else                 ! no slip
         dat(n_r_max-1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)        =0.0_cp
            dat(2,nR_out)        =0.0_cp
            dat(n_r_max-1,nR_out)=0.0_cp
            dat(n_r_max,nR_out)  =0.0_cp
         end do
      end if

      !----- Bulk points:
      do nR_out=1,n_r_max
         do nR=3,n_r_max-2
            dat(nR,nR_out)=rscheme_oc%rnorm* (                              &
            &          -O_dt*dLh*or2(nR)*orho1(nR)*(                        &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out) -   &
            &                    dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )   &
            &    + alpha*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (            &
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

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,n_r_max
         wMat_fac(nR,1)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*wMat_fac(nR,1)
      end do

      ! also compute the rowsum of each column
      do nR=1,n_r_max
         wMat_fac(nR,2)=one/maxval(abs(dat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,n_r_max
         dat(:,nR) = dat(:,nR)*wMat_fac(nR,2)
      end do

      !-- Array copy
      call wMat%set_data(dat)

#ifdef TOTO
      block
         use radial_functions, only: r
         integer :: fileHandle

         if ( l == 4 ) then
            open(newunit=fileHandle, file='wMat', form='unformatted', access='stream')

            write(fileHandle) r
            write(fileHandle) dat
            write(fileHandle) wMat%dat

            close(fileHandle)
         end if

      end block
#endif

      call wMat%prepare(info)

      if ( info /= 0 ) call abortRun('Singular matrix wMat!')

   end subroutine get_wMat
!-----------------------------------------------------------------------------
   subroutine get_p0Mat(pMat)

      !-- Output variables:
      class(type_realmat), intent(inout) :: pMat

      !-- Local variables:
      real(cp) :: dat(n_r_max,n_r_max), delr
      integer :: info, nCheb, nR_out, nR, nCheb_in

      !-- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max
            dat(nR,nR_out)= rscheme_oc%rnorm * (              &
            &                    rscheme_oc%drMat(nR,nR_out)- &
            &            beta(nR)*rscheme_oc%rMat(nR,nR_out) )
         end do
      end do

      !-- Boundary condition for spherically-symmetric pressure
      !-- The integral of rho' r^2 dr vanishes
      if ( ThExpNb*ViscHeatFac /= 0 .and. ktopp==1 ) then

         work(:) = ThExpNb*ViscHeatFac*ogrun(:)*alpha0(:)*r(:)*r(:)
         call rscheme_oc%costf1(work)
         work(:)      =work(:)*rscheme_oc%rnorm
         work(1)      =rscheme_oc%boundary_fac*work(1)
         work(n_r_max)=rscheme_oc%boundary_fac*work(n_r_max)

         if ( rscheme_oc%version == 'cheb' ) then

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

         else

            !!-- In the finite differences case, we restrict the integral boundary
            !!-- condition to a trapezoidal rule of integration
            !do nR_out=2,rscheme_oc%n_max-1
            !   dat(1,nR_out)=half*work(nR_out)*( r(nR_out+1)-r(nR_out-1) )
            !end do
            !dat(1,1)=half*work(1)*(r(2)-r(1))
            dat(1,:)=rscheme_oc%rMat(1,:)
            delr = r(n_r_max)-r(n_r_max-1)

            !-- First-order on the last point
            dat(n_r_max,n_r_max)    =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)  =-one/delr
            dat(n_r_max,1:n_r_max-2)=0.0_cp
         end if

      else

         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
         if ( rscheme_oc%version == 'fd' ) then
            delr = r(n_r_max)-r(n_r_max-1)
            !-- First-order on the last point
            dat(n_r_max,n_r_max)  =one/delr-beta(n_r_max)
            dat(n_r_max,n_r_max-1)=-one/delr
         end if

      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      !-- Array copy
      call pMat%set_data(dat)

#ifdef TOTO
      block 

         integer :: fileHandle 

         open(newunit=fileHandle, file='p0Mat', form='unformatted', access='stream')

         write(fileHandle) r
         write(fileHandle) dat
         write(fileHandle) pMat%dat
         close(fileHandle)

         stop

      end block
#endif

      !---- LU decomposition:
      call pMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix p0Mat!')

   end subroutine get_p0Mat
!-----------------------------------------------------------------------------
end module updateWP_mod
