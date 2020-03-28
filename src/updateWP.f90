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
   use physical_parameters, only: kbotv, ktopv, ra, BuoFac, ChemFac,   &
       &                          ViscHeatFac, ThExpNb, ktopp
   use num_param, only: dct_counter, solve_counter
   use blocking, only: lo_sub_map, lo_map, st_sub_map, llm, ulm
   use horizontal_data, only: hdif_V
   use logic, only: l_update_v, l_chemical_conv, l_RMS, l_double_curl, &
       &            l_fluxProfs, l_finite_diff, l_full_sphere, l_heat
   use RMS, only: DifPol2hInt, DifPolLMr
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize, rank, n_procs, get_openmp_blocks
   use RMS_helpers, only:  hInt2Pol
   use radial_der, only: get_dddr, get_ddr, get_dr
   use integration, only: rInt_R
   use fields, only: work_LMloc
   use constants, only: zero, one, two, three, four, third, half
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
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
   real(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: rhs0(:,:,:)
   real(cp), allocatable :: wpMat_fac(:,:,:)
   class(type_realmat), allocatable :: wpMat(:), p0Mat, pMat(:)
   logical, public, allocatable :: lWPmat(:)
   logical, allocatable :: lPmat(:)
   integer :: maxThreads, size_rhs1

   public :: initialize_updateWP, finalize_updateWP, updateWP, assemble_pol, &
   &         finish_exp_pol, get_pol_rhs_imp, updateP

contains

   subroutine initialize_updateWP(tscheme)

      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer, pointer :: nLMBs2(:)
      integer :: ll, n_bands

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

         if ( tscheme%l_assembly ) then
            allocate( type_densemat :: pMat(nLMBs2(1+rank)) )
            do ll=1,nLMBs2(1+rank)
               call pMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do
            allocate( lPmat(0:l_max) )
            lPmat(:) = .false.
            allocate( rhs0(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
            rhs0(:,:,:)=zero
         end if
      end if

      allocate( lWPmat(0:l_max) )
      bytes_allocated=bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      if ( l_double_curl ) then
         allocate( ddddw(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         if ( l_RMS .or. l_FluxProfs ) then
            allocate( dwold(llm:ulm,n_r_max) )
            bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
            dwold(:,:)=zero
         end if
      end if

      allocate( work(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL

      allocate( Dif(llm:ulm) )
      allocate( Pre(llm:ulm) )
      allocate( Buo(llm:ulm) )
      bytes_allocated = bytes_allocated+3*(ulm-llm+1)*SIZEOF_DEF_COMPLEX

      if ( l_double_curl ) then
         size_rhs1 = n_r_max
         allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated=bytes_allocated+n_r_max*maxThreads* &
         &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
      else
         size_rhs1 = 2*n_r_max
         allocate( rhs1(2*n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
         &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_updateWP
!-----------------------------------------------------------------------------
   subroutine finalize_updateWP(tscheme)

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variables:
      integer, pointer :: nLMBs2(:)
      integer :: ll

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

      do ll=1,nLMBs2(1+rank)
         call wpMat(ll)%finalize()
      end do
      call p0Mat%finalize()

      if ( .not. l_finite_diff .and. tscheme%l_assembly ) then
         deallocate( lPmat ) 
         do ll=1,nLMBs2(1+rank)
            call pMat(ll)%finalize()
         end do
         deallocate( rhs0 )
      end if

      deallocate( wpMat_fac,lWPmat, rhs1, work )
      deallocate( Dif, Pre, Buo )
      if ( l_double_curl ) then
         deallocate( ddddw )
         if ( l_RMS .or. l_FluxProfs ) then
            deallocate( dwold )
         end if
      end if

   end subroutine finalize_updateWP
!-----------------------------------------------------------------------------
   subroutine updateWP(s, xi, w, dw, ddw, dwdt, p, dp, dpdt, tscheme, &
              &        lRmsNext, lPressNext)
      !
      !  updates the poloidal velocity potential w, the pressure p,  and
      !  their derivatives
      !  adds explicit part to time derivatives of w and p
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

      complex(cp),       intent(out) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm,lmB     ! position of (l,m) in array
      integer :: lmStart_00     ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR             ! counts radial grid points
      integer :: n_r_out         ! counts cheb modes
      real(cp) :: rhs(n_r_max)  ! real RHS for l=m=0
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

      nLMB       =1+rank
      lmStart_00 =max(2,llm)

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dwdt, llm, ulm, n_r_max)
      if ( .not. l_double_curl ) then
         call tscheme%set_imex_rhs(ddw, dpdt, llm, ulm, n_r_max)
      end if

      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      !PERFON('upWP_ssol')
      !$OMP SINGLE
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)

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
                  call get_wMat(tscheme,l1,hdif_V(lm2(l1,0)),    &
                       &        wpMat(nLMB2),wpMat_fac(:,:,nLMB2))
               else
                  call get_wpMat(tscheme,l1,hdif_V(lm2(l1,0)),   &
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
                           &        r(nR)*real(s(lm2(0,0),nR))
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
                        rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm2(0,0),nR))+   &
                        &       rho0(nR)*ChemFac*rgrav(nR)*real(xi(lm2(0,0),nR))+ &
                        &       real(dwdt%expl(lm2(0,0),nR,tscheme%istage))
                     end do
                  else
                     do nR=2,n_r_max
                        rhs(nR)=rho0(nR)*BuoFac*rgrav(nR)*real(s(lm2(0,0),nR))+  &
                        &       real(dwdt%expl(lm2(0,0),nR,tscheme%istage))
                     end do
                  end if

                  call p0Mat%solve(rhs)

               else ! l1 /= 0
                  lmB=lmB+1
                  rhs1(1,2*lmB-1,threadid)      =0.0_cp
                  rhs1(1,2*lmB,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lmB-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lmB,threadid)  =0.0_cp
                  if ( l_double_curl ) then
                     rhs1(2,2*lmB-1,threadid)        =0.0_cp
                     rhs1(2,2*lmB,threadid)          =0.0_cp
                     rhs1(n_r_max-1,2*lmB-1,threadid)=0.0_cp
                     rhs1(n_r_max-1,2*lmB,threadid)  =0.0_cp
                     do nR=3,n_r_max-2
                        rhs1(nR,2*lmB-1,threadid)= real(work_LMloc(lm1,nR))
                        rhs1(nR,2*lmB,threadid)  =aimag(work_LMloc(lm1,nR))
                     end do

                     if ( l_heat ) then
                        do nR=3,n_r_max-2
                           rhs1(nR,2*lmB-1,threadid)=rhs1(nR,2*lmB-1,threadid)+ &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *      &
                           &      or2(nR)*BuoFac*rgrav(nR)*real(s(lm1,nR))
                           rhs1(nR,2*lmB,threadid)  =rhs1(nR,2*lmB,threadid)+   &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) *      &
                           &      or2(nR)*BuoFac*rgrav(nR)*aimag(s(lm1,nR))
                        end do
                     end if

                     if ( l_chemical_conv ) then
                        do nR=3,n_r_max-2
                           rhs1(nR,2*lmB-1,threadid)=rhs1(nR,2*lmB-1,threadid)+ &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) * &
                           &      or2(nR)*ChemFac*rgrav(nR)*real(xi(lm1,nR))
                           rhs1(nR,2*lmB,threadid)  =rhs1(nR,2*lmB,threadid)+   &
                           &      tscheme%wimp_lin(1)*real(l1*(l1+1),cp) * &
                           &      or2(nR)*ChemFac*rgrav(nR)*aimag(xi(lm1,nR))
                        end do
                     end if
                  else
                     rhs1(n_r_max+1,2*lmB-1,threadid)=0.0_cp
                     rhs1(n_r_max+1,2*lmB,threadid)  =0.0_cp
                     rhs1(2*n_r_max,2*lmB-1,threadid)=0.0_cp
                     rhs1(2*n_r_max,2*lmB,threadid)  =0.0_cp
                     do nR=2,n_r_max-1
                        rhs1(nR,2*lmB-1,threadid)        = real(work_LMloc(lm1,nR))
                        rhs1(nR,2*lmB,threadid)          =aimag(work_LMloc(lm1,nR))
                        rhs1(nR+n_r_max,2*lmB-1,threadid)= real(ddw(lm1,nR)) ! ddw is a work array
                        rhs1(nR+n_r_max,2*lmB,threadid)  =aimag(ddw(lm1,nR))
                     end do

                     if ( l_heat ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,2*lmB-1,threadid)=rhs1(nR,2*lmB-1,threadid)+ &
                           &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*       &
                           &                      rgrav(nR)*real(s(lm1,nR))
                           rhs1(nR,2*lmB,threadid)  =rhs1(nR,2*lmB,threadid)+   &
                           &         tscheme%wimp_lin(1)*rho0(nR)*BuoFac*       &
                           &                      rgrav(nR)*aimag(s(lm1,nR))
                        end do
                     end if

                     if ( l_chemical_conv ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,2*lmB-1,threadid)=rhs1(nR,2*lmB-1,threadid)+ &
                           &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*      &
                           &                      rgrav(nR)*real(xi(lm1,nR))
                           rhs1(nR,2*lmB,threadid)  =rhs1(nR,2*lmB,threadid)+   &
                           &         tscheme%wimp_lin(1)*rho0(nR)*ChemFac*      &
                           &                      rgrav(nR)*aimag(xi(lm1,nR))
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
                     rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)* &
                     &                        wpMat_fac(nR,1,nLMB2)
                     rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)* &
                     &                        wpMat_fac(nR,1,nLMB2)
                  end do
               end do
               call wpMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*lmB,threadid), &
                                       2*(lmB-lmB0))
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,size_rhs1
                     rhs1(nR,2*lm-1,threadid)=rhs1(nR,2*lm-1,threadid)* &
                     &                        wpMat_fac(nR,2,nLMB2)
                     rhs1(nR,2*lm,threadid)  =rhs1(nR,2*lm,threadid)* &
                     &                        wpMat_fac(nR,2,nLMB2)
                  end do
               end do
            end if
            !PERFOFF

            if ( l_double_curl .and. lPressNext .and. tscheme%istage == 1) then
               ! Store old dw
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
                           w(lm1,n_r_out)  =cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                           &                      rhs1(n_r_out,2*lmB,threadid),cp)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  = cmplx(rhs1(n_r_out,2*lmB-1,threadid),&
                           &                       0.0_cp,kind=cp)
                        end do
                     end if
                  else
                     if ( m1 > 0 ) then
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                           &                    rhs1(n_r_out,2*lmB,threadid),cp)
                           p(lm1,n_r_out)=cmplx(rhs1(n_r_max+n_r_out,2*lmB-1,   &
                           &                    threadid),rhs1(n_r_max+n_r_out, &
                           &                    2*lmB,threadid),cp)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lmB-1,threadid), &
                           &                    0.0_cp,kind=cp)
                           p(lm1,n_r_out)= cmplx(rhs1(n_r_max+n_r_out,2*lmB-1, &
                           &                    threadid),0.0_cp,kind=cp)
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
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            w(lm1,n_r_out)=zero
            p(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dwdt, llm, ulm, n_r_max)
      if ( .not. l_double_curl ) call tscheme%rotate_imex(dpdt, llm, ulm, n_r_max)

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
   subroutine updateP(s, xi, w, dw, ddw, dwdt, p, dp, dpdt, tscheme)
      !
      !  updates the poloidal velocity potential w, the pressure p,  and
      !  their derivatives
      !  adds explicit part to time derivatives of w and p
      !

      !-- Input/output of scalar fields:
      class(type_tscheme), intent(in) :: tscheme

      type(type_tarray), intent(inout) :: dpdt
      type(type_tarray), intent(inout) :: dwdt
      complex(cp),       intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(in) :: xi(llm:ulm,n_r_max)
      complex(cp),       intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(in) :: ddw(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)

      complex(cp),       intent(out) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm,lmB     ! position of (l,m) in array
      integer :: lmStart_00     ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR             ! counts radial grid points
      integer :: n_r_out         ! counts cheb modes
      real(cp) :: dL
      integer :: nLMB
      complex(cp) :: dpoldr(llm:ulm,n_r_max)
      complex(cp) :: top_bc(llm:ulm)
      complex(cp) :: bot_bc(llm:ulm)

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

      !-- Compute the right hand side
      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=lmStart_00; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp do private(nR,lm)
      !-- Assemble the terms that will need to be derived
      do nR=1,n_r_max
         do lm=lmStart_00,ulm
            work_LMloc(lm,nR)=r(nR)*r(nR)*dwdt%expl(lm,nR,1)
            if ( l_heat ) work_LMloc(lm,nR)=work_LMloc(lm,nR) + &
            &                      r(nR)*r(nR)*BuoFac*rho0(nR)*rgrav(nR)*s(lm,nR)
            if ( l_chemical_conv ) work_LMloc(lm,nR)=work_LMloc(lm,nR) + &
            &                    r(nR)*r(nR)*ChemFac*rho0(nR)*rgrav(nR)*xi(lm,nR)
         end do
      end do
      !$omp end do

      !-- Radial derivative
      call get_dr( work_LMloc, dpoldr, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true.)
      !$omp barrier

      !-- Now add viscosity + horizontal divergence + take 1/r^2 * d/dr
      !$omp do private(nR,lm,l1,dL)
      do nR=1,n_r_max
         do lm=lmStart_00,ulm
            l1 = lm2l(lm)
            dL = real(l1*(l1+1),cp)
            work_LMloc(lm,nR)=or2(nR)*dpoldr(lm,nR)+dpdt%expl(lm,nR,1)  &
            &                 +hdif_V(lm)*visc(nR)*dL*or2(nR)* (        & 
            &         (dLvisc(nR)+two*third*beta(nR)) *      ddw(lm,nR)-&
            &      (two*(third*dbeta(nR)-ddLvisc(nR))+third*dLvisc(nR)* &
            &       beta(nR)+two*or1(nR)*(dLvisc(nR)-third*             &
            &                              beta(nR)))*        dw(lm,nR)-&
            &      (four*third*(ddbeta(nR)+ddLVisc(nR)*beta(nR) +       &
            &       dLvisc(nR)*dbeta(nR)+or1(nR)*(three*(               &
            &       ddLvisc(nR)-dLvisc(nR)*or1(nR))+dbeta(nR)-          &
            &       beta(nR)*or1(nR)))+dL*or2(nR)*(dLvisc(nR)+          &
            &       two*third*beta(nR)))*                      w(lm,nR))
         end do
      end do
      !$omp end do

      !-- Boundary conditions for pressure
      do lm=lmStart_00,ulm
         l1 = lm2l(lm)
         dL = real(l1*(l1+1),cp)
         top_bc(lm)=dwdt%expl(lm,1,1)+hdif_V(lm)*visc(1)*dL*or2(1)*(   &
         &                                               ddw(lm,1)   &
         &        +(two*dLvisc(1)-third*beta(1))*     dw(lm,1)   &
         &        -( dL*or2(1)+four*third*( dbeta(1)+dLvisc(1)*  &
         &           beta(1)+(three*dLvisc(1)+beta(1))*or1(1)))&
         &                                         *       w(lm,1) )
         if ( l_heat ) top_bc(lm)=top_bc(lm)+BuoFac*rho0(1)*rgrav(1)*s(lm,1)
         if ( l_chemical_conv ) top_bc(lm)=top_bc(lm)+ChemFac*rho0(1)*rgrav(1)*xi(lm,1)

         bot_bc(lm)=dwdt%expl(lm,n_r_max,1)+hdif_V(lm)*visc(n_r_max)*dL*or2(n_r_max)*(   &
         &                                               ddw(lm,n_r_max)   &
         &        +(two*dLvisc(n_r_max)-third*beta(n_r_max))*     dw(lm,n_r_max)   &
         &        -( dL*or2(n_r_max)+four*third*( dbeta(n_r_max)+dLvisc(n_r_max)*  &
         &           beta(n_r_max)+(three*dLvisc(n_r_max)+beta(n_r_max))*or1(n_r_max)))&
         &                                         *       w(lm,n_r_max) )
         if ( l_heat ) bot_bc(lm)=bot_bc(lm)+BuoFac*rho0(n_r_max)* &
         &                             rgrav(n_r_max)*s(lm,n_r_max)
         if ( l_chemical_conv ) bot_bc(lm)=bot_bc(lm)+ChemFac*     &
         &                             rho0(n_r_max)*rgrav(n_r_max)*xi(lm,n_r_max)

         !top_bc(lm)=dwdt%expl(lm,1,1)+hdif_V(lm)*visc(1)*dL*or2(1)* (            & 
         !&                        (dLvisc(1)+two*third*beta(1)) *   ddw(lm,1)-   &
         !&          (two*(third*dbeta(1)-ddLvisc(1))+third*dLvisc(1)*beta(1)+    &
         !&           two*or1(1)*(dLvisc(1)-third*beta(1)))*          dw(lm,1)-   &
         !&          (four*third*(ddbeta(1)+ddLVisc(1)*beta(1)+dLvisc(1)*dbeta(1) &
         !&           +or1(1)*(three*(ddLvisc(1)-dLvisc(1)*or1(1))+dbeta(1)-      &
         !&           beta(1)*or1(1)))+dL*or2(1)*(dLvisc(1)+two*third*beta(1)))*  &
         !&                                                            w(lm,1) )
         !if ( l_heat ) top_bc(lm)=top_bc(lm)+BuoFac*rho0(1)*rgrav(1)*s(lm,1)
         !if ( l_chemical_conv ) top_bc(lm)=top_bc(lm)+ChemFac*rho0(1)*rgrav(1)* &
         !&                                 xi(lm,1)
      end do


      !$omp single
      call solve_counter%start_count()
      !$omp end single
      !PERFON('upWP_ssol')
      !$OMP SINGLE
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)

         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,iChunk,nChunks,size_of_last_chunk,threadid) &
         !$OMP shared(dwold,nLMB,nLMBs2,rhs0)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)

         if ( l1 > 0 .and. .not. lPmat(l1) ) then
            call get_pMat(l1, pMat(nLMB2))
            lPmat(l1) = .true.
         end if

         do iChunk=1,nChunks
            !$OMP TASK if (nChunks>1) default(shared) &
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
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 /= 0 ) then
                  lmB=lmB+1
                  rhs0(1,2*lmB-1,threadid)      = real(top_bc(lm1))
                  rhs0(1,2*lmB,threadid)        =aimag(top_bc(lm1))
                  rhs0(n_r_max,2*lmB-1,threadid)= real(bot_bc(lm1))
                  rhs0(n_r_max,2*lmB,threadid)  =aimag(bot_bc(lm1))
                  do nR=2,n_r_max-1
                     rhs0(nR,2*lmB-1,threadid)= real(work_LMloc(lm1,nR))
                     rhs0(nR,2*lmB,threadid)  =aimag(work_LMloc(lm1,nR))
                  end do
               end if
            end do

            if ( lmB > lmB0 ) then
               call pMat(nLMB2)%solve(rhs0(:,2*(lmB0+1)-1:2*lmB,threadid), &
                                       2*(lmB-lmB0))
            end if

            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 /= 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        p(lm1,n_r_out)  =cmplx(rhs0(n_r_out,2*lmB-1,threadid), &
                        &                      rhs0(n_r_out,2*lmB,threadid),cp)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        p(lm1,n_r_out)  = cmplx(rhs0(n_r_out,2*lmB-1,threadid),&
                        &                       0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do   ! end of loop over l1 subblocks
      !$OMP END SINGLE
      !$omp single
      call solve_counter%stop_count()
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            p(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !-- Construct the radial derivative of pressure
      call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
           &       n_r_max, rscheme_oc, l_dct_in=.false.)
      call rscheme_oc%costf1(p, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1)
      !$omp barrier

      !-- Now finally correct the implicit term required in pol/pre equations
      if ( tscheme%l_imp_calc_rhs(1) ) then
         !$omp do private(nR,lm,l1,dL)
         do nR=1,n_r_max
            do lm=lmStart_00,ulm
               l1 = lm2l(lm)
               dL = real(l1*(l1+1),cp)
               dwdt%impl(lm,nR,1)=dwdt%impl(lm,nR,1)-dp(lm,nR)+beta(nR)*p(lm,nR)
               dpdt%impl(lm,nR,1)=dpdt%impl(lm,nR,1)+dL*or2(nR)*p(lm,nR)
            end do
         end do
         !$omp end do
      end if

      !$omp end parallel

   end subroutine updateP
!------------------------------------------------------------------------------
   subroutine finish_exp_pol(dVxVhLM, dw_exp_last)

      !-- Input variables
      complex(cp), intent(inout) :: dVxVhLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dw_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      integer :: n_r, start_lm, stop_lm

      !$omp parallel default(shared) private(start_lm,stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVxVhLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
           &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do
      do n_r=1,n_r_max
         dw_exp_last(:,n_r)= dw_exp_last(:,n_r)+or2(n_r)*work_LMloc(:,n_r)
      end do
      !$omp end do
      !$omp end parallel

   end subroutine finish_exp_pol
!------------------------------------------------------------------------------
   subroutine get_pol_rhs_imp(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, tscheme,     &
              &               istage, l_calc_lin, lPressNext, lRmsNext, dp_expl, &
              &               l_in_cheb_space)

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

      !-- Output variable
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dp(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)

      !-- Local variables 
      logical :: l_in_cheb
      integer :: n_r_top, n_r_bot, l1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:),lm2m(:)
      real(cp) :: dL

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
      lmStart_00 =max(2,llm)

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
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
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      if ( istage == 1 ) then
         if ( l_double_curl ) then
            !$omp do private(n_r,lm,l1,dL)
            do n_r=2,n_r_max-1
               do lm=lmStart_00,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dwdt%old(lm,n_r,istage)=dL*or2(n_r)* ( -orho1(n_r)*(  &
                  &                   ddw(lm,n_r)-beta(n_r)*dw(lm,n_r)- &
                  &                             dL*or2(n_r)* w(lm,n_r) ) )
               end do
            end do
            !$omp end do
         else
            !$omp do private(n_r,lm,l1,dL) 
            do n_r=2,n_r_max-1
               do lm=lmStart_00,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dwdt%old(lm,n_r,istage)= dL*or2(n_r)*w(lm,n_r)
                  dpdt%old(lm,n_r,istage)=-dL*or2(n_r)*dw(lm,n_r)
               end do
            end do
            !$omp end do
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

         !PERFON('upWP_ex')
         !-- Calculate explicit time step part:
         if ( l_double_curl ) then

            if ( lPressNext ) then
               n_r_top=n_r_cmb
               n_r_bot=n_r_icb
            end if

            !$omp do private(n_r,lm,l1,Dif,Buo,dL)
            do n_r=n_r_top,n_r_bot
               do lm=lmStart_00,ulm
                  l1=lm2l(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif(lm)=-hdif_V(lm)*dL*or2(n_r)*visc(n_r)*orho1(n_r)*      (      &
                  &                                                  ddddw(lm,n_r)  &
                  &            +two*( dLvisc(n_r)-beta(n_r) ) * work_LMloc(lm,n_r)  &
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
                  &          beta(n_r)-three*or1(n_r))+dL*or2(n_r) ) *   w(lm,n_r) )

                  Buo(lm) = zero
                  if ( l_heat ) Buo(lm) = BuoFac*dL*or2(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*dL*or2(n_r)*&
                  &                                rgrav(n_r)*xi(lm,n_r)

                  dwdt%impl(lm,n_r,istage)=Dif(lm)+Buo(lm)

                  if ( l1 /= 0 .and. lPressNext .and. &
                  &    tscheme%istage==tscheme%nstages) then
                     ! In the double curl formulation, we can estimate the pressure
                     ! if required.
                     p(lm,n_r)=-r(n_r)*r(n_r)/dL*                 dp_expl(lm,n_r)  &
                     &            -one/tscheme%dt(1)*(dw(lm,n_r)-dwold(lm,n_r))+   &
                     &              hdif_V(lm)*visc(n_r)* ( work_LMloc(lm,n_r)     &
                     &                       - (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                     &            - ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                     &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                     &                                              ) * dw(lm,n_r) &
                     &             + dL*or2(n_r)*(two*or1(n_r)+two*third*beta(n_r) &
                     &                     +dLvisc(n_r) )   *            w(lm,n_r) )
                  end if

                  if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                     !-- In case RMS force balance is required, one needs to also
                     !-- compute the classical diffusivity that is used in the non
                     !-- double-curl version
                     Dif(lm) = hdif_V(lm)*dL*or2(n_r)*visc(n_r) *  ( ddw(lm,n_r)   &
                     &        +(two*dLvisc(n_r)-third*beta(n_r))*     dw(lm,n_r)   &
                     &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*  &
                     &           beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)))&
                     &                                         *       w(lm,n_r) )
                  end if
               end do
               if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                  call hInt2Pol(Dif,llm,ulm,n_r,lmStart_00,ulm, &
                       &        DifPolLMr(llm:ulm,n_r),         &
                       &        DifPol2hInt(:,n_r),lo_map)
               end if
            end do
            !$omp end do

         else

            !$omp do private(n_r,lm,l1,Dif,Buo,Pre,dL)
            do n_r=n_r_top,n_r_bot
               do lm=lmStart_00,ulm
                  l1=lm2l(lm)
                  dL=real(l1*(l1+1),cp)

                  Dif(lm) = hdif_V(lm)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)  & 
                  &        +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
                  &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
                  &          beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
                  &                                                   w(lm,n_r)  )
                  Pre(lm) = -dp(lm,n_r)+beta(n_r)*p(lm,n_r)
                  Buo(lm) = zero
                  if ( l_heat )  Buo(lm) = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
                  if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*rho0(n_r)* &
                  &                                rgrav(n_r)*xi(lm,n_r)
                  dwdt%impl(lm,n_r,istage)=Pre(lm)+Dif(lm)+Buo(lm)
                  dpdt%impl(lm,n_r,istage)=               dL*or2(n_r)*p(lm,n_r) &
                  &            + hdif_V(lm)*visc(n_r)*dL*or2(n_r)               &
                  &                                     * ( -work_LMloc(lm,n_r) &
                  &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
                  &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
                  &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
                  &                                           ) *    dw(lm,n_r) &
                  &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
                  &                     +dLvisc(n_r) )   *           w(lm,n_r)  )
               end do
               if ( lRmsNext .and. tscheme%istage==tscheme%nstages ) then
                  call hInt2Pol(Dif,llm,ulm,n_r,lmStart_00,ulm, &
                       &        DifPolLMr(llm:ulm,n_r),         &
                       &        DifPol2hInt(:,n_r),lo_map)
               end if
            end do
            !$omp end do

         end if

      end if

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext .and. l_double_curl ) then
         call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, rscheme_oc)
         !$omp barrier
      end if

      !$omp end parallel

   end subroutine get_pol_rhs_imp
!------------------------------------------------------------------------------
   subroutine assemble_pol(s, xi, w, dw, ddw, dwdt, dpdt, tscheme, lRmsNext)

      !-- Input variables
      complex(cp),         intent(in) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: xi(llm:ulm,n_r_max)
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lRmsNext

      !-- Output variable
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      complex(cp),       intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ddw(llm:ulm,n_r_max)

      !-- Local variables 
      real(cp) :: fac_top, fac_bot
      integer :: n_r_top, n_r_bot, l1, lmStart_00
      integer :: n_r, lm, start_lm, stop_lm
      integer, pointer :: lm2l(:)
      real(cp) :: dL

      lm2l(1:lm_max) => lo_map%lm2l
      lmStart_00 =max(2,llm)

      call tscheme%assemble_imex(work_LMloc, dwdt, llm, ulm, n_r_max)
      call tscheme%assemble_imex(ddw, dpdt, llm, ulm, n_r_max) ! Use ddw as a work array

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !-- Now get the poloidal from the assembly
      !$omp do private(n_r,lm,l1,dL)
      do n_r=2,n_r_max-1
         do lm=lmStart_00,ulm
            l1 = lm2l(lm)
            dL = real(l1*(l1+1),cp)
            w(lm,n_r) = r(n_r)*r(n_r)/dL*work_LMloc(lm,n_r)
            dw(lm,n_r)=-r(n_r)*r(n_r)/dL*ddw(lm,n_r)
         end do
      end do
      !$omp end do

      !-- Non-penetration: u_r=0 -> w_lm=0 on both boundaries
      do lm=lmStart_00,ulm
         w(lm,1)      =zero
         w(lm,n_r_max)=zero
      end do

      !-- Other boundary condition: stress-free or rigid
      if ( l_full_sphere ) then
         if ( ktopv == 1 ) then ! Stress-free
            fac_top=-two*or1(1)-beta(1)
            do lm=lmStart_00,ulm
               l1 = lm2l(lm)
               if ( l1 == 1 ) then
                  call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
               else
                  call rscheme_oc%robin_bc(one, fac_top, zero, one, 0.0_cp, zero, dw(lm,:))
               end if
            end do
         else
            do lm=lmStart_00,ulm
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
            do lm=lmStart_00,ulm
               dw(lm,1)      =zero
               dw(lm,n_r_max)=zero
            end do
         else if ( ktopv /= 1 .and. kbotv == 1 ) then ! Rigid top/Stress-free bottom
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(0.0_cp, one, zero, one, fac_bot, zero, dw(lm,:))
            end do
         else if ( ktopv == 1 .and. kbotv /= 1 ) then ! Rigid bottom/Stress-free top
            fac_top=-two*or1(1)-beta(1)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(one, fac_top, zero, 0.0_cp, one, zero, dw(lm,:))
            end do
         else if ( ktopv == 1 .and. kbotv == 1 ) then ! Stress-free at both boundaries
            fac_bot=-two*or1(n_r_max)-beta(n_r_max)
            fac_top=-two*or1(1)-beta(1)
            do lm=lmStart_00,ulm
               call rscheme_oc%robin_bc(one, fac_top, zero, one, fac_bot, zero, dw(lm,:))
            end do
         end if
      end if

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr( dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1, &
           &         stop_lm-llm+1, n_r_max, rscheme_oc)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count()
      !$omp end single

      !$omp do private(n_r,lm,l1,dL) 
      do n_r=2,n_r_max-1
         do lm=lmStart_00,ulm
            l1 = lm2l(lm)
            dL = real(l1*(l1+1),cp)
            dwdt%old(lm,n_r,1)= dL*or2(n_r)*w(lm,n_r)
            dpdt%old(lm,n_r,1)=-dL*or2(n_r)*dw(lm,n_r)
         end do
      end do
      !$omp end do

      if ( tscheme%l_imp_calc_rhs(1) .or. lRmsNext ) then
         if ( lRmsNext ) then
            n_r_top=n_r_cmb
            n_r_bot=n_r_icb
         else
            n_r_top=n_r_cmb+1
            n_r_bot=n_r_icb-1
         end if

         !$omp do private(n_r,lm,l1,Dif,Buo,Pre,dL)
         do n_r=n_r_top,n_r_bot
            do lm=lmStart_00,ulm
               l1=lm2l(lm)
               dL=real(l1*(l1+1),cp)

               Dif(lm) = hdif_V(lm)*dL*or2(n_r)*visc(n_r)*(       ddw(lm,n_r)  & 
               &        +(two*dLvisc(n_r)-third*beta(n_r))*        dw(lm,n_r)  &
               &        -( dL*or2(n_r)+four*third*( dbeta(n_r)+dLvisc(n_r)*    &
               &          beta(n_r)+(three*dLvisc(n_r)+beta(n_r))*or1(n_r)) )* &
               &                                                   w(lm,n_r)  )
               Buo(lm) = zero
               if ( l_heat )  Buo(lm) = BuoFac*rho0(n_r)*rgrav(n_r)*s(lm,n_r)
               if ( l_chemical_conv ) Buo(lm) = Buo(lm)+ChemFac*rho0(n_r)* &
               &                                rgrav(n_r)*xi(lm,n_r)
               dwdt%impl(lm,n_r,1)=Dif(lm)+Buo(lm)
               dpdt%impl(lm,n_r,1)=hdif_V(lm)*visc(n_r)*dL*or2(n_r)*         & 
               &                                       ( -work_LMloc(lm,n_r) &
               &                       + (beta(n_r)-dLvisc(n_r))*ddw(lm,n_r) &
               &            + ( dL*or2(n_r)+dLvisc(n_r)*beta(n_r)+dbeta(n_r) &
               &                  + two*(dLvisc(n_r)+beta(n_r))*or1(n_r)     &
               &                                           ) *    dw(lm,n_r) &
               &            - dL*or2(n_r)* ( two*or1(n_r)+two*third*beta(n_r)&
               &                     +dLvisc(n_r) )   *            w(lm,n_r) )
            end do

            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,n_r,lmStart_00,ulm, &
                    &        DifPolLMr(llm:ulm,n_r),         &
                    &        DifPol2hInt(:,n_r),lo_map)
            end if
         end do
         !$omp end do
      end if

      !$omp end parallel

   end subroutine assemble_pol
!------------------------------------------------------------------------------
   subroutine get_wpMat(tscheme,l,hdif,wpMat,wpMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables:
      class(type_realmat), intent(inout) :: wpMat
      real(cp), intent(out) :: wpMat_fac(2*n_r_max,2)

      !-- local variables:
      integer :: nR,nR_out,nR_p,nR_out_p
      integer :: info
      real(cp) :: dLh

#ifdef MATRIX_CHECK
      integer ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
      real(cp) :: work(8*n_r_max),anorm,linesum,rcond
      real(cp) :: temp_wpMat(2*n_r_max,2*n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
      logical :: first_run=.true.
#endif

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
            ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
            nR_p=nR+n_r_max
            wpMat%dat(nR,nR_out)= rscheme_oc%rnorm *  (                     &
            &               dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)          &
            &         - tscheme%wimp_lin(1)*hdif*visc(nR)*dLh*or2(nR) * (   &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &        +(two*dLvisc(nR)-third*beta(nR))*                      &
            &                               rscheme_oc%drMat(nR,nR_out)     &
            &        -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)        &
            &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
            &          )                    *rscheme_oc%rMat(nR,nR_out)  )  )

            wpMat%dat(nR,nR_out_p)= rscheme_oc%rnorm*tscheme%wimp_lin(1)*(  &
            &                            rscheme_oc%drMat(nR,nR_out)        &
            &                  -beta(nR)* rscheme_oc%rMat(nR,nR_out))
            ! the following part gives sometimes very large
            ! matrix entries
            wpMat%dat(nR_p,nR_out)=rscheme_oc%rnorm * (                       &
            &                  -dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)       &
            &         -tscheme%wimp_lin(1)*hdif*visc(nR)*dLh*or2(nR)      *(  &
            &                                 - rscheme_oc%d3rMat(nR,nR_out)  &
            &          +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out)  &
            &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
            &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
            &                                    rscheme_oc%drMat(nR,nR_out)  &
            &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
            &          +two*third*beta(nR)   )*   rscheme_oc%rMat(nR,nR_out) ) )

            wpMat%dat(nR_p,nR_out_p)= -rscheme_oc%rnorm*tscheme%wimp_lin(1)*  &
            &                          dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)
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
!------------------------------------------------------------------------------
   subroutine get_pMat(l, pMat)
      !
      !  Purpose of this subroutine is to contruct the matrix needed
      !  for the pressure Poisson solver.
      !

      !-- Input variables:
      integer, intent(in) :: l       ! degree

      !-- Output variables:
      class(type_realmat), intent(inout) :: pMat

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh =real(l*(l+1),kind=cp)

      !----- Boundary conditions: ensure pressure gradient
      dat(1,:)=rscheme_oc%rnorm*( rscheme_oc%drMat(1,:) &
      &                        -beta(1)*rscheme_oc%rMat(1,:) )
      dat(n_r_max,:)=rscheme_oc%rnorm*(                   &
      &                       rscheme_oc%drMat(n_r_max,:) &
      &         -beta(n_r_max)*rscheme_oc%rMat(n_r_max,:) )

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Other points:
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm *  (                           &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &       +(two*or1(nR)-beta(nR))*rscheme_oc%drMat(nR,nR_out)     &
            &       -(dLh*or2(nR)+two*beta(nR)*or1(nR)+dbeta(nR))*          &
            &                                rscheme_oc%rMat(nR,nR_out) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)          =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)    =rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

      !-- Array copy
      call pMat%set_data(dat)

      !-- Lu factorisation
      call pMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix pMat!')

   end subroutine get_pMat
!-----------------------------------------------------------------------------
   subroutine get_wMat(tscheme,l,hdif,wMat,wMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix
      !  wpmat  for the NS equation. This matrix corresponds here to the
      !  radial component of the double-curl of the Navier-Stokes equation.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables:
      class(type_realmat), intent(inout) :: wMat
      real(cp), intent(out) :: wMat_fac(n_r_max,2)

      !-- local variables:
      integer :: nR, nR_out
      integer :: info
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

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
            &          -dLh*or2(nR)*orho1(nR)*(                             &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out) -   &
            &                    dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )   &
            &  + tscheme%wimp_lin(1)*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (&
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

      !---- LU decomposition:
      call pMat%prepare(info)
      if ( info /= 0 ) call abortRun('! Singular matrix p0Mat!')

   end subroutine get_p0Mat
!-----------------------------------------------------------------------------
end module updateWP_mod
