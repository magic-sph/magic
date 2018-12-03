#include "perflib_preproc.cpp"
module updateWP_mod
   
   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, l_max
   use radial_data, only: n_r_cmb,n_r_icb
   use radial_functions, only: or1, or2, rho0, rgrav, visc, dLvisc, r, &
       &                       alpha0, temp0, beta, dbeta, ogrun,      &
       &                       rscheme_oc, ddLvisc, ddbeta, orho1
   use physical_parameters, only: kbotv, ktopv, ra, BuoFac, ChemFac,    &
       &                          ViscHeatFac, ThExpNb, ktopp
   use num_param, only: alpha
   use blocking, only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map, &
       &               lmStartB,lmStopB
   use horizontal_data, only: hdif_V, dLh
   use logic, only: l_update_v, l_chemical_conv, l_RMS, l_double_curl, &
       &            l_fluxProfs
   use RMS, only: DifPol2hInt, dtVPolLMr, dtVPol2hInt, DifPolLMr
   use algebra, only: prepare_mat, solve_mat
   use LMLoop_data, only: llm, ulm
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize
   use RMS_helpers, only:  hInt2Pol
   use radial_der, only: get_dddr, get_ddr, get_dr
   use integration, only: rInt_R
   use fields, only: work_LMloc
   use constants, only: zero, one, two, three, four, third, half
   use useful, only: abortRun

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: workB(:,:), ddddw(:,:)
   complex(cp), allocatable :: dwold(:,:)
   real(cp), allocatable :: work(:)
   complex(cp), allocatable :: Dif(:),Pre(:),Buo(:),dtV(:)
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: wpMat(:,:,:), wpMat_fac(:,:,:)
   real(cp), allocatable :: p0Mat(:,:)
   integer, allocatable :: wpPivot(:,:), p0Pivot(:)
   logical, public, allocatable :: lWPmat(:)
   integer :: maxThreads

   public :: initialize_updateWP, finalize_updateWP, updateWP

contains

   subroutine initialize_updateWP

      allocate( wpMat(2*n_r_max,2*n_r_max,l_max), p0Mat(n_r_max,n_r_max) )
      allocate( wpMat_fac(2*n_r_max,2,l_max) )
      allocate( wpPivot(2*n_r_max,l_max), p0Pivot(n_r_max) )
      allocate( lWPmat(0:l_max) )
      bytes_allocated=bytes_allocated+((4*n_r_max+4)*(l_max)+n_r_max)*n_r_max* &
      &               SIZEOF_DEF_REAL+(2*n_r_max*l_max+n_r_max)*SIZEOF_INTEGER+&
      &               (l_max+1)*SIZEOF_LOGICAL

      if ( l_RMS ) then
         allocate( workB(llm:ulm,n_r_max) )
         bytes_allocated = bytes_allocated+(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      end if

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
      allocate( dtV(llm:ulm) )
      bytes_allocated = bytes_allocated+4*(ulm-llm+1)*SIZEOF_DEF_COMPLEX

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate( rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
      &               lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX


   end subroutine initialize_updateWP
!-----------------------------------------------------------------------------
   subroutine finalize_updateWP

      deallocate( wpMat, wpMat_fac, wpPivot, lWPmat )
      deallocate( p0Mat, p0Pivot )
      deallocate( rhs1, work )
      deallocate( Dif, Pre, Buo, dtV )
      if ( l_RMS ) deallocate( workB )
      if ( l_double_curl ) then
         deallocate( ddddw )
         if ( l_RMS .or. l_FluxProfs ) then
            deallocate( dwold )
         end if
      end if

   end subroutine finalize_updateWP
!-----------------------------------------------------------------------------
   subroutine updateWP(w,dw,ddw,dVxVhLM,dwdt,dwdtLast,p,dp,dpdt,dpdtLast,s,xi, &
        &              w1,coex,dt,nLMB,lRmsNext,lPressNext)
      !
      !  updates the poloidal velocity potential w, the pressure p,  and
      !  their derivatives
      !  adds explicit part to time derivatives of w and p
      !

      !-- Input/output of scalar fields:
      real(cp),    intent(in) :: w1       ! weight for time step !
      real(cp),    intent(in) :: coex     ! factor depending on alpha
      real(cp),    intent(in) :: dt       ! time step
      integer,     intent(in) :: nLMB     ! block number
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
      integer :: lmStart,lmStop ! max and min number of orders m
      integer :: lmStart_00     ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR             ! counts radial grid points
      integer :: n_r_out         ! counts cheb modes
      real(cp) :: rhs(n_r_max)  ! real RHS for l=m=0
      integer :: n_r_top, n_r_bot

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: iThread,start_lm,stop_lm,all_lms,per_thread,nThreads
      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid

      if ( .not. l_update_v ) return

      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !allocate(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,nLMBs2(nLMB)))

      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)
      lmStart_00  =max(2,lmStart)

      w2  =one-w1
      O_dt=one/dt

      if ( l_double_curl ) then
         !PERFON('upS_fin')
         !$OMP PARALLEL  &
         !$OMP private(iThread,start_lm,stop_lm,nR,lm) &
         !$OMP shared(all_lms,per_thread) &
         !$OMP shared(dVxVhLM,rscheme_oc,dwdt) &
         !$OMP shared(or2,lmStart,lmStop) &
         !$OMP shared(n_r_max,work_LMloc,nThreads,llm,ulm)
         !$OMP SINGLE
#ifdef WITHOMP
         nThreads=omp_get_num_threads()
#else
         nThreads=1
#endif
         !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
         all_lms=lmStop-lmStart+1
         per_thread=all_lms/nThreads
         !$OMP END SINGLE
         !$OMP BARRIER
         !$OMP DO
         do iThread=0,nThreads-1
            start_lm=lmStart+iThread*per_thread
            stop_lm = start_lm+per_thread-1
            if (iThread == nThreads-1) stop_lm=lmStop

            !--- Finish calculation of dsdt:
            call get_dr( dVxVhLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
                 &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
         end do
         !$OMP end do

         !$OMP DO
         do nR=1,n_r_max
            do lm=lmStart,lmStop
               dwdt(lm,nR)= dwdt(lm,nR)+or2(nR)*work_LMloc(lm,nR)
            end do
         end do
         !$OMP end do
         !$OMP END PARALLEL
         !PERFOFF
      end if

      !PERFON('upWP_ssol')
      !$OMP PARALLEL default(shared) &
      !$OMP private(nLMB2,lm,lm1,l1,m1,lmB)
      !write(*,"(I3,A)") omp_get_thread_num(),": before SINGLE"
      !$OMP SINGLE
      ! each of the nLMBs2(nLMB) subblocks have one l value
      do nLMB2=1,nLMBs2(nLMB)
         !write(*,"(2(A,I3))") "Constructing next task for ",nLMB2,"/",nLMBs2(nLMB)

         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,iChunk,nChunks,size_of_last_chunk,threadid) &
         !$OMP shared(workB,dwold,nLMB,nLMBs2,rhs1)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk=chunksize+(sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 == 0 ) then
            if ( .not. lWPmat(l1) ) then
               call get_p0Mat(p0Mat,p0Pivot)
               lWPmat(l1)=.true.
            end if
         else
            if ( .not. lWPmat(l1) ) then
               !PERFON('upWP_mat')
               if ( l_double_curl ) then
                  call get_wMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                       &        wpMat(:,:,l1),wpPivot(:,l1),wpMat_fac(:,:,l1))
               else
                  call get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                       &         wpMat(:,:,l1),wpPivot(:,l1),wpMat_fac(:,:,l1))
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
                     do nR=1,n_r_max
                        work(nR)=ThExpNb*alpha0(nR)*temp0(nR)*rho0(nR)*r(nR)*&
                        &        r(nR)*real(s(st_map%lm2(0,0),nR))
                     end do
                     rhs(1)=rInt_R(work,r,rscheme_oc)
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

                  call solve_mat(p0Mat,n_r_max,n_r_max,p0Pivot,rhs)

               else ! l1 /= 0
                  lmB=lmB+1
                  rhs1(1,lmB,threadid)        =0.0_cp
                  rhs1(n_r_max,lmB,threadid)  =0.0_cp
                  rhs1(n_r_max+1,lmB,threadid)=0.0_cp
                  rhs1(2*n_r_max,lmB,threadid)=0.0_cp
                  if ( l_double_curl ) then
                     if ( l_chemical_conv ) then
                        do nR=2,n_r_max-1
                           rhs1(nR,lmB,threadid)=dLh(st_map%lm2(l1,m1))*or2(nR)* (   &
                           &                     -orho1(nR)*O_dt*(    ddw(lm1,nR)    &
                           &                     -beta(nR)*dw(lm1,nR)-               &
                           &                     dLh(st_map%lm2(l1,m1))*or2(nR)*     &
                           &                                w(lm1,nR) ) +            &
                           &                     alpha*BuoFac *rgrav(nR)* s(lm1,nR)+ &
                           &                     alpha*ChemFac*rgrav(nR)*xi(lm1,nR) )&
                           &                     +w1*dwdt(lm1,nR) +                  &
                           &                     w2*dwdtLast(lm1,nR)
                           rhs1(nR+n_r_max,lmB,threadid)=0.0_cp
                        end do
                     else
                        do nR=2,n_r_max-1
                           rhs1(nR,lmB,threadid)=dLh(st_map%lm2(l1,m1))*or2(nR)* (   &
                           &                     -orho1(nR)*O_dt*(    ddw(lm1,nR)    &
                           &                     -beta(nR)*dw(lm1,nR)-               &
                           &                     dLh(st_map%lm2(l1,m1))*or2(nR)*     &
                           &                                w(lm1,nR) ) +            &
                           &                     alpha*BuoFac *rgrav(nR)* s(lm1,nR) )&
                           &                     +w1*dwdt(lm1,nR) +                  &
                           &                     w2*dwdtLast(lm1,nR)
                           rhs1(nR+n_r_max,lmB,threadid)=0.0_cp
                        end do
                     end if
                  else
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
                  do nR=1,2*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,1,l1)
                  end do
               end do
               call solve_mat(wpMat(:,:,l1),2*n_r_max,2*n_r_max,    &
                    &         wpPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,2*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,2,l1)
                  end do
               end do
            end if
            !PERFOFF

            if ( lRmsNext ) then ! Store old w
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     workB(lm1,nR)=w(lm1,nR)
                  end do
               end do
            end if

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
                           ddw(lm1,n_r_out)=rhs1(n_r_max+n_r_out,lmB,threadid)
                        end do
                     else
                        do n_r_out=1,rscheme_oc%n_max
                           w(lm1,n_r_out)  = cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                           &                    0.0_cp,kind=cp)
                           ddw(lm1,n_r_out)= cmplx(real(rhs1(n_r_max+n_r_out,lmB, &
                           &                    threadid)),0.0_cp,kind=cp)
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
      !$OMP END PARALLEL
      !PERFOFF
      !write(*,"(A,I3,4ES22.12)") "w,p after: ",nLMB,get_global_SUM(w),get_global_SUM(p)

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=lmStart,lmStop
            w(lm1,n_r_out)=zero
            p(lm1,n_r_out)=zero
            if ( l_double_curl ) then
               ddw(lm1,n_r_out)=zero
            end if
         end do
      end do


      !PERFON('upWP_drv')
      all_lms=lmStop-lmStart+1
#ifdef WITHOMP
      if (all_lms < omp_get_max_threads()) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL  &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStart_00,lmStop) &
      !$OMP shared(w,dw,ddw,p,dp) &
      !$OMP shared(rscheme_oc,n_r_max,nThreads,ddddw,work_LMloc,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads = 1
#endif
      !$OMP END SINGLE
      !$OMP BARRIER
      per_thread=all_lms/nThreads
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop
         !write(*,"(2(A,I3),2(A,I5))") "iThread=",iThread," on thread ", &
         !     & omp_get_thread_num()," lm = ",start_lm,":",stop_lm

         !-- Transform to radial space and get radial derivatives
         !   using dwdtLast, dpdtLast as work arrays:

         if ( l_double_curl ) then
            call get_dr( w, dw, ulm-llm+1, start_lm-llm+1,  &
                 &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
            call get_ddr( ddw, work_LMloc, ddddw, ulm-llm+1,                  &
                 &        start_lm-llm+1, stop_lm-llm+1, n_r_max, rscheme_oc, &
                 &        l_dct_in=.false. )
            call rscheme_oc%costf1(ddw,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
         else
            call get_dddr( w, dw, ddw, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
                 &         stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
            call get_dr( p, dp, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
                 &       n_r_max, rscheme_oc, l_dct_in=.false. )
         end if
         call rscheme_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
         call rscheme_oc%costf1(p,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(omp_get_max_threads())
#endif
      !PERFOFF

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

         do nR=n_r_top,n_r_bot
            do lm1=lmStart_00,lmStop
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
                  dtV(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                  &             ( w(lm1,nR)-workB(lm1,nR) )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr(llm:,nR), &
                    &        DifPol2hInt(:,nR,1),lo_map)
               call hInt2Pol(dtV,llm,ulm,nR,lmStart_00,lmStop, &
                    &        dtVPolLMr(llm:,nR),dtVPol2hInt(:,nR,1),lo_map)
            end if
         end do

      else

         do nR=n_r_top,n_r_bot
            do lm1=lmStart_00,lmStop
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
               if ( lRmsNext ) then
                  dtV(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                  &             ( w(lm1,nR)-workB(lm1,nR) )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr(llm:,nR), &
                    &        DifPol2hInt(:,nR,1),lo_map)
               call hInt2Pol(dtV,llm,ulm,nR,lmStart_00,lmStop, &
                    &        dtVPolLMr(llm:,nR),dtVPol2hInt(:,nR,1),lo_map)
            end if
         end do

      end if
      !PERFOFF

      ! In case pressure is needed in the double curl formulation
      ! we also have to compute the radial derivative of p
      if ( lPressNext .and. l_double_curl ) then
         !PERFON('upWP_drv')
         all_lms=lmStop-lmStart+1
#ifdef WITHOMP
         if (all_lms < omp_get_max_threads()) then
            call omp_set_num_threads(all_lms)
         end if
#endif
         !$OMP PARALLEL  &
         !$OMP private(iThread,start_lm,stop_lm) &
         !$OMP shared(all_lms,per_thread,lmStop) &
         !$OMP shared(p,dp,rscheme_oc,n_r_max,nThreads,llm,ulm)
         !$OMP SINGLE
#ifdef WITHOMP
         nThreads=omp_get_num_threads()
#else
         nThreads = 1
#endif
         !$OMP END SINGLE
         !$OMP BARRIER
         per_thread=all_lms/nThreads
         !$OMP DO
         do iThread=0,nThreads-1
            start_lm=lmStart+iThread*per_thread
            stop_lm = start_lm+per_thread-1
            if (iThread == nThreads-1) stop_lm=lmStop

            call get_dr( p, dp, ulm-llm+1, start_lm-llm+1,  &
                 &         stop_lm-llm+1, n_r_max, rscheme_oc)
         end do
         !$OMP end do
         !$OMP END PARALLEL
#ifdef WITHOMP
         call omp_set_num_threads(omp_get_max_threads())
#endif
      end if

   end subroutine updateWP
!------------------------------------------------------------------------------
   subroutine get_wpMat(dt,l,hdif,wpMat,wpPivot,wpMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix  
      !  wpmat  for the NS equation.                                    
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: wpMat(2*n_r_max,2*n_r_max)
      real(cp), intent(out) :: wpMat_fac(2*n_r_max,2)
      integer,  intent(out) :: wpPivot(2*n_r_max)

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
    
         wpMat(1,nR_out)        =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         wpMat(1,nR_out_p)      =0.0_cp
         wpMat(n_r_max,nR_out)  =rscheme_oc%rnorm* &
         &                       rscheme_oc%rMat(n_r_max,nR_out)
         wpMat(n_r_max,nR_out_p)=0.0_cp
    
         if ( ktopv == 1 ) then  ! free slip !
            wpMat(n_r_max+1,nR_out)= rscheme_oc%rnorm * (           &
            &                        rscheme_oc%d2rMat(1,nR_out) -  &
            &    (two*or1(1)+beta(1))*rscheme_oc%drMat(1,nR_out) )
         else                    ! no slip, note exception for l=1,m=0
            wpMat(n_r_max+1,nR_out)=rscheme_oc%rnorm* &
            &                       rscheme_oc%drMat(1,nR_out)
         end if
         wpMat(n_r_max+1,nR_out_p)=0.0_cp
    
         if ( kbotv == 1 ) then  ! free slip !
            wpMat(2*n_r_max,nR_out)=rscheme_oc%rnorm * (           &
            &                  rscheme_oc%d2rMat(n_r_max,nR_out) - &
            &      ( two*or1(n_r_max)+beta(n_r_max))*              &
            &                  rscheme_oc%drMat(n_r_max,nR_out))
         else                 ! no slip, note exception for l=1,m=0
            wpMat(2*n_r_max,nR_out)=rscheme_oc%rnorm * &
            &                       rscheme_oc%drMat(n_r_max,nR_out)
         end if
         wpMat(2*n_r_max,nR_out_p)=0.0_cp
    
      end do   !  loop over nR_out
    
      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_p=nR_out+n_r_max
            wpMat(1,nR_out)          =0.0_cp
            wpMat(n_r_max,nR_out)    =0.0_cp
            wpMat(n_r_max+1,nR_out)  =0.0_cp
            wpMat(2*n_r_max,nR_out)  =0.0_cp
            wpMat(1,nR_out_p)        =0.0_cp
            wpMat(n_r_max,nR_out_p)  =0.0_cp
            wpMat(n_r_max+1,nR_out_p)=0.0_cp
            wpMat(2*n_r_max,nR_out_p)=0.0_cp
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
            wpMat(nR,nR_out)= rscheme_oc%rnorm *  (                         &
            &               O_dt*dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out)     &
            &            - alpha*hdif*visc(nR)*dLh*or2(nR) * (              &
            &                              rscheme_oc%d2rMat(nR,nR_out)     &
            &        +(two*dLvisc(nR)-third*beta(nR))*                      &
            &                               rscheme_oc%drMat(nR,nR_out)     &
            &        -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)        &
            &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )     &
            &          )                    *rscheme_oc%rMat(nR,nR_out)  )  )
    
            wpMat(nR,nR_out_p)= rscheme_oc%rnorm*alpha*(             &
            &                            rscheme_oc%drMat(nR,nR_out) &
            &                  -beta(nR)* rscheme_oc%rMat(nR,nR_out))
            ! the following part gives sometimes very large 
            ! matrix entries
            wpMat(nR_p,nR_out)=rscheme_oc%rnorm * (                           &
            &                  -O_dt*dLh*or2(nR)*rscheme_oc%drMat(nR,nR_out)  &
            &         -alpha*hdif*visc(nR)*dLh*or2(nR)      *(                &
            &                                 - rscheme_oc%d3rMat(nR,nR_out)  &
            &          +( beta(nR)-dLvisc(nR) )*rscheme_oc%d2rMat(nR,nR_out)  &
            &          +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
            &          +two*(dLvisc(nR)+beta(nR))*or1(nR) )*                  &
            &                                    rscheme_oc%drMat(nR,nR_out)  &
            &          -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
            &          +two*third*beta(nR)   )*   rscheme_oc%rMat(nR,nR_out) ) )
    
            wpMat(nR_p,nR_out_p)= -rscheme_oc%rnorm*alpha*dLh*or2(nR)* &
            &                      rscheme_oc%rMat(nR,nR_out)
         end do
      end do
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         wpMat(nR,1)          =rscheme_oc%boundary_fac*wpMat(nR,1)
         wpMat(nR,n_r_max)    =rscheme_oc%boundary_fac*wpMat(nR,n_r_max)
         wpMat(nR,n_r_max+1)  =rscheme_oc%boundary_fac*wpMat(nR,n_r_max+1)
         wpMat(nR,2*n_r_max)  =rscheme_oc%boundary_fac*wpMat(nR,2*n_r_max)
         wpMat(nR_p,1)        =rscheme_oc%boundary_fac*wpMat(nR_p,1)
         wpMat(nR_p,n_r_max)  =rscheme_oc%boundary_fac*wpMat(nR_p,n_r_max)
         wpMat(nR_p,n_r_max+1)=rscheme_oc%boundary_fac*wpMat(nR_p,n_r_max+1)
         wpMat(nR_p,2*n_r_max)=rscheme_oc%boundary_fac*wpMat(nR_p,2*n_r_max)
      end do
    
      ! compute the linesum of each line
      do nR=1,2*n_r_max
         wpMat_fac(nR,1)=one/maxval(abs(wpMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         wpMat(nR,:) = wpMat(nR,:)*wpMat_fac(nR,1)
      end do
    
      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         wpMat_fac(nR,2)=one/maxval(abs(wpMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         wpMat(:,nR) = wpMat(:,nR)*wpMat_fac(nR,2)
      end do

#ifdef MATRIX_CHECK
      ! copy the wpMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1
      
      do i=1,2*n_r_max
         do j=1,2*n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") wpMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_wpMat=wpMat
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

      call prepare_mat(wpMat,2*n_r_max,2*n_r_max,wpPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix wpMat!')
      end if

   end subroutine get_wpMat
!-----------------------------------------------------------------------------
   subroutine get_wMat(dt,l,hdif,wMat,wPivot,wMat_fac)
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
      real(cp), intent(out) :: wMat(2*n_r_max,2*n_r_max)
      real(cp), intent(out) :: wMat_fac(2*n_r_max,2)
      integer,  intent(out) :: wPivot(2*n_r_max)

      !-- local variables:
      integer :: nR,nR_out,nR_ddw,nR_out_ddw
      integer :: info
      real(cp) :: O_dt,dLh

      O_dt=one/dt
      dLh =real(l*(l+1),kind=cp)
    
      !-- Now mode l>0
    
      !----- Boundary conditions, see above:
      do nR_out=1,rscheme_oc%n_max
         nR_out_ddw=nR_out+n_r_max
    
         wMat(1,nR_out)      =rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         wMat(1,nR_out_ddw)  =0.0_cp
         wMat(n_r_max,nR_out)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR_out)
         wMat(n_r_max,nR_out_ddw)=0.0_cp
    
         if ( ktopv == 1 ) then  ! free slip !
            wMat(n_r_max+1,nR_out)= -rscheme_oc%rnorm *         &
            &                       (two*or1(1)+beta(1))*rscheme_oc%drMat(1,nR_out) 
            wMat(n_r_max+1,nR_out_ddw)= rscheme_oc%rnorm * rscheme_oc%rMat(1,nR_out)
         else                    ! no slip, note exception for l=1,m=0
            wMat(n_r_max+1,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(1,nR_out)
            wMat(n_r_max+1,nR_out_ddw)=0.0_cp
         end if
    
         if ( kbotv == 1 ) then  ! free slip !
            wMat(2*n_r_max,nR_out)= -rscheme_oc%rnorm *              &
            &                      (two*or1(n_r_max)+beta(n_r_max))* &
            &                      rscheme_oc%drMat(n_r_max,nR_out)
            wMat(2*n_r_max,nR_out_ddw)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,nR_out)
         else                 ! no slip, note exception for l=1,m=0
            wMat(2*n_r_max,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,nR_out)
            wMat(2*n_r_max,nR_out_ddw)=0.0_cp
         end if
    
      end do   !  loop over nR_out
    
      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            nR_out_ddw=nR_out+n_r_max
            wMat(1,nR_out)            =0.0_cp
            wMat(n_r_max,nR_out)      =0.0_cp
            wMat(n_r_max+1,nR_out)    =0.0_cp
            wMat(2*n_r_max,nR_out)    =0.0_cp
            wMat(1,nR_out_ddw)        =0.0_cp
            wMat(n_r_max,nR_out_ddw)  =0.0_cp
            wMat(n_r_max+1,nR_out_ddw)=0.0_cp
            wMat(2*n_r_max,nR_out_ddw)=0.0_cp
         end do
      end if
    
      !----- Other points:
      do nR_out=1,n_r_max
         nR_out_ddw=nR_out+n_r_max
         do nR=2,n_r_max-1

            !-- Here instead of solving one single matrix that would require
            !-- the 4th derivative, we instead solve a matrix of size 
            !-- (2*n_r_max, 2*n_r_max) that solve for w and ddw

            nR_ddw=nR+n_r_max
            wMat(nR,nR_out)= rscheme_oc%rnorm *  (                          &
            &             -O_dt*dLh*or2(nR)*orho1(nR)*(                     &
            &                     -beta(nR)*rscheme_oc%drMat(nR,nR_out) -   &
            &                    dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )   &
            &    + alpha*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (            &
            &     ( -ddbeta(nR)-dbeta(nR)*(two*dLvisc(nR)-beta(nR)+         &
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

            wMat(nR,nR_out_ddw)= rscheme_oc%rnorm *  (                      &
            &              -O_dt*dLh*or2(nR)*orho1(nR)*                     &
            &                                 rscheme_oc%rMat(nR,nR_out) +  &
            &      alpha*orho1(nR)*hdif*visc(nR)*dLh*or2(nR) * (            &
            &                               rscheme_oc%d2rMat(nR,nR_out)    &
            &    +two*(dLvisc(nR)-beta(nR))* rscheme_oc%drMat(nR,nR_out)    &
            &    +( ddLvisc(nR)-two*dbeta(nR)+dLvisc(nR)*dLvisc(nR)+        &
            &       beta(nR)*beta(nR)-three*dLvisc(nR)*beta(nR)-            &
            &       two*or1(nR)*(dLvisc(nR)+beta(nR))-two*dLh*or2(nR) ) *   &
            &                                 rscheme_oc%rMat(nR,nR_out) ) )

            wMat(nR_ddw,nR_out)    =-rscheme_oc%rnorm*rscheme_oc%d2rMat(nR,nR_out)
            wMat(nR_ddw,nR_out_ddw)=   rscheme_oc%rnorm*rscheme_oc%rMat(nR,nR_out)
         end do
      end do
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_ddw=nR+n_r_max
         wMat(nR,1)            =rscheme_oc%boundary_fac*wMat(nR,1)
         wMat(nR,n_r_max)      =rscheme_oc%boundary_fac*wMat(nR,n_r_max)
         wMat(nR,n_r_max+1)    =rscheme_oc%boundary_fac*wMat(nR,n_r_max+1)
         wMat(nR,2*n_r_max)    =rscheme_oc%boundary_fac*wMat(nR,2*n_r_max)
         wMat(nR_ddw,1)        =rscheme_oc%boundary_fac*wMat(nR_ddw,1)
         wMat(nR_ddw,n_r_max)  =rscheme_oc%boundary_fac*wMat(nR_ddw,n_r_max)
         wMat(nR_ddw,n_r_max+1)=rscheme_oc%boundary_fac*wMat(nR_ddw,n_r_max+1)
         wMat(nR_ddw,2*n_r_max)=rscheme_oc%boundary_fac*wMat(nR_ddw,2*n_r_max)
      end do
    
      ! compute the linesum of each line
      do nR=1,2*n_r_max
         wMat_fac(nR,1)=one/maxval(abs(wMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         wMat(nR,:) = wMat(nR,:)*wMat_fac(nR,1)
      end do
    
      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         wMat_fac(nR,2)=one/maxval(abs(wMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         wMat(:,nR) = wMat(:,nR)*wMat_fac(nR,2)
      end do

      call prepare_mat(wMat,2*n_r_max,2*n_r_max,wPivot,info)

      if ( info /= 0 ) then
         call abortRun('Singular matrix wMat!')
      end if

   end subroutine get_wMat
!-----------------------------------------------------------------------------
   subroutine get_p0Mat(pMat,pPivot)

      !-- Output variables:
      real(cp), intent(out) :: pMat(n_r_max,n_r_max)
      integer,  intent(out) :: pPivot(n_r_max)

      !-- Local variables:
      integer :: info, nCheb, nR_out, nR, nCheb_in


      do nR_out=1,n_r_max
         do nR=2,n_r_max
            pMat(nR,nR_out)= rscheme_oc%rnorm * ( rscheme_oc%drMat(nR,nR_out)- &
            &                              beta(nR)*rscheme_oc%rMat(nR,nR_out) )
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
               pMat(1,nCheb)=0.0_cp
               do nCheb_in=1,rscheme_oc%n_max
                  if ( mod(nCheb+nCheb_in-2,2)==0 ) then
                     pMat(1,nCheb)=pMat(1,nCheb)+ &
                     &             ( one/(one-real(nCheb_in-nCheb,cp)**2)    + &
                     &               one/(one-real(nCheb_in+nCheb-2,cp)**2) )* &
                     &               work(nCheb_in)*half*rscheme_oc%rnorm
                  end if
               end do
            end do

         else

            !-- In the finite differences case, we restrict the integral boundary
            !-- condition to a trapezoidal rule of integration
            do nR_out=2,rscheme_oc%n_max-1
               pMat(1,nR_out)=half*work(nR_out)*( r(nR_out+1)-r(nR_out-1) )
            end do
            pMat(1,1)=half*work(1)*(r(2)-r(1))

         end if

      else

         do nR_out=1,rscheme_oc%n_max
            pMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         end do

      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros
         do nR_out=rscheme_oc%n_max+1,n_r_max
            pMat(1,nR_out)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         pMat(nR,1)      =rscheme_oc%boundary_fac*pMat(nR,1)
         pMat(nR,n_r_max)=rscheme_oc%boundary_fac*pMat(nR,n_r_max)
      end do

      !---- LU decomposition:
      call prepare_mat(pMat,n_r_max,n_r_max,pPivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular matrix p0Mat!')
      end if

   end subroutine get_p0Mat
!-----------------------------------------------------------------------------
end module updateWP_mod
