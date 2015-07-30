!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
module updateWP_mod
   use omp_lib
   use truncation, only: lm_max, n_cheb_max, n_r_max
   use radial_data, only: n_r_cmb,n_r_icb
   use radial_functions, only: drx,ddrx,dddrx,or1,or2,rho0,agrav,rgrav, &
                             & i_costf_init,d_costf_init,visc,dlvisc,   &
                             & beta,dbeta,cheb,dcheb,d2cheb,d3cheb,     &
                             & cheb_norm
   use physical_parameters, only: kbotv, ktopv, ra
   use num_param, only: alpha
   use blocking, only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map, &
                     & lmStartB,lmStopB
   use horizontal_data, only: hdif_V, dLh
   use logic, only: l_update_v, l_RMStest
   use matrices, only: wpMat, wpPivot, lWPmat, wpMat_fac
   use RMS, only: DifPol2hInt, DifPolAs2hInt, dtVPolLMr, dtVPol2hInt, &
                  dtVPolAs2hInt, DifPolLMr
#ifdef WITH_MKL_LU
   use lapack95, only: getrs, getrf
#else
   use algebra, only: cgeslML, sgefa
#endif
   use LMLoop_data, only:llm,ulm, llm_real,ulm_real
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize
   use RMS_helpers, only:  hInt2Pol

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(kind=8), allocatable :: workA(:,:),workB(:,:)
   complex(kind=8), allocatable :: Dif(:),Pre(:),Buo(:)
   complex(kind=8), allocatable :: rhs1(:,:,:)
   integer :: maxThreads

   public :: initialize_updateWP, updateWP

contains

   subroutine initialize_updateWP

      allocate( workA(llm:ulm,n_r_max) )
      allocate( workB(llm:ulm,n_r_max) )
      allocate( Dif(llm:ulm) )
      allocate( Pre(llm:ulm) )
      allocate( Buo(llm:ulm) )
#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif

      allocate( rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )

   end subroutine initialize_updateWP

   subroutine finalize_updateWP

      deallocate( workA )
      deallocate( workB )
      deallocate( Dif )
      deallocate( Pre )
      deallocate( Buo )
      deallocate( rhs1 )

   end subroutine finalize_updateWP

   subroutine updateWP(w,dw,ddw,dwdt,dwdtLast,p,dp,dpdt,dpdtLast,s, &
        &              w1,coex,dt,nLMB,lRmsNext)
      !-----------------------------------------------------------------------

      !  updates the poloidal velocity potential w, the pressure p,  and
      !  their derivatives
      !  adds explicit part to time derivatives of w and p

      !-------------------------------------------------------------------------

      !-- Input/output of scalar fields:
      real(kind=8),    intent(in) :: w1        ! weight for time step !
      real(kind=8),    intent(in) :: coex      ! factor depending on alpha
      real(kind=8),    intent(in) :: dt        ! time step
      integer,         intent(in) :: nLMB     ! block number
      logical,         intent(in) :: lRmsNext
      complex(kind=8), intent(in) :: dwdt(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: dpdt(llm:ulm,n_r_max)
      complex(kind=8), intent(in) :: s(llm:ulm,n_r_max)

      complex(kind=8), intent(inout) :: w(llm:ulm,n_r_max)
      complex(kind=8), intent(inout) :: dw(llm:ulm,n_r_max)
      complex(kind=8), intent(inout) :: dwdtLast(llm:ulm,n_r_max)
      complex(kind=8), intent(inout) :: p(llm:ulm,n_r_max)
      complex(kind=8), intent(inout) :: dpdtLast(llm:ulm,n_r_max)

      complex(kind=8), intent(out) :: ddw(llm:ulm,n_r_max)
      complex(kind=8), intent(out) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      real(kind=8) :: w2                  ! weight of second time step
      real(kind=8) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lm,lmB         ! position of (l,m) in array
      integer :: lmStart,lmStop ! max and min number of orders m
      integer :: lmStart_real      ! range of lm for real array
      integer :: lmStop_real       !
      integer :: lmStart_00        ! excluding l=0,m=0
      integer :: nLMB2
      integer :: nR                ! counts radial grid points
      integer :: n_cheb             ! counts cheb modes

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
      lmStart_real=2*lmStart_00-1
      lmStop_real =2*lmStop

      w2  =1.D0-w1
      O_dt=1.D0/dt
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
         !$OMP shared(workB,nLMB,nLMBs2,rhs1)

         ! determine the number of chunks of m
         ! total number for l1 is sizeLMB2(nLMB2,nLMB)
         ! chunksize is given
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         l1=lm22l(1,nLMB2,nLMB)
         if ( l1 > 0 ) then
            if ( .not. lWPmat(l1) ) then
               !PERFON('upWP_mat')
               call get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                    wpMat(1,1,l1),wpPivot(1,l1),wpMat_fac(1,1,l1))
               lWPmat(l1)=.TRUE.
               !PERFOFF
            end if
            !write(*,"(2(A,I3))") "Running ",nChunks," for task with l1=",l1
            do iChunk=1,nChunks
               !$OMP TASK if (nChunks>1) default(shared) &
               !$OMP firstprivate(iChunk) &
               !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_cheb) &
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

                  lmB=lmB+1
                  rhs1(1,lmB,threadid)        =0.D0
                  rhs1(n_r_max,lmB,threadid)  =0.D0
                  rhs1(n_r_max+1,lmB,threadid)=0.D0
                  rhs1(2*n_r_max,lmB,threadid)=0.D0
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=                         &
                          & O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*w(lm1,nR) + &
                          & rho0(nR)*agrav(nR)*s(lm1,nR) + &
                          & w1*dwdt(lm1,nR) + &
                          & w2*dwdtLast(lm1,nR)
                     rhs1(nR+n_r_max,lmB,threadid)=                 &
                          -O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*dw(lm1,nR) + &
                          w1*dpdt(lm1,nR) + &
                          w2*dpdtLast(lm1,nR)
                  end do
               end do
               !PERFOFF
               !PERFON('upWP_sol')
               !if ( lmB > 0 ) then

               ! use the mat_fac(:,1) to scale the rhs
               do lm=lmB0+1,lmB
                  do nR=1,2*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,1,l1)
                  end do
               end do
#ifdef WITH_MKL_LU
               call getrs(cmplx(wpMat(:,:,l1),0.D0,kind=kind(0.D0)), &
                    &       wpPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid))
#else
               call cgeslML(wpMat(:,:,l1),2*n_r_max,2*n_r_max,    &
                    &       wpPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),2*n_r_max,lmB-lmB0)
#endif
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,2*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,2,l1)
                  end do
               end do
            !end if
               !PERFOFF

               if ( lRmsNext ) then ! Store old w
                  do nR=1,n_r_max
                     do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                        lm1=lm22lm(lm,nLMB2,nLMB)
                        workB(lm1,nR)=w(lm1,nR)
                     end do
                  end do
               end if

               !PERFON('upWP_aft')
               lmB=lmB0
               do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                  lm1=lm22lm(lm,nLMB2,nLMB)
                  !l1 =lm22l(lm,nLMB2,nLMB)
                  m1 =lm22m(lm,nLMB2,nLMB)
                  !if ( l1 > 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                        p(lm1,n_cheb)=rhs1(n_r_max+n_cheb,lmB,threadid)
                     end do
                  else
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)= cmplx(real(rhs1(n_cheb,lmB, &
                                       threadid)),0.D0,kind=kind(0d0))
                        p(lm1,n_cheb)= cmplx(real(rhs1(n_r_max+n_cheb, & 
                                       lmB,threadid)),0.D0,kind=kind(0d0))
                     end do
                  end if
               end do
               !PERFOFF
               !$OMP END TASK
            end do
         end if
         !write(*,"(3(A,I3))") "End of task ",nLMB2,"/",nLMBs2(nLMB)," on thread ",omp_get_thread_num()
         !$OMP END TASK
      end do   ! end of loop over l1 subblocks
      !$OMP END SINGLE
      !$OMP END PARALLEL
      !PERFOFF
      !write(*,"(A,I3,4ES22.12)") "w,p after: ",nLMB,get_global_SUM(w),get_global_SUM(p)

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      do n_cheb=n_cheb_max+1,n_r_max
         do lm1=lmStart_00,lmStop
            w(lm1,n_cheb)=cmplx(0.D0,0.D0,kind=kind(0d0))
            p(lm1,n_cheb)=cmplx(0.D0,0.D0,kind=kind(0d0))
         end do
      end do


      !PERFON('upWP_drv')
      all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
      if (all_lms < omp_get_max_threads()) then
         call omp_set_num_threads(all_lms)
      end if
#endif
      !$OMP PARALLEL default(none) &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(all_lms,per_thread,lmStart_real,lmStop_real) &
      !$OMP shared(w,dw,ddw,p,dp,dwdtLast,dpdtLast) &
      !$OMP shared(i_costf_init,d_costf_init,drx,ddrx,dddrx) &
      !$OMP shared(n_r_max,n_cheb_max,nThreads,llm_real,ulm_real,workA)
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
         start_lm=lmStart_real+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop_real
         !write(*,"(2(A,I3),2(A,I5))") "iThread=",iThread," on thread ",omp_get_thread_num(),&
         !     & " lm = ",start_lm,":",stop_lm

         !-- Transform to radial space and get radial derivatives
         !   using dwdtLast, dpdtLast as work arrays:
         call costf1( w, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
              &       dwdtLast,i_costf_init,d_costf_init)
         call get_dddr( w, dw, ddw, workA, &
              &         ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
              &         n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
              &         i_costf_init,d_costf_init,drx,ddrx,dddrx)
         call costf1( p, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
              &       dwdtLast,i_costf_init,d_costf_init)
         call get_dr( p, dp, &
              &       ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
              &       n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
              &       i_costf_init,d_costf_init,drx)
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(omp_get_max_threads())
#endif
      !PERFOFF

      !PERFON('upWP_ex')
      !-- Calculate explicit time step part:
      if ( ra /= 0.D0 ) then
         do nR=n_r_cmb+1,n_r_icb-1
            do lm1=lmStart_00,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)

               Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
                          or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                    &   +(2.d0*dLvisc(nR)-beta(nR)/3.d0)*        dw(lm1,nR) &
                    &   -( dLh(st_map%lm2(l1,m1))*or2(nR)+4.d0/3.d0* (      &
                    &        dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                    &        +(3.d0*dLvisc(nR)+beta(nR))*or1(nR) )   )*     &
                    &                                            w(lm1,nR)  )
               Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
               Buo(lm1) = rho0(nR)*rgrav(nR)*s(lm1,nR)
               dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
               dpdtLast(lm1,nR)= dpdt(lm1,nR) - coex*(                    &
                    &            dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
                    &          + hdif_V(st_map%lm2(l1,m1))*               &
                    &            visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR)  &
                    &                                  * ( -workA(lm1,nR) &
                    &                  + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                    &          + ( dLh(st_map%lm2(l1,m1))*or2(nR)         &
                    &             + dLvisc(nR)*beta(nR)+ dbeta(nR)        &
                    &             + 2.d0*(dLvisc(nR)+beta(nR))*or1(nR)    &
                    &                                      ) * dw(lm1,nR) &
                    &          - dLh(st_map%lm2(l1,m1))*or2(nR)           &
                    &             * ( 2.d0*or1(nR)+2.d0/3.d0*beta(nR)     &
                                     +dLvisc(nR) )   *         w(lm1,nR)  &
                    &                                    ) )
               if ( lRmsNext ) then
                  workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                       &        ( w(lm1,nR)-workB(lm1,nR) )
                  if ( l_RMStest ) workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                    DifPol2hInt(nR,1),DifPolAs2hInt(nR,1),lo_map)
               !write(*,"(A,I4,3ES22.14)") "upWP, work=",nR,SUM(workB(:,nR)),dtVPol2hInt(nR,nTh)
               call hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                    dtVPolLMr,dtVPol2hInt(nR,1),dtVPolAs2hInt(nR,1),lo_map)
               !write(*,"(A,2I4,ES22.14)") "upWP: ",nR,nTh,dtVPol2hInt(nR,nTh)
            end if
         end do

      else  ! no s-contribution !

         do nR=n_r_cmb+1,n_r_icb-1
            do lm1=lmStart_00,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)
               Dif(lm1)=hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
                                                        or2(nR)*visc(nR)* &
                                 (                          ddw(lm1,nR) + &
                             (2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) - &
                    (dLh(st_map%lm2(l1,m1))*or2(nR)+4.d0/3.d0*(dbeta(nR)+ &
                                dLvisc(nR)*beta(nR)+                      &
                    (3.d0*dLvisc(nR)+beta(nR))*or1(nR))) *                &
                                                              w(lm1,nR) )
               Pre(lm1)=-dp(lm1,nR)+beta(nR)*p(lm1,nR)
               dwdtLast(lm1,nR) = dwdt(lm1,nR) - coex*(Pre(lm1)+Dif(lm1))
               dpdtLast(lm1,nR)=        dpdt(lm1,nR) - coex*(  &
                    dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) + &
                    hdif_V(st_map%lm2(l1,m1))*visc(nR)*        &
                    dLh(st_map%lm2(l1,m1))*or2(nR) * (         &
                                       -workA(lm1,nR)        + &
                    (beta(nR)-dLvisc(nR))*ddw(lm1,nR)        + &
                    ( dLh(st_map%lm2(l1,m1))*or2(nR)+          &
                      dLvisc(nR)*beta(nR) + dbeta(nR)+         &
                     2.d0*(dLvisc(nR)+beta(nR))*or1(nR))    *  &
                                           dw(lm1,nR)        - &
                    dLh(st_map%lm2(l1,m1))*or2(nR)           * &
                    (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+          &
                    dLvisc(nR))*            w(lm1,nR)          )  )
               if ( lRmsNext ) then
                  workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                       ( w(lm1,nR)-workB(lm1,nR) )
                  if ( l_RMStest ) workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                    DifPol2hInt(nR,1),DifPolAs2hInt(nR,1),lo_map)
               call hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                    dtVPolLMr, dtVPol2hInt(nR,1),dtVpolAs2hInt(nR,1),lo_map)
            end if
         end do
         
      end if
      !PERFOFF

      !deallocate(rhs1)

      !  Note: workA=dddw not needed beyond this point!

   end subroutine updateWP
   !------------------------------------------------------------------------------
   subroutine get_wpMat(dt,l,hdif,wpMat,wpPivot,wpMat_fac)
      !  +-------------+----------------+------------------------------------+
      !  |                                                                   |
      !  |  Purpose of this subroutine is to contruct the time step matrix   | 
      !  |  wpmat  for the NS equation.                                      |
      !  |                                                                   |
      !  +-------------------------------------------------------------------+

      !-- Input variables:
      real(kind=8), intent(in) :: dt
      real(kind=8), intent(in) :: hdif
      integer,      intent(in) :: l

      !-- Output variables:
      real(kind=8), intent(out) :: wpMat(2*n_r_max,2*n_r_max)
      real(kind=8), intent(out) :: wpMat_fac(2*n_r_max,2)
      integer,      intent(out) :: wpPivot(2*n_r_max)

      !-- local variables:
      integer :: nR,nCheb,nR_p,nCheb_p
      integer :: info
      real(kind=8) :: O_dt,dLh

#ifdef MATRIX_CHECK
      integer ::ipiv(2*n_r_max),iwork(2*n_r_max),i,j
      real(kind=8) :: work(8*n_r_max),anorm,linesum,rcond
      real(kind=8) :: temp_wpMat(2*n_r_max,2*n_r_max)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
      logical :: first_run=.true.
#endif

#if 0
      if (first_run) then
         open(NEWUNIT=filehandle,file="cheb.dat")
         do nR=1,n_r_max
            do nCheb=1,n_r_max
               write(filehandle,"(ES20.12)",advance='no') cheb(nCheb,nR)
            end do
            write(filehandle,"(A)") ""
         end do
         close(filehandle)
         open(NEWUNIT=filehandle,file="dcheb.dat")
         do nR=1,n_r_max
            do nCheb=1,n_r_max
               write(filehandle,"(ES20.12)",advance='no') dcheb(nCheb,nR)
            end do
            write(filehandle,"(A)") ""
         end do
         close(filehandle)
         open(NEWUNIT=filehandle,file="d2cheb.dat")
         do nR=1,n_r_max
            do nCheb=1,n_r_max
               write(filehandle,"(ES20.12)",advance='no') d2cheb(nCheb,nR)
            end do
            write(filehandle,"(A)") ""
         end do
         close(filehandle)
         open(NEWUNIT=filehandle,file="d3cheb.dat")
         do nR=1,n_r_max
            do nCheb=1,n_r_max
               write(filehandle,"(ES20.12)",advance='no') d3cheb(nCheb,nR)
            end do
            write(filehandle,"(A)") ""
         end do
         close(filehandle)
         first_run=.false.
      end if
#endif

      O_dt=1.D0/dt
      dLh =dble(l*(l+1))
    
      !-- Now mode l>0
    
      !----- Boundary conditions, see above:
      do nCheb=1,n_cheb_max
         nCheb_p=nCheb+n_r_max
    
         wpMat(1,nCheb)        =cheb_norm*cheb(nCheb,1)
         wpMat(1,nCheb_p)      =0.D0
         wpMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
         wpMat(n_r_max,nCheb_p)=0.D0
    
         if ( ktopv == 1 ) then  ! free slip !
            wpMat(n_r_max+1,nCheb)=   cheb_norm * ( &
                 d2cheb(nCheb,1) - (2.d0*or1(1)+beta(1))*dcheb(nCheb,1) )
         else                    ! no slip, note exception for l=1,m=0
            wpMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,1)
         end if
         wpMat(n_r_max+1,nCheb_p)=0.D0
    
         if ( kbotv == 1 ) then  ! free slip !
            wpMat(2*n_r_max,nCheb)=        cheb_norm * ( &
                 d2cheb(nCheb,n_r_max) - &
                 (2.d0*or1(n_r_max)+beta(n_r_max))*dcheb(nCheb,n_r_max))
         else                 ! no slip, note exception for l=1,m=0
            wpMat(2*n_r_max,nCheb)=cheb_norm * dcheb(nCheb,n_r_max)
         end if
         wpMat(2*n_r_max,nCheb_p)=0.D0
    
      end do   !  loop over nCheb
    
      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            wpMat(1,nCheb)          =0.D0
            wpMat(n_r_max,nCheb)    =0.D0
            wpMat(n_r_max+1,nCheb)  =0.D0
            wpMat(2*n_r_max,nCheb)  =0.D0
            wpMat(1,nCheb_p)        =0.D0
            wpMat(n_r_max,nCheb_p)  =0.D0
            wpMat(n_r_max+1,nCheb_p)=0.D0
            wpMat(2*n_r_max,nCheb_p)=0.D0
         end do
      end if
    
      !----- Other points:
      do nCheb=1,n_r_max
         nCheb_p=nCheb+n_r_max
         do nR=2,n_r_max-1
            !write(*,"(I3,A,6ES11.3)") nR,", visc,beta,dLvisc,dbeta = ",&
            !     & visc(nR),beta(nR),dLvisc(nR),dbeta(nR),hdif,alpha
            ! in the BM2 case: visc=1.0,beta=0.0,dLvisc=0.0,dbeta=0.0
            nR_p=nR+n_r_max
            wpMat(nR,nCheb)= cheb_norm *  ( O_dt*dLh*or2(nR)*cheb(nCheb,nR) &
                 &   - alpha*hdif*visc(nR)*dLh*or2(nR) * ( d2cheb(nCheb,nR) &
                 &          +(2.*dLvisc(nR)-beta(nR)/3.D0)*dcheb(nCheb,nR)  &
                 &         -( dLh*or2(nR)+4.D0/3.D0*( dLvisc(nR)*beta(nR)   &
                 &          +(3.d0*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) ) &
                 &          )                               *cheb(nCheb,nR) &
                 &                                       )  )
    
            wpMat(nR,nCheb_p)= cheb_norm*alpha*(dcheb(nCheb,nR) &
                               -beta(nR)* cheb(nCheb,nR))
            ! the following part gives sometimes very large 
            ! matrix entries
            wpMat(nR_p,nCheb)= cheb_norm * ( -O_dt*dLh*or2(nR)*dcheb(nCheb,nR) &
                 &  -alpha*hdif*visc(nR)*dLh*or2(nR) &    *(- d3cheb(nCheb,nR) &
                 &                   +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
                 &      +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
                 &       +2.D0*(dLvisc(nR)+beta(nR))*or1(nR) )*dcheb(nCheb,nR) &
                 &      -dLh*or2(nR)*( 2.d0*or1(nR)+dLvisc(nR)                 &
                 &                     +2.d0/3.d0*beta(nR)   )* cheb(nCheb,nR) &
                 &                                        ) )
    
            wpMat(nR_p,nCheb_p)= -cheb_norm*alpha*dLh*or2(nR)*cheb(nCheb,nR)
         end do
      end do
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         wpMat(nR,1)          =0.5D0*wpMat(nR,1)
         wpMat(nR,n_r_max)    =0.5D0*wpMat(nR,n_r_max)
         wpMat(nR,n_r_max+1)  =0.5D0*wpMat(nR,n_r_max+1)
         wpMat(nR,2*n_r_max)  =0.5D0*wpMat(nR,2*n_r_max)
         wpMat(nR_p,1)        =0.5D0*wpMat(nR_p,1)
         wpMat(nR_p,n_r_max)  =0.5D0*wpMat(nR_p,n_r_max)
         wpMat(nR_p,n_r_max+1)=0.5D0*wpMat(nR_p,n_r_max+1)
         wpMat(nR_p,2*n_r_max)=0.5D0*wpMat(nR_p,2*n_r_max)
      end do
    
      ! compute the linesum of each line
      do nR=1,2*n_r_max
         wpMat_fac(nR,1)=1.0D0/maxval(abs(wpMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         wpMat(nR,:) = wpMat(nR,:)*wpMat_fac(nR,1)
      end do
    
      ! also compute the rowsum of each column
      do nR=1,2*n_r_max
         wpMat_fac(nR,2)=1.0D0/maxval(abs(wpMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,2*n_r_max
         wpMat(:,nR) = wpMat(:,nR)*wpMat_fac(nR,2)
      end do

#ifdef MATRIX_CHECK
      ! copy the wpMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "wpMat_",l,"_",counter,".dat"
      open(NEWUNIT=filehandle,file=trim(filename))
      counter= counter+1
      
      do i=1,2*n_r_max
         do j=1,2*n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") wpMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_wpMat=wpMat
      anorm = 0.0D0
      do i=1,2*n_r_max
         linesum = 0.0D0
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

#ifdef WITH_MKL_LU
      call getrf(wpMat,wpPivot,info)
#else
      call sgefa(wpMat,2*n_r_max,2*n_r_max,wpPivot,info)
#endif
      if ( info /= 0 ) then
         write(*,*) 'Singular matrix wpmat!'
         stop '35'
      end if

   end subroutine get_wpMat
!-----------------------------------------------------------------------------
end module updateWP_mod
