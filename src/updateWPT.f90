#include "perflib_preproc.cpp"
module updateWPT_mod
   
   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_cheb_max, n_r_max, l_max
   use radial_data, only: n_r_cmb,n_r_icb
   use radial_functions, only: drx,ddrx,dddrx,or1,or2,rho0,rgrav,       &
       &                       chebt_oc,visc,dLvisc,orho1,              &
       &                       beta,dbeta,cheb,dcheb,d2cheb,d3cheb,     &
       &                       cheb_norm, dLkappa, dLtemp0, ddLtemp0,   &
       &                       alpha0,dLalpha0, ddLalpha0, otemp1,      &
       &                       kappa, orho2, dentropy0, temp0, r
   use physical_parameters, only: kbotv, ktopv, ktops, kbots, ra, opr, &
       &                          ViscHeatFac, ThExpNb, ogrun, BuoFac, &
       &                          CorFac, ktopp
   use num_param, only: alpha
   use init_fields, only: tops, bots
   use blocking, only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map, &
       &               lmStartB,lmStopB
   use horizontal_data, only: hdif_V, hdif_S, dLh
   use logic, only: l_update_v, l_temperature_diff
   use RMS, only: DifPol2hInt, dtVPolLMr, dtVPol2hInt, DifPolLMr
   use RMS_helpers, only:  hInt2Pol
   use algebra, only: cgeslML, sgefa, sgesl
   use LMLoop_data, only: llm, ulm
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize, rank
   use cosine_transform_odd
   use radial_der, only: get_dddr, get_ddr, get_dr
   use integration, only: rInt_R
   use constants, only: zero, one, two, three, four, third, half, pi, osq4pi

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: workA(:,:),workB(:,:), workC(:,:),workD(:,:)
   complex(cp), allocatable :: Dif(:),Pre(:),Buo(:),dtV(:)
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: pt0Mat(:,:), pt0Mat_fac(:)
   integer, allocatable :: pt0Pivot(:)
   real(cp), allocatable :: wptMat(:,:,:)
   integer, allocatable :: wptPivot(:,:)
   real(cp), allocatable :: wptMat_fac(:,:,:)
   logical, public, allocatable :: lWPTmat(:)
   real(cp) :: Cor00_fac
   integer :: maxThreads

   public :: initialize_updateWPT, finalize_updateWPT, updateWPT

contains

   subroutine initialize_updateWPT

      allocate( pt0Mat(2*n_r_max,2*n_r_max) )
      allocate( pt0Mat_fac(2*n_r_max) )
      allocate( pt0Pivot(2*n_r_max) )
      bytes_allocated = bytes_allocated+(4*n_r_max+2)*n_r_max*SIZEOF_DEF_REAL &
      &                 +2*n_r_max*SIZEOF_INTEGER
      allocate( wptMat(3*n_r_max,3*n_r_max,l_max) )
      allocate(wptMat_fac(3*n_r_max,2,l_max))
      allocate ( wptPivot(3*n_r_max,l_max) )
      bytes_allocated = bytes_allocated+(9*n_r_max*l_max+6*n_r_max*l_max)*&
      &                 SIZEOF_DEF_REAL+3*n_r_max*l_max*SIZEOF_INTEGER
      allocate( lWPTmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL


      allocate( workA(llm:ulm,n_r_max) )
      allocate( workB(llm:ulm,n_r_max) )
      allocate( workC(llm:ulm,n_r_max) )
      allocate( workD(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+4*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

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

      allocate( rhs1(3*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated=bytes_allocated+2*n_r_max*maxThreads* &
                      lo_sub_map%sizeLMB2max*SIZEOF_DEF_COMPLEX

      Cor00_fac=four/sqrt(three)


   end subroutine initialize_updateWPT
!-----------------------------------------------------------------------------
   subroutine finalize_updateWPT

      deallocate( wptMat, wptMat_fac, wptPivot )
      deallocate( pt0Mat, pt0Mat_fac, pt0Pivot, lWPTmat )
      deallocate( workA, workB, workC, rhs1 )
      deallocate( Dif, Pre, Buo, dtV )

   end subroutine finalize_updateWPT
!-----------------------------------------------------------------------------
   subroutine updateWPT(w,dw,ddw,z10,dwdt,dwdtLast,p,dp,dpdt,dpdtLast,tt,  &
        &               dtt,dVTrLM,dVPrLM,dttdt,dttdtLast,w1,coex,dt,nLMB, &
        &               lRmsNext)
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
      complex(cp), intent(in) :: dwdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdt(llm:ulm,n_r_max)
      real(cp),    intent(in) :: z10(n_r_max)

      complex(cp), intent(inout) :: dttdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dwdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dpdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: tt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dtt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dttdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVTrLM(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVPrLM(llm:ulm,n_r_max)

      complex(cp), intent(out) :: ddw(llm:ulm,n_r_max)
      complex(cp), intent(out) :: dp(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm,lmB     ! position of (l,m) in array
      integer :: lmStart,lmStop ! max and min number of orders m
      integer :: nLMB2
      integer :: nR             ! counts radial grid points
      integer :: n_cheb         ! counts cheb modes
      real(cp) :: rhs(2*n_r_max)  ! real RHS for l=m=0
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

      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)

      w2  =one-w1
      O_dt=one/dt

      !$OMP PARALLEL  &
      !$OMP private(iThread,start_lm,stop_lm,nR,lm) &
      !$OMP shared(all_lms,per_thread,lmStart,lmStop) &
      !$OMP shared(dVTrLM,dVPrLM,chebt_oc,drx,dttdt,orho1,or2) &
      !$OMP shared(alpha0,temp0,dLalpha0) &
      !$OMP shared(n_r_max,n_cheb_max,workA,workB,workC,workD) &
      !$OMP shared(nThreads,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of tt: workA,dttdtLast used as work arrays
      all_lms=lmStop-lmStart+1
      per_thread=all_lms/nThreads
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop

         !--- Finish calculation of dttdt:
         call get_dr( dVTrLM,workA,ulm-llm+1,start_lm-llm+1,  &
              &       stop_lm-llm+1,n_r_max,n_cheb_max,workB, &
              &       workC,chebt_oc,drx)

         call get_dr( dVPrLM,workB,ulm-llm+1,start_lm-llm+1,  &
              &       stop_lm-llm+1,n_r_max,n_cheb_max,workC, &
              &       workD,chebt_oc,drx)
      end do
      !$OMP end do

      !$OMP DO
      do nR=1,n_r_max
         do lm=lmStart,lmStop
            dttdt(lm,nR)=          orho1(nR)*dttdt(lm,nR) -             &
            &            or2(nR)*orho1(nR)*workA(lm,nR) +               &
            &            ViscHeatFac*ThExpNb*or2(nR)*alpha0(nR)*        &
            &            temp0(nR)*orho2(nR)*workB(lm,nR) +             &
            &            or2(nR)*orho1(nR)*dLtemp0(nR)*dVTrLM(lm,nR)+   &
            &            ViscHeatFac*ThExpNb*or2(nR)*alpha0(nR)*        &
            &            temp0(nR)*orho2(nR)*( dLalpha0(nR)-beta(nR) )* &   
            &                                           dVPrLM(lm,nR)
         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL


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
         if ( l1 == 0 ) then
            if ( .not. lWPTmat(l1) ) then
               call get_pt0Mat(dt,pt0Mat,pt0Pivot,pt0Mat_fac)
               lWPTmat(l1)=.true.
            end if
         else
            if ( .not. lWPTmat(l1) ) then
               call get_wptMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                    &          hdif_S(st_map%lm2(l1,0)),       &
                    &          wptMat(1,1,l1),wptPivot(1,l1),  &
                    &          wptMat_fac(1,1,l1))
               lWPTmat(l1)=.true.
            end if
         end if

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

               if ( l1 == 0 ) then

                  do nR=1,n_r_max
                     rhs(nR)=real(tt(lm1,nR))*O_dt-ViscHeatFac*ThExpNb*&
                     &       alpha0(nR)*temp0(nR)*orho1(nR)*           &
                     &       real(p(lm1,nR))*O_dt+                     &
                     &       w1*real(dttdt(lm1,nR)) +                  &
                     &       w2*real(dttdtLast(lm1,nR))
                     rhs(nR+n_r_max)=real(dwdt(lm1,nR))+Cor00_fac*&
                     &               CorFac*or1(nR)*z10(nR)
                  end do
                  rhs(1)        =real(tops(0,0))
                  rhs(n_r_max)  =real(bots(0,0))
                  rhs(n_r_max+1)=0.0_cp

                  rhs = pt0Mat_fac*rhs

                  call sgesl(pt0Mat,2*n_r_max,2*n_r_max,pt0Pivot,rhs)

               else ! l1 /= 0
                  lmB=lmB+1
                  rhs1(1,lmB,threadid)          =0.0_cp
                  rhs1(n_r_max,lmB,threadid)    =0.0_cp
                  rhs1(n_r_max+1,lmB,threadid)  =0.0_cp
                  rhs1(2*n_r_max,lmB,threadid)  =0.0_cp
                  rhs1(2*n_r_max+1,lmB,threadid)=tops(l1,m1)
                  rhs1(3*n_r_max,lmB,threadid)  =bots(l1,m1)
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=O_dt*dLh(st_map%lm2(l1,m1))*  &
                     &                     or2(nR)*w(lm1,nR) +           &
                     &                     w1*dwdt(lm1,nR) +             &
                     &                     w2*dwdtLast(lm1,nR)
                     rhs1(nR+n_r_max,lmB,threadid)=-O_dt*                 &
                     &                             dLh(st_map%lm2(l1,m1))*&
                     &                             or2(nR)*dw(lm1,nR) +   &
                     &                             w1*dpdt(lm1,nR) +      &
                     &                             w2*dpdtLast(lm1,nR)
                     rhs1(nR+2*n_r_max,lmB,threadid)=tt(lm1,nR)*O_dt -       &
                     &                               ViscHeatFac*ThExpNb*    &
                     &                               alpha0(nR)*temp0(nR)*   &
                     &                               orho1(nR)*p(lm1,nR)*O_dt&
                     &                               + w1*dttdt(lm1,nR) +    &
                     &                               w2*dttdtLast(lm1,nR)
                  end do
               end if
            end do
            !PERFOFF

            !PERFON('upWP_sol')
            if ( lmB > 0 ) then

               ! use the mat_fac(:,1) to scale the rhs
               do lm=lmB0+1,lmB
                  do nR=1,3*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wptMat_fac(nR,1,l1)
                  end do
               end do
               call cgeslML(wptMat(:,:,l1),3*n_r_max,3*n_r_max,        &
                    &       wptPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),&
                    &       3*n_r_max,lmB-lmB0)
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,3*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wptMat_fac(nR,2,l1)
                  end do
               end do
            end if
            !PERFOFF

            if ( lRmsNext ) then ! Store old w
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     workD(lm1,nR)=w(lm1,nR)
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
                  do n_cheb=1,n_cheb_max
                     tt(lm1,n_cheb)=rhs(n_cheb)
                     p(lm1,n_cheb) =rhs(n_cheb+n_r_max)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb) =rhs1(n_cheb,lmB,threadid)
                        p(lm1,n_cheb) =rhs1(n_r_max+n_cheb,lmB,threadid)
                        tt(lm1,n_cheb)=rhs1(2*n_r_max+n_cheb,lmB,threadid)
                     end do
                  else
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)= cmplx(real(rhs1(n_cheb,lmB,threadid)), &
                                       &     0.0_cp,kind=cp)
                        p(lm1,n_cheb)= cmplx(real(rhs1(n_r_max+n_cheb,lmB,threadid)), &
                                       &     0.0_cp,kind=cp)
                        tt(lm1,n_cheb)= cmplx(real(rhs1(2*n_r_max+n_cheb,lmB,threadid)), &
                                       &     0.0_cp,kind=cp)
                     end do
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

      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      do n_cheb=n_cheb_max+1,n_r_max
         do lm1=lmStart,lmStop
            w(lm1,n_cheb)=zero
            p(lm1,n_cheb)=zero
            tt(lm1,n_cheb)=zero
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
      !$OMP shared(all_lms,per_thread,lmStop) &
      !$OMP shared(w,dw,ddw,p,dp,tt,dtt,dwdtLast,dpdtLast,dttdtLast) &
      !$OMP shared(chebt_oc,drx,ddrx,dddrx) &
      !$OMP shared(n_r_max,n_cheb_max,nThreads,workA,workB,workC,llm,ulm)
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

         call chebt_oc%costf1(w,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,dwdtLast)
         call get_dddr( w, dw, ddw, workA, ulm-llm+1, start_lm-llm+1,  &
              &         stop_lm-llm+1, n_r_max,n_cheb_max,dwdtLast,    &
              &         dpdtLast,chebt_oc,drx,ddrx,dddrx)

         call chebt_oc%costf1(p,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,dwdtLast)
         call get_ddr( p, dp, workC, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max,n_cheb_max,dwdtLast,dpdtLast,                    &
              &       chebt_oc,drx,ddrx)

         call chebt_oc%costf1(tt,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,dttdtLast)
         call get_ddr(tt, dtt, workB, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, n_cheb_max, dwdtLast, dttdtLast,                &
              &       chebt_oc,drx,ddrx)

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

      !-- Calculate explicit time step part:
      if ( l_temperature_diff ) then
         do nR=n_r_top,n_r_bot
            do lm1=lmStart,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)

               Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
                    &     or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                    &   +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
                    &   -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
                    &        dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                    &        +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
                    &                                             w(lm1,nR) )
               Pre(lm1) = -dp(lm1,nR)-BuoFac*ViscHeatFac*(        &
                    &     ThExpNb*alpha0(nR)*temp0(nR)+ogrun)*    &
                    &     alpha0(nR)*rgrav(nR)*p(lm1,nR)
               Buo(lm1) = BuoFac*rho0(nR)*rgrav(nR)*alpha0(nR)*tt(lm1,nR)
               dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
               dpdtLast(lm1,nR)= dpdt(lm1,nR) - coex*(                    &
                    &            dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
                    &          + hdif_V(st_map%lm2(l1,m1))*               &
                    &            visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR)  &
                    &                                  * ( -workA(lm1,nR) &
                    &                  + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                    &          + ( dLh(st_map%lm2(l1,m1))*or2(nR)         &
                    &             + dLvisc(nR)*beta(nR)+ dbeta(nR)        &
                    &             + two*(dLvisc(nR)+beta(nR))*or1(nR)     &
                    &                                      ) * dw(lm1,nR) &
                    &          - dLh(st_map%lm2(l1,m1))*or2(nR)           &
                    &             * ( two*or1(nR)+two*third*beta(nR)      &
                                     +dLvisc(nR) )   *         w(lm1,nR)  &
                    &                                    ) )
               dttdtLast(lm1,nR)=dttdt(lm1,nR) &
                    & - coex*opr*hdif_S(st_map%lm2(l1,m1)) * kappa(nR) *  &
                    &   (             workB(lm1,nR)                       &
                    &     + ( beta(nR)+two*or1(nR)+dLkappa(nR) ) *        &
                    &                    dtt(lm1,nR) -                    &
                    &       dLh(st_map%lm2(l1,m1))*or2(nR)  *             &
                    &                     tt(lm1,nR)  )+                  &
                    &   coex*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR) &
                    &   *temp0(nR)*orho1(nR)*dentropy0(nR)*w(lm1,nR)
               if ( lRmsNext ) then
                  dtV(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                  &        ( w(lm1,nR)-workD(lm1,nR) )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart,lmStop,DifPolLMr(1,nR), &
                    &        DifPol2hInt(:,nR,1),lo_map)
               call hInt2Pol(dtV,llm,ulm,nR,lmStart,lmStop, &
                    &        dtVPolLMr(1,nR),dtVPol2hInt(:,nR,1),lo_map)
            end if
         end do

      else ! entropy diffusion

         do nR=n_r_top,n_r_bot
            do lm1=lmStart,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)

               Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))* &
                    &     or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                    &   +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
                    &   -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
                    &        dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                    &        +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
                    &                                            w(lm1,nR)  )
               Pre(lm1) = -dp(lm1,nR)-BuoFac*ViscHeatFac*(        &
                    &     ThExpNb*alpha0(nR)*temp0(nR)+ogrun)*    &
                    &     alpha0(nR)*rgrav(nR)*p(lm1,nR)
               Buo(lm1) = BuoFac*rho0(nR)*rgrav(nR)*alpha0(nR)*tt(lm1,nR)
               dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
               dpdtLast(lm1,nR)= dpdt(lm1,nR) - coex*(                    &
                    &            dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
                    &          + hdif_V(st_map%lm2(l1,m1))*               &
                    &            visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR)  &
                    &                                  * ( -workA(lm1,nR) &
                    &                  + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                    &          + ( dLh(st_map%lm2(l1,m1))*or2(nR)         &
                    &             + dLvisc(nR)*beta(nR)+ dbeta(nR)        &
                    &             + two*(dLvisc(nR)+beta(nR))*or1(nR)     &
                    &                                      ) * dw(lm1,nR) &
                    &          - dLh(st_map%lm2(l1,m1))*or2(nR)           &
                    &             * ( two*or1(nR)+two*third*beta(nR)      &
                                     +dLvisc(nR) )   *         w(lm1,nR)  &
                    &                                    ) )
               dttdtLast(lm1,nR)=dttdt(lm1,nR) &
                    &               - coex*opr*hdif_S(st_map%lm2(l1,m1))* &
                    &                                         kappa(nR) * &
                    &   (                             workB(lm1,nR)       &
                    &     + ( beta(nR) - dLtemp0(nR) +                    &
                    &       two*or1(nR) + dLkappa(nR) ) * dtt(lm1,nR)     &
                    &     - ( ddLtemp0(nR)+dLtemp0(nR)*(dLkappa(nR)+      &
                    &         beta(nR)+two*or1(nR) ) +                    &
                    &             dLh(st_map%lm2(l1,m1))*or2(nR) )        &
                    &                                  *  tt(lm1,nR)    -  &
                    &  ViscHeatFac*ThExpNb*alpha0(nR)*orho1(nR)*temp0(nR)*& 
                    &   (                             workC(lm1,nR) +     &
                    &   ( dLkappa(nR)+dLtemp0(nR)+two*dLalpha0(nR)        &
                    &     +two*or1(nR)-beta(nR))*        dp(lm1,nR) +     &
                    &   ( (dLalpha0(nR)-beta(nR))*(two*or1(nR)+           &
                    &      dLalpha0(nR)+dLkappa(nR)+dLtemp0(nR))+         &
                    &      ddLalpha0(nR)-dbeta(nR) -                      &
                    &      dLh(st_map%lm2(l1,m1))*or2(nR) )*              &
                    &                                     p(lm1,nR) ) ) + &
                    &   coex*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR) &
                    &   *orho1(nR)*dentropy0(nR)*w(lm1,nR)
               if ( lRmsNext ) then
                  dtV(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                  &        ( w(lm1,nR)-workD(lm1,nR) )
               end if
            end do
            if ( lRmsNext ) then
               call hInt2Pol(Dif,llm,ulm,nR,lmStart,lmStop,DifPolLMr(1,nR), &
                    &        DifPol2hInt(:,nR,1),lo_map)
               call hInt2Pol(dtV,llm,ulm,nR,lmStart,lmStop, &
                    &        dtVPolLMr(1,nR),dtVPol2hInt(:,nR,1),lo_map)
            end if
         end do
      end if

   end subroutine updateWPT
!------------------------------------------------------------------------------
   subroutine get_wptMat(dt,l,hdif_vel,hdif_t,wptMat,wptPivot,wptMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix  
      !  wpmat  for the NS equation.                                    
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif_vel
      real(cp), intent(in) :: hdif_t
      integer,  intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: wptMat(3*n_r_max,3*n_r_max)
      real(cp), intent(out) :: wptMat_fac(3*n_r_max,2)
      integer,  intent(out) :: wptPivot(3*n_r_max)

      !-- local variables:
      integer :: nR,nCheb,nR_p,nR_t,nCheb_p,nCheb_t
      integer :: info
      real(cp) :: O_dt,dLh

      O_dt=one/dt
      dLh =real(l*(l+1),kind=cp)
    
      !-- Now mode l>0
    
      !----- Boundary conditions, see above:
      do nCheb=1,n_cheb_max
         nCheb_p=nCheb+n_r_max
         nCheb_t=nCheb+2*n_r_max
    
         wptMat(1,nCheb)        =cheb_norm*cheb(nCheb,1)
         wptMat(1,nCheb_p)      =0.0_cp
         wptMat(1,nCheb_t)      =0.0_cp
         wptMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
         wptMat(n_r_max,nCheb_p)=0.0_cp
         wptMat(n_r_max,nCheb_t)=0.0_cp
    
         if ( ktopv == 1 ) then  ! free slip !
            wptMat(n_r_max+1,nCheb)=   cheb_norm * ( &
                 d2cheb(nCheb,1) - (two*or1(1)+beta(1))*dcheb(nCheb,1) )
         else                    ! no slip, note exception for l=1,m=0
            wptMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,1)
         end if
         wptMat(n_r_max+1,nCheb_p)=0.0_cp
         wptMat(n_r_max+1,nCheb_t)=0.0_cp

         if ( kbotv == 1 ) then  ! free slip !
            wptMat(2*n_r_max,nCheb)=        cheb_norm * ( &
                 d2cheb(nCheb,n_r_max) - &
                 (two*or1(n_r_max)+beta(n_r_max))*dcheb(nCheb,n_r_max))
         else                 ! no slip, note exception for l=1,m=0
            wptMat(2*n_r_max,nCheb)=cheb_norm * dcheb(nCheb,n_r_max)
         end if
         wptMat(2*n_r_max,nCheb_p)=0.0_cp
         wptMat(2*n_r_max,nCheb_t)=0.0_cp

         if ( ktops == 1 ) then ! fixed entropy
            wptMat(2*n_r_max+1,nCheb_t)=cheb_norm*otemp1(1)
            wptMat(2*n_r_max+1,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb* &
            &                           alpha0(1)*orho1(1)
         else if ( ktops == 2 ) then ! fixed entropy flux
            wptMat(2*n_r_max+1,nCheb_t)=cheb_norm*otemp1(1)*( dcheb(nCheb,1)- &
            &                           dLtemp0(1)*   cheb(nCheb,1) )
            wptMat(2*n_r_max+1,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*alpha0(1)*&
            &                            orho1(1)*(            dcheb(nCheb,1)+   &
            &                           (dLalpha0(1)-beta(1))* cheb(nCheb,1) )
         else if ( ktops == 3 ) then ! fixed temperature
            wptMat(2*n_r_max+1,nCheb_t)=cheb_norm
            wptMat(2*n_r_max+1,nCheb_p)=0.0_cp
         else if ( ktops == 4 ) then ! fixed temperature flux
            wptMat(2*n_r_max+1,nCheb_t)=cheb_norm*dcheb(nCheb,1)
            wptMat(2*n_r_max+1,nCheb_p)=0.0_cp
         end if
         wptMat(2*n_r_max+1,nCheb)  =0.0_cp

         if ( kbots == 1 ) then ! fixed entropy
            wptMat(3*n_r_max,nCheb_t)=cheb_norm*cheb(nCheb,n_r_max)* &
            &                         otemp1(n_r_max) 
            wptMat(3*n_r_max,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*  &
            &                         alpha0(n_r_max)*orho1(n_r_max)* &
            &                         cheb(nCheb,n_r_max)
         else if ( kbots == 2) then ! fixed entropy flux
            wptMat(3*n_r_max,nCheb_t)=cheb_norm*otemp1(n_r_max)*(         &
            &                                       dcheb(nCheb,n_r_max)- &
            &                       dLtemp0(n_r_max)*cheb(nCheb,n_r_max) )
            wptMat(3*n_r_max,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*     &
            &                       alpha0(n_r_max)*orho1(n_r_max)*(      &
            &                                       dcheb(nCheb,n_r_max)+ &
            &                       (dLalpha0(n_r_max)-beta(n_r_max))*    &
            &                                        cheb(nCheb,n_r_max) )
         else if ( kbots == 3) then ! fixed temperature
            wptMat(3*n_r_max,nCheb_t)=cheb_norm*cheb(nCheb,n_r_max)
            wptMat(3*n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 4) then ! fixed temperature flux
            wptMat(3*n_r_max,nCheb_t)=cheb_norm*dcheb(nCheb,n_r_max)
            wptMat(3*n_r_max,nCheb_p)=0.0_cp
         end if
         wptMat(3*n_r_max,nCheb)  =0.0_cp

    
      end do   !  loop over nCheb
    
      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_t=nCheb+2*n_r_max
            wptMat(1,nCheb)            =0.0_cp
            wptMat(n_r_max,nCheb)      =0.0_cp
            wptMat(n_r_max+1,nCheb)    =0.0_cp
            wptMat(2*n_r_max,nCheb)    =0.0_cp
            wptMat(2*n_r_max+1,nCheb)  =0.0_cp
            wptMat(3*n_r_max,nCheb)    =0.0_cp
            wptMat(1,nCheb_p)          =0.0_cp
            wptMat(n_r_max,nCheb_p)    =0.0_cp
            wptMat(n_r_max+1,nCheb_p)  =0.0_cp
            wptMat(2*n_r_max,nCheb_p)  =0.0_cp
            wptMat(2*n_r_max+1,nCheb_p)=0.0_cp
            wptMat(3*n_r_max,nCheb_p)  =0.0_cp
            wptMat(1,nCheb_t)          =0.0_cp
            wptMat(n_r_max,nCheb_t)    =0.0_cp
            wptMat(n_r_max+1,nCheb_t)  =0.0_cp
            wptMat(2*n_r_max,nCheb_t)  =0.0_cp
            wptMat(2*n_r_max+1,nCheb_t)=0.0_cp
            wptMat(3*n_r_max,nCheb_t)  =0.0_cp
         end do
      end if
    
      if ( l_temperature_diff ) then ! temperature diffusion

         do nCheb=1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_t=nCheb+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_t=nR+2*n_r_max

               ! W equation
               wptMat(nR,nCheb)= cheb_norm *  ( O_dt*dLh*or2(nR)*cheb(nCheb,nR)  &
                    &   - alpha*hdif_vel*visc(nR)*dLh*or2(nR) * ( d2cheb(nCheb,nR)   &
                    &          +(two*dLvisc(nR)-third*beta(nR))*dcheb(nCheb,nR)  &
                    &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)    &
                    &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )  &
                    &          )                               *cheb(nCheb,nR)   &
                    &                                       )  )

               ! Buoyancy
               wptMat(nR,nCheb_t)=-cheb_norm*alpha*BuoFac*rgrav(nR)*rho0(nR)* &
                    &                     alpha0(nR)*cheb(nCheb,nR)
       
               ! Pressure gradient
               wptMat(nR,nCheb_p)= cheb_norm*alpha*(     dcheb(nCheb,nR) &
                    &              +BuoFac*ViscHeatFac*(                 &
                    &              ThExpNb*alpha0(nR)*temp0(nR)+ogrun)*  &
                    &              alpha0(nR)*rgrav(nR)*  cheb(nCheb,nR))

               ! P equation
               wptMat(nR_p,nCheb)= cheb_norm * ( -O_dt*dLh*or2(nR)*dcheb(nCheb,nR)&
                    &  -alpha*hdif_vel*visc(nR)*dLh*or2(nR)      *(- d3cheb(nCheb,nR) &
                    &                   +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
                    &      +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
                    &       +two*(dLvisc(nR)+beta(nR))*or1(nR) )*dcheb(nCheb,nR)  &
                    &      -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
                    &                     +two*third*beta(nR)   )* cheb(nCheb,nR) &
                    &                                        ) )
       
               wptMat(nR_p,nCheb_p)= -cheb_norm*alpha*dLh*or2(nR)*cheb(nCheb,nR)

               wptMat(nR_p,nCheb_t)=0.0_cp

               ! T equation
               wptMat(nR_t,nCheb_t)= cheb_norm * (                      &
               &                               O_dt*cheb(nCheb,nR) -    &
               &      alpha*opr*hdif_t*kappa(nR)*(  d2cheb(nCheb,nR) +  &
               &      ( beta(nR)+two*or1(nR)+dLkappa(nR) )*             &
               &                                   dcheb(nCheb,nR) -    &
               &           dLh*or2(nR)  *           cheb(nCheb,nR) ) )

               wptMat(nR_t,nCheb_p)= -cheb_norm*ViscHeatFac*ThExpNb*    &
               &                     alpha0(nR)*temp0(nR)*orho1(nR)*    &
               &                     O_dt*cheb(nCheb,nR)


               !Advection of the background entropy u_r * dso/dr
               wptMat(nR_t,nCheb)=cheb_norm*alpha*dLh*or2(nR)*temp0(nR)* &
               &                  dentropy0(nR)*orho1(nR)*cheb(nCheb,nR)

            end do
         end do

      else ! entropy diffusion

         do nCheb=1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_t=nCheb+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_t=nR+2*n_r_max

               ! W equation
               wptMat(nR,nCheb)= cheb_norm *  ( O_dt*dLh*or2(nR)*cheb(nCheb,nR)  &
                    &   - alpha*hdif_vel*visc(nR)*dLh*or2(nR) * ( d2cheb(nCheb,nR)   &
                    &          +(two*dLvisc(nR)-third*beta(nR))*dcheb(nCheb,nR)  &
                    &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)    &
                    &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )  &
                    &          )                               *cheb(nCheb,nR)   &
                    &                                       )  )

               ! Buoyancy
               wptMat(nR,nCheb_t)=-cheb_norm*alpha*BuoFac*rgrav(nR)*rho0(nR)* &
               &                                 alpha0(nR)*cheb(nCheb,nR)
       
               ! Pressure gradient
               wptMat(nR,nCheb_p)= cheb_norm*alpha*(     dcheb(nCheb,nR) &
                    &             +BuoFac*ViscHeatFac*(                  &
                    &              ThExpNb*alpha0(nR)*temp0(nR)+ogrun)*  &
                    &              alpha0(nR)*rgrav(nR)*  cheb(nCheb,nR))

               ! P equation
               wptMat(nR_p,nCheb)= cheb_norm * ( -O_dt*dLh*or2(nR)*dcheb(nCheb,nR)&
                    &  -alpha*hdif_vel*visc(nR)*dLh*or2(nR)      *(- d3cheb(nCheb,nR) &
                    &                   +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
                    &      +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
                    &       +two*(dLvisc(nR)+beta(nR))*or1(nR) )*dcheb(nCheb,nR)  &
                    &      -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
                    &                     +two*third*beta(nR)   )* cheb(nCheb,nR) &
                    &                                        ) )
       
               wptMat(nR_p,nCheb_p)= -cheb_norm*alpha*dLh*or2(nR)*cheb(nCheb,nR)

               wptMat(nR_p,nCheb_t)=0.0_cp

               ! T equation
               wptMat(nR_t,nCheb_t)= cheb_norm * (                      &
               &                               O_dt*cheb(nCheb,nR) -    &
               &      alpha*opr*kappa(nR)*hdif_t*(  d2cheb(nCheb,nR) +  &
               &      ( beta(nR)-dLtemp0(nR)+two*or1(nR)+dLkappa(nR) )* &
               &                              dcheb(nCheb,nR) -         &
               &      (ddLtemp0(nR)+dLtemp0(nR)*(dLkappa(nR)+beta(nR)   &
               &       +two*or1(nR))+dLh*or2(nR))*                      &
               &                               cheb(nCheb,nR) ) )

               wptMat(nR_t,nCheb_p)= -cheb_norm * ViscHeatFac*ThExpNb*  &
               &                  alpha0(nR)*temp0(nR)*orho1(nR)* (     &
               &                               O_dt*cheb(nCheb,nR) -    &
               &      alpha*opr*kappa(nR)*hdif_t*(  d2cheb(nCheb,nR) +  &
               &      ( dLkappa(nR)+dLtemp0(nR)+two*or1(nR)+            &
               &        two*dLalpha0(nR)-beta(nR) ) *                   &
               &                              dcheb(nCheb,nR) +         &
               &      ((dLalpha0(nR)-beta(nR))*( two*or1(nR)+           &
               &        dLalpha0(nR)+dLkappa(nR)+dLtemp0(nR) )          &
               &        +ddLalpha0(nR)-dbeta(nR)-dLh*or2(nR) ) *        &
               &                               cheb(nCheb,nR) ) )

               !Advection of the background entropy u_r * dso/dr
               wptMat(nR_t,nCheb)=cheb_norm*alpha*dLh*or2(nR)*temp0(nR)* &
               &                  dentropy0(nR)*orho1(nR)*cheb(nCheb,nR)

            end do
         end do
      end if
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         nR_t=nR+2*n_r_max
         wptMat(nR,1)            =half*wptMat(nR,1)
         wptMat(nR,n_r_max)      =half*wptMat(nR,n_r_max)
         wptMat(nR,n_r_max+1)    =half*wptMat(nR,n_r_max+1)
         wptMat(nR,2*n_r_max)    =half*wptMat(nR,2*n_r_max)
         wptMat(nR,2*n_r_max+1)  =half*wptMat(nR,2*n_r_max+1)
         wptMat(nR,3*n_r_max)    =half*wptMat(nR,3*n_r_max)
         wptMat(nR_p,1)          =half*wptMat(nR_p,1)
         wptMat(nR_p,n_r_max)    =half*wptMat(nR_p,n_r_max)
         wptMat(nR_p,n_r_max+1)  =half*wptMat(nR_p,n_r_max+1)
         wptMat(nR_p,2*n_r_max)  =half*wptMat(nR_p,2*n_r_max)
         wptMat(nR_p,2*n_r_max+1)=half*wptMat(nR_p,2*n_r_max+1)
         wptMat(nR_p,3*n_r_max)  =half*wptMat(nR_p,3*n_r_max)
         wptMat(nR_t,1)          =half*wptMat(nR_t,1)
         wptMat(nR_t,n_r_max)    =half*wptMat(nR_t,n_r_max)
         wptMat(nR_t,n_r_max+1)  =half*wptMat(nR_t,n_r_max+1)
         wptMat(nR_t,2*n_r_max)  =half*wptMat(nR_t,2*n_r_max)
         wptMat(nR_t,2*n_r_max+1)=half*wptMat(nR_t,2*n_r_max+1)
         wptMat(nR_t,3*n_r_max)  =half*wptMat(nR_t,3*n_r_max)
      end do
    
      ! compute the linesum of each line
      do nR=1,3*n_r_max
         wptMat_fac(nR,1)=one/maxval(abs(wptMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,3*n_r_max
         wptMat(nR,:) = wptMat(nR,:)*wptMat_fac(nR,1)
      end do
    
      ! also compute the rowsum of each column
      do nR=1,3*n_r_max
         wptMat_fac(nR,2)=one/maxval(abs(wptMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,3*n_r_max
         wptMat(:,nR) = wptMat(:,nR)*wptMat_fac(nR,2)
      end do

      call sgefa(wptMat,3*n_r_max,3*n_r_max,wptPivot,info)
      if ( info /= 0 ) then
         write(*,*) 'Singular matrix wptMat!'
         stop '35'
      end if

   end subroutine get_wptMat
!-----------------------------------------------------------------------------
   subroutine get_pt0Mat(dt,ptMat,ptPivot,ptMat_fac)

      !-- Input variables
      real(cp), intent(in) :: dt

      !-- Output variables:
      real(cp), intent(out) :: ptMat(2*n_r_max,2*n_r_max)
      integer,  intent(out) :: ptPivot(2*n_r_max)
      real(cp), intent(out) :: ptMat_fac(2*n_r_max)

      !-- Local variables:
      integer :: info,nCheb,nCheb_p,nR,nR_p,n_cheb_in
      real(cp) :: work(n_r_max),work2(n_r_max),work3(n_r_max)
      real(cp) :: O_dt

      O_dt=one/dt

      if ( l_temperature_diff ) then ! temperature diffusion

         do nCheb=1,n_r_max
            nCheb_p=nCheb+n_r_max
            do nR=1,n_r_max
               nR_p=nR+n_r_max

               ptMat(nR,nCheb)= cheb_norm * (                           &
               &                          O_dt*cheb(nCheb,nR) -         &
               &      alpha*opr*kappa(nR)*(  d2cheb(nCheb,nR) +         &
               &      ( beta(nR)+ two*or1(nR)+dLkappa(nR) )*            &
               &                              dcheb(nCheb,nR) ) )

               ptMat(nR,nCheb_p)= -cheb_norm * ViscHeatFac*ThExpNb*     &
               &                  alpha0(nR)*temp0(nR)*orho1(nR)*       &
               &                               O_dt*cheb(nCheb,nR)

               ptMat(nR_p,nCheb)  = -cheb_norm*rho0(nR)*alpha0(nR)*  &
               &                     BuoFac*rgrav(nR)*cheb(nCheb,nR)
               ptMat(nR_p,nCheb_p)= cheb_norm*(          dcheb(nCheb,nR) &
               &                   +BuoFac*ViscHeatFac*(                 &
               &                   ThExpNb*alpha0(nR)*temp0(nR)+ogrun )* &
               &                   alpha0(nR)*rgrav(nR)*  cheb(nCheb,nR))
            end do
         end do

      else ! entropy diffusion

         do nCheb=1,n_r_max
           nCheb_p=nCheb+n_r_max
            do nR=1,n_r_max
               nR_p=nR+n_r_max

               ptMat(nR,nCheb)= cheb_norm * (                           &
               &                          O_dt*cheb(nCheb,nR) -         &
               &      alpha*opr*kappa(nR)*(  d2cheb(nCheb,nR) +         &
               &      ( beta(nR)-dLtemp0(nR)+two*or1(nR)+dLkappa(nR) )* &
               &                              dcheb(nCheb,nR) -         &
               &      (ddLtemp0(nR)+dLtemp0(nR)*(dLkappa(nR)+beta(nR)   &
               &          +two*or1(nR)))*     cheb(nCheb,nR) ) )

               ptMat(nR,nCheb_p)= -cheb_norm * ViscHeatFac*ThExpNb*     &
               &                  alpha0(nR)*temp0(nR)*orho1(nR)* (     &
               &                          O_dt*cheb(nCheb,nR) -         &
               &      alpha*opr*kappa(nR)*(  d2cheb(nCheb,nR) +         &
               &      ( dLkappa(nR)+dLtemp0(nR)+two*or1(nR)+            &
               &        two*dLalpha0(nR)-beta(nR) ) *                   &
               &                              dcheb(nCheb,nR) +         &
               &      ((dLalpha0(nR)-beta(nR))*( two*or1(nR)+           &
               &        dLalpha0(nR)+dLkappa(nR)+dLtemp0(nR) )          &
               &                    +ddLalpha0(nR)-dbeta(nR) ) *        &
               &                               cheb(nCheb,nR)))

               ptMat(nR_p,nCheb)  = -cheb_norm*rho0(nR)*alpha0(nR)*&
               &                     BuoFac*rgrav(nR)*cheb(nCheb,nR)
               ptMat(nR_p,nCheb_p)= cheb_norm*(         dcheb(nCheb,nR)  &
               &                   +BuoFac*ViscHeatFac*(                 &
               &                   ThExpNb*alpha0(nR)*temp0(nR)+ogrun )* &
               &                   alpha0(nR)*rgrav(nR)*  cheb(nCheb,nR))
            end do
         end do

      end if


      !----- Boundary condition:
      do nCheb=1,n_cheb_max
         nCheb_p=nCheb+n_r_max

         if ( ktops == 1 ) then
            !--------- Constant entropy at outer boundary:
            ptMat(1,nCheb)  =otemp1(1)*cheb_norm
            ptMat(1,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*alpha0(1)*&
            &                orho1(1)
         else if ( ktops == 2) then
            !--------- Constant entropy flux at outer boundary:
            ptMat(1,nCheb)  =cheb_norm*otemp1(1)*( dcheb(nCheb,1)- &
            &                         dLtemp0(1)*   cheb(nCheb,1) )
            ptMat(1,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*alpha0(1)*&
            &                 orho1(1)*(            dcheb(nCheb,1)+   &
            &                 (dLalpha0(1)-beta(1))* cheb(nCheb,1) )
         else if ( ktops == 3) then
            !--------- Constant temperature at outer boundary:
            ptMat(1,nCheb)  =cheb_norm
            ptMat(1,nCheb_p)=0.0_cp
         else if ( ktops == 4) then
            !--------- Constant temperature flux at outer boundary:
            ptMat(1,nCheb)  =cheb_norm*dcheb(nCheb,1)
            ptMat(1,nCheb_p)=0.0_cp
         end if

         if ( kbots == 1 ) then
            !--------- Constant entropy at inner boundary:
            ptMat(n_r_max,nCheb)  =cheb(nCheb,n_r_max)*cheb_norm*otemp1(n_r_max)
            ptMat(n_r_max,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*  &
            &                       alpha0(n_r_max)*orho1(n_r_max)* &
            &                       cheb(nCheb,n_r_max)
         else if ( kbots == 2) then
            !--------- Constant entropy flux at inner boundary:
            ptMat(n_r_max,nCheb)  =cheb_norm*otemp1(n_r_max)*(            &
            &                                       dcheb(nCheb,n_r_max)- &
            &                       dLtemp0(n_r_max)*cheb(nCheb,n_r_max) )
            ptMat(n_r_max,nCheb_p)=-cheb_norm*ViscHeatFac*ThExpNb*         &
            &                       alpha0(n_r_max)*orho1(n_r_max)*(       &
            &                                       dcheb(nCheb,n_r_max)+  &
            &                       (dLalpha0(n_r_max)-beta(n_r_max))*     &
            &                                        cheb(nCheb,n_r_max) )
         else if ( kbots == 3) then
            !--------- Constant temperature at inner boundary:
            ptMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
            ptMat(n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 4) then
            !--------- Constant temperature flux at inner boundary:
            ptMat(n_r_max,nCheb)  =cheb_norm*dcheb(nCheb,n_r_max)
            ptMat(n_r_max,nCheb_p)=0.0_cp
         end if

         ptMat(2*n_r_max,nCheb)  =0.0_cp
         ptMat(2*n_r_max,nCheb_p)=0.0_cp
      end do

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes
      if ( ViscHeatFac*ThExpNb /= 0.0_cp .and. ktopp==1 ) then

         work(:)=ViscHeatFac*alpha0(:)*(ThExpNb*alpha0(:)*temp0(:)+ogrun)*r(:)*r(:)
         call chebt_oc%costf1(work,work2)
         work         =work*cheb_norm
         work(1)      =half*work(1)
         work(n_r_max)=half*work(n_r_max)

         work2(:)=-alpha0(:)*rho0(:)*r(:)*r(:)
         call chebt_oc%costf1(work2,work3)
         work2         =work2*cheb_norm
         work2(1)      =half*work2(1)
         work2(n_r_max)=half*work2(n_r_max)

         do nCheb=1,n_cheb_max
            nCheb_p=nCheb+n_r_max
            ptMat(n_r_max+1,nCheb_p)=0.0_cp
            ptMat(n_r_max+1,nCheb)  =0.0_cp
            do n_cheb_in=1,n_cheb_max
               if (mod(nCheb+n_cheb_in-2,2)==0) then
                  ptMat(n_r_max+1,nCheb_p)=ptMat(n_r_max+1,nCheb_p)+            &
                  &                    (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                  &                    one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                  &                       work(n_cheb_in)*half*cheb_norm
                  ptMat(n_r_max+1,nCheb)  =ptMat(n_r_max+1,nCheb)+              &
                  &                    (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                  &                    one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                  &                    work2(n_cheb_in)*half*cheb_norm
               end if
            end do
         end do

      else

         do nCheb=1,n_cheb_max
            nCheb_p=nCheb+n_r_max
            ptMat(n_r_max+1,nCheb)  =0.0_cp
            ptMat(n_r_max+1,nCheb_p)=cheb_norm
         end do

      end if


      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            ptMat(1,nCheb)          =0.0_cp
            ptMat(n_r_max,nCheb)    =0.0_cp
            ptMat(n_r_max+1,nCheb)  =0.0_cp
            ptMat(2*n_r_max,nCheb)  =0.0_cp
            ptMat(1,nCheb_p)        =0.0_cp
            ptMat(n_r_max,nCheb_p)  =0.0_cp
            ptMat(n_r_max+1,nCheb_p)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         ptMat(nR,1)          =half*ptMat(nR,1)
         ptMat(nR,n_r_max)    =half*ptMat(nR,n_r_max)
         ptMat(nR,n_r_max+1)  =half*ptMat(nR,n_r_max+1)
         ptMat(nR,2*n_r_max)  =half*ptMat(nR,2*n_r_max)
         ptMat(nR_p,1)        =half*ptMat(nR_p,1)
         ptMat(nR_p,n_r_max)  =half*ptMat(nR_p,n_r_max)
         ptMat(nR_p,n_r_max+1)=half*ptMat(nR_p,n_r_max+1)
         ptMat(nR_p,2*n_r_max)=half*ptMat(nR_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         ptMat_fac(nR)=one/maxval(abs(ptMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         ptMat(nR,:) = ptMat(nR,:)*ptMat_fac(nR)
      end do

      !---- LU decomposition:
      call sgefa(ptMat,2*n_r_max,2*n_r_max,ptPivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular matrix pt0Mat!'
         stop '29'
      end if

   end subroutine get_pt0Mat
!-----------------------------------------------------------------------------
end module updateWPT_mod
