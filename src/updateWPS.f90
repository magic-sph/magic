#include "perflib_preproc.cpp"
module updateWPS_mod
   
   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_cheb_max, n_r_max, l_max
   use radial_data, only: n_r_cmb,n_r_icb
   use radial_functions, only: drx,ddrx,dddrx,or1,or2,rho0,rgrav,       &
       &                       chebt_oc,visc,dLvisc,                    &
       &                       beta,dbeta,cheb,dcheb,d2cheb,d3cheb,     &
       &                       cheb_norm, dLkappa, dLtemp0, ddLtemp0,   &
       &                       alpha0,dLalpha0, ddLalpha0,              &
       &                       kappa, orho1, dentropy0, temp0, r
   use physical_parameters, only: kbotv, ktopv, ktops, kbots, ra, opr, &
       &                          ViscHeatFac, ThExpNb, ogrun, BuoFac, &
       &                          CorFac
   use num_param, only: alpha
   use init_fields, only: tops, bots
   use blocking, only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map, &
       &               lmStartB,lmStopB
   use horizontal_data, only: hdif_V, hdif_S, dLh
   use logic, only: l_update_v, l_temperature_diff, l_RMS
   use RMS, only: DifPol2hInt, dtVPolLMr, dtVPol2hInt, DifPolLMr
   use RMS_helpers, only:  hInt2Pol
   use algebra, only: cgeslML, sgefa, sgesl
   use LMLoop_data, only: llm, ulm
   use communications, only: get_global_sum
   use parallel_mod, only: chunksize, rank
   use cosine_transform_odd
   use radial_der, only: get_dddr, get_ddr, get_dr, get_drNS
   use integration, only: rInt_R
   use constants, only: zero, one, two, three, four, third, half, pi, osq4pi

   implicit none

   private

   !-- Input of recycled work arrays:
   complex(cp), allocatable :: workA(:,:),workB(:,:), workC(:,:), workD(:,:)
   complex(cp), allocatable :: Dif(:),Pre(:),Buo(:),dtV(:)
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: ps0Mat(:,:), ps0Mat_fac(:)
   integer, allocatable :: ps0Pivot(:)
   real(cp), allocatable :: wpsMat(:,:,:)
   integer, allocatable :: wpsPivot(:,:)
   real(cp), allocatable :: wpsMat_fac(:,:,:)
   real(cp) :: Cor00_fac
   logical, public, allocatable :: lWPSmat(:)

   integer :: maxThreads

   public :: initialize_updateWPS, finalize_updateWPS, updateWPS

contains

   subroutine initialize_updateWPS

      allocate( ps0Mat(2*n_r_max,2*n_r_max) )
      allocate( ps0Mat_fac(2*n_r_max) )
      allocate( ps0Pivot(2*n_r_max) )
      bytes_allocated = bytes_allocated+(4*n_r_max+2)*n_r_max*SIZEOF_DEF_REAL &
      &                 +2*n_r_max*SIZEOF_INTEGER
      allocate( wpsMat(3*n_r_max,3*n_r_max,l_max) )
      allocate(wpsMat_fac(3*n_r_max,2,l_max))
      allocate ( wpsPivot(3*n_r_max,l_max) )
      bytes_allocated = bytes_allocated+(9*n_r_max*l_max+6*n_r_max*l_max)*&
      &                 SIZEOF_DEF_REAL+3*n_r_max*l_max*SIZEOF_INTEGER
      allocate( lWPSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

      allocate( workA(llm:ulm,n_r_max) )
      allocate( workB(llm:ulm,n_r_max) )
      allocate( workC(llm:ulm,n_r_max) )
      bytes_allocated = bytes_allocated+3*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

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

      if ( l_RMS ) then
         allocate( workD(llm:ulm, n_r_max) )
         bytes_allocated=bytes_allocated+n_r_max*(ulm-llm+1)*SIZEOF_DEF_REAL
      end if

      Cor00_fac=four/sqrt(three)

   end subroutine initialize_updateWPS
!-----------------------------------------------------------------------------
   subroutine finalize_updateWPS

      deallocate( ps0Mat, ps0Mat_fac, ps0Pivot )
      deallocate( wpsMat, wpsMat_fac, wpsPivot, lWPSmat )
      deallocate( workA, workB, workC, rhs1)
      deallocate( Dif, Pre, Buo, dtV )
      if ( l_RMS ) deallocate( workD )

   end subroutine finalize_updateWPS
!-----------------------------------------------------------------------------
   subroutine updateWPS(w,dw,ddw,z10,dwdt,dwdtLast,p,dp,dpdt,dpdtLast,s, &
              &         ds,dVSrLM,dsdt,dsdtLast,w1,coex,dt,nLMB,lRmsNext)
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

      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dwdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dpdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: ds(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdtLast(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

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

      !allocate(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,nLMBs2(nLMB)))

      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)

      w2  =one-w1
      O_dt=one/dt

      !$OMP PARALLEL  &
      !$OMP private(iThread,start_lm,stop_lm,nR,lm) &
      !$OMP shared(all_lms,per_thread,lmStart,lmStop) &
      !$OMP shared(dVSrLM,chebt_oc,drx,dsdt,orho1,or2) &
      !$OMP shared(n_r_max,n_cheb_max,workA,workB,nThreads,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
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
         call get_drNS( dVSrLM,workA,ulm-llm+1,start_lm-llm+1,  &
              &         stop_lm-llm+1,n_r_max,n_cheb_max,workB, &
              &         chebt_oc,drx)
      end do
      !$OMP end do

      !$OMP DO
      do nR=1,n_r_max
         do lm=lmStart,lmStop
            dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*workA(lm,nR))
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
            if ( .not. lWPSmat(l1) ) then
               call get_ps0Mat(dt,ps0Mat,ps0Pivot,ps0Mat_fac)
               lWPSmat(l1)=.true.
            end if
         else
            if ( .not. lWPSmat(l1) ) then
               call get_wpsMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                    &          hdif_S(st_map%lm2(l1,0)),       &
                    &          wpsMat(1,1,l1),wpsPivot(1,l1),  &
                    &          wpsMat_fac(1,1,l1))
               lWPSmat(l1)=.true.
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
                     rhs(nR)        =real(s(lm1,nR))*O_dt+ &
                     &             w1*real(dsdt(lm1,nR)) + &
                     &             w2*real(dsdtLast(lm1,nR))
                     !rhs(nR+n_r_max)=w1*real(dwdt(lm1,nR))+Cor00_fac*&
                     !&               alpha*CorFac*or1(nR)*z10(nR)+   &
                     !&               w2*real(dwdtLast(lm1,nR))
                     rhs(nR+n_r_max)=real(dwdt(lm1,nR))+Cor00_fac*&
                     &               CorFac*or1(nR)*z10(nR)
                  end do
                  rhs(1)        =real(tops(0,0))
                  rhs(n_r_max)  =real(bots(0,0))
                  rhs(n_r_max+1)=0.0_cp

                  rhs = ps0Mat_fac*rhs

                  call sgesl(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,rhs)

               else ! l1 /= 0
                  lmB=lmB+1
                  rhs1(1,lmB,threadid)          =0.0_cp
                  rhs1(n_r_max,lmB,threadid)    =0.0_cp
                  rhs1(n_r_max+1,lmB,threadid)  =0.0_cp
                  rhs1(2*n_r_max,lmB,threadid)  =0.0_cp
                  rhs1(2*n_r_max+1,lmB,threadid)=tops(l1,m1)
                  rhs1(3*n_r_max,lmB,threadid)  =bots(l1,m1)
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=O_dt*dLh(st_map%lm2(l1,m1))* &
                     &                     or2(nR)*w(lm1,nR) +          &
                     &                     w1*dwdt(lm1,nR) +            &
                     &                     w2*dwdtLast(lm1,nR)
                     rhs1(nR+n_r_max,lmB,threadid)=-O_dt*                 &
                     &                             dLh(st_map%lm2(l1,m1))*&
                     &                             or2(nR)*dw(lm1,nR) +   &
                     &                             w1*dpdt(lm1,nR) +      &
                     &                             w2*dpdtLast(lm1,nR)
                     rhs1(nR+2*n_r_max,lmB,threadid)=s(lm1,nR)*O_dt +  &
                     &                               w1*dsdt(lm1,nR) + &
                     &                               w2*dsdtLast(lm1,nR)
                  end do
               end if
            end do
            !PERFOFF

            !PERFON('upWP_sol')
            if ( lmB > 0 ) then

               ! use the mat_fac(:,1) to scale the rhs
               do lm=lmB0+1,lmB
                  do nR=1,3*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpsMat_fac(nR,1,l1)
                  end do
               end do
               call cgeslML(wpsMat(:,:,l1),3*n_r_max,3*n_r_max,        &
                    &       wpsPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),&
                    &       3*n_r_max,lmB-lmB0)
               ! rescale the solution with mat_fac(:,2)
               do lm=lmB0+1,lmB
                  do nR=1,3*n_r_max
                     rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpsMat_fac(nR,2,l1)
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
                     s(lm1,n_cheb)=rhs(n_cheb)
                     p(lm1,n_cheb)=rhs(n_cheb+n_r_max)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                        p(lm1,n_cheb)=rhs1(n_r_max+n_cheb,lmB,threadid)
                        s(lm1,n_cheb)=rhs1(2*n_r_max+n_cheb,lmB,threadid)
                     end do
                  else
                     do n_cheb=1,n_cheb_max
                        w(lm1,n_cheb)= cmplx(real(rhs1(n_cheb,lmB,threadid)), &
                                       &     0.0_cp,kind=cp)
                        p(lm1,n_cheb)= cmplx(real(rhs1(n_r_max+n_cheb,lmB,threadid)), &
                                       &     0.0_cp,kind=cp)
                        s(lm1,n_cheb)= cmplx(real(rhs1(2*n_r_max+n_cheb,lmB,threadid)), &
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
            s(lm1,n_cheb)=zero
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
      !$OMP shared(w,dw,ddw,p,dp,s,ds,dwdtLast,dpdtLast,dsdtLast) &
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

         call get_ddr( p, dp, workC,ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max,n_cheb_max,dwdtLast,dpdtLast,            &
              &       chebt_oc,drx,ddrx)

         call chebt_oc%costf1(s,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,dsdtLast)
         call get_ddr(s, ds, workB, ulm-llm+1, start_lm-llm+1, stop_lm-llm+1, &
              &       n_r_max, n_cheb_max, dwdtLast, dsdtLast,                &
              &       chebt_oc,drx,ddrx)

      end do
      !$OMP end do
      !$OMP END PARALLEL


      !do nR=1,n_r_max
      !   rhoprime(nR)=osq4pi*ThExpNb*alpha0(nR)*( -rho0(nR)*temp0(nR)* &
      !                real(s(1,nR))+ViscHeatFac*ogrun*real(p(1,nR)) )
      !end do
      !mass = rInt_R(rhoprime*r**2,n_r_max,n_r_max,drx,chebt_oc)

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
                          or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                    &   +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
                    &   -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
                    &        dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                    &        +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
                    &                                            w(lm1,nR)  )
               Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
               Buo(lm1) = BuoFac*rho0(nR)*rgrav(nR)*s(lm1,nR)
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
               dsdtLast(lm1,nR)=dsdt(lm1,nR) &
                    & - coex*opr*hdif_S(st_map%lm2(l1,m1)) * kappa(nR) *    &
                    &   (             workB(lm1,nR)                         &
                    &     + ( beta(nR)+two*dLtemp0(nR)+two*or1(nR)+         &
                    &    dLkappa(nR) ) * ds(lm1,nR) +                       &
                    &   ( ddLtemp0(nR)+ dLtemp0(nR)*(                       &
                    &      two*or1(nR)+dLkappa(nR)+dLtemp0(nR)+beta(nR)) -  &
                    &       dLh(st_map%lm2(l1,m1))*or2(nR) ) *  &
                    &                     s(lm1,nR)  +                      &
                    &   alpha0(nR)*orho1(nR)*ViscHeatFac*ThExpNb*(          &
                    &                 workC(lm1,nR)   +                     &
                    &     ( dLkappa(nR)+two*(dLtemp0(nR)+dLalpha0(nR))  +   &
                    &       two*or1(nR)-beta(nR) ) * dp(lm1,nR) +           &
                    & ( (dLkappa(nR)+dLtemp0(nR)+dLalpha0(nR)+two*or1(nR))* &
                    &     (dLalpha0(nR)+dLtemp0(nR)-beta(nR)) +             &
                    &     ddLtemp0(nR)+ddLalpha0(nR)-dbeta(nR) -            &
                    &     dLh(st_map%lm2(l1,m1)) * or2(nR)                  &
                    &   )*                p(lm1,nR) ) ) +                   &
                    &   coex*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)   &
                    &   *orho1(nR)*dentropy0(nR)*w(lm1,nR)
               if ( lRmsNext ) then
                  dtV(lm1)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                  &        ( w(lm1,nR)-workD(lm1,nR) )
               end if
               !if ( l1 == 0 ) then
               !   dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1) + &
               !   &                 Cor00_fac*CorFac*or1(nR)*z10(nR))
               !end if
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
                          or2(nR)*visc(nR) *                  ( ddw(lm1,nR) &
                    &   +(two*dLvisc(nR)-third*beta(nR))*        dw(lm1,nR) &
                    &   -( dLh(st_map%lm2(l1,m1))*or2(nR)+four*third* (     &
                    &        dbeta(nR)+dLvisc(nR)*beta(nR)                  &
                    &        +(three*dLvisc(nR)+beta(nR))*or1(nR) )   )*    &
                    &                                            w(lm1,nR)  )
               Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
               Buo(lm1) = BuoFac*rho0(nR)*rgrav(nR)*s(lm1,nR)
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
               dsdtLast(lm1,nR)=dsdt(lm1,nR) &
                    & - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1)))* &
                    &                                         kappa(nR) * &
                    &   ( workB(lm1,nR)                                   &
                    &     + ( beta(nR) + dLtemp0(nR) +                    &
                    &       two*or1(nR) + dLkappa(nR) ) * ds(lm1,nR)      &
                    &     - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)  &
                    &     *  s(lm1,nR)    ) +                             &
                    &   coex*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR) &
                    &   *orho1(nR)*dentropy0(nR)*w(lm1,nR)

               !if ( l1 == 0 ) then
               !   dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1) + &
               !   &                 Cor00_fac*CorFac*or1(nR)*z10(nR))
               !end if
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

   end subroutine updateWPS
   !------------------------------------------------------------------------------
   subroutine get_wpsMat(dt,l,hdif_vel,hdif_s,wpsMat,wpsPivot,wpsMat_fac)
      !
      !  Purpose of this subroutine is to contruct the time step matrix  
      !  wpmat  for the NS equation.                                    
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif_vel
      real(cp), intent(in) :: hdif_s
      integer,  intent(in) :: l

      !-- Output variables:
      real(cp), intent(out) :: wpsMat(3*n_r_max,3*n_r_max)
      real(cp), intent(out) :: wpsMat_fac(3*n_r_max,2)
      integer,  intent(out) :: wpsPivot(3*n_r_max)

      !-- local variables:
      integer :: nR,nCheb,nR_p,nR_s,nCheb_p,nCheb_s
      integer :: info
      real(cp) :: O_dt,dLh

      O_dt=one/dt
      dLh =real(l*(l+1),kind=cp)
    
      !-- Now mode l>0
    
      !----- Boundary conditions, see above:
      do nCheb=1,n_cheb_max
         nCheb_p=nCheb+n_r_max
         nCheb_s=nCheb+2*n_r_max
    
         wpsMat(1,nCheb)        =cheb_norm*cheb(nCheb,1)
         wpsMat(1,nCheb_p)      =0.0_cp
         wpsMat(1,nCheb_s)      =0.0_cp
         wpsMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)
         wpsMat(n_r_max,nCheb_p)=0.0_cp
         wpsMat(n_r_max,nCheb_s)=0.0_cp
    
         if ( ktopv == 1 ) then  ! free slip !
            wpsMat(n_r_max+1,nCheb)=   cheb_norm * ( &
                 d2cheb(nCheb,1) - (two*or1(1)+beta(1))*dcheb(nCheb,1) )
         else                    ! no slip, note exception for l=1,m=0
            wpsMat(n_r_max+1,nCheb)=cheb_norm*dcheb(nCheb,1)
         end if
         wpsMat(n_r_max+1,nCheb_p)=0.0_cp
         wpsMat(n_r_max+1,nCheb_s)=0.0_cp

         if ( kbotv == 1 ) then  ! free slip !
            wpsMat(2*n_r_max,nCheb)=        cheb_norm * ( &
                 d2cheb(nCheb,n_r_max) - &
                 (two*or1(n_r_max)+beta(n_r_max))*dcheb(nCheb,n_r_max))
         else                 ! no slip, note exception for l=1,m=0
            wpsMat(2*n_r_max,nCheb)=cheb_norm * dcheb(nCheb,n_r_max)
         end if
         wpsMat(2*n_r_max,nCheb_p)=0.0_cp
         wpsMat(2*n_r_max,nCheb_s)=0.0_cp

         if ( ktops == 1 ) then ! fixed entropy
            wpsMat(2*n_r_max+1,nCheb_s)=cheb_norm
            wpsMat(2*n_r_max+1,nCheb_p)=0.0_cp
         else if ( ktops == 2 ) then ! fixed entropy flux
            wpsMat(2*n_r_max+1,nCheb_s)=cheb_norm*dcheb(nCheb,1)
            wpsMat(2*n_r_max+1,nCheb_p)=0.0_cp
         else if ( ktops == 3 ) then ! fixed temperature
            wpsMat(2*n_r_max+1,nCheb_s)=cheb_norm*temp0(1)
            wpsMat(2*n_r_max+1,nCheb_p)=cheb_norm*orho1(1)*alpha0(1)* &
                                        temp0(1)*ViscHeatFac*ThExpNb
         else if ( ktops == 4 ) then ! fixed temperature flux
            wpsMat(2*n_r_max+1,nCheb_s)=cheb_norm*temp0(1)*( dcheb(nCheb,1)+ &
              &                                    dLtemp0(1)*cheb(nCheb,1) )
            wpsMat(2*n_r_max+1,nCheb_p)=cheb_norm*orho1(1)*alpha0(1)*   &
              &                         temp0(1)*ViscHeatFac*ThExpNb*(  &
              &                         dcheb(nCheb,1)+(dLalpha0(1)+    &
              &                         dLtemp0(1)-beta(1))*cheb(nCheb,1) )
         end if
         wpsMat(2*n_r_max+1,nCheb)  =0.0_cp

         if ( kbots == 1 ) then ! fixed entropy
            wpsMat(3*n_r_max,nCheb_s)=cheb_norm*cheb(nCheb,n_r_max)
            wpsMat(3*n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 2) then ! fixed entropy flux
            wpsMat(3*n_r_max,nCheb_s)=cheb_norm*dcheb(nCheb,n_r_max)
            wpsMat(3*n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 3) then ! fixed temperature
            wpsMat(3*n_r_max,nCheb_s)=cheb_norm*cheb(nCheb,n_r_max)*temp0(n_r_max)
            wpsMat(3*n_r_max,nCheb_p)=cheb_norm*cheb(nCheb,n_r_max)* &
              &                      orho1(n_r_max)*alpha0(n_r_max)* &
              &                      temp0(n_r_max)*ViscHeatFac*ThExpNb
         else if ( kbots == 4) then ! fixed temperature flux
            wpsMat(3*n_r_max,nCheb_s)=cheb_norm*temp0(n_r_max)*(            &
              &                                       dcheb(nCheb,n_r_max)+ &
              &                       dLtemp0(n_r_max)*cheb(nCheb,n_r_max) )
            wpsMat(3*n_r_max,nCheb_p)=cheb_norm*orho1(n_r_max)*alpha0(n_r_max)* &
              &                         temp0(n_r_max)*ViscHeatFac*ThExpNb*(    &
              &                                       dcheb(nCheb,n_r_max)+     &
              &                         (dLalpha0(n_r_max)+dLtemp0(n_r_max)-    &
              &                         beta(n_r_max))*cheb(nCheb,n_r_max) )
         end if
         wpsMat(3*n_r_max,nCheb)  =0.0_cp

    
      end do   !  loop over nCheb
    
      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_s=nCheb+2*n_r_max
            wpsMat(1,nCheb)            =0.0_cp
            wpsMat(n_r_max,nCheb)      =0.0_cp
            wpsMat(n_r_max+1,nCheb)    =0.0_cp
            wpsMat(2*n_r_max,nCheb)    =0.0_cp
            wpsMat(2*n_r_max+1,nCheb)  =0.0_cp
            wpsMat(3*n_r_max,nCheb)    =0.0_cp
            wpsMat(1,nCheb_p)          =0.0_cp
            wpsMat(n_r_max,nCheb_p)    =0.0_cp
            wpsMat(n_r_max+1,nCheb_p)  =0.0_cp
            wpsMat(2*n_r_max,nCheb_p)  =0.0_cp
            wpsMat(2*n_r_max+1,nCheb_p)=0.0_cp
            wpsMat(3*n_r_max,nCheb_p)  =0.0_cp
            wpsMat(1,nCheb_s)          =0.0_cp
            wpsMat(n_r_max,nCheb_s)    =0.0_cp
            wpsMat(n_r_max+1,nCheb_s)  =0.0_cp
            wpsMat(2*n_r_max,nCheb_s)  =0.0_cp
            wpsMat(2*n_r_max+1,nCheb_s)=0.0_cp
            wpsMat(3*n_r_max,nCheb_s)  =0.0_cp
         end do
      end if
    
      if ( l_temperature_diff ) then ! temperature diffusion

         do nCheb=1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_s=nCheb+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_s=nR+2*n_r_max

               ! W equation
               wpsMat(nR,nCheb)= cheb_norm *  ( O_dt*dLh*or2(nR)*cheb(nCheb,nR)  &
                    &   - alpha*hdif_vel*visc(nR)*dLh*or2(nR) * ( d2cheb(nCheb,nR)   &
                    &          +(two*dLvisc(nR)-third*beta(nR))*dcheb(nCheb,nR)  &
                    &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)    &
                    &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )  &
                    &          )                               *cheb(nCheb,nR)   &
                    &                                       )  )

               ! Buoyancy
               wpsMat(nR,nCheb_s)=-cheb_norm*alpha*BuoFac*rgrav(nR)*rho0(nR)* &
                                                      cheb(nCheb,nR)
       
               ! Pressure gradient
               wpsMat(nR,nCheb_p)= cheb_norm*alpha*(dcheb(nCheb,nR) &
                                  -beta(nR)* cheb(nCheb,nR))

               ! P equation
               wpsMat(nR_p,nCheb)= cheb_norm * ( -O_dt*dLh*or2(nR)*dcheb(nCheb,nR)&
                    &  -alpha*hdif_vel*visc(nR)*dLh*or2(nR)      *(- d3cheb(nCheb,nR) &
                    &                   +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
                    &      +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
                    &       +two*(dLvisc(nR)+beta(nR))*or1(nR) )*dcheb(nCheb,nR)  &
                    &      -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
                    &                     +two*third*beta(nR)   )* cheb(nCheb,nR) &
                    &                                        ) )
       
               wpsMat(nR_p,nCheb_p)= -cheb_norm*alpha*dLh*or2(nR)*cheb(nCheb,nR)

               wpsMat(nR_p,nCheb_s)=0.0_cp

               ! S equation
               wpsMat(nR_s,nCheb_s)= cheb_norm * (                      &
               &                               O_dt*cheb(nCheb,nR) -    &
               &      alpha*opr*hdif_s*kappa(nR)*(  d2cheb(nCheb,nR) +  &
               &      ( beta(nR)+two*dLtemp0(nR)+                       &
               &        two*or1(nR)+dLkappa(nR) )* dcheb(nCheb,nR) +    &
               &      ( ddLtemp0(nR)+dLtemp0(nR)*(                      &
               &  two*or1(nR)+dLkappa(nR)+dLtemp0(nR)+beta(nR) )   -    &   
               &           dLh*or2(nR) )*           cheb(nCheb,nR) ) )

               wpsMat(nR_s,nCheb_p)= -alpha*cheb_norm*hdif_s*kappa(nR)* &
               &      opr*alpha0(nR)*orho1(nR)*ViscHeatFac*ThExpNb*(    &
               &                                  d2cheb(nCheb,nR) +    &
               &      ( dLkappa(nR)+two*(dLalpha0(nR)+dLtemp0(nR)) -    &
               &        beta(nR) +two*or1(nR) ) *  dcheb(nCheb,nR) +    &
               & ( (dLkappa(nR)+dLalpha0(nR)+dLtemp0(nR)+two*or1(nR)) * &
               &        (dLalpha0(nR)+dLtemp0(nR)-beta(nR)) +           &
               &        ddLalpha0(nR)+ddLtemp0(nR)-dbeta(nR)-           &
               &        dLh*or2(nR) ) *             cheb(nCheb,nR) )


               !Advection of the background entropy u_r * dso/dr
               wpsMat(nR_s,nCheb)=cheb_norm*alpha*dLh*or2(nR)*dentropy0(nR)* &
                                 orho1(nR)*cheb(nCheb,nR)

            end do
         end do

      else ! entropy diffusion

         do nCheb=1,n_r_max
            nCheb_p=nCheb+n_r_max
            nCheb_s=nCheb+2*n_r_max
            do nR=2,n_r_max-1
               nR_p=nR+n_r_max
               nR_s=nR+2*n_r_max

               ! W equation
               wpsMat(nR,nCheb)= cheb_norm *  ( O_dt*dLh*or2(nR)*cheb(nCheb,nR)  &
                    &   - alpha*hdif_vel*visc(nR)*dLh*or2(nR) * ( d2cheb(nCheb,nR)   &
                    &          +(two*dLvisc(nR)-third*beta(nR))*dcheb(nCheb,nR)  &
                    &         -( dLh*or2(nR)+four*third*( dLvisc(nR)*beta(nR)    &
                    &          +(three*dLvisc(nR)+beta(nR))*or1(nR)+dbeta(nR) )  &
                    &          )                               *cheb(nCheb,nR)   &
                    &                                       )  )

               ! Buoyancy
               wpsMat(nR,nCheb_s)=-cheb_norm*alpha*BuoFac*rgrav(nR)*rho0(nR)* &
               &                                      cheb(nCheb,nR)
       
               ! Pressure gradient
               wpsMat(nR,nCheb_p)= cheb_norm*alpha*(dcheb(nCheb,nR) &
               &                  -beta(nR)* cheb(nCheb,nR))

               ! P equation
               wpsMat(nR_p,nCheb)= cheb_norm * ( -O_dt*dLh*or2(nR)*dcheb(nCheb,nR)&
                    &  -alpha*hdif_vel*visc(nR)*dLh*or2(nR)  *(- d3cheb(nCheb,nR) &
                    &                   +( beta(nR)-dLvisc(nR) )*d2cheb(nCheb,nR) &
                    &      +( dLh*or2(nR)+dbeta(nR)+dLvisc(nR)*beta(nR)           &
                    &       +two*(dLvisc(nR)+beta(nR))*or1(nR) )*dcheb(nCheb,nR)  &
                    &      -dLh*or2(nR)*( two*or1(nR)+dLvisc(nR)                  &
                    &                     +two*third*beta(nR)   )* cheb(nCheb,nR) &
                    &                                        ) )
       
               wpsMat(nR_p,nCheb_p)= -cheb_norm*alpha*dLh*or2(nR)*cheb(nCheb,nR)

               wpsMat(nR_p,nCheb_s)=0.0_cp

               ! S equation
               wpsMat(nR_s,nCheb_s)= cheb_norm * (                   &
               &                               O_dt*cheb(nCheb,nR) - &
               &     alpha*opr*hdif_s*kappa(nR)*(  d2cheb(nCheb,nR)+ &
               &      ( beta(nR)+dLtemp0(nR)+                        &
               &        two*or1(nR)+dLkappa(nR) )*dcheb(nCheb,nR) -  &
               &           dLh*or2(nR)*            cheb(nCheb,nR) ) )

               wpsMat(nR_s,nCheb_p)=0.0_cp ! temperature diffusion

               !Advection of the background entropy u_r * dso/dr
               wpsMat(nR_s,nCheb)=cheb_norm*alpha*dLh*or2(nR)*dentropy0(nR)* &
                                 orho1(nR)*cheb(nCheb,nR)

            end do
         end do
      end if
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         nR_s=nR+2*n_r_max
         wpsMat(nR,1)            =half*wpsMat(nR,1)
         wpsMat(nR,n_r_max)      =half*wpsMat(nR,n_r_max)
         wpsMat(nR,n_r_max+1)    =half*wpsMat(nR,n_r_max+1)
         wpsMat(nR,2*n_r_max)    =half*wpsMat(nR,2*n_r_max)
         wpsMat(nR,2*n_r_max+1)  =half*wpsMat(nR,2*n_r_max+1)
         wpsMat(nR,3*n_r_max)    =half*wpsMat(nR,3*n_r_max)
         wpsMat(nR_p,1)          =half*wpsMat(nR_p,1)
         wpsMat(nR_p,n_r_max)    =half*wpsMat(nR_p,n_r_max)
         wpsMat(nR_p,n_r_max+1)  =half*wpsMat(nR_p,n_r_max+1)
         wpsMat(nR_p,2*n_r_max)  =half*wpsMat(nR_p,2*n_r_max)
         wpsMat(nR_p,2*n_r_max+1)=half*wpsMat(nR_p,2*n_r_max+1)
         wpsMat(nR_p,3*n_r_max)  =half*wpsMat(nR_p,3*n_r_max)
         wpsMat(nR_s,1)          =half*wpsMat(nR_s,1)
         wpsMat(nR_s,n_r_max)    =half*wpsMat(nR_s,n_r_max)
         wpsMat(nR_s,n_r_max+1)  =half*wpsMat(nR_s,n_r_max+1)
         wpsMat(nR_s,2*n_r_max)  =half*wpsMat(nR_s,2*n_r_max)
         wpsMat(nR_s,2*n_r_max+1)=half*wpsMat(nR_s,2*n_r_max+1)
         wpsMat(nR_s,3*n_r_max)  =half*wpsMat(nR_s,3*n_r_max)
      end do
    
      ! compute the linesum of each line
      do nR=1,3*n_r_max
         wpsMat_fac(nR,1)=one/maxval(abs(wpsMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,3*n_r_max
         wpsMat(nR,:) = wpsMat(nR,:)*wpsMat_fac(nR,1)
      end do
    
      ! also compute the rowsum of each column
      do nR=1,3*n_r_max
         wpsMat_fac(nR,2)=one/maxval(abs(wpsMat(:,nR)))
      end do
      ! now divide each row by the rowsum
      do nR=1,3*n_r_max
         wpsMat(:,nR) = wpsMat(:,nR)*wpsMat_fac(nR,2)
      end do

      call sgefa(wpsMat,3*n_r_max,3*n_r_max,wpsPivot,info)
      if ( info /= 0 ) then
         write(*,*) 'Singular matrix wpsMat!'
         stop '35'
      end if

   end subroutine get_wpsMat
!-----------------------------------------------------------------------------
   subroutine get_ps0Mat(dt,psMat,psPivot,psMat_fac)

      !-- Input variables
      real(cp), intent(in) :: dt

      !-- Output variables:
      real(cp), intent(out) :: psMat(2*n_r_max,2*n_r_max)
      integer,  intent(out) :: psPivot(2*n_r_max)
      real(cp), intent(out) :: psMat_fac(2*n_r_max)

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

               psMat(nR,nCheb)= cheb_norm * (                           &
               &                               O_dt*cheb(nCheb,nR) -    &
               &      alpha*opr*kappa(nR)*(  d2cheb(nCheb,nR) +         &
               &      ( beta(nR)+two*dLtemp0(nR)+                       &
               &        two*or1(nR)+dLkappa(nR) )* dcheb(nCheb,nR) +    &
               &      ( ddLtemp0(nR)+dLtemp0(nR)*(                      &
               &  two*or1(nR)+dLkappa(nR)+dLtemp0(nR)+beta(nR) ) ) *    &   
               &                                    cheb(nCheb,nR) ) )

               psMat(nR,nCheb_p)= -alpha*cheb_norm*kappa(nR)*opr*       &
               &     alpha0(nR)*orho1(nR)*ViscHeatFac*ThExpNb*(         &
               &                                  d2cheb(nCheb,nR) +    &
               &      ( dLkappa(nR)+two*(dLalpha0(nR)+dLtemp0(nR)) -    &
               &        beta(nR) +two*or1(nR) ) *  dcheb(nCheb,nR) +    &
               & ( (dLkappa(nR)+dLalpha0(nR)+dLtemp0(nR)+two*or1(nR)) * &
               &        (dLalpha0(nR)+dLtemp0(nR)-beta(nR)) +           &
               &        ddLalpha0(nR)+ddLtemp0(nR)-dbeta(nR) ) *        &
               &                                    cheb(nCheb,nR) )

               psMat(nR_p,nCheb)  = -cheb_norm*rho0(nR)*          &
               &                     BuoFac*rgrav(nR)*cheb(nCheb,nR)
               psMat(nR_p,nCheb_p)= cheb_norm *( dcheb(nCheb,nR)- &
               &                         beta(nR)*cheb(nCheb,nR) )
            end do
         end do

      else ! entropy diffusion

         do nCheb=1,n_r_max
           nCheb_p=nCheb+n_r_max
            do nR=1,n_r_max
               nR_p=nR+n_r_max

               psMat(nR,nCheb)    = cheb_norm * (                    &
               &                               O_dt*cheb(nCheb,nR) - &
               &         alpha*opr*kappa(nR)*(    d2cheb(nCheb,nR) + &
               &    (beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR))*  &
               &                                   dcheb(nCheb,nR) ) )
               psMat(nR,nCheb_p)  =0.0_cp ! entropy diffusion

               psMat(nR_p,nCheb)  = -cheb_norm*BuoFac*rho0(nR)* &
               &                     rgrav(nR)*cheb(nCheb,nR)
               psMat(nR_p,nCheb_p)= cheb_norm*( dcheb(nCheb,nR)- &
               &                        beta(nR)*cheb(nCheb,nR) )
            end do
         end do

      end if


      !----- Boundary condition:
      do nCheb=1,n_cheb_max
         nCheb_p=nCheb+n_r_max

         if ( ktops == 1 ) then
            !--------- Constant entropy at CMB:
            psMat(1,nCheb)=cheb_norm
            psMat(1,nCheb_p)=0.0_cp
         else if ( ktops == 2) then
            !--------- Constant entropy flux at CMB:
            psMat(1,nCheb)=cheb_norm*dcheb(nCheb,1)
            psMat(1,nCheb_p)=0.0_cp
         else if ( ktops == 3) then
            !--------- Constant temperature at CMB:
            psMat(1,nCheb)  =cheb_norm*temp0(1)
            psMat(1,nCheb_p)=cheb_norm*orho1(1)*alpha0(1)*temp0(1)* &
            &                ViscHeatFac*ThExpNb
         else if ( ktops == 4) then
            !--------- Constant temperature flux at CMB:
            psMat(1,nCheb)  =cheb_norm*temp0(1)*( dcheb(nCheb,1)+ &
              &                         dLtemp0(1)*cheb(nCheb,1) )
            psMat(1,nCheb_p)=cheb_norm*orho1(1)*alpha0(1)*   &
              &              temp0(1)*ViscHeatFac*ThExpNb*(  &
              &              dcheb(nCheb,1)+(dLalpha0(1)+    &
              &              dLtemp0(1)-beta(1))*cheb(nCheb,1) )
         end if

         if ( kbots == 1 ) then
            !--------- Constant entropy at ICB:
            psMat(n_r_max,nCheb)=cheb_norm*cheb(nCheb,n_r_max)
            psMat(n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 2) then
            !--------- Constant entropy flux at ICB:
            psMat(n_r_max,nCheb)=cheb_norm*dcheb(nCheb,n_r_max)
            psMat(n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 3) then
            !--------- Constant temperature at ICB:
            psMat(n_r_max,nCheb)  =cheb_norm*cheb(nCheb,n_r_max)*temp0(n_r_max)
            psMat(n_r_max,nCheb_p)=cheb_norm*cheb(nCheb,n_r_max)*   &
              &                    alpha0(n_r_max)*temp0(n_r_max)*  &
              &                    orho1(n_r_max)*ViscHeatFac*ThExpNb
         else if ( kbots == 4) then
            !--------- Constant temperature flux at ICB:
            psMat(n_r_max,nCheb)  =cheb_norm*temp0(n_r_max)*(            &
              &                                    dcheb(nCheb,n_r_max)+ &
              &                     dLtemp0(n_r_max)*cheb(nCheb,n_r_max) )
            psMat(n_r_max,nCheb_p)=cheb_norm*orho1(n_r_max)*alpha0(n_r_max)* &
              &                      temp0(n_r_max)*ViscHeatFac*ThExpNb*(    &
              &                                    dcheb(nCheb,n_r_max)+     &
              &                      (dLalpha0(n_r_max)+dLtemp0(n_r_max)-    &
              &                       beta(n_r_max))*cheb(nCheb,n_r_max) )
         end if

         if ( ViscHeatFac*ThExpNb == 0.0_cp ) then ! No feedback of density on pressure
            psMat(n_r_max+1,nCheb_p)=cheb_norm
         else
            psMat(n_r_max+1,nCheb_p)=0.0_cp
         end if
         psMat(n_r_max+1,nCheb)  =0.0_cp
         psMat(2*n_r_max,nCheb)  =0.0_cp
         psMat(2*n_r_max,nCheb_p)=0.0_cp
      end do

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes
      if ( ViscHeatFac*ThExpNb /= 0.0_cp ) then

         work(:)=ThExpNb*ViscHeatFac*ogrun*alpha0(:)*r(:)*r(:)
         call chebt_oc%costf1(work,work2)
         work         =work*cheb_norm
         work(1)      =half*work(1)
         work(n_r_max)=half*work(n_r_max)

         work2(:)=-ThExpNb*alpha0(:)*temp0(:)*rho0(:)*r(:)*r(:)
         call chebt_oc%costf1(work2,work3)
         work2         =work2*cheb_norm
         work2(1)      =half*work2(1)
         work2(n_r_max)=half*work2(n_r_max)

         do nCheb=1,n_cheb_max
            nCheb_p=nCheb+n_r_max
            psMat(n_r_max+1,nCheb_p)=0.0_cp
            psMat(n_r_max+1,nCheb)  =0.0_cp
            do n_cheb_in=1,n_cheb_max
               if (mod(nCheb+n_cheb_in-2,2)==0) then
                  psMat(n_r_max+1,nCheb_p)=psMat(n_r_max+1,nCheb_p)+             &
                  &                     (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                  &                     one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                  &                       work(n_cheb_in)*half*cheb_norm
                  psMat(n_r_max+1,nCheb)  =psMat(n_r_max+1,nCheb)+               &
                  &                     (one/(one-real(n_cheb_in-nCheb,cp)**2)+  &
                  &                     one/(one-real(n_cheb_in+nCheb-2,cp)**2))*&
                  &                     work2(n_cheb_in)*half*cheb_norm
               end if
            end do
         end do

      end if


      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            nCheb_p=nCheb+n_r_max
            psMat(1,nCheb)          =0.0_cp
            psMat(n_r_max,nCheb)    =0.0_cp
            psMat(n_r_max+1,nCheb)  =0.0_cp
            psMat(2*n_r_max,nCheb)  =0.0_cp
            psMat(1,nCheb_p)        =0.0_cp
            psMat(n_r_max,nCheb_p)  =0.0_cp
            psMat(n_r_max+1,nCheb_p)=0.0_cp
         end do
      end if

      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         nR_p=nR+n_r_max
         psMat(nR,1)          =half*psMat(nR,1)
         psMat(nR,n_r_max)    =half*psMat(nR,n_r_max)
         psMat(nR,n_r_max+1)  =half*psMat(nR,n_r_max+1)
         psMat(nR,2*n_r_max)  =half*psMat(nR,2*n_r_max)
         psMat(nR_p,1)        =half*psMat(nR_p,1)
         psMat(nR_p,n_r_max)  =half*psMat(nR_p,n_r_max)
         psMat(nR_p,n_r_max+1)=half*psMat(nR_p,n_r_max+1)
         psMat(nR_p,2*n_r_max)=half*psMat(nR_p,2*n_r_max)
      end do

      ! compute the linesum of each line
      do nR=1,2*n_r_max
         psMat_fac(nR)=one/maxval(abs(psMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,2*n_r_max
         psMat(nR,:) = psMat(nR,:)*psMat_fac(nR)
      end do

      !---- LU decomposition:
      call sgefa(psMat,2*n_r_max,2*n_r_max,psPivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular matrix ps0Mat!'
         stop '29'
      end if

   end subroutine get_ps0Mat
!-----------------------------------------------------------------------------
end module updateWPS_mod
