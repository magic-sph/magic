#include "perflib_preproc.cpp"
module updateZ_mod

   use init_fields
   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, n_cheb_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: visc, or1, or2, cheb, dcheb, d2cheb, &
                               cheb_norm, dLvisc, beta, rho0, r_icb,&
                               r_cmb,  drx, ddrx, r, beta, dbeta,   &
                               chebt_oc
   use physical_parameters, only: kbotv, ktopv, LFfac
   use num_param, only: alpha, AMstart
   use torsional_oscillations, only: ddzASL
   use blocking, only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map, &
                     & lmStartB,lmStopB
   use horizontal_data, only: dLh, hdif_V
   use logic, only: l_rot_ma, l_rot_ic, l_SRMA, l_SRIC, l_z10mat, &
                    l_RMStest, l_correct_AMe, l_correct_AMz,      &
                    l_update_v, l_TO
   use matrices, only: z10Mat, z10Pivot, lZmat, lZ10mat, &
#ifdef WITH_PRECOND_Z
                     & zMat_fac, &
#endif
#ifdef WITH_PRECOND_Z10
                     & z10Mat_fac, &
#endif
                     & zMat, zPivot
   use RMS, only: DifTor2hInt, DifTorAs2hInt, dtVTor2hInt, dtVTorAs2hInt
   use constants, only: c_lorentz_ma, c_lorentz_ic, c_dt_z10_ma, c_dt_z10_ic, &
                    c_moi_ma, c_moi_ic, c_z10_omega_ma, c_z10_omega_ic,   &
                    c_moi_oc, y10_norm, y11_norm, zero, one, two, four,   &
                    half
   use parallel_mod
   use algebra, only: cgeslML, cgesl, sgefa
   use LMLoop_data, only: llm,ulm
   use communications, only:get_global_sum
   use outRot, only: get_angular_moment
   use RMS_helpers, only: hInt2Pol, hInt2Tor
   use radial_der, only: get_ddr
   use cosine_transform_odd
 
   implicit none
 
   private
 
   !-- Input of recycled work arrays:
   complex(cp), allocatable :: workA(:,:)  ! Work array
   complex(cp), allocatable :: workB(:,:)  ! Work array
   complex(cp), allocatable :: workC(:,:)  ! Work array
   complex(cp), allocatable :: rhs1(:,:,:) ! RHS for other modes
   integer :: maxThreads
   
   public :: updateZ,initialize_updateZ

contains

   subroutine initialize_updateZ

      allocate(workA(llm:ulm,n_r_max))
      allocate(workB(llm:ulm,n_r_max))
      allocate(workC(llm:ulm,n_r_max))

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate(rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))

   end subroutine initialize_updateZ

   subroutine updateZ(z,dz,dzdt,dzdtLast,time, &
     &             omega_ma,d_omega_ma_dtLast, &
     &             omega_ic,d_omega_ic_dtLast, &
     &             lorentz_torque_ma,lorentz_torque_maLast, &
     &             lorentz_torque_ic,lorentz_torque_icLast, &
     &             w1,coex,dt,lRmsNext)
      !
      !  updates the toroidal potential z and its radial derivatives
      !  adds explicit part to time derivatives of z
      !
    
      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: z(llm:ulm,n_r_max)        ! Toroidal velocity potential z
      complex(cp), intent(in)    :: dzdt(llm:ulm,n_r_max)     ! Time derivative of z
      complex(cp), intent(inout) :: dzdtLast(llm:ulm,n_r_max) ! Time derivative of z of previous step
      real(cp),    intent(inout) :: d_omega_ma_dtLast         ! Time derivative of OC rotation of previous step
      real(cp),    intent(inout) :: d_omega_ic_dtLast         ! Time derivative of IC rotation of previous step
      real(cp),    intent(in) :: lorentz_torque_ma            ! Lorentz torque (for OC rotation)
      real(cp),    intent(in) :: lorentz_torque_maLast        ! Lorentz torque (for OC rotation) of previous step
      real(cp),    intent(in) :: lorentz_torque_ic            ! Lorentz torque (for IC rotation)
      real(cp),    intent(in) :: lorentz_torque_icLast        ! Lorentz torque (for IC rotation) of previous step
    
      !-- Input of other variables:
      real(cp),    intent(in) :: time       ! Current time
      real(cp),    intent(in) :: w1         ! Weight for time step
      real(cp),    intent(in) :: coex       ! Factor depending on alpha
      real(cp),    intent(in) :: dt         ! Time step interval
      logical,     intent(in) :: lRmsNext   ! Logical for storing update if (l_RMS.and.l_logNext)

      !-- Output variables
      complex(cp), intent(out) :: dz(llm:ulm,n_r_max)   ! Radial derivative of z
      real(cp),    intent(out) :: omega_ma              ! Calculated OC rotation
      real(cp),    intent(out) :: omega_ic              ! Calculated IC rotation
    
      !-- local variables:
      real(cp) :: w2                  ! weight of second time step
      real(cp) :: O_dt
      real(cp) :: d_omega_ic_dt,d_omega_ma_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lm,lmB         ! position of (l,m) in array
      integer :: lmStart_00         ! excluding l=0,m=0
      integer :: lmStart,lmStop ! max and min number of orders m
      integer :: nLMB2
      integer :: nR                 ! counts radial grid points
      integer :: n_cheb             ! counts cheb modes
      complex(cp) :: rhs(n_r_max)   ! RHS of matrix multiplication
      !complex(cp) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! RHS for other modes
      complex(cp) :: z10(n_r_max),z11(n_r_max) ! toroidal flow scalar components
      real(cp) :: angular_moment(3)   ! total angular momentum
      real(cp) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
      real(cp) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
      real(cp) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
      complex(cp) :: corr_l1m0      ! correction factor for z(l=1,m=0)
      complex(cp) :: corr_l1m1      ! correction factor for z(l=1,m=1)
      real(cp) :: r_E_2             ! =r**2
      real(cp) :: nomi              ! nominator for Z10 AM correction
      integer :: l1m0,l1m1          ! position of (l=1,m=0) and (l=1,m=1) in lm.
      integer :: i                  ! counter
      logical :: l10
      integer :: nLMB
      real(cp) :: ddzASL_loc(l_max+1)
    
      complex(cp) :: Dif(lm_max)
    
      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)
    
      logical :: DEBUG_OUTPUT=.false.
      integer :: nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
      integer :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid
      complex(cp) :: rhs_sum
    
      !call mpi_barrier(MPI_COMM_WORLD,ierr)
      !write(*,"(3(A,2ES20.12))") "begin upZ: dzdt = ",get_global_sum( dzdt ),&
      !     &", z = ",get_global_sum( z ),&
      !     &", dzdtLast = ",get_global_sum( dzdtLast )
      !call mpi_barrier(MPI_COMM_WORLD,ierr)
    
      if ( .not. l_update_v ) return
    
      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m
    
    
      nLMB = 1+rank
      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)
      lmStart_00  =max(2,lmStart)
      l1m0        =lm2(1,0)
    
      w2  =one-w1
      O_dt=one/dt
    
      l10=.false.
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lmB,lm,lm1,l1,m1,n_cheb,nR) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk) &
         !$OMP private(tOmega_ma1,tOmega_ma2,threadid) &
         !$OMP private(tOmega_ic1,tOmega_ic2,rhs_sum)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)
    
         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=", &
         !     &   nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."
    
         if ( l1 /= 0 ) then
            if ( .not. lZmat(l1) ) then
#ifdef WITH_PRECOND_Z
               call get_zMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                   &        zMat(1,1,l1),zPivot(1,l1),zMat_fac(1,l1))
#else
               call get_zMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                   &        zMat(1,1,l1),zPivot(1,l1))
#endif
               lZmat(l1)=.TRUE.
            !write(*,"(A,I3,A,2ES20.12)") "zMat(",l1,") = ",SUM(zMat(:,:,l1))
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_cheb) &
            !$OMP private(tOmega_ma1,tOmega_ma2) &
            !$OMP private(tOmega_ic1,tOmega_ic2,rhs_sum) &
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
               if ( lm1 == l1m0 ) l10= .true. 
    
               if ( l_z10mat .and. lm1 == l1m0 ) then
                  !write(*,"(A,3I3)") "l_z10mat and lm1=",lm1,l1,m1
                  !PERFON('upZ_z10')
                  !----- Special treatment of z10 component if ic or mantle
                  !      are allowed to rotate about z-axis (l_z10mat=.true.) and
                  !      we use no slip boundary condition (ktopv=1,kbotv=1):
                  !      Lorentz torque is the explicit part of this time integration
                  !      at the boundaries!
                  !      Note: no angular momentum correction necessary for this case !
                  if ( .not. lZ10mat ) then
#ifdef WITH_PRECOND_Z10
                     call get_z10Mat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                                     z10Mat,z10Pivot,z10Mat_fac)
#else
                     call get_z10Mat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                                     z10Mat,z10Pivot)
#endif
                     lZ10mat=.true.
                  end if

                  if ( l_SRMA ) then
                     tOmega_ma1=time+tShift_ma1
                     tOmega_ma2=time+tShift_ma2
                     omega_ma= omega_ma1*cos(omegaOsz_ma1*tOmega_ma1) + &
                          &    omega_ma2*cos(omegaOsz_ma2*tOmega_ma2)
                     rhs(1)=omega_ma
                  else if ( ktopv == 2 .and. l_rot_ma ) then  ! time integration
                     d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_ma
                     rhs(1)=O_dt*c_dt_z10_ma*z(lm1,1) + &
                            w1*d_omega_ma_dt + w2*d_omega_ma_dtLast
                  else
                     rhs(1)=0.0_cp
                  end if
    
                  if ( l_SRIC ) then
                     tOmega_ic1=time+tShift_ic1
                     tOmega_ic2=time+tShift_ic2
                     omega_ic= omega_ic1*cos(omegaOsz_ic1*tOmega_ic1) + &
                          &    omega_ic2*cos(omegaOsz_ic2*tOmega_ic2)
                     rhs(n_r_max)=omega_ic
                  else if ( kbotv == 2 .and. l_rot_ic ) then  ! time integration
                     d_omega_ic_dt = LFfac*c_lorentz_ic*lorentz_torque_ic
                     rhs(n_r_max)=O_dt*c_dt_z10_ic*z(lm1,n_r_max) + &
                                  w1*d_omega_ic_dt + w2*d_omega_ic_dtLast
                  else
                     rhs(n_r_max)=0.0_cp
                  end if
    
                  !----- This is the normal RHS for the other radial grid points:
                  do nR=2,n_r_max-1
                     rhs(nR)=O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))* &
                             or2(nR)*z(lm1,nR)+ w1*dzdt(lm1,nR)+ w2*dzdtLast(lm1,nR)
                  end do

#ifdef WITH_PRECOND_Z10
                  do nR=1,n_r_max
                     rhs(nR) = z10Mat_fac(nR)*rhs(nR)
                  end do
#endif
                  if ( DEBUG_OUTPUT ) then
                     rhs_sum=sum(rhs)
                     write(*,"(2I3,A,2(I4,F20.16))") nLMB2,lm1,             &
                          & ":rhs_sum (z10) before = ",                     &
                          & exponent(real(rhs_sum)),fraction(real(rhs_sum)),&
                          & exponent(aimag(rhs_sum)),fraction(aimag(rhs_sum))
                     !do nR=1,n_r_max
                     !   write(*,"(3I4,A,2(I4,F20.16))")                        &
                     !        &nLMB2,lm1,nR,":rhs (z10) before = ",             &
                     !        & EXPONENT(real(rhs(nR))),FRACTION(real(rhs(nR))),&
                     !        & EXPONENT(AIMAG(rhs(nR))),FRACTION(AIMAG(rhs(nR)))
                     !end do
                  end if
                  call cgesl(z10Mat,n_r_max,n_r_max,z10Pivot,rhs)
                  if ( DEBUG_OUTPUT ) then
                     !do nR=1,n_r_max
                     !   write(*,"(3I4,A,2(I4,F20.16))")                        &
                     !        & nLMB2,lm1,nR,":rhs (z10) after = ",             &
                     !        & EXPONENT(real(rhs(nR))),FRACTION(real(rhs(nR))),&
                     !        & EXPONENT(AIMAG(rhs(nR))),FRACTION(AIMAG(rhs(nR)))
                     !end do
                     rhs_sum=sum(rhs)
                     write(*,"(2I3,A,2(I4,F20.16))") nLMB2,lm1,             &
                          & ":rhs_sum (z10) after = ",                      &
                          & exponent(real(rhs_sum)),fraction(real(rhs_sum)),&
                          & exponent(aimag(rhs_sum)),fraction(aimag(rhs_sum))
                  end if
    
    
               else if ( l1 /= 0 ) then
                  !PERFON('upZ_ln0')
                  lmB=lmB+1
                  rhs1(1,lmB,threadid)      =0.0_cp
                  rhs1(n_r_max,lmB,threadid)=0.0_cp
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)= O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))* &
                          & or2(nR)*z(lm1,nR) + w1*dzdt(lm1,nR) + w2*dzdtLast(lm1,nR)
#ifdef WITH_PRECOND_Z
                     rhs1(nR,lmB,threadid) = zMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                  end do
                  !PERFOFF
               end if
            end do

            !PERFON('upZ_sol')
            if ( lmB > lmB0 ) then
               call cgeslML(zMat(:,:,l1),n_r_max,n_r_max, &
                    &       zPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),n_r_max,lmB-lmB0)
            end if
            !PERFOFF
            if ( lRmsNext ) then ! Store old z
               do nR=1,n_r_max
                  do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                     lm1=lm22lm(lm,nLMB2,nLMB)
                     workB(lm1,nR)=z(lm1,nR)
                  end do
               end do
            end if
    
            lmB=lmB0
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               
               if ( l_z10mat .and. lm1 == l1m0 ) then
                  do n_cheb=1,n_cheb_max
                     z(lm1,n_cheb)=real(rhs(n_cheb))
                  end do
               else if ( l1 /= 0 ) then
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_cheb=1,n_cheb_max
                        z(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                     end do
                  else
                     do n_cheb=1,n_cheb_max
                        z(lm1,n_cheb)= &
                             cmplx(real(rhs1(n_cheb,lmB,threadid)),0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do       ! end of loop over lm blocks
      !$OMP END SINGLE
      !$OMP END PARALLEL
    
    
      !-- set cheb modes > n_cheb_max to zero (dealiazing)
      do n_cheb=n_cheb_max+1,n_r_max
         do lm1=lmStart,lmStop
            z(lm1,n_cheb)=zero
         end do
      end do
    
      !PERFON('upZ_drv')
      all_lms=lmStop-lmStart_00+1
#ifdef WITHOMP
      if (all_lms < omp_get_max_threads()) then
         call omp_set_num_threads(all_lms)
         per_thread=1
      else
         per_thread=all_lms/omp_get_max_threads()
      end if
#else
      per_thread=all_lms
#endif
      !$OMP PARALLEL &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(per_thread,lmStart_00,lmStop,nThreads) &
      !$OMP shared(z,dz,dzdtLast,chebt_oc,drx,ddrx) &
      !$OMP shared(n_r_max,n_cheb_max,workA,workC,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm = lmStart_00+iThread*per_thread
         stop_lm  = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop
         !write(*,"(3(A,I5))") "thread ",omp_get_thread_num()," from ",start_lm," to ",stop_lm
         !-- Get derivatives:
         call chebt_oc%costf1(z,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1,dzdtLast)
         call get_ddr(z, dz, workA, ulm-llm+1, start_lm-llm+1,     &
                      stop_lm-llm+1,n_r_max, n_cheb_max, dzdtLast, &
                      workC,chebt_oc,drx,ddrx)
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(omp_get_max_threads())
#endif
      !PERFOFF
      
      !PERFON('upZ_icma')
      !--- Update of inner core and mantle rotation:
      if ( l10 ) then
         if ( l_rot_ma .and. .not. l_SRMA ) then
            if ( ktopv == 1 ) then  ! free slip, explicit time stepping of omega !
               omega_ma=O_dt*omega_ma + LFfac/c_moi_ma * &
                    (w1*lorentz_torque_ma+w2*lorentz_torque_maLast)
               omega_ma=dt*omega_ma
            else if ( ktopv == 2 ) then ! no slip, omega given by z10
               omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
            end if
            omega_ma1=omega_ma
         end if
         if ( l_rot_ic .and. .not. l_SRIC ) then
            if ( kbotv == 1 ) then  ! free slip, explicit time stepping of omega !
               omega_ic=O_dt*omega_ic + LFfac/c_moi_ic * &
                    (w1*lorentz_torque_ic+w2*lorentz_torque_icLast)
               omega_ic=dt*omega_ic
            else if ( kbotv == 2 ) then ! no slip, omega given by z10
               omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
            end if
            omega_ic1=omega_ic
            !write(*,"(A,I4,A,ES20.13)") "after ic update, nLMB = ",nLMB,", omega_ic = ",omega_ic
         end if
      end if  ! l=1,m=0 contained in block ?
      !PERFOFF
      !PERFON('upZ_ang')
      !--- We correct so that the angular moment about axis in the equatorial plane
      !    vanish and the angular moment about the (planetary) rotation axis
      !    is kept constant.
      l1m1=lm2(1,1)
      if ( l_correct_AMz .and.  l1m0 > 0 .and. &
           lmStart_00 <= l1m0 .and. lmStop >= l1m0 ) then
    
         do nR=1,n_r_max
            z10(nR)=z(l1m0,nR)
         end do
         call get_angular_moment(z10,z11,omega_ic,omega_ma, &
                                         angular_moment_oc, &
                                 angular_moment_ic,angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
                              angular_moment_ma(i)
         end do
         if ( ( ktopv == 2 .and. l_rot_ma ) .and. ( kbotv == 2 .and. l_rot_ic ) ) then
            nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb + &
                 c_moi_ic*c_z10_omega_ic*r_icb*r_icb + &
                 c_moi_oc*y10_norm
         else if ( ktopv == 2 .and. l_rot_ma ) then
            nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb+c_moi_oc*y10_norm
         else if ( kbotv == 2 .and. l_rot_ic ) then
            nomi=c_moi_ic*c_z10_omega_ic*r_icb*r_icb+c_moi_oc*y10_norm
         else
            nomi=c_moi_oc*y10_norm
         end if
         corr_l1m0=cmplx(angular_moment(3)-AMstart,0.0_cp,kind=cp)/nomi
    
         !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
         !         derivatives:
         !$OMP PARALLEL do default(shared) &
         !$OMP PRIVATE(r_E_2,nR)
         do nR=1,n_r_max
            r_E_2=r(nR)*r(nR)
            z(l1m0,nR)  =z(l1m0,nR)  - rho0(nR)*r_E_2*corr_l1m0
            dz(l1m0,nR) =dz(l1m0,nR) - rho0(nR)*( &
                 two*r(nR)+r_E_2*beta(nR))*corr_l1m0
            workA(l1m0,nR)=workA(l1m0,nR)-rho0(nR)*( &
                 two+four*beta(nR)*r(nR) + &
                 dbeta(nR)*r_E_2 + &
                 beta(nR)*beta(nR)*r_E_2 )*corr_l1m0
         end do
         !$OMP END PARALLEL DO
         if ( ktopv == 2 .and. l_rot_ma ) &
              omega_ma=c_z10_omega_ma*real(z(l1m0,n_r_cmb))
         if ( kbotv == 2 .and. l_rot_ic ) &
              omega_ic=c_z10_omega_ic*real(z(l1m0,n_r_icb))
         omega_ic1=omega_ic
         omega_ma1=omega_ma
    
      end if ! l=1,m=0 contained in lm-block ?
    
      if ( l_correct_AMe .and.  l1m1 > 0 .and. &
           lmStart_00 <= l1m1 .and. lmStop >= l1m1 ) then
    
         do nR=1,n_r_max
            z11(nR)=z(l1m1,nR)
         end do
         call get_angular_moment(z10,z11,omega_ic,omega_ma, &
                                         angular_moment_oc, &
                                 angular_moment_ic,angular_moment_ma)
         do i=1,3
            angular_moment(i)=angular_moment_oc(i) + angular_moment_ic(i) + &
                              angular_moment_ma(i)
         end do
         corr_l1m1=cmplx(angular_moment(1),-angular_moment(2),kind=cp) / &
              (two*y11_norm*c_moi_oc)
    
         !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
         !         derivatives:
         !$OMP PARALLEL do default(shared) &
         !$OMP private(nR,r_E_2)
         do nR=1,n_r_max
            r_E_2=r(nR)*r(nR)
            z(l1m1,nR)  =z(l1m1,nR)  -  rho0(nR)*r_E_2*corr_l1m1
            dz(l1m1,nR) =dz(l1m1,nR) -  rho0(nR)*( &
                         two*r(nR)+r_E_2*beta(nR))*corr_l1m1
            workA(l1m1,nR)=workA(l1m1,nR)-rho0(nR)*( &
                          two+four*beta(nR)*r(nR) + &
                                     dbeta(nR)*r_E_2 + &
                             beta(nR)*beta(nR)*r_E_2 )*corr_l1m1
         end do
         !$OMP END PARALLEL DO
      end if ! l=1,m=1 contained in lm-block ?
    
      !if (DEBUG_OUTPUT) then
      !   do nR=1,n_r_max
      !      if ((nR == 1).or.(nR == 5)) then
      !         do lm=llm,ulm
      !            write(*,"(4X,A,2I3,4ES22.14)") "upZ_new: ",nR,lm,z(lm,nR),dz(lm,nR)
      !         end do
      !      end if
      !      write(*,"(A,I3,4ES22.14)") "upZ_new: ",nR,get_global_SUM( z(:,nR) ),get_global_SUM( dz(:,nR) )
      !   end do
      !end if
      !-- Calculate explicit time step part:
      !$OMP PARALLEL default(shared) &
      !$OMP private(nR,lm1,Dif)
      !$OMP DO
      do nR=n_r_cmb+1,n_r_icb-1
         !write(*,"(A,I4,5ES20.12)") "r-dependent : ",nR,dLvisc(nR),beta(nR),or1(nR),or2(nR),dbeta(nR)
         do lm1=lmStart_00,lmStop
            Dif(lm1)=hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1)))*                &
                     dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*visc(nR)*  &
                 & ( workA(lm1,nR)   +(dLvisc(nR)-beta(nR))  *dz(lm1,nR)     &
                 &   -( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)           &
                 &      + dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)       &
                 &      + dbeta(nR)+ two*beta(nR)*or1(nR) ) * z(lm1,nR) )
    
    !        if (nR == 2) then
    !           write(*,"(2I4,8ES20.12)") nR,lm1,workA(lm1,nR),Dif(lm1),z(lm1,nR),dz(lm1,nR)
    !        end if
            dzdtLast(lm1,nR)=dzdt(lm1,nR)-coex*Dif(lm1)
            if ( lRmsNext ) then
               workB(lm1,nR)= O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))* &
                              or2(nR)*(z(lm1,nR)-workB(lm1,nR))
               if ( l_RMStest ) workB(lm1,nR)= workB(lm1,nR)-Dif(lm1)
            end if
         end do
         if ( lRmsNext ) then
            call hInt2Tor(Dif,1,lm_max,nR,lmStart_00,lmStop, &
                          DifTor2hInt(nR,1),DifTorAs2hInt(nR,1),lo_map)
            call hInt2Tor(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                          dtVTor2hInt(nR,1),dtVTorAs2hInt(nR,1),lo_map)
         end if
      end do
      !$OMP end do
      !--- Note: from ddz=workA only the axisymmetric contributions are needed
      !    beyond this point for the TO calculation.
      !    Parallization note: Very likely, all axisymmetric modes m=0 are
      !    located on the first processor #0.
      if ( l_TO ) then
         !$OMP do private(nR,lm1,l1,m1)
         do nR=1,n_r_max
            ddzASL_loc(:)=0.0_cp
            do lm1=lmStart_00,lmStop
               l1=lm2l(lm1)
               m1=lm2m(lm1)
               if ( m1 == 0 ) ddzASL_loc(l1+1)=real(workA(lm1,nR))
            end do
#ifdef WITH_MPI
            call MPI_Allreduce(ddzASL_loc, ddzASL(:,nR), l_max+1, MPI_DEF_REAL, &
                 &             MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            ddzASL(:,nR)=ddzASL_loc(:)
#endif
         end do
         !$OMP end do
      end if
      !$OMP END PARALLEL
    
      !----- Special thing for l=1,m=0 for rigid boundaries and
      !      if IC or mantle are allowed to rotate:
      if ( l10 .and. l_z10mat ) then ! z10 term !
         lm1=lm2(1,0)
         !----- NOTE opposite sign of visouse torque on ICB and CMB:
         if ( .not. l_SRMA .and. ktopv == 2 .and. l_rot_ma ) then
            d_omega_ma_dtLast=d_omega_ma_dt -            &
                 coex * ( two*or1(1)*real(z(lm1,1))-real(dz(lm1,1)) )
         end if
         if ( .not. l_SRIC .and. kbotv == 2 .and. l_rot_ic ) THEn
            d_omega_ic_dtLast=d_omega_ic_dt +                     &
                 coex * ( two*or1(n_r_max)*real(z(lm1,n_r_max))- &
                 real(dz(lm1,n_r_max)) )
         end if
      end if
      !PERFOFF
   end subroutine updateZ
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_Z10
   subroutine get_z10Mat(dt,l,hdif,zMat,zPivot,zMat_fac)
#else
   subroutine get_z10Mat(dt,l,hdif,zMat,zPivot)
#endif
      !
      !  Purpose of this subroutine is to construct and LU-decompose the  
      !  inversion matrix z10mat for the implicit time step of the       
      !  toroidal velocity potential z of degree l=1 and order m=0.       
      !  This differs from the the normal zmat only if either the ICB or  
      !  CMB have no-slip boundary condition and inner core or mantle are 
      !  chosen to rotate freely (either kbotv=1 and/or ktopv=1).         
      !
      
      real(cp), intent(in) :: dt      ! Time step internal
      real(cp), intent(in) :: hdif    ! Value of hyperdiffusivity in zMat terms
      integer,  intent(in) :: l       ! Variable to loop over l's

      !-- Output: z10Mat and pivot_z10
      real(cp), intent(out) :: zMat(n_r_max,n_r_max) ! LHS matrix to calculate z
      integer,  intent(out) :: zPivot(n_r_max)       ! Pivot to invert zMat
#ifdef WITH_PRECOND_Z10
      real(cp), intent(out) :: zMat_fac(n_r_max)     ! Inverse of max(zMat) for inversion
#endif

      !-- local variables:
      integer :: nR,nCheb,info
      real(cp) :: O_dt,dLh

      O_dt=one/dt
      dLh=real(l*(l+1),kind=cp)

      !-- Boundary conditions:
      do nCheb=1,n_cheb_max
         !----- CMB condition:
         !        Note opposite sign of viscous torques (-dz+2 z /r )
         !        for CMB and ICB!

         if ( ktopv == 1 ) then  ! free slip
            zMat(1,nCheb)=                   cheb_norm * &
                ( (two*or1(1)+beta(1))*cheb(nCheb,1) - &
                                       dcheb(nCheb,1) )
         else if ( ktopv == 2 ) then ! no slip
            if ( l_SRMA ) then
               zMat(1,nCheb)= cheb_norm * c_z10_omega_ma*cheb(nCheb,1)
            else if ( l_rot_ma ) then
               zMat(1,nCheb)= cheb_norm *               ( &
                         c_dt_z10_ma*O_dt*cheb(nCheb,1) - &
                      alpha*( two*or1(1)*cheb(nCheb,1) - &
                                         dcheb(nCheb,1) ) )
            else
               zMat(1,nCheb)= cheb_norm*cheb(nCheb,1)
            end if
         end if

         !----- ICB condition:
         if ( kbotv == 1 ) then  ! free slip
            zMat(n_r_max,nCheb)=                          cheb_norm * &
            ( (two*or1(n_r_max)+beta(n_r_max))*cheb(nCheb,n_r_max) - &
                                               dcheb(nCheb,n_r_max) )
         else if ( kbotv == 2 ) then ! no slip
            if ( l_SRIC ) then
               zMat(n_r_max,nCheb)= cheb_norm * &
                                   c_z10_omega_ic*cheb(nCheb,n_r_max)
            else if ( l_rot_ic ) then     !  time integration of z10
               zMat(n_r_max,nCheb)= cheb_norm *             (     &
                           c_dt_z10_ic*O_dt*cheb(nCheb,n_r_max) + &
                  alpha*( two*or1(n_r_max)*cheb(nCheb,n_r_max) - &
                                           dcheb(nCheb,n_r_max) ) )
            else
               zMat(n_r_max,nCheb)= cheb_norm * cheb(nCheb,n_r_max)
            end if
         end if

      end do

      !----- Other points: (same as zMat)
      do nCheb=1,n_r_max
         do nR=2,n_r_max-1
            zMat(nR,nCheb)=                   cheb_norm * ( &
                          O_dt*dLh*or2(nR)*cheb(nCheb,nR) - &
                        alpha*hdif*dLh*visc(nR)*or2(nR) * ( &
                                         d2cheb(nCheb,nR) + &
             (dLvisc(nR)- beta(nR))*      dcheb(nCheb,nR) - &
            (dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)  + &
              dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR))* &
                                           cheb(nCheb,nR) ) )

         end do
      end do

      !-- Normalisation
      do nR=1,n_r_max
         zMat(nR,1)      =half*zMat(nR,1)
         zMat(nR,n_r_max)=half*zMat(nR,n_r_max)
      end do

      !-- Fill up with zeros:
      do nCheb=n_cheb_max+1,n_r_max
         zMat(1,nCheb)      =0.0_cp
         zMat(n_r_max,nCheb)=0.0_cp
      end do

#ifdef WITH_PRECOND_Z10
      ! compute the linesum of each line
      do nR=1,n_r_max
         zMat_fac(nR)=one/maxval(abs(zMat(nR,:)))
         zMat(nR,:) = zMat(nR,:)*zMat_fac(nR)
      end do
#endif

!-- LU-decomposition of z10mat:
      call sgefa(zMat,n_r_max,n_r_max,zPivot,info)

      if ( info /= 0 ) then
         write(*,*) 'ERROR MESSAGE FROM subroutine GET_z10MAT:'
         write(*,*) 'singular matrix z10Mat!'
         stop
      end if

   end subroutine get_z10Mat
!-----------------------------------------------------------------------
#ifdef WITH_PRECOND_Z
   subroutine get_zMat(dt,l,hdif,zMat,zPivot,zMat_fac)
#else
   subroutine get_zMat(dt,l,hdif,zMat,zPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
      !  zmat(i,j) for the NS equation.                                   
      !
    
      !-- Input variables:
      real(cp), intent(in) :: dt                     ! Time interval
      integer,  intent(in) :: l                      ! Variable to loop over degrees
      real(cp), intent(in) :: hdif                   ! Hyperdiffusivity
    
      !-- Output variables:
      real(cp), intent(out) :: zMat(n_r_max,n_r_max) ! Matrix with LHS of z equation
      integer,  intent(out) :: zPivot(n_r_max)       ! Pivot for zMat inversion
#ifdef WITH_PRECOND_Z
      real(cp), intent(out) :: zMat_fac(n_r_max)     !  Inverse of max(zMat) for the inversion
#endif

      !-- local variables:
      integer :: nR,nCheb
      integer :: info
      real(cp) :: O_dt,dLh

#ifdef MATRIX_CHECK
      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer, save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

      O_dt=one/dt
      dLh=real(l*(l+1),kind=cp)
    
      !----- Boundary conditions, see above:
      do nCheb=1,n_cheb_max
         if ( ktopv == 1 ) then  ! free slip !
            zMat(1,nCheb)=cheb_norm *             ( &
                 dcheb(nCheb,1) - (two*or1(1)+beta(1))*cheb(nCheb,1) )
         else                    ! no slip, note exception for l=1,m=0
            zMat(1,nCheb)=cheb_norm*cheb(nCheb,1)
         end if
    
         if ( kbotv == 1 ) then  ! free slip !
            zMat(n_r_max,nCheb)= cheb_norm *            ( &
                 dcheb(nCheb,n_r_max) - &
                 (two*or1(n_r_max)+beta(n_r_max))*cheb(nCheb,n_r_max) )
         else                    ! no slip, note exception for l=1,m=0
            zMat(n_r_max,nCheb)= cheb_norm * cheb(nCheb,n_r_max)
         end if
    
      end do  !  loop over nCheb
    
      if ( n_cheb_max < n_r_max ) then ! fill with zeros !
         do nCheb=n_cheb_max+1,n_r_max
            zMat(1,nCheb)      =0.0_cp
            zMat(n_r_max,nCheb)=0.0_cp
         end do
      end if
    
      !----- Other points:
      do nCheb=1,n_r_max
         do nR=2,n_r_max-1
            zMat(nR,nCheb)= cheb_norm * ( O_dt*dLh*or2(nR)* cheb(nCheb,nR) &
                 &   -alpha*hdif*dLh*visc(nR)*or2(nR) * ( d2cheb(nCheb,nR) &
                 &      + (dLvisc(nR)- beta(nR)) *         dcheb(nCheb,nR) &
                 &      - ( dLvisc(nR)*beta(nR)+two*dLvisc(nR)*or1(nR)    &
                 &          +dLh*or2(nR)+dbeta(nR)+two*beta(nR)*or1(nR)   &
                 &                                      ) * cheb(nCheb,nR) ) )
         end do
      end do
    
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         zMat(nR,1)      =half*zMat(nR,1)
         zMat(nR,n_r_max)=half*zMat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_Z
      ! compute the linesum of each line
      do nR=1,n_r_max
         zMat_fac(nR)=one/maxval(abs(zMat(nR,:)))
         zMat(nR,:) = zMat(nR,:)*zMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      ! copy the zMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "zMat_",l,"_",counter,".dat"
      open(NEWUNIT=filehandle,file=trim(filename))
      counter= counter+1
 
      do i=1,n_r_max
         do j=1,n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") zMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=zMat
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
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of zMat for l=",l," is ",rcond
#endif

  !----- LU decomposition:
      call sgefa(zMat,n_r_max,n_r_max,zPivot,info)

      if ( info /= 0 ) then
         write(*,*) 'Singular matrix zmat for l=',l,", info = ",info
         stop '34'
      end if

   end subroutine get_zMat
!-----------------------------------------------------------------------------
end module updateZ_mod
