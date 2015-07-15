!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
MODULE updateZ_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE torsional_oscillations
  USE init_fields
  USE blocking,ONLY: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map,&
       &lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE matrices
  USE RMS
  USE const
  use parallel_mod
#ifdef WITH_MKL_LU
  USE lapack95, ONLY: getrs
#else
  USE algebra, ONLY: cgeslML,cgesl
#endif
  USE LMLoop_data,ONLY: llm,ulm,llm_real,ulm_real
  USE communications, only:get_global_sum
  IMPLICIT NONE

  PRIVATE
  !-- Input of recycled work arrays:
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: workA,workB,workC
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: rhs1 ! RHS for other modes
  INTEGER :: maxThreads
  
  PUBLIC :: updateZ,initialize_updateZ

CONTAINS
  SUBROUTINE initialize_updateZ
    allocate(workA(llm:ulm,n_r_max))
    allocate(workB(llm:ulm,n_r_max))
    allocate(workC(llm:ulm,n_r_max))

#ifdef WITHOMP
    maxThreads=omp_get_max_threads()
#else
    maxThreads=1
#endif
    ALLOCATE(rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
  END SUBROUTINE initialize_updateZ

  SUBROUTINE updateZ(z,dz,dzdt,dzdtLast,time, &
     &             omega_ma,d_omega_ma_dtLast, &
     &             omega_ic,d_omega_ic_dtLast, &
     &             lorentz_torque_ma,lorentz_torque_maLast, &
     &             lorentz_torque_ic,lorentz_torque_icLast, &
     &             w1,coex,dt,lRmsNext)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------------------------------------------------------------+
  !  |                                                                   |
  !  |     Changed by J.Wicht 20.07.2000                                 |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+

  !  its radial derivatives
  !  adds explicit part to time derivatives of z

  !  Input:  w1 - weight for dbdt-contribution from current time step
  !               (w2=1-w1: weight for contrib. from previous step)
  !          coex - factor depending on weighting alpha of implicit contribution
  !          m1,m2- range of mca-indices in which field is updated
  !                 (harmonic order m=(mca-1)*minc)

  !--------------------------------------------------------------------


  !-- Input/output of scalar fields:
  COMPLEX(kind=8),INTENT(INOUT) :: z(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(OUT) :: dz(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(IN)    :: dzdt(llm:ulm,n_r_max)
  COMPLEX(kind=8),INTENT(INOUT) :: dzdtLast(llm:ulm,n_r_max)
  REAL(kind=8),intent(INOUT) :: d_omega_ma_dtLast,d_omega_ic_dtLast
  REAL(kind=8),intent(OUT) :: omega_ma,omega_ic
  REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_maLast
  REAL(kind=8),intent(IN) :: lorentz_torque_ic,lorentz_torque_icLast
  !-- Output:
  !       Updated z,dz,dzdtLast,d_omega_maLast,d_omega_icLast,ddzAS

  !-- Input of other variables:
  REAL(kind=8),intent(IN) :: time
  REAL(kind=8),intent(IN) :: w1    ! weight for time step !
  REAL(kind=8),intent(IN) :: coex  ! factor depending on alpha
  REAL(kind=8),intent(IN) :: dt
  LOGICAL,intent(IN) :: lRmsNext


  !-- local variables:
  REAL(kind=8) :: w2                  ! weight of second time step
  REAL(kind=8) :: O_dt
  REAL(kind=8) :: d_omega_ic_dt,d_omega_ma_dt
  INTEGER :: l1,m1              ! degree and order
  INTEGER :: lm1,lm,lmB         ! position of (l,m) in array
  INTEGER :: lmStart_real       ! range of lm for real array
  INTEGER :: lmStop_real        !
  INTEGER :: lmStart_00         ! excluding l=0,m=0
  INTEGER :: lmStart,lmStop ! max and min number of orders m
  INTEGER :: nLMB2
  INTEGER :: nR                 ! counts radial grid points
  INTEGER :: n_cheb             ! counts cheb modes
  COMPLEX(kind=8) :: rhs(n_r_max)   ! RHS of matrix multiplication
  !COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! RHS for other modes
  COMPLEX(kind=8) :: z10(n_r_max),z11(n_r_max) ! toroidal flow scalar components
  REAL(kind=8) :: angular_moment(3)   ! total angular momentum
  REAL(kind=8) :: angular_moment_oc(3)! x,y,z component of outer core angular mom.
  REAL(kind=8) :: angular_moment_ic(3)! x,y,z component of inner core angular mom.
  REAL(kind=8) :: angular_moment_ma(3)! x,y,z component of mantle angular mom.
  COMPLEX(kind=8) :: corr_l1m0       ! correction factor for z(l=1,m=0)
  COMPLEX(kind=8) :: corr_l1m1       ! correction factor for z(l=1,m=1)
  REAL(kind=8) :: r_E_2               ! =r**2
  REAL(kind=8) :: nomi                ! nominator for Z10 AM correction
  INTEGER :: l1m0,l1m1          ! position of (l=1,m=0) and (l=1,m=1) in lm.
  INTEGER :: i                  ! counter
  LOGICAL :: l10
  integer :: nLMB

  COMPLEX(kind=8) :: Dif(lm_max)

  INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
  INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
  INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

  logical :: DEBUG_OUTPUT=.false.
  INTEGER :: nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
  INTEGER :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid
  complex(kind=8) :: rhs_sum

  !-- end of declaration
  !-----------------------------------------------------------------------

  !CALL mpi_barrier(MPI_COMM_WORLD,ierr)
  !WRITE(*,"(3(A,2ES20.12))") "begin upZ: dzdt = ",get_global_sum( dzdt ),&
  !     &", z = ",get_global_sum( z ),&
  !     &", dzdtLast = ",get_global_sum( dzdtLast )
  !CALL mpi_barrier(MPI_COMM_WORLD,ierr)

  IF ( .NOT. l_update_v ) RETURN

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
  lmStart_00  =MAX(2,lmStart)
  lmStart_real=2*lmStart_00-1
  lmStop_real =2*lmStop
  l1m0        =lm2(1,0)

  w2  =1.D0-w1
  O_dt=1.D0/dt

  l10=.FALSE.
  !$OMP PARALLEL default(shared)
  !$OMP SINGLE
  DO nLMB2=1,nLMBs2(nLMB)
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
     !WRITE(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

     IF ( l1 .NE. 0 ) THEN
        IF ( .NOT. lZmat(l1) ) THEN
#ifdef WITH_PRECOND_Z
           CALL get_zMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                &        zMat(1,1,l1),zPivot(1,l1),zMat_fac(1,l1))
#else
           CALL get_zMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                &        zMat(1,1,l1),zPivot(1,l1))
#endif
           lZmat(l1)=.TRUE.
           !WRITE(*,"(A,I3,A,2ES20.12)") "zMat(",l1,") = ",SUM(zMat(:,:,l1))
        END IF
     END IF

     DO iChunk=1,nChunks
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

        DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
           !DO lm=1,sizeLMB2(nLMB2,nLMB)
           lm1=lm22lm(lm,nLMB2,nLMB)
           !l1 =lm22l(lm,nLMB2,nLMB)
           m1 =lm22m(lm,nLMB2,nLMB)
           IF ( lm1 == l1m0 ) l10= .TRUE. 

           IF ( l_z10mat .AND. lm1 == l1m0 ) THEN
              !WRITE(*,"(A,3I3)") "l_z10mat and lm1=",lm1,l1,m1
              !PERFON('upZ_z10')
              !----- Special treatment of z10 component if ic or mantle
              !      are allowed to rotate about z-axis (l_z10mat=.true.) and
              !      we use no slip boundary condition (ktopv=1,kbotv=1):
              !      Lorentz torque is the explicit part of this time integration
              !      at the boundaries!
              !      Note: no angular momentum correction necessary for this case !
              IF ( .NOT. lZ10mat ) THEN
#ifdef WITH_PRECOND_Z10
                 CALL get_z10Mat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                      z10Mat,z10Pivot,z10Mat_fac)
#else
                 CALL get_z10Mat(dt,l1,hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1))), &
                      z10Mat,z10Pivot)
#endif
                 lZ10mat=.TRUE.
              END IF

              IF ( l_SRMA ) THEN
                 tOmega_ma1=time+tShift_ma1
                 tOmega_ma2=time+tShift_ma2
                 omega_ma= omega_ma1*DCOS(omegaOsz_ma1*tOmega_ma1) + &
                      &    omega_ma2*DCOS(omegaOsz_ma2*tOmega_ma2)
                 rhs(1)=omega_ma
              ELSE IF ( ktopv == 2 .AND. l_rot_ma ) THEN  ! time integration
                 d_omega_ma_dt=LFfac*c_lorentz_ma*lorentz_torque_ma
                 rhs(1)=O_dt*c_dt_z10_ma*z(lm1,1) + &
                      w1*d_omega_ma_dt + &
                      w2*d_omega_ma_dtLast
              ELSE
                 rhs(1)=0.d0
              END IF

              IF ( l_SRIC ) THEN
                 tOmega_ic1=time+tShift_ic1
                 tOmega_ic2=time+tShift_ic2
                 omega_ic= omega_ic1*COS(omegaOsz_ic1*tOmega_ic1) + &
                      &    omega_ic2*COS(omegaOsz_ic2*tOmega_ic2)
                 rhs(n_r_max)=omega_ic
              ELSE if ( kbotv == 2 .AND. l_rot_ic ) then  ! time integration
                 d_omega_ic_dt = LFfac*c_lorentz_ic*lorentz_torque_ic
                 rhs(n_r_max)=O_dt*c_dt_z10_ic*z(lm1,n_r_max) + &
                      w1*d_omega_ic_dt + &
                      w2*d_omega_ic_dtLast
              ELSE
                 rhs(n_r_max)=0.d0
              END IF

              !----- This is the normal RHS for the other radial grid points:
              DO nR=2,n_r_max-1
                 rhs(nR)=O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*z(lm1,nR)+ &
                      w1*dzdt(lm1,nR)+ &
                      w2*dzdtLast(lm1,nR)
              END DO

#ifdef WITH_PRECOND_Z10
              DO nR=1,n_r_max
                 rhs(nR) = z10Mat_fac(nR)*rhs(nR)
              END DO
#endif
              IF (DEBUG_OUTPUT) THEN
                 rhs_sum=SUM(rhs)
                 WRITE(*,"(2I3,A,2(I4,F20.16))") nLMB2,lm1,":rhs_sum (z10) before = ",&
                      & EXPONENT(REAL(rhs_sum)),FRACTION(REAL(rhs_sum)),&
                      & EXPONENT(AIMAG(rhs_sum)),FRACTION(AIMAG(rhs_sum))
                 !DO nR=1,n_r_max
                 !   WRITE(*,"(3I4,A,2(I4,F20.16))") nLMB2,lm1,nR,":rhs (z10) before = ",&
                 !        & EXPONENT(REAL(rhs(nR))),FRACTION(REAL(rhs(nR))),&
                 !        & EXPONENT(AIMAG(rhs(nR))),FRACTION(AIMAG(rhs(nR)))
                 !END DO
              END IF
#ifdef WITH_MKL_LU
              CALL getrs(CMPLX(z10Mat,0.D0,KIND=KIND(0.D0)),z10Pivot,rhs)
#else
              CALL cgesl(z10Mat,n_r_max,n_r_max,z10Pivot,rhs)
#endif
              IF (DEBUG_OUTPUT) THEN
                 !DO nR=1,n_r_max
                 !   WRITE(*,"(3I4,A,2(I4,F20.16))") nLMB2,lm1,nR,":rhs (z10) after = ",&
                 !        & EXPONENT(REAL(rhs(nR))),FRACTION(REAL(rhs(nR))),&
                 !        & EXPONENT(AIMAG(rhs(nR))),FRACTION(AIMAG(rhs(nR)))
                 !END DO
                 rhs_sum=SUM(rhs)
                 WRITE(*,"(2I3,A,2(I4,F20.16))") nLMB2,lm1,":rhs_sum (z10) after = ",&
                      & EXPONENT(REAL(rhs_sum)),FRACTION(REAL(rhs_sum)),&
                      & EXPONENT(AIMAG(rhs_sum)),FRACTION(AIMAG(rhs_sum))
              END IF


           ELSE IF ( l1 /= 0 ) THEN
              !PERFON('upZ_ln0')
              lmB=lmB+1
              rhs1(1,lmB,threadid)      =0.D0
              rhs1(n_r_max,lmB,threadid)=0.D0
              DO nR=2,n_r_max-1
                 rhs1(nR,lmB,threadid)=&
                      & O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*z(lm1,nR)&
                      & + w1*dzdt(lm1,nR) &
                      & + w2*dzdtLast(lm1,nR)
#ifdef WITH_PRECOND_Z
                 rhs1(nR,lmB,threadid) = zMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
              END DO
              !PERFOFF
           END IF
        END DO

        !PERFON('upZ_sol')
        IF ( lmB > lmB0 ) THEN
#ifdef WITH_MKL_LU
           CALL getrs(CMPLX(zMat(:,:,l1),0.D0,KIND=KIND(0.D0)), &
                &       zPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid))
#else
           CALL cgeslML(zMat(:,:,l1),n_r_max,n_r_max, &
                &       zPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),n_r_max,lmB-lmB0)
#endif
        END IF
        !PERFOFF
        IF ( lRmsNext ) THEN ! Store old z
           DO nR=1,n_r_max
              DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                 lm1=lm22lm(lm,nLMB2,nLMB)
                 workB(lm1,nR)=z(lm1,nR)
              END DO
           END DO
        END IF

        lmB=lmB0
        DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
           !DO lm=1,sizeLMB2(nLMB2,nLMB)
           lm1=lm22lm(lm,nLMB2,nLMB)
           !l1 =lm22l(lm,nLMB2,nLMB)
           m1 =lm22m(lm,nLMB2,nLMB)
           
           IF ( l_z10mat .AND. lm1 == l1m0 ) THEN
              DO n_cheb=1,n_cheb_max
                 z(lm1,n_cheb)=REAL(rhs(n_cheb))
              END DO
           ELSE IF ( l1 /= 0 ) THEN
              lmB=lmB+1
              IF ( m1 > 0 ) THEN
                 DO n_cheb=1,n_cheb_max
                    z(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                 END DO
              ELSE
                 DO n_cheb=1,n_cheb_max
                    z(lm1,n_cheb)= &
                         CMPLX(REAL(rhs1(n_cheb,lmB,threadid)),0.D0,KIND=KIND(0d0))
                 END DO
              END IF
           END IF
        END DO
        !PERFOFF
        !$OMP END TASK
     end DO
     !$OMP END TASK
  END DO       ! end of loop over lm blocks
  !$OMP END SINGLE
  !$OMP END PARALLEL


  !-- set cheb modes > n_cheb_max to zero (dealiazing)
  DO n_cheb=n_cheb_max+1,n_r_max
     DO lm1=lmStart,lmStop
        z(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
     END DO
  END DO

  !PERFON('upZ_drv')
  all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
  IF (all_lms < omp_get_max_threads()) THEN
     CALL omp_set_num_threads(all_lms)
     per_thread=1
  ELSE
     per_thread=all_lms/omp_get_max_threads()
  END IF
#else
  per_thread=all_lms
#endif
  !$OMP PARALLEL default(none) &
  !$OMP private(iThread,start_lm,stop_lm) &
  !$OMP shared(per_thread,lmStart_real,lmStop_real,nThreads) &
  !$OMP shared(z,dz,dzdtLast,i_costf_init,d_costf_init,drx,ddrx) &
  !$OMP shared(n_r_max,n_cheb_max,workA,workC,llm_real,ulm_real)
  !$OMP SINGLE
#ifdef WITHOMP
  nThreads=omp_get_num_threads()
#else
  nThreads=1
#endif
  !$OMP END SINGLE
  !$OMP BARRIER
  !$OMP DO
  DO iThread=0,nThreads-1
     start_lm = lmStart_real+iThread*per_thread
     stop_lm  = start_lm+per_thread-1
     IF (iThread == nThreads-1) stop_lm=lmStop_real
     !WRITE(*,"(3(A,I5))") "thread ",omp_get_thread_num()," from ",start_lm," to ",stop_lm
     !-- Get derivatives:
     CALL costf1(z, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
          &      dzdtLast, i_costf_init, d_costf_init)
     CALL get_ddr(z, dz, workA, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
          &       n_r_max, n_cheb_max, dzdtLast, workC, &
          &       i_costf_init,d_costf_init,drx,ddrx)
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
#ifdef WITHOMP
  CALL omp_set_num_threads(omp_get_max_threads())
#endif
  !PERFOFF
  
  !PERFON('upZ_icma')
  !--- Update of inner core and mantle rotation:
  IF ( l10 ) THEN
     IF ( l_rot_ma .AND. .NOT. l_SRMA ) THEN
        IF ( ktopv == 1 ) THEN  ! free slip, explicit time stepping of omega !
           omega_ma=O_dt*omega_ma + LFfac/c_moi_ma * &
                (w1*lorentz_torque_ma+w2*lorentz_torque_maLast)
           omega_ma=dt*omega_ma
        ELSE IF ( ktopv == 2 ) THEN ! no slip, omega given by z10
           omega_ma=c_z10_omega_ma*REAL(z(l1m0,n_r_cmb))
        END IF
        omega_ma1=omega_ma
     END IF
     IF ( l_rot_ic .AND. .NOT. l_SRIC ) THEN
        IF ( kbotv == 1 ) THEN  ! free slip, explicit time stepping of omega !
           omega_ic=O_dt*omega_ic + LFfac/c_moi_ic * &
                (w1*lorentz_torque_ic+w2*lorentz_torque_icLast)
           omega_ic=dt*omega_ic
        ELSE IF ( kbotv == 2 ) THEN ! no slip, omega given by z10
           omega_ic=c_z10_omega_ic*REAL(z(l1m0,n_r_icb))
        END IF
        omega_ic1=omega_ic
        !write(*,"(A,I4,A,ES20.13)") "after ic update, nLMB = ",nLMB,", omega_ic = ",omega_ic
     END IF
  END IF  ! l=1,m=0 contained in block ?
  !PERFOFF
  !PERFON('upZ_ang')
  !--- We correct so that the angular moment about axis in the equatorial plane
  !    vanish and the angular moment about the (planetary) rotation axis
  !    is kept constant.
  l1m1=lm2(1,1)
  IF ( l_correct_AMz .AND.  l1m0 > 0 .AND. &
       lmStart_00 <= l1m0 .AND. lmStop >= l1m0 ) THEN

     DO nR=1,n_r_max
        z10(nR)=z(l1m0,nR)
     END DO
     CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
          angular_moment_oc, &
          angular_moment_ic,angular_moment_ma)
     DO i=1,3
        angular_moment(i)=angular_moment_oc(i) + &
             angular_moment_ic(i) + &
             angular_moment_ma(i)
     END DO
     IF ( ( ktopv == 2 .AND. l_rot_ma ) .AND. &
          ( kbotv == 2 .AND. l_rot_ic ) ) THEN
        nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb + &
             c_moi_ic*c_z10_omega_ic*r_icb*r_icb + &
             c_moi_oc*y10_norm
     ELSE IF ( ktopv == 2 .AND. l_rot_ma ) THEN
        nomi=c_moi_ma*c_z10_omega_ma*r_cmb*r_cmb+c_moi_oc*y10_norm
     ELSE IF ( kbotv == 2 .AND. l_rot_ic ) THEN
        nomi=c_moi_ic*c_z10_omega_ic*r_icb*r_icb+c_moi_oc*y10_norm
     ELSE
        nomi=c_moi_oc*y10_norm
     END IF
     corr_l1m0=CMPLX(angular_moment(3)-AMstart,0.d0,KIND=KIND(0d0))/nomi

     !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
     !         derivatives:
     !$OMP PARALLEL DO default(shared) &
     !$OMP PRIVATE(r_E_2,nR)
     DO nR=1,n_r_max
        r_E_2=r(nR)*r(nR)
        z(l1m0,nR)  =z(l1m0,nR)  - rho0(nR)*r_E_2*corr_l1m0
        dz(l1m0,nR) =dz(l1m0,nR) - rho0(nR)*( &
             2.d0*r(nR)+r_E_2*beta(nR))*corr_l1m0
        workA(l1m0,nR)=workA(l1m0,nR)-rho0(nR)*( &
             2.d0+4.d0*beta(nR)*r(nR) + &
             dbeta(nR)*r_E_2 + &
             beta(nR)*beta(nR)*r_E_2 )*corr_l1m0
     END DO
     !$OMP END PARALLEL DO
     IF ( ktopv == 2 .AND. l_rot_ma ) &
          omega_ma=c_z10_omega_ma*REAL(z(l1m0,n_r_cmb))
     IF ( kbotv == 2 .AND. l_rot_ic ) &
          omega_ic=c_z10_omega_ic*REAL(z(l1m0,n_r_icb))
     omega_ic1=omega_ic
     omega_ma1=omega_ma

  END IF ! l=1,m=0 contained in lm-block ?

  IF ( l_correct_AMe .AND.  l1m1 > 0 .AND. &
       lmStart_00 <= l1m1 .AND. lmStop >= l1m1 ) THEN

     DO nR=1,n_r_max
        z11(nR)=z(l1m1,nR)
     END DO
     CALL get_angular_moment(z10,z11,omega_ic,omega_ma, &
          angular_moment_oc, &
          angular_moment_ic,angular_moment_ma)
     DO i=1,3
        angular_moment(i)=angular_moment_oc(i) + &
             angular_moment_ic(i) + &
             angular_moment_ma(i)
     END DO
     corr_l1m1=CMPLX(angular_moment(1),-angular_moment(2),KIND=KIND(0d0)) / &
          (2.d0*y11_norm*c_moi_oc)

     !-------- Correct z(2,nR) and z(l_max+2,nR) plus the respective
     !         derivatives:
     !$OMP PARALLEL DO default(shared) &
     !$OMP private(nR,r_E_2)
     DO nR=1,n_r_max
        r_E_2=r(nR)*r(nR)
        z(l1m1,nR)  =z(l1m1,nR)  -  rho0(nR)*r_E_2*corr_l1m1
        dz(l1m1,nR) =dz(l1m1,nR) -  rho0(nR)*( &
             2.d0*r(nR)+r_E_2*beta(nR))*corr_l1m1
        workA(l1m1,nR)=workA(l1m1,nR)-rho0(nR)*( &
             2.d0+4.d0*beta(nR)*r(nR) + &
             dbeta(nR)*r_E_2 + &
             beta(nR)*beta(nR)*r_E_2 )*corr_l1m1
     END DO
     !$OMP END PARALLEL DO
  END IF ! l=1,m=1 contained in lm-block ?

  !IF (DEBUG_OUTPUT) THEN
  !   DO nR=1,n_r_max
  !      IF ((nR == 1).OR.(nR == 5)) THEN
  !         DO lm=llm,ulm
  !            WRITE(*,"(4X,A,2I3,4ES22.14)") "upZ_new: ",nR,lm,z(lm,nR),dz(lm,nR)
  !         END DO
  !      END IF
  !      WRITE(*,"(A,I3,4ES22.14)") "upZ_new: ",nR,get_global_SUM( z(:,nR) ),get_global_SUM( dz(:,nR) )
  !   END DO
  !END IF
  !-- Calculate explicit time step part:
  !$OMP PARALLEL default(shared) &
  !$OMP private(nR,lm1,Dif)
  !$OMP DO
  DO nR=n_r_cmb+1,n_r_icb-1
     !WRITE(*,"(A,I4,5ES20.12)") "r-dependent : ",nR,dLvisc(nR),beta(nR),or1(nR),or2(nR),dbeta(nR)
     DO lm1=lmStart_00,lmStop
        Dif(lm1)=hdif_V(st_map%lm2(lm2l(lm1),lm2m(lm1)))*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*visc(nR)* &
             & ( workA(lm1,nR) &
             &   +(dLvisc(nR)-beta(nR))*dz(lm1,nR) &
             &   -( dLvisc(nR)*beta(nR) &
             &      + 2.d0*dLvisc(nR)*or1(nR) &
             &      + dLh(st_map%lm2(lm2l(lm1),lm2m(lm1))) * or2(nR)&
             &      + dbeta(nR)&
             &      + 2.d0*beta(nR)*or1(nR) &
             &    ) * z(lm1,nR) &
             & )

!        IF (nR == 2) THEN
!           WRITE(*,"(2I4,8ES20.12)") nR,lm1,workA(lm1,nR),Dif(lm1),z(lm1,nR),dz(lm1,nR)
!        END IF
        dzdtLast(lm1,nR)=dzdt(lm1,nR)-coex*Dif(lm1)
        IF ( lRmsNext ) THEN
           workB(lm1,nR)= &
                O_dt*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*(z(lm1,nR)-workB(lm1,nR))
           IF ( l_RMStest ) &
                workB(lm1,nR)= workB(lm1,nR)-Dif(lm1)
        END IF
     END DO
     IF ( lRmsNext ) THEN
        CALL hInt2Tor(Dif,1,lm_max,nR,lmStart_00,lmStop, &
             DifTor2hInt(nR,1),DifTorAs2hInt(nR,1),lo_map)
        CALL hInt2Tor(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
             dtVTor2hInt(nR,1),dtVTorAs2hInt(nR,1),lo_map)
     END IF
  END DO
  !$OMP END DO
  !--- Note: from ddz=workA only the axisymmetric contributions are needed
  !    beyond this point for the TO calculation.
  !    Parallization note: Very likely, all axisymmetric modes m=0 are
  !    located on the first processor #0.
  IF ( l_TO ) THEN
     !$OMP DO private(nR,lm1,l1,m1)
     DO nR=1,n_r_max
        DO lm1=lmStart_00,lmStop
           l1=lm2l(lm1)
           m1=lm2m(lm1)
           IF ( m1 == 0 ) ddzASL(l1+1,nR)=REAL(workA(lm1,nR))
        END DO
     END DO
     !$OMP END DO
  END IF
  !$OMP END PARALLEL

  !----- Special thing for l=1,m=0 for rigid boundaries and
  !      if IC or mantle are allowed to rotate:
  IF ( l10 .AND. l_z10mat ) THEN ! z10 term !
     lm1=lm2(1,0)
     !----- NOTE opposite sign of visouse torque on ICB and CMB:
     IF ( .NOT. l_SRMA .AND. ktopv == 2 .AND. l_rot_ma ) THEN
        d_omega_ma_dtLast=d_omega_ma_dt -            &
             coex * ( 2.D0*or1(1)*REAL(z(lm1,1))  - &
             REAL(dz(lm1,1)) )
     END IF
     IF ( .NOT. l_SRIC .AND. kbotv == 2 .AND. l_rot_ic ) THEn
        d_omega_ic_dtLast=d_omega_ic_dt +                        &
             coex * ( 2.D0*or1(n_r_max)*REAL(z(lm1,n_r_max))  - &
             REAL(dz(lm1,n_r_max)) )
     END IF
  END IF
  !PERFOFF
end SUBROUTINE updateZ

END MODULE updateZ_mod
