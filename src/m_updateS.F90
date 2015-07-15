!$Id$
!*************************************************************************
#include "perflib_preproc.cpp"
MODULE updateS_mod
  use omp_lib
  USE truncation
  USE radial_functions,ONLY: i_costf_init,d_costf_init,orho1,or1,or2, &
                            &           n_r_cmb,n_r_icb,beta,drx,ddrx,&
                            &           kappa,dlkappa,dtemp0,otemp1,  &
                            &           temp0,dentropy0
  USE physical_parameters,ONLY: opr
  USE init_fields,ONLY: tops,bots
  USE blocking,ONLY: nLMBs,st_map,lo_map,lo_sub_map,lmStartB,lmStopB
  USE horizontal_data,ONLY: dLh,hdif_S
  USE logic,only: l_update_s
  USE matrices,only: lSmat,s0Mat,s0Pivot,&
#ifdef WITH_PRECOND_S
       & sMat_fac, &
#endif
#ifdef WITH_PRECOND_S0
       & s0Mat_fac, &
#endif
       & sMat,sPivot

  USE LMLoop_data, ONLY: llm,ulm,llm_real,ulm_real
  USE parallel_mod,ONLY: rank,chunksize
#ifdef WITH_MKL_LU
  USE lapack95, ONLY: getrs
#else
  USE algebra, ONLY: cgeslML,sgesl
#endif
  IMPLICIT NONE

  PRIVATE
  !-- Local work arrays
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: workA,workB
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: rhs1

  integer :: maxThreads
  PUBLIC :: initialize_updateS,updateS,updateS_ala

CONTAINS
  SUBROUTINE initialize_updateS
    ALLOCATE(workA(llm:ulm,n_r_max))
    ALLOCATE(workB(llm:ulm,n_r_max))

    !COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! complex RHS for l>0
#ifdef WITHOMP
    maxThreads=omp_get_max_threads()
#else
    maxThreads=1
#endif
    ALLOCATE(rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
  END SUBROUTINE initialize_updateS
  
  SUBROUTINE updateS(s,ds,dVSrLM,dsdt,dsdtLast, &
       &             w1,coex,dt,nLMB)
    !*************************************************************************

    !-------------------------------------------------------------------------

    !  updates the entropy field s and its radial derivatives
    !  adds explicit part to time derivatives of s

    !-------------------------------------------------------------------------


    !-- Input of variables:
    REAL(kind=8),intent(IN) :: w1        ! weight for time step !
    REAL(kind=8),intent(IN) :: coex      ! factor depending on alpha
    REAL(kind=8),intent(IN) :: dt        ! time step
    INTEGER,intent(IN) :: nLMB

    !-- Input/output of scalar fields:
    COMPLEX(kind=8),intent(INOUT) :: s(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(OUT) :: ds(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dVSrLM(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dsdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dsdtLast(llm:ulm,n_r_max)
    !-- Output: udpated s,ds,dsdtLast

    !-- Local variables:
    REAL(kind=8) :: w2            ! weight of second time step
    REAL(kind=8) :: O_dt
    INTEGER :: l1,m1              ! degree and order
    INTEGER :: lm1,lmB,lm         ! position of (l,m) in array
    INTEGER :: lmStart_real       ! range of lm for real array
    INTEGER :: lmStop_real        !
    INTEGER :: lmStart,lmStop
    INTEGER :: nLMB2
    INTEGER :: nR                 ! counts radial grid points
    INTEGER :: n_cheb             ! counts cheb modes
    REAL(kind=8) ::  rhs(n_r_max) ! real RHS for l=m=0
    !COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! complex RHS for l>0

    INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
    INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
    INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

    INTEGER :: threadid,nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
    INTEGER :: iChunk,nChunks,size_of_last_chunk,lmB0
    !-- end of declaration
    !---------------------------------------------------------------------

    IF ( .NOT. l_update_s ) RETURN

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
    lmStart_real=2*lmStart-1
    lmStop_real =2*lmStop
    w2  =1.-w1
    O_dt=1.D0/dt


    !PERFON('upS_fin')
    !$OMP PARALLEL default(none) &
    !$OMP private(iThread,start_lm,stop_lm,nR,lm) &
    !$OMP shared(all_lms,per_thread,lmStart_real,lmStop_real) &
    !$OMP shared(dVSrLM,i_costf_init,d_costf_init,drx,dsdt,orho1,or2,lmStart,lmStop) &
    !$OMP shared(n_r_max,n_cheb_max,workA,workB,nThreads,llm_real,ulm_real)
    !$OMP SINGLE
#ifdef WITHOMP
    nThreads=omp_get_num_threads()
#else
    nThreads=1
#endif
    !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
    all_lms=lmStop_real-lmStart_real+1
    per_thread=all_lms/nThreads
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP DO
    DO iThread=0,nThreads-1
       start_lm=lmStart_real+iThread*per_thread
       stop_lm = start_lm+per_thread-1
       IF (iThread == nThreads-1) stop_lm=lmStop_real

       !--- Finish calculation of dsdt:
       CALL get_drNS( dVSrLM,workA, &
            &         ulm_real-llm_real+1,start_lm-llm_real+1,stop_lm-llm_real+1, &
            &         n_r_max,n_cheb_max,workB, &
            &         i_costf_init,d_costf_init,drx)
    END DO
    !$OMP END DO

    !$OMP DO
    DO nR=1,n_r_max
       DO lm=lmStart,lmStop
          dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*workA(lm,nR))
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    !PERFOFF

    ! one subblock is linked to one l value and needs therefore once the matrix
    !$OMP PARALLEL default(shared)
    !$OMP SINGLE
    DO nLMB2=1,nLMBs2(nLMB)
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
       !WRITE(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

       IF ( l1 == 0 ) THEN
          IF ( .NOT. lSmat(l1) ) THEN
#ifdef WITH_PRECOND_S0
             CALL get_s0Mat(dt,s0Mat,s0Pivot,s0Mat_fac)
#else
             CALL get_s0Mat(dt,s0Mat,s0Pivot)
#endif
             lSmat(l1)=.TRUE.
          END IF
       ELSE
          IF ( .NOT. lSmat(l1) ) THEN
#ifdef WITH_PRECOND_S
             CALL get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                  sMat(1,1,l1),sPivot(1,l1),sMat_fac(1,l1))
#else
             CALL get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                  sMat(1,1,l1),sPivot(1,l1))
#endif
             lSmat(l1)=.TRUE.
             !WRITE(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
          END IF
       END IF

       DO iChunk=1,nChunks
          !$OMP TASK default(shared) &
          !$OMP firstprivate(iChunk) &
          !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_cheb) &
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

             IF ( l1 == 0 ) THEN
                rhs(1)=      REAL(tops(0,0))
                rhs(n_r_max)=REAL(bots(0,0))
                DO nR=2,n_r_max-1
                   rhs(nR)=REAL(s(lm1,nR))*O_dt+ &
                        w1*REAL(dsdt(lm1,nR)) + &
                        w2*REAL(dsdtLast(lm1,nR))
                END DO

#ifdef WITH_PRECOND_S0
                rhs = s0Mat_fac*rhs
#endif

#ifdef WITH_MKL_LU
                CALL getrs(s0Mat,s0Pivot,rhs)
#else
                CALL sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
#endif

             ELSE ! l1 .ne. 0
                lmB=lmB+1

                rhs1(1,lmB,threadid)=      tops(l1,m1)
                rhs1(n_r_max,lmB,threadid)=bots(l1,m1)
#ifdef WITH_PRECOND_S
                rhs1(1,lmB,threadid)=      sMat_fac(1,l1)*rhs1(1,lmB,threadid)
                rhs1(n_r_max,lmB,threadid)=sMat_fac(1,l1)*rhs1(n_r_max,lmB,threadid)
#endif
                DO nR=2,n_r_max-1
                   rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt + &
                        w1*dsdt(lm1,nR) + &
                        w2*dsdtLast(lm1,nR)
#ifdef WITH_PRECOND_S
                   rhs1(nR,lmB,threadid) = sMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                END DO
             END IF
          END DO
          !PERFOFF

          !PERFON('upS_sol')
          IF ( lmB  >  lmB0 ) THEN
#ifdef WITH_MKL_LU
             CALL getrs(CMPLX(sMat(:,:,l1),0.D0,KIND=KIND(0.d0)), &
                  &       sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid))
#else
             CALL cgeslML(sMat(:,:,l1),n_r_max,n_r_max, &
                  &       sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),n_r_max,lmB-lmB0)
#endif
          END IF
          !PERFOFF

          lmB=lmB0
          !PERFON('upS_af')
          DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
           !DO lm=1,sizeLMB2(nLMB2,nLMB)
             lm1=lm22lm(lm,nLMB2,nLMB)
             !l1 =lm22l(lm,nLMB2,nLMB)
             m1 =lm22m(lm,nLMB2,nLMB)
             IF ( l1 == 0 ) THEN
                DO n_cheb=1,n_cheb_max
                   s(lm1,n_cheb)=rhs(n_cheb)
                END DO
             ELSE
                lmB=lmB+1
                IF ( m1 > 0 ) THEN
                   DO n_cheb=1,n_cheb_max
                      s(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                   END DO
                ELSE
                   DO n_cheb=1,n_cheb_max
                      s(lm1,n_cheb)= CMPLX(REAL(rhs1(n_cheb,lmB,threadid)),0.D0,KIND=KIND(0d0))
                   END DO
                END IF
             END IF
          END DO
          !PERFOFF
          !$OMP END TASK
       END DO
       !$OMP END TASK
    END DO     ! loop over lm blocks
    !$OMP END SINGLE
    !$OMP END PARALLEL

    !WRITE(*,"(A,2ES22.12)") "s after = ",SUM(s)
    !-- set cheb modes > n_cheb_max to zero (dealiazing)
    DO n_cheb=n_cheb_max+1,n_r_max
       DO lm1=lmStart,lmStop
          s(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       END DO
    END DO

    !PERFON('upS_drv')
    all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
    IF (all_lms < maxThreads) THEN
       call omp_set_num_threads(all_lms)
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
    !$OMP shared(s,ds,dsdtLast,i_costf_init,d_costf_init,drx,ddrx) &
    !$OMP shared(n_r_max,n_cheb_max,workA,workB,llm_real,ulm_real) &
    !$OMP shared(n_r_cmb,n_r_icb,lmStart,lmStop,dsdt,coex,opr,hdif_S) &
    !$OMP shared(st_map,lm2l,lm2m,kappa,beta,otemp1,dtemp0,or1,dLkappa,dLh,or2)
    !$OMP DO
    DO iThread=0,nThreads-1
       start_lm=lmStart_real+iThread*per_thread
       stop_lm = start_lm+per_thread-1
       IF (iThread == nThreads-1) stop_lm=lmStop_real
       CALL costf1(s, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &      dsdtLast, i_costf_init, d_costf_init)
       CALL get_ddr(s, ds, workA, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &       n_r_max, n_cheb_max, workB, dsdtLast, &
            &       i_costf_init,d_costf_init,drx,ddrx)
    END DO
    !$OMP END DO

    !-- Calculate explicit time step part:
    !$OMP DO private(nR,lm1)
    DO nR=n_r_cmb+1,n_r_icb-1
       DO lm1=lmStart,lmStop
          dsdtLast(lm1,nR)=dsdt(lm1,nR) &
               & - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1))) * kappa(nR) * &
               &   ( workA(lm1,nR) &
               &     + ( beta(nR) + otemp1(nR)*dtemp0(nR) + &
               &       2.D0*or1(nR) + dLkappa(nR) ) * ds(lm1,nR) &
               &     - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1))) * or2(nR)   *  s(lm1,nR) &
               &   )
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
#ifdef WITHOMP
    call omp_set_num_threads(maxThreads)
#endif
    !PERFOFF
    !-- workA=dds not needed further after this point, used as work array later


    RETURN
  end SUBROUTINE updateS
  !*************************************************************************
  SUBROUTINE updateS_ala(s,ds,w,dVSrLM,dsdt,dsdtLast,w1,coex,dt,nLMB)
  !*************************************************************************

    !-------------------------------------------------------------------------

    !  updates the entropy field s and its radial derivatives
    !  adds explicit part to time derivatives of s

    !-------------------------------------------------------------------------


    USE num_param, ONLY: alpha

    IMPLICIT NONE

    !-- Input of variables:
    REAL(kind=8),intent(IN) :: w1        ! weight for time step !
    REAL(kind=8),intent(IN) :: coex      ! factor depending on alpha
    REAL(kind=8),intent(IN) :: dt        ! time step
    INTEGER,intent(IN) :: nLMB

    !-- Input/output of scalar fields:
    COMPLEX(kind=8),intent(INOUT) :: s(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(OUT) :: ds(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dVSrLM(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dsdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dsdtLast(llm:ulm,n_r_max)
    !-- Output: udpated s,ds,dsdtLast

    !-- Local variables:
    REAL(kind=8) :: w2            ! weight of second time step
    REAL(kind=8) :: O_dt
    INTEGER :: l1,m1              ! degree and order
    INTEGER :: lm1,lmB,lm         ! position of (l,m) in array
    INTEGER :: lmStart_real       ! range of lm for real array
    INTEGER :: lmStop_real        !
    INTEGER :: lmStart,lmStop
    INTEGER :: nLMB2
    INTEGER :: nR                 ! counts radial grid points
    INTEGER :: n_cheb             ! counts cheb modes
    REAL(kind=8) ::  rhs(n_r_max) ! real RHS for l=m=0
    !COMPLEX(kind=8) :: rhs1(n_r_max,lo_sub_map%sizeLMB2max) ! complex RHS for l>0

    INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
    INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
    INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

    INTEGER :: threadid,nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
    INTEGER :: iChunk,nChunks,size_of_last_chunk,lmB0
    !-- end of declaration
    !---------------------------------------------------------------------

    IF ( .NOT. l_update_s ) RETURN

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
    lmStart_real=2*lmStart-1
    lmStop_real =2*lmStop
    w2  =1.-w1
    O_dt=1.D0/dt


    !PERFON('upS_fin')
    !$OMP PARALLEL default(none) &
    !$OMP private(iThread,start_lm,stop_lm,nR,lm) &
    !$OMP shared(all_lms,per_thread,lmStart_real,lmStop_real) &
    !$OMP shared(dVSrLM,i_costf_init,d_costf_init,drx,dsdt,orho1) &
    !$OMP shared(otemp1,dtemp0,or2,lmStart,lmStop) &
    !$OMP shared(n_r_max,n_cheb_max,workA,workB,nThreads,llm_real,ulm_real)
    !$OMP SINGLE
#ifdef WITHOMP
    nThreads=omp_get_num_threads()
#else
    nThreads=1
#endif
    !-- Get radial derivatives of s: workA,dsdtLast used as work arrays
    all_lms=lmStop_real-lmStart_real+1
    per_thread=all_lms/nThreads
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP DO
    DO iThread=0,nThreads-1
       start_lm=lmStart_real+iThread*per_thread
       stop_lm = start_lm+per_thread-1
       IF (iThread == nThreads-1) stop_lm=lmStop_real

       !--- Finish calculation of dsdt:
       CALL get_drNS( dVSrLM,workA, &
            &         ulm_real-llm_real+1,start_lm-llm_real+1,stop_lm-llm_real+1, &
            &         n_r_max,n_cheb_max,workB, &
            &         i_costf_init,d_costf_init,drx)
    END DO
    !$OMP END DO

    !$OMP DO
    DO nR=1,n_r_max
       DO lm=lmStart,lmStop
          dsdt(lm,nR)=          orho1(nR)*dsdt(lm,nR)  - & 
              &         or2(nR)*orho1(nR)*workA(lm,nR) + &
              &       orho1(nR)*otemp1(nR)*dtemp0(nR)*dVSrLM(lm,nR)
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    !PERFOFF

    ! one subblock is linked to one l value and needs therefore once the matrix
    !$OMP PARALLEL default(shared)
    !$OMP SINGLE
    DO nLMB2=1,nLMBs2(nLMB)
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
       !WRITE(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

       IF ( l1 == 0 ) THEN
          IF ( .NOT. lSmat(l1) ) THEN
#ifdef WITH_PRECOND_S0
             CALL get_s0Mat(dt,s0Mat,s0Pivot,s0Mat_fac)
#else
             CALL get_s0Mat(dt,s0Mat,s0Pivot)
#endif
             lSmat(l1)=.TRUE.
          END IF
       ELSE
          IF ( .NOT. lSmat(l1) ) THEN
#ifdef WITH_PRECOND_S
             CALL get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                  sMat(1,1,l1),sPivot(1,l1),sMat_fac(1,l1))
#else
             CALL get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                  sMat(1,1,l1),sPivot(1,l1))
#endif
             lSmat(l1)=.TRUE.
             !WRITE(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
          END IF
       END IF

       DO iChunk=1,nChunks
          !$OMP TASK default(shared) &
          !$OMP firstprivate(iChunk) &
          !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_cheb) &
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

             IF ( l1 == 0 ) THEN
                rhs(1)=      REAL(tops(0,0))
                rhs(n_r_max)=REAL(bots(0,0))
                DO nR=2,n_r_max-1
                   rhs(nR)=REAL(s(lm1,nR))*O_dt+ &
                        w1*REAL(dsdt(lm1,nR)) + &
                        w2*REAL(dsdtLast(lm1,nR))
                END DO

#ifdef WITH_PRECOND_S0
                rhs = s0Mat_fac*rhs
#endif

#ifdef WITH_MKL_LU
                CALL getrs(s0Mat,s0Pivot,rhs)
#else
                CALL sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
#endif

             ELSE ! l1 .ne. 0
                lmB=lmB+1

                rhs1(1,lmB,threadid)=      tops(l1,m1)
                rhs1(n_r_max,lmB,threadid)=bots(l1,m1)
#ifdef WITH_PRECOND_S
                rhs1(1,lmB,threadid)=      sMat_fac(1,l1)*rhs1(1,lmB,threadid)
                rhs1(n_r_max,lmB,threadid)=sMat_fac(1,l1)*rhs1(n_r_max,lmB,threadid)
#endif
                DO nR=2,n_r_max-1
                   rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt +            &
                       &                w1*dsdt(lm1,nR) +            &
                       &            w2*dsdtLast(lm1,nR) -            &
                       &  alpha*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1))) &
                       &  *or2(nR)*orho1(nR)*temp0(nR)*              &
                       &        dentropy0(nR)*w(lm1,nR)
#ifdef WITH_PRECOND_S
                   rhs1(nR,lmB,threadid) = sMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                END DO
             END IF
          END DO
          !PERFOFF

          !PERFON('upS_sol')
          IF ( lmB  >  lmB0 ) THEN
#ifdef WITH_MKL_LU
             CALL getrs(CMPLX(sMat(:,:,l1),0.D0,KIND=KIND(0.d0)),sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid))
#else
             CALL cgeslML(sMat(:,:,l1),n_r_max,n_r_max, &
                  &       sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),n_r_max,lmB-lmB0)
#endif
          END IF
          !PERFOFF

          lmB=lmB0
          !PERFON('upS_af')
          DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
           !DO lm=1,sizeLMB2(nLMB2,nLMB)
             lm1=lm22lm(lm,nLMB2,nLMB)
             !l1 =lm22l(lm,nLMB2,nLMB)
             m1 =lm22m(lm,nLMB2,nLMB)
             IF ( l1 == 0 ) THEN
                DO n_cheb=1,n_cheb_max
                   s(lm1,n_cheb)=rhs(n_cheb)
                END DO
             ELSE
                lmB=lmB+1
                IF ( m1 > 0 ) THEN
                   DO n_cheb=1,n_cheb_max
                      s(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                   END DO
                ELSE
                   DO n_cheb=1,n_cheb_max
                      s(lm1,n_cheb)= CMPLX(REAL(rhs1(n_cheb,lmB,threadid)),0.D0,KIND=KIND(0d0))
                   END DO
                END IF
             END IF
          END DO
          !PERFOFF
          !$OMP END TASK
       END DO
       !$OMP END TASK
    END DO     ! loop over lm blocks
    !$OMP END SINGLE
    !$OMP END PARALLEL

    !WRITE(*,"(A,2ES22.12)") "s after = ",SUM(s)
    !-- set cheb modes > n_cheb_max to zero (dealiazing)
    DO n_cheb=n_cheb_max+1,n_r_max
       DO lm1=lmStart,lmStop
          s(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       END DO
    END DO

    !PERFON('upS_drv')
    all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
    IF (all_lms < maxThreads) THEN
       call omp_set_num_threads(all_lms)
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
    !$OMP shared(s,ds,w,dsdtLast,i_costf_init,d_costf_init,drx,ddrx) &
    !$OMP shared(n_r_max,n_cheb_max,workA,workB,llm_real,ulm_real,temp0) &
    !$OMP shared(n_r_cmb,n_r_icb,lmStart,lmStop,dsdt,coex,opr,hdif_S,dentropy0) &
    !$OMP shared(st_map,lm2l,lm2m,kappa,beta,otemp1,dtemp0,or1,dLkappa,dLh,or2) &
    !$OMP shared(orho1)
    !$OMP DO
    DO iThread=0,nThreads-1
       start_lm=lmStart_real+iThread*per_thread
       stop_lm = start_lm+per_thread-1
       IF (iThread == nThreads-1) stop_lm=lmStop_real
       CALL costf1(s, ulm_real-llm_real+1, start_lm-llm_real+1, &
            &      stop_lm-llm_real+1,dsdtLast, i_costf_init,   &
            &      d_costf_init)
       CALL get_ddr(s, ds, workA, ulm_real-llm_real+1,       &
            &       start_lm-llm_real+1, stop_lm-llm_real+1, &
            &       n_r_max, n_cheb_max, workB, dsdtLast,    &
            &       i_costf_init,d_costf_init,drx,ddrx)
    END DO
    !$OMP END DO

    !-- Calculate explicit time step part:
    !$OMP DO private(nR,lm1)
    DO nR=n_r_cmb+1,n_r_icb-1
       DO lm1=lmStart,lmStop
         dsdtLast(lm1,nR)=dsdt(lm1,nR) &
              & - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1)))*kappa(nR) * &
              &   (                                              workA(lm1,nR) &
              &           + ( beta(nR)+2.D0*or1(nR)+dLkappa(nR) ) * ds(lm1,nR) &
              &     - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*  s(lm1,nR) &
              &   ) + coex*dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)        &
              &           *orho1(nR)*temp0(nR)*dentropy0(nR)*        w(lm1,nR)
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
#ifdef WITHOMP
    call omp_set_num_threads(maxThreads)
#endif
    !PERFOFF
    !-- workA=dds not needed further after this point, used as work array later


    RETURN
  end SUBROUTINE updateS_ala
  !-------------------------------------------------------------------------------
END MODULE updateS_mod
!-------------------------------------------------------------------------------
