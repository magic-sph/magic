!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"
MODULE updateWP_mod
  use omp_lib
  USE truncation
  USE radial_data,ONLY: n_r_cmb,n_r_icb
  USE radial_functions, ONLY: drx,ddrx,dddrx,or1,or2,rho0,agrav,rgrav,&
       &i_costf_init,d_costf_init,&
       &visc,dlvisc,beta,dbeta
  USE physical_parameters
  USE num_param
  USE blocking,only: nLMBs,lo_sub_map,lo_map,st_map,st_sub_map,&
       &lmStartB,lmStopB
  USE horizontal_data
  USE logic
  USE matrices
  USE RMS
#ifdef WITH_MKL_LU
  USE lapack95, ONLY: getrs
#else
  USE algebra, ONLY: cgeslML
#endif
  USE LMLoop_data, ONLY:llm,ulm, llm_real,ulm_real
  USE communications, only: get_global_sum
  USE parallel_mod,only: chunksize
  IMPLICIT NONE

  PRIVATE
  !-- Input of recycled work arrays:
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:) :: workA,workB
  !COMPLEX(kind=8),DIMENSION(:,:,:),ALLOCATABLE :: rhs1
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:) :: Dif,Pre,Buo
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:) :: rhs1  
  INTEGER :: maxThreads
  
  PUBLIC :: initialize_updateWP, updateWP

contains
  SUBROUTINE initialize_updateWP
    ALLOCATE(workA(llm:ulm,n_r_max))
    ALLOCATE(workB(llm:ulm,n_r_max))
    ALLOCATE(Dif(llm:ulm))
    ALLOCATE(Pre(llm:ulm))
    ALLOCATE(Buo(llm:ulm))
#ifdef WITHOMP
    maxThreads=omp_get_max_threads()
#else
    maxThreads=1
#endif

    ALLOCATE(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1))
  END SUBROUTINE initialize_updateWP

  SUBROUTINE finalize_updateWP
    DEALLOCATE(workA)
    DEALLOCATE(workB)
    DEALLOCATE(Dif)
    DEALLOCATE(Pre)
    DEALLOCATE(Buo)
    DEALLOCATE(rhs1)
  END SUBROUTINE finalize_updateWP

  SUBROUTINE updateWP(w,dw,ddw,dwdt,dwdtLast, &
       &              p,dp,dpdt,dpdtLast,s, &
       &              w1,coex,dt,nLMB,lRmsNext)
    !***********************************************************************

    !-----------------------------------------------------------------------

    !  updates the poloidal velocity potential w, the pressure p,  and
    !  their derivatives
    !  adds explicit part to time derivatives of w and p

    !-------------------------------------------------------------------------


    !-- Input/output of scalar fields:
    COMPLEX(kind=8),INTENT(IN) :: dwdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dpdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),INTENT(IN) :: s(llm:ulm,n_r_max)

    COMPLEX(kind=8),INTENT(INOUT) :: w(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dw(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(OUT) :: ddw(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dwdtLast(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: p(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(OUT) :: dp(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(INOUT) :: dpdtLast(llm:ulm,n_r_max)
    !-- Output: updated w,dw,ddw,p,dp,dwdtLast,dpdtLast

    !-- Input of other variables:
    REAL(kind=8),intent(IN) :: w1        ! weight for time step !
    REAL(kind=8),intent(IN) :: coex      ! factor depending on alpha
    REAL(kind=8),intent(IN) :: dt        ! time step
    INTEGER,intent(IN) :: nLMB     ! block number
    LOGICAL,intent(IN) :: lRmsNext

    !-- Local variables:
    REAL(kind=8) :: w2                  ! weight of second time step
    REAL(kind=8) :: O_dt
    INTEGER :: l1,m1              ! degree and order
    INTEGER :: lm1,lm,lmB         ! position of (l,m) in array
    INTEGER :: lmStart,lmStop ! max and min number of orders m
    INTEGER :: lmStart_real      ! range of lm for real array
    INTEGER :: lmStop_real       !
    INTEGER :: lmStart_00        ! excluding l=0,m=0
    INTEGER :: nLMB2
    INTEGER :: nR                ! counts radial grid points
    INTEGER :: n_cheb             ! counts cheb modes

    !COMPLEX(kind=8) :: Dif(llm:ulm),Pre(llm:ulm),Buo(llm:ulm)

    !COMPLEX(kind=8) :: rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,lo_sub_map%nLMBs2(nLMB))
    !COMPLEX(kind=8),DIMENSION(:,:,:),allocatable :: rhs1
    !COMPLEX(kind=8) :: rhs1_sum,rhs2_sum

    INTEGER, DIMENSION(:),POINTER :: nLMBs2,lm2l,lm2m
    INTEGER, DIMENSION(:,:),POINTER :: sizeLMB2,lm2
    INTEGER, DIMENSION(:,:,:),POINTER :: lm22lm,lm22l,lm22m

    INTEGER :: iThread,start_lm,stop_lm,all_lms,per_thread,nThreads
    INTEGER :: nChunks,iChunk,lmB0,size_of_last_chunk,threadid
    !-- end of declaration
    !---------------------------------------------------------------------

    IF ( .NOT. l_update_v ) RETURN

    
    nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
    sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
    lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
    lm22l(1:,1:,1:) => lo_sub_map%lm22l
    lm22m(1:,1:,1:) => lo_sub_map%lm22m
    lm2(0:,0:) => lo_map%lm2
    lm2l(1:lm_max) => lo_map%lm2l
    lm2m(1:lm_max) => lo_map%lm2m

    !ALLOCATE(rhs1(2*n_r_max,lo_sub_map%sizeLMB2max,nLMBs2(nLMB)))

    lmStart     =lmStartB(nLMB)
    lmStop      =lmStopB(nLMB)
    lmStart_00  =MAX(2,lmStart)
    lmStart_real=2*lmStart_00-1
    lmStop_real =2*lmStop

    w2  =1.D0-w1
    O_dt=1.D0/dt
    !PERFON('upWP_ssol')
    !$OMP PARALLEL default(shared) &
    !$OMP private(nLMB2,lm,lm1,l1,m1,lmB)
    !WRITE(*,"(I3,A)") omp_get_thread_num(),": before SINGLE"
    !$OMP SINGLE
    ! each of the nLMBs2(nLMB) subblocks have one l value
    DO nLMB2=1,nLMBs2(nLMB)
       !WRITE(*,"(2(A,I3))") "Constructing next task for ",nLMB2,"/",nLMBs2(nLMB)

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
       IF ( l1 > 0 ) THEN
          IF ( .NOT. lWPmat(l1) ) THEN
             !PERFON('upWP_mat')
             CALL get_wpMat(dt,l1,hdif_V(st_map%lm2(l1,0)), &
                  wpMat(1,1,l1),wpPivot(1,l1),wpMat_fac(1,1,l1))
             lWPmat(l1)=.TRUE.
             !PERFOFF
          END IF
          !WRITE(*,"(2(A,I3))") "Running ",nChunks," for task with l1=",l1
          DO iChunk=1,nChunks
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
             DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
             !DO lm=1,sizeLMB2(nLMB2,nLMB)
                lm1=lm22lm(lm,nLMB2,nLMB)
                m1 =lm22m(lm,nLMB2,nLMB)

                lmB=lmB+1
                rhs1(1,lmB,threadid)        =0.D0
                rhs1(n_r_max,lmB,threadid)  =0.D0
                rhs1(n_r_max+1,lmB,threadid)=0.D0
                rhs1(2*n_r_max,lmB,threadid)=0.D0
                DO nR=2,n_r_max-1
                   rhs1(nR,lmB,threadid)=                         &
                        & O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*w(lm1,nR) + &
                        & rho0(nR)*agrav(nR)*s(lm1,nR) + &
                        & w1*dwdt(lm1,nR) + &
                        & w2*dwdtLast(lm1,nR)
                   rhs1(nR+n_r_max,lmB,threadid)=                 &
                        -O_dt*dLh(st_map%lm2(l1,m1))*or2(nR)*dw(lm1,nR) + &
                        w1*dpdt(lm1,nR) + &
                        w2*dpdtLast(lm1,nR)
                END DO
             END DO
             !PERFOFF
             !PERFON('upWP_sol')
             !IF ( lmB > 0 ) THEN

             ! use the mat_fac(:,1) to scale the rhs
             DO lm=lmB0+1,lmB
                DO nR=1,2*n_r_max
                   rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,1,l1)
                END DO
             END DO
#ifdef WITH_MKL_LU
             CALL getrs(CMPLX(wpMat(:,:,l1),0.D0,KIND=KIND(0.D0)), &
                  &       wpPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid))
#else
             CALL cgeslML(wpMat(:,:,l1),2*n_r_max,2*n_r_max,    &
                  &       wpPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),2*n_r_max,lmB-lmB0)
#endif
             ! rescale the solution with mat_fac(:,2)
             DO lm=lmB0+1,lmB
                DO nR=1,2*n_r_max
                   rhs1(nR,lm,threadid)=rhs1(nR,lm,threadid)*wpMat_fac(nR,2,l1)
                END DO
             END DO
          !END IF
             !PERFOFF

             IF ( lRmsNext ) THEN ! Store old w
                DO nR=1,n_r_max
                   DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                      lm1=lm22lm(lm,nLMB2,nLMB)
                      workB(lm1,nR)=w(lm1,nR)
                   END DO
                END DO
             END IF

             !PERFON('upWP_aft')
             lmB=lmB0
             DO lm=lmB0+1,MIN(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
                lm1=lm22lm(lm,nLMB2,nLMB)
                !l1 =lm22l(lm,nLMB2,nLMB)
                m1 =lm22m(lm,nLMB2,nLMB)
                !IF ( l1 > 0 ) THEN
                lmB=lmB+1
                IF ( m1 > 0 ) THEN
                   DO n_cheb=1,n_cheb_max
                      w(lm1,n_cheb)=rhs1(n_cheb,lmB,threadid)
                      p(lm1,n_cheb)=rhs1(n_r_max+n_cheb,lmB,threadid)
                   END DO
                ELSE
                   DO n_cheb=1,n_cheb_max
                      w(lm1,n_cheb)= &
                           CMPLX(REAL(rhs1(n_cheb,lmB,threadid)),0.D0,KIND=KIND(0d0))
                      p(lm1,n_cheb)= &
                           CMPLX(REAL(rhs1(n_r_max+n_cheb,lmB,threadid)),0.D0,KIND=KIND(0d0))
                   END DO
                END IF
             END DO
             !PERFOFF
             !$OMP END TASK
          END DO
       END IF
       !WRITE(*,"(3(A,I3))") "End of task ",nLMB2,"/",nLMBs2(nLMB)," on thread ",omp_get_thread_num()
       !$OMP END TASK
    END DO   ! end of loop over l1 subblocks
    !$OMP END SINGLE
    !$OMP END PARALLEL
    !PERFOFF
    !WRITE(*,"(A,I3,4ES22.12)") "w,p after: ",nLMB,get_global_SUM(w),get_global_SUM(p)

    !-- set cheb modes > n_cheb_max to zero (dealiazing)
    DO n_cheb=n_cheb_max+1,n_r_max
       DO lm1=lmStart_00,lmStop
          w(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
          p(lm1,n_cheb)=CMPLX(0.D0,0.D0,KIND=KIND(0d0))
       END DO
    END DO


    !PERFON('upWP_drv')
    all_lms=lmStop_real-lmStart_real+1
#ifdef WITHOMP
    IF (all_lms < omp_get_max_threads()) THEN
       call omp_set_num_threads(all_lms)
    END IF
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
    DO iThread=0,nThreads-1
       start_lm=lmStart_real+iThread*per_thread
       stop_lm = start_lm+per_thread-1
       IF (iThread == nThreads-1) stop_lm=lmStop_real
       !WRITE(*,"(2(A,I3),2(A,I5))") "iThread=",iThread," on thread ",omp_get_thread_num(),&
       !     & " lm = ",start_lm,":",stop_lm

       !-- Transform to radial space and get radial derivatives
       !   using dwdtLast, dpdtLast as work arrays:
       CALL costf1( w, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &       dwdtLast,i_costf_init,d_costf_init)
       CALL get_dddr( w, dw, ddw, workA, &
            &         ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &         n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
            &         i_costf_init,d_costf_init,drx,ddrx,dddrx)
       CALL costf1( p, ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &       dwdtLast,i_costf_init,d_costf_init)
       CALL get_dr( p, dp, &
            &       ulm_real-llm_real+1, start_lm-llm_real+1, stop_lm-llm_real+1, &
            &       n_r_max,n_cheb_max,dwdtLast,dpdtLast, &
            &       i_costf_init,d_costf_init,drx)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
#ifdef WITHOMP
    call omp_set_num_threads(omp_get_max_threads())
#endif
    !PERFOFF

    !PERFON('upWP_ex')
    !-- Calculate explicit time step part:
    IF ( ra /= 0.D0 ) THEN
       DO nR=n_r_cmb+1,n_r_icb-1
          DO lm1=lmStart_00,lmStop
             l1=lm2l(lm1)
             m1=lm2m(lm1)

             Dif(lm1) = hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR)*visc(nR) * &
                  & ( ddw(lm1,nR) &
                  &   +(2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) &
                  &   -( dLh(st_map%lm2(l1,m1))*or2(nR)&
                  &      +4.d0/3.d0*( dbeta(nR)&
                  &                   +dLvisc(nR)*beta(nR)&
                  &                   +(3.d0*dLvisc(nR)+beta(nR))*or1(nR)&
                  &                 )&
                  &    ) * w(lm1,nR) &
                  & )
             Pre(lm1) = -dp(lm1,nR)+beta(nR)*p(lm1,nR)
             Buo(lm1) = rho0(nR)*rgrav(nR)*s(lm1,nR)
             dwdtLast(lm1,nR)=dwdt(lm1,nR) - coex*(Pre(lm1)+Buo(lm1)+Dif(lm1))
             dpdtLast(lm1,nR)=&
                  & dpdt(lm1,nR) &
                  & - coex*( dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) &
                  &          + hdif_V(st_map%lm2(l1,m1))*visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR) &
                  &            * ( -workA(lm1,nR) &
                  &                + (beta(nR)-dLvisc(nR))*ddw(lm1,nR)&
                  &                + ( dLh(st_map%lm2(l1,m1))*or2(nR)&
                  &                    + dLvisc(nR)*beta(nR) &
                  &                    + dbeta(nR) &
                  &                    +2.d0*(dLvisc(nR)+beta(nR))*or1(nR)&
                  &                  ) * dw(lm1,nR) &
                  &                - dLh(st_map%lm2(l1,m1))*or2(nR) &
                  &                  * (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))&
                  &                  * w(lm1,nR)    &
                  &              )  &
                  &        )
             IF ( lRmsNext ) THEN
                workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                     &        ( w(lm1,nR)-workB(lm1,nR) )
                IF ( l_RMStest ) workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
             END IF
          END DO
          IF ( lRmsNext ) THEN
             CALL hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                  DifPol2hInt(nR,1),DifPolAs2hInt(nR,1),lo_map)
             !WRITE(*,"(A,I4,3ES22.14)") "upWP, work=",nR,SUM(workB(:,nR)),dtVPol2hInt(nR,nTh)
             CALL hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                  dtVPolLMr,dtVPol2hInt(nR,1),dtVPolAs2hInt(nR,1),lo_map)
             !WRITE(*,"(A,2I4,ES22.14)") "upWP: ",nR,nTh,dtVPol2hInt(nR,nTh)
          END IF
       END DO

    ELSE  ! no s-contribution !

       DO nR=n_r_cmb+1,n_r_icb-1
          DO lm1=lmStart_00,lmStop
             l1=lm2l(lm1)
             m1=lm2m(lm1)
             Dif(lm1)=     hdif_V(st_map%lm2(l1,m1))*dLh(st_map%lm2(l1,m1))*or2(nR)*visc(nR)* &
                  (                   ddw(lm1,nR)        + &
                  (2.d0*dLvisc(nR)-beta(nR)/3.d0)*dw(lm1,nR) - &
                  (dLh(st_map%lm2(l1,m1))*or2(nR)+4.d0/3.d0*(dbeta(nR)           + &
                  dLvisc(nR)*beta(nR)+ &
                  (3.d0*dLvisc(nR)+beta(nR))*or1(nR)))              * &
                  w(lm1,nR) )
             Pre(lm1)=-dp(lm1,nR)+beta(nR)*p(lm1,nR)
             dwdtLast(lm1,nR) = dwdt(lm1,nR) - coex*(Pre(lm1)+Dif(lm1))
             dpdtLast(lm1,nR)=        dpdt(lm1,nR) - coex*( &
                  dLh(st_map%lm2(l1,m1))*or2(nR)*p(lm1,nR) + &
                  hdif_V(st_map%lm2(l1,m1))*visc(nR)*dLh(st_map%lm2(l1,m1))*or2(nR) * ( &
                  -workA(lm1,nR)          + &
                  (beta(nR)-dLvisc(nR))*ddw(lm1,nR)          + &
                  ( dLh(st_map%lm2(l1,m1))*or2(nR)+dLvisc(nR)*beta(nR)    + &
                  dbeta(nR) &
                  +2.d0*(dLvisc(nR)+beta(nR))*or1(nR))    * &
                  dw(lm1,nR)          - &
                  dLh(st_map%lm2(l1,m1))*or2(nR)    * &
                  (2.d0*or1(nR)+2.d0/3.d0*beta(nR)+dLvisc(nR))* &
                  w(lm1,nR)          )  )
             IF ( lRmsNext ) THEN
                workB(lm1,nR)=O_dt*dLh(st_map%lm2(l1,m1))*or2(nR) * &
                     ( w(lm1,nR)-workB(lm1,nR) )
                IF ( l_RMStest ) &
                     workB(lm1,nR)=workB(lm1,nR)-Dif(lm1)
             END IF
          END DO
          IF ( lRmsNext ) THEN
             CALL hInt2Pol(Dif,llm,ulm,nR,lmStart_00,lmStop,DifPolLMr, &
                  DifPol2hInt(nR,1),DifPolAs2hInt(nR,1),lo_map)
             CALL hInt2Pol(workB(llm,nR),llm,ulm,nR,lmStart_00,lmStop, &
                  dtVPolLMr, dtVPol2hInt(nR,1),dtVpolAs2hInt(nR,1),lo_map)
          END IF
       END DO
       
    END IF
    !PERFOFF

    !DEALLOCATE(rhs1)

    !  Note: workA=dddw not needed beyond this point!

  end SUBROUTINE updateWP

  !------------------------------------------------------------------------------
END MODULE updateWP_mod
