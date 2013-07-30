!$Id$
!***********************************************************************
#include "perflib_preproc.cpp"

    SUBROUTINE LMLoop(w1,coex,time,dt,lMat,lRmsNext, &
       dVxBhLM,dVSrLM,dsdt,dwdt,dzdt,dpdt,dbdt,djdt, &
                lorentz_torque_ma,lorentz_torque_ic, &
                       b_nl_cmb,aj_nl_cmb,aj_nl_icb)
!***********************************************************************

!    !------------ This is release 2 level 10  --------------!
!    !------------ Created on 2/5/02  by JW. -----------

!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  This subroutine performs the actual time-stepping.               |
!  |                                                                   |
!  +-------------------------------------------------------------------+

    USE truncation
    USE radial_functions
    USE physical_parameters
    USE blocking
    USE horizontal_data
    USE logic
    USE matrices
    USE fields
    USE fieldsLast
    use omp_lib
    USE output_data
    USE parallel_mod
    USE timing, only: wallTime,subTime,writeTime

    IMPLICIT NONE

!-- Input of variables:
    REAL(kind=8),intent(IN) :: w1,coex
    REAL(kind=8),intent(IN) :: dt,time
    LOGICAL,intent(IN) :: lMat
    LOGICAL,intent(IN) :: lRmsNext

!--- Input from radialLoop:
!    These fields are provided in the R-distributed space!
    COMPLEX(kind=8),intent(IN) :: dVxBhLM(lm_maxMag,n_r_maxMag) ! for djdt in update_b
    COMPLEX(kind=8),intent(IN) :: dVSrLM(lm_max,n_r_max)        ! for dsdt in updata_s
    COMPLEX(kind=8),intent(IN) :: dsdt(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dwdt(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dzdt(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dpdt(lm_max,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dbdt(lm_maxMag,n_r_maxMag)
    COMPLEX(kind=8),intent(IN) :: djdt(lm_maxMag,n_r_maxMag)
    REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_ic
    COMPLEX(kind=8),intent(IN) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
    COMPLEX(kind=8),intent(IN)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
    COMPLEX(kind=8),intent(IN)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

!--- Local counter
    INTEGER :: nLMB
    INTEGER :: lmStart,lmStop
    INTEGER :: nTh
    INTEGER :: l
    INTEGER :: tStart(4),tStop(4),tPassed(4)

!--- Inner core rotation from last time step
    REAL(kind=8) :: omega_icLast
    SAVE omega_icLast
    COMPLEX(kind=8) :: workA(lm_max,n_r_max)


!-- end of declaration

    PERFON('LMloop')
    IF ( lVerbose ) CALL safeOpen(nLF,log_file)

    omega_icLast=omega_ic

!$OMP PARALLEL PRIVATE(l)
    IF ( lMat ) THEN ! update matricies:
    !---- The following logicals tell whether the respective inversion
    !     matricies have been updated. lMat=.TRUE. when a general
    !     update is necessary. These logicals are THREADPRIVATE and
    !     stored in the common block contained in c_mat.f:
        lZ10mat=.FALSE.
        DO l=0,l_max
            lSmat(l) =.FALSE.
            lZmat(l) =.FALSE.
            lWPmat(l)=.FALSE.
            lBmat(l) =.FALSE.
        END DO
    END IF
!$OMP END PARALLEL


!$OMP PARALLEL &
!$OMP PRIVATE(lmStart,lmStop,nTh,l,tStart,tStop,tPassed)

#ifdef WITHOMP
    !$OMP MASTER
    nThreadsLMmax=omp_get_num_threads()
    !$OMP END MASTER
#endif
    !$OMP DO SCHEDULE(STATIC,1) 
    DO nLMB=1,nLMBs ! Blocking of loop over all (l,m)
#ifdef WITHOMP
        nTh=OMP_GET_THREAD_NUM()+1
#else
        nTH=1
#endif
        !IF ( nTh > nThreadsLMmax ) nThreadsLMmax=nTh
        IF ( lVerbose ) WRITE(*,'(/," ! lm block no:",i3)') nLMB
        IF ( lVerbose ) WRITE(*,'(/," !   thread no:",i3)') nTh

        IF ( lVerbose ) CALL wallTime(tStart)

        lmStart=lmStartB(nLMB)
        lmStop =lmStopB(nLMB)

        IF ( l_heat ) & ! dp,workA usead as work arrays
            CALL updateS(s,ds,dVSrLM,dsdt,dsdtLast, &
                         dp,workA,w1,coex,dt,nLMB)
        IF ( l_conv ) THEN
            ! dp, dVSrLM, workA used as work arrays
            CALL updateZ(     z,dz,dzdt,dzdtLast,time, &
                           omega_ma,d_omega_ma_dtLast, &
                           omega_ic,d_omega_ic_dtLast, &
              lorentz_torque_ma,lorentz_torque_maLast, &
              lorentz_torque_ic,lorentz_torque_icLast, &
                                      dp,dVSrLM,workA, &
                         w1,coex,dt,nLMB,lRmsNext,nTh)
            ! dVSrLM, workA used as work arrays
            CALL updateWP(     w,dw,ddw,dwdt,dwdtLast, &
                                 p,dp,dpdt,dpdtLast,s, &
            dVSrLM,workA,w1,coex,dt,nLMB,lRmsNext,nTh)
        END IF
        IF ( l_mag )  & ! dwdt,dpdt used as work arrays
           CALL updateB(      b,db,ddb,aj,dj,ddj,dVxBhLM, &
                             dbdt,dbdtLast,djdt,djdtLast, &
                    b_ic,db_ic,ddb_ic,aj_ic,dj_ic,ddj_ic, &
                                 dbdt_icLast,djdt_icLast, &
               b_nl_cmb,aj_nl_cmb,aj_nl_icb,omega_icLast, &
             dwdt,dpdt,w1,coex,dt,time,nLMB,lRmsNext,nTh)

        IF ( lVerbose ) THEN
            CALL wallTime(tStop)
            CALL subTime(tStart,tStop,tPassed)
            WRITE(nLF,*) '! Thread no:',nTh
            WRITE(nLF,*) 'lmStart,lmStop:',lmStart,lmStop
            CALL writeTime(nLF,'! Time for thread:',tPassed)
        END IF


    END DO

    !$OMP END DO
!$OMP END PARALLEL    ! END OF SMP PARALLEL LOOP OVER LM blocks !


    lorentz_torque_maLast=lorentz_torque_ma
    lorentz_torque_icLast=lorentz_torque_ic

    IF ( lVerbose ) CALL safeClose(nLF)

    PERFOFF
    RETURN
    end SUBROUTINE LMLoop

!--------------------------------------------------------------------------------
