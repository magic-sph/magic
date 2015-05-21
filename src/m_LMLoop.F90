!$Id: s_LMLoop.F90 449 2013-02-27 12:39:23Z dannert $
!***********************************************************************
#include "perflib_preproc.cpp"

MODULE LMLoop_mod
  USE truncation
  USE radial_functions
  USE physical_parameters
  USE blocking
  USE horizontal_data
  USE logic
  USE matrices
  USE fields
  USE fieldsLast
  USE omp_lib
  USE output_data
  USE parallel_mod
  USE timing, ONLY: wallTime,subTime,writeTime
  USE LMLoop_data
  USE debugging, ONLY: debug_write
  USE communications,ONLY: get_global_sum, lo2r_redist_start,&
       & lo2r_s, lo2r_z, lo2r_p,&
       & lo2r_b, lo2r_aj, &
       & lo2r_w
  USE updateS_mod,ONLY: initialize_updateS,updateS,updateS_ala
  USE updateZ_mod,ONLY: initialize_updateZ,updateZ
  USE updateWP_mod,ONLY: initialize_updateWP,updateWP
  USE updateB_mod,ONLY: initialize_updateB,updateB
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LMLoop,initialize_LMLoop

CONTAINS
  SUBROUTINE initialize_LMLoop
    call initialize_updateS
    call initialize_updateZ
    call initialize_updateWP
    call initialize_updateB
  END SUBROUTINE initialize_LMLoop

  SUBROUTINE LMLoop(w1,coex,time,dt,lMat,lRmsNext, &
       &            dVxBhLM,dVSrLM,dsdt,dwdt,dzdt,dpdt,dbdt,djdt, &
       &            lorentz_torque_ma,lorentz_torque_ic, &
       &            b_nl_cmb,aj_nl_cmb,aj_nl_icb,n_time_step)
    !***********************************************************************
    
    !    !------------ This is release 2 level 10  --------------!
    !    !------------ Created on 2/5/02  by JW. -----------
    
    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine performs the actual time-stepping.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+
!-- Input of variables:
    REAL(kind=8),intent(IN) :: w1,coex
    REAL(kind=8),intent(IN) :: dt,time
    LOGICAL,intent(IN) :: lMat
    LOGICAL,intent(IN) :: lRmsNext

!--- Input from radialLoop:
!    These fields are provided in the R-distributed space!
    COMPLEX(kind=8),INTENT(IN) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag) ! for djdt in update_b
    COMPLEX(kind=8),intent(IN) :: dVSrLM(llm:ulm,n_r_max)        ! for dsdt in updata_s
    INTEGER, intent(IN) :: n_time_step

    !--- Input from radialLoop and then redistributed:
    COMPLEX(kind=8),intent(INOUT) :: dsdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dwdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dzdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dpdt(llm:ulm,n_r_max)
    COMPLEX(kind=8),intent(IN) :: dbdt(llmMag:ulmMag,n_r_maxMag)
    COMPLEX(kind=8),intent(INOUT) :: djdt(llmMag:ulmMag,n_r_maxMag)
    REAL(kind=8),intent(IN) :: lorentz_torque_ma,lorentz_torque_ic
    COMPLEX(kind=8),intent(IN) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
    COMPLEX(kind=8),intent(IN)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
    COMPLEX(kind=8),intent(IN)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

!--- Local counter
    INTEGER :: nLMB
    INTEGER :: lmStart,lmStop
    !INTEGER :: nTh,lm
    INTEGER :: l
    INTEGER :: tStart(4),tStop(4),tPassed(4)

    LOGICAL,PARAMETER :: DEBUG_OUTPUT=.false.
!--- Inner core rotation from last time step
    REAL(kind=8) :: omega_icLast
    SAVE omega_icLast
    complex(kind=8) :: sum_dwdt

!-- end of declaration

    PERFON('LMloop')
    !LIKWID_ON('LMloop')
    IF ( lVerbose ) CALL safeOpen(nLF,log_file)

    omega_icLast=omega_ic

    IF ( lMat ) THEN ! update matrices:
    !---- The following logicals tell whether the respective inversion
    !     matrices have been updated. lMat=.TRUE. when a general
    !     update is necessary. These logicals are THREADPRIVATE and
    !     stored in the MODULE matrices in m_mat.F90:
        lZ10mat=.FALSE.
        DO l=0,l_max
            lSmat(l) =.FALSE.
            lZmat(l) =.FALSE.
            lWPmat(l)=.FALSE.
            lBmat(l) =.FALSE.
        END DO
    END IF

    !nThreadsLMmax = 1
    nLMB=1+rank
    !nTh=1
    IF ( lVerbose ) THEN
       WRITE(*,'(/," ! lm block no:",i3)') nLMB
       CALL wallTime(tStart)
    END IF

    IF (DEBUG_OUTPUT) THEN
       lmStart=lmStartB(nLMB)
       lmStop =lmStopB(nLMB)
    END IF
    !CALL debug_write(dwdt(lmStart:lmStop,:),lmStop-lmStart+1,n_r_max,"dwdt",n_time_step*1000+nLMB*100,"E")

    IF ( l_heat ) THEN ! dp,workA usead as work arrays
       IF (DEBUG_OUTPUT) THEN
          WRITE(*,"(A,I2,8ES20.12)") "s_before: ",nLMB,GET_GLOBAL_SUM( s_LMloc(:,:) ),&
               & GET_GLOBAL_SUM( ds_LMloc(:,:) ), GET_GLOBAL_SUM( dsdt(:,:) ), &
               & GET_GLOBAL_SUM( dsdtLast_LMloc(:,:) )
       end if
       !CALL debug_write(dsdt,ulm-llm+1,n_r_max,"dsdt_LMloc",n_time_step*1000+nLMB*100,"E")
       PERFON('up_S')
       IF ( l_anelastic_liquid ) THEN
          CALL updateS_ala(s_LMloc,ds_LMloc,w_LMloc,dVSrLM,dsdt,    & 
               &       dsdtLast_LMloc,w1,coex,dt,nLMB)
       ELSE
          CALL updateS(s_LMloc,ds_LMloc,dVSrLM,dsdt,dsdtLast_LMloc, &
               &       w1,coex,dt,nLMB)

       END IF
       PERFOFF
       ! Here one could start the redistribution of s_LMloc,ds_LMloc etc. with a 
       ! nonblocking send
       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       !PERFON('rdstSst')
       CALL lo2r_redist_start(lo2r_s,s_LMloc_container,s_Rloc_container)
       !PERFOFF

       IF (DEBUG_OUTPUT) THEN
          WRITE(*,"(A,I2,4ES20.12)") "s_after : ",nLMB,GET_GLOBAL_SUM( s_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( ds_LMloc(:,:) )
          WRITE(*,"(A,I2,8ES22.14)") "s_after(bnd_r): ",nLMB,&
               &GET_GLOBAL_SUM( s_LMloc(:,n_r_icb) ),GET_GLOBAL_SUM( s_LMloc(:,n_r_cmb) ),&
               &GET_GLOBAL_SUM( ds_LMloc(:,n_r_icb) ),GET_GLOBAL_SUM( ds_LMloc(:,n_r_cmb) )
       END IF
    END IF
    IF ( l_conv ) THEN
       if (DEBUG_OUTPUT) then
          WRITE(*,"(A,I2,6ES20.12)") "z_before: ",nLMB,GET_GLOBAL_SUM( z_LMloc(:,:) ),&
               & GET_GLOBAL_SUM( dz_LMloc(:,:) ),GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
       end if
       PERFON('up_Z')
       ! dp, dVSrLM, workA used as work arrays
       !CALL updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_LMloc, time, &
       CALL updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_lo, time, &
            &        omega_ma,d_omega_ma_dtLast, &
            &        omega_ic,d_omega_ic_dtLast, &
            &        lorentz_torque_ma,lorentz_torque_maLast, &
            &        lorentz_torque_ic,lorentz_torque_icLast, &
            &        w1,coex,dt,lRmsNext)
       PERFOFF

       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       !PERFON('rdstZst')
       CALL lo2r_redist_start(lo2r_z,z_LMloc_container,z_Rloc_container)
       !PERFOFF

       if (DEBUG_OUTPUT) then
          !DO lm=lmStart,lmStop
          !   WRITE(*,"(A,I4,6ES20.12)") "z_after : ",lm,SUM( z_LMloc(lm,:) ),&
          !        & SUM( dz_LMloc(lm,:) ),SUM( dzdtLast_lo(lm,:) )
          !END DO
          WRITE(*,"(A,I2,6ES20.12)") "z_after: ",nLMB,GET_GLOBAL_SUM( z_LMloc(:,:) ),&
               & GET_GLOBAL_SUM( dz_LMloc(:,:) ),GET_GLOBAL_SUM( dzdtLast_lo(:,:) )
       end if
       ! dVSrLM, workA used as work arrays
       !CALL debug_write(dwdt(:,:),lmStop-lmStart+1,n_r_max,"dwdt",n_time_step*1000+nLMB*100,"E")
       IF (DEBUG_OUTPUT) THEN
          sum_dwdt=GET_GLOBAL_SUM( dwdt(:,:) )
          WRITE(*,"(A,I2,8ES22.14,4(I3,F19.16))") "wp_before: ",nLMB,&
               &GET_GLOBAL_SUM( w_LMloc(:,:) ),GET_GLOBAL_SUM( p_LMloc(:,:) ),&
               & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ), GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),&
               & EXPONENT(REAL(sum_dwdt)),FRACTION(REAL(sum_dwdt)),&
               & EXPONENT(AIMAG(sum_dwdt)),FRACTION(AIMAG(sum_dwdt))
       end if
       PERFON('up_WP')
       CALL updateWP( w_LMloc, dw_LMloc, ddw_LMloc, dwdt, dwdtLast_LMloc, &
            &         p_LMloc, dp_LMloc, dpdt, dpdtLast_LMloc, s_LMloc, &
            &         w1,coex,dt,nLMB,lRmsNext)
       PERFOFF

       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       !PERFON('rdstWPst')
       CALL lo2r_redist_start(lo2r_w,w_LMloc_container,w_Rloc_container)

       CALL lo2r_redist_start(lo2r_p,p_LMloc_container,p_Rloc_container)
       !PERFOFF

       IF (DEBUG_OUTPUT) THEN
          WRITE(*,"(A,I2,12ES22.14)") "wp_after: ",nLMB,&
               &GET_GLOBAL_SUM( w_LMloc(:,:) ),GET_GLOBAL_SUM( p_LMloc(:,:) ),&
               & GET_GLOBAL_SUM( dwdtLast_LMloc(:,:) ), GET_GLOBAL_SUM( dpdtLast_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( dw_LMloc(:,:) )
          WRITE(*,"(A,I2,4ES22.14)") "wp_after(bnd_r): ",nLMB,&
               &GET_GLOBAL_SUM( w_LMloc(:,n_r_icb) ),GET_GLOBAL_SUM( w_LMloc(:,n_r_cmb) )
       end if
    END IF
    IF ( l_mag ) THEN ! dwdt,dpdt used as work arrays
       IF (DEBUG_OUTPUT) THEN
          write(*,"(A,I2,14ES20.12)") "b_before: ",nLMB,&
               &GET_GLOBAL_SUM(  b_LMloc(:,:) ),& 
               &GET_GLOBAL_SUM( aj_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( dVxBhLM(:,:) )
       END IF
       !LIKWID_ON('up_B')
       PERFON('up_B')
       CALL updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc,dVxBhLM, &
            &        dbdt,dbdtLast_LMloc,djdt,djdtLast_LMloc, &
            &        b_ic_LMloc,db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,dj_ic_LMloc,ddj_ic_LMloc, &
            &        dbdt_icLast_LMloc, djdt_icLast_LMloc, &
            &        b_nl_cmb,aj_nl_cmb,aj_nl_icb,omega_icLast, &
            &        w1,coex,dt,time,nLMB,lRmsNext)
       PERFOFF
       !LIKWID_OFF('up_B')
       CALL lo2r_redist_start(lo2r_b,  b_LMloc_container,b_Rloc_container)
       
       CALL lo2r_redist_start(lo2r_aj, aj_LMloc_container,aj_Rloc_container)

       if (DEBUG_OUTPUT) then
          WRITE(*,"(A,I2,8ES20.12)") "b_after: ",nLMB,&
               &GET_GLOBAL_SUM(  b_LMloc(:,:) ),& 
               &GET_GLOBAL_SUM( aj_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( dbdtLast_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( djdtLast_LMloc(:,:) )
          WRITE(*,"(A,I2,8ES20.12)") "b_ic_after: ",nLMB,&
               &GET_GLOBAL_SUM( b_ic_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( aj_ic_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( dbdt_icLast_LMloc(:,:) ),&
               &GET_GLOBAL_SUM( djdt_icLast_LMloc(:,:) )
       end if
    END IF

    IF ( lVerbose ) THEN
       CALL wallTime(tStop)
       CALL subTime(tStart,tStop,tPassed)
       !WRITE(nLF,*) '! Thread no:',nTh
       WRITE(nLF,*) 'lmStart,lmStop:',lmStartB(nLMB),lmStopB(nLMB)
       CALL writeTime(nLF,'! Time for thread:',tPassed)
    END IF


    lorentz_torque_maLast=lorentz_torque_ma
    lorentz_torque_icLast=lorentz_torque_ic

    IF ( lVerbose ) CALL safeClose(nLF)

    !LIKWID_OFF('LMloop')
    PERFOFF
    RETURN
  END SUBROUTINE LMLoop
  !--------------------------------------------------------------------------------

END MODULE LMLoop_mod
