!$Id: s_step_time.F90 436 2013-02-20 11:17:48Z gastine $
#include "perflib_preproc.cpp"
MODULE step_time_mod
    USE truncation
    USE radial_functions
    USE physical_parameters
    USE num_param
    USE init_fields
    USE blocking
    USE horizontal_data
    USE logic
    USE movie_data
    use radialLoop
    USE LMLoop_data, only: llm,ulm,llmMag,ulmMag,lm_per_rank,lm_on_last_rank
    USE LMLoop_mod,ONLY: LMLoop
    USE output_data
    USE output_mod, only: output
    USE const
    USE timing
    use parallel_mod
    USE fields,ONLY: w_LMloc,dw_LMloc,ddw_LMloc,s_LMloc,ds_LMloc,p_LMloc,b_LMloc,aj_LMloc,&
         &s_Rloc,ds_Rloc,z_Rloc,dz_Rloc,w_Rloc,dw_Rloc,ddw_Rloc,p_Rloc,dp_Rloc,b_Rloc,db_Rloc,ddb_Rloc,&
         &aj_Rloc,dj_Rloc,&
         &dp,dp_LMloc,db_LMloc,ddb_LMloc,dj_LMloc,ddj_LMloc,b_ic_LMloc,db_ic_LMloc,ddb_ic_LMloc,aj_ic_LMloc,dj_ic_LMloc,&
         &ddj_ic_LMloc,s,ds,z,dz,dw,ddw,p,b,db,ddb,aj,dj,ddj,b_ic,aj_ic,db_ic,dj_ic,ddb_ic,ddj_ic,omega_ic,omega_ma,w,&
         & z_LMloc,dz_LMloc
    USE fieldsLast
    USE charmanip, only: capitalize,dble2str
    USE usefull, only: l_correct_step
    USE communications,ONLY: get_global_sum,& !lm2r_redist,lo2r_redist,r2lm_redist,&
         & r2lo_redist,&
         & lm2r_type,&
         & lo2r_redist_start,lo2r_redist_wait,&
         &lo2r_s,lo2r_ds,lo2r_z,lo2r_dz, lo2r_w,lo2r_dw,lo2r_ddw,lo2r_p,lo2r_dp,&
         & lo2r_b, lo2r_db, lo2r_ddb, lo2r_aj, lo2r_dj
    IMPLICIT NONE 

  private

  public :: initialize_step_time,step_time
contains
  SUBROUTINE initialize_step_time
  END SUBROUTINE initialize_step_time

  !***********************************************************************
  SUBROUTINE step_time(time,dt,dtNew,n_time_step)
    !***********************************************************************

    !    !------------ This is release 2 level 10  --------------!
    !    !------------ Created on 2/5/02  by JW. -----------

    !  +-------------+----------------+------------------------------------+
    !  |                                                                   |
    !  |  This subroutine performs the actual time-stepping.               |
    !  |                                                                   |
    !  +-------------------------------------------------------------------+

    !-- Input from initialization:
    !   time and n_time_step updated and returned to magic.f
    REAL(kind=8) :: time
    REAL(kind=8) :: dt,dtNew
    INTEGER :: n_time_step

    !--- Local variables:

    !--- Logicals controlling output/calculation:
    LOGICAL :: l_graph          ! 
    LOGICAL :: l_spectrum
    LOGICAL :: l_cour           ! Check Courant criteria
    LOGICAL :: lCourChecking    ! Ongoing Courant criteria check
    LOGICAL :: l_store          ! Store output in restart file 
    !LOGICAL :: l_storeSched     ! Storage of rst-file scheduled
    LOGICAL :: l_new_rst_file   ! Use new rst file
    LOGICAL :: l_log            ! Log output
    LOGICAL :: l_stop_time      ! Stop time stepping
    LOGICAL :: l_frame          ! Movie frame output
    LOGICAL :: lTOframe         ! TO movie frame output
    LOGICAL :: l_cmb            ! Store set of b at CMB
    LOGICAL :: lHelCalc         ! Calculate helicity for output
    LOGICAL :: lTOCalc          ! Calculate TO stuff
    LOGICAL :: lTONext,lTONext2 ! TO stuff for next steps
    LOGICAL :: lTOframeNext,lTOframeNext2
    LOGICAL :: lTOZhelp,lTOZwrite
    LOGICAL :: l_logNext,l_logNext2
    LOGICAL :: l_Bpot,l_Vpot,l_Tpot
    LOGICAL :: lRmsCalc,lRmsNext
    LOGICAL :: lMat             ! update matricies

    !--- Counter:
    INTEGER :: n                ! Counter
    INTEGER :: n_graph          ! No. of graphic file
    INTEGER :: n_frame          ! No. of movie frames
    INTEGER :: n_cmb_sets       ! No. of stored sets of b at CMB
    !INTEGER :: nR,lm

    !--- Stuff needed to construct output files:
    CHARACTER(len=20) :: string
    !INTEGER :: length

    CHARACTER(len=255) :: message
    CHARACTER(len=76) :: SIG    

    !--- Courant criteria/diagnosis:
    REAL(kind=8) :: dtr,dth
    !REAL(kind=8) :: dtrkc(n_r_max),dthkc(n_r_max)   ! Saves values for time step
    REAL(kind=8),dimension(nRstart:nRstop) :: dtrkc_Rloc,dthkc_Rloc   ! Saves values for time step

    !--- Explicit part of time stepping, calculated in s_radialLoopG.f and
    !    passed to s_LMLoop.f where the time step is preformed.
    !    Note that the respective arrays for the changes in inner-core
    !    magnetic field are calculated in s_updateB.f and are only
    !    needed there.
    !--- dVSrLM and dVxBhLM are help arrays for calculating dsdt and djdt (see s_updateS.f,s_updateB.f):
    COMPLEX(kind=8),DIMENSION(lm_max,nRstart:nRstop) :: dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dsdt_Rloc,&
         & dVSrLM_Rloc
    COMPLEX(kind=8),DIMENSION(lm_maxMag,nRstartMag:nRstopMag) :: dbdt_Rloc,djdt_Rloc,dVxBhLM_Rloc
    TARGET :: dbdt_Rloc

    ! The same arrays, but now the LM local part
    COMPLEX(kind=8),DIMENSION(llm:ulm,n_r_max) :: dwdt_LMloc,dzdt_LMloc,dpdt_LMloc,dsdt_LMloc,&
         & dVSrLM_LMloc
    COMPLEX(kind=8),DIMENSION(llm:ulm,n_r_max) :: dzdt_lo
    COMPLEX(kind=8),DIMENSION(llmMag:ulmMag,n_r_maxMag) :: dbdt_LMloc,djdt_LMloc,dVxBhLM_LMloc

    !--- Lorentz torques:
    REAL(kind=8) :: lorentz_torque_ma,lorentz_torque_ic

    !--- Help arrays for calculating axsymmetric helicity in s_outMisc.f:
    REAL(kind=8),DIMENSION(l_max+1,nRstart:nRstop) :: HelLMr_Rloc,Hel2LMr_Rloc,&
         & HelnaLMr_Rloc,Helna2LMr_Rloc
    REAL(kind=8),DIMENSION(l_max+1,nRstart:nRstop) :: uhLMr_Rloc,duhLMr_Rloc

    !REAL(kind=8),DIMENSION(l_max+1,n_r_max) :: HelLMr,Hel2LMr,HelnaLMr,Helna2LMr
    !REAL(kind=8),DIMENSION(l_max+1,n_r_max) :: uhLMr,duhLMr

    !--- Nonlinear magnetic boundary conditions needed in s_updateB.f :
    COMPLEX(kind=8) :: br_vt_lm_cmb(lmP_max)    ! product br*vt at CMB
    COMPLEX(kind=8) :: br_vp_lm_cmb(lmP_max)    ! product br*vp at CMB
    COMPLEX(kind=8) :: br_vt_lm_icb(lmP_max)    ! product br*vt at ICB
    COMPLEX(kind=8) :: br_vp_lm_icb(lmP_max)    ! product br*vp at ICB
    COMPLEX(kind=8) :: b_nl_cmb(lm_max)         ! nonlinear bc for b at CMB
    COMPLEX(kind=8) :: aj_nl_cmb(lm_max)        ! nonlinear bc for aj at CMB
    COMPLEX(kind=8) :: aj_nl_icb(lm_max)        ! nonlinear bc for dr aj at ICB


    !--- Various stuff for time control:
    REAL(kind=8) :: timeLast
    REAL(kind=8) :: dtLast
    REAL(kind=8) :: w1,w2,w2New,coex
    INTEGER :: n_time_steps_go,n_time_cour
    LOGICAL :: l_new_dt        ! causes call of matbuild !
    LOGICAL :: l_new_dtNext    ! causes call of matbuild !
    LOGICAL :: l_new_dtHit     ! causes call of matbuild !
    INTEGER :: n_dt_changed    
    INTEGER :: n_dt_check        
    REAL(kind=8) :: timeScaled        ! Scaled time for output.
    INTEGER :: nPercent         ! percentage of finished time stepping
    REAL(kind=8) :: tenth_n_time_steps

    !-- Interupt procedure:
    INTEGER :: n_stop_signal     ! =1 causes run to stop
    INTEGER :: n_graph_signal    ! =1 causes output of graphic file
    INTEGER :: n_spec_signal     ! =1 causes output of a spec file
    INTEGER :: n_rst_signal      ! =1 causes output of rst file

    !--- Timing
    INTEGER :: runTimePassed(4)
    INTEGER :: runTimeRstart(4),runTimeRstop(4)
    INTEGER :: runTimeTstart(4),runTimeTstop(4)
    INTEGER :: runTimeR(4),runTimeLM(4),runTimeT(4)
    INTEGER :: runTimeTL(4),runTimeTM(4)
    INTEGER :: nTimeT,nTimeTL,nTimeTM,nTimeR,nTimeLM

    logical,parameter :: DEBUG_OUTPUT=.false.
    INTEGER :: lmStart,lmStop,lmStart_00
    INTEGER :: nR_i1,nR_i2
    INTEGER :: lm,l,m
#ifdef WITH_MPI
    ! MPI related variables
    INTEGER :: i,j,info,irank, send_pe, recv_pe
    INTEGER :: sendcount
    !integer :: lmStart_on_rank,lmStop_on_rank,nR
    INTEGER,ALLOCATABLE,DIMENSION(:) :: recvcounts,displs
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: error_string
    INTEGER :: length_of_error,nR,nLMB
#endif
    COMPLEX(KIND=8),POINTER,DIMENSION(:) :: ptr_dbdt_CMB
    !REAL(kind=8) :: start_time, end_time
    !-- end of declaration


    IF ( lVerbose ) WRITE(*,'(/,'' ! STARTING STEP_TIME !'')')
    
#ifdef WITH_MPI
    ! allocate the buffers for MPI gathering
    ALLOCATE(recvcounts(0:n_procs-1),displs(0:n_procs-1))
    CALL MPI_INFO_CREATE(info,ierr)
#endif

    l_log       =.FALSE.
    l_stop_time =.FALSE. 
    l_new_dt    =.TRUE.   ! Invokes calculation of t-step matricies
    l_new_dtNext=.TRUE.  
    w2New       =-0.5D0*dtNew/dt
    n_dt_changed=0        ! No of time steps since dt changed
    n_dt_check  =4        ! No of courant checks after dt changed

    tenth_n_time_steps=DBLE(n_time_steps)/10.D0
    nPercent=9

    !---- Set Lorentz torques to zero:
    lorentz_torque_ic=0.D0
    lorentz_torque_ma=0.D0

    !---- Counter for output files/sets:
    n_graph         =0    ! No. of graphic file
    n_frame         =0    ! No. of movie frames
    n_cmb_sets      =0    ! No. of store dt_b sets at CMB

    !---- Prepare signalling via file signal
    n_stop_signal =0     ! Stop signal returned to calling program
    n_graph_signal=0     ! Graph signal returned to calling program
    n_spec_signal=0      ! Spec signal
    n_rst_signal=0
    IF (rank.EQ.0) THEN
       message='signal'//'.'//tag
       OPEN(19,FILE=TRIM(message),STATUS='unknown')
       WRITE(19,'(A3)') 'NOT'
       CLOSE(19)
    END IF

    !-- STARTING THE TIME STEPPING LOOP:
    IF (rank.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '! Starting time integration!'
    END IF
    nTimeT =0
    nTimeTL=0
    nTimeTM=0
    nTimeR =0
    nTimeLM=0
    DO n=1,4
       runTime(n)  =0
       runTimeT(n) =0
       runTimeTM(n)=0
       runTimeTL(n)=0
       runTimeR(n) =0
       runTimeLM(n)=0
    END DO

!!!!! Time loop starts !!!!!!
    n_time_cour=-2 ! Causes a Courant check after first update
    IF ( n_time_steps.EQ.1 ) THEN
       n_time_steps_go=1 ! Output only, for example G-file/movie etc.
    ELSE IF ( n_time_steps.EQ.2 ) THEN
       n_time_steps_go=2 ! 
    ELSE
       n_time_steps_go=n_time_steps+1  ! Last time step for output only !
    END IF

    CALL mpi_barrier(MPI_COMM_WORLD,ierr)

    PERFON('tloop')
    DO n_time_step=1,n_time_steps_go 
       n_time_cour=n_time_cour+1
       
       IF ( lVerbose .or. DEBUG_OUTPUT ) THEN 
          WRITE(*,*)
          WRITE(*,*) '! Starting time step ',n_time_step
       END IF

       ! =================================== BARRIER ======================
       PERFON('barr_0')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       ! ==================================================================

       PERFON('signals')
       IF (rank.EQ.0) THEN
          !----- Signalling via file signal:
          message='signal'//'.'//tag
          OPEN(19,FILE=TRIM(message),STATUS='old')
          READ(19,*) SIG
          CLOSE(19)
          IF ( LEN(TRIM(SIG)).GT.0 ) THEN ! Non blank string ?
             CALL capitalize(SIG)
             IF ( INDEX(SIG,'END')/=0 ) n_stop_signal=1
             n_graph_signal=0
             IF ( INDEX(SIG,'GRA')/=0 ) THEN 
                n_graph_signal=1
                OPEN(19,FILE=TRIM(message),STATUS='unknown')
                WRITE(19,'(A3)') 'NOT'
                CLOSE(19)
             END IF
             n_rst_signal=0
             IF ( INDEX(SIG,'RST')/=0 ) THEN
                n_rst_signal=1
                OPEN(19,FILE=TRIM(message),STATUS='unknown')
                WRITE(19,'(A3)') 'NOT'
                CLOSE(19)
             END IF
             n_spec_signal=0
             IF ( INDEX(SIG,'SPE')/=0 ) THEN
                n_spec_signal=1
                OPEN(19,FILE=TRIM(message),STATUS='unknown')
                WRITE(19,'(A3)') 'NOT'
                CLOSE(19)
             END IF
          END IF
       END IF
       ! Broadcast the results from the signal file to all processes
       ! =======> THIS IS A GLOBAL SYNCHRONIZATION POINT <==========
       CALL MPI_Bcast(n_stop_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Bcast(n_graph_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Bcast(n_spec_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       CALL MPI_Bcast(n_rst_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       PERFOFF
       CALL wallTime(runTimeTstart)

       !PERFON('chk_stop')
       !--- Various reasons to stop the time integration:
       IF ( l_runTimeLimit ) THEN
          call MPI_Allreduce(MPI_IN_PLACE,runTime,4,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
          IF ( lTimeLimit(runTime,runTimeLimit) ) THEN
             WRITE(message,'("! Run time limit exeeded !")')
             CALL logWrite(message)
             l_stop_time=.TRUE.
          END IF
       END IF
       IF ( n_stop_signal.GT.0 .OR.                                 &
            &          n_time_step.EQ.n_time_steps_go )                        &
            &        l_stop_time=.TRUE.   ! last time step !

       !--- Another reasons to stop the time integration:
       IF ( time.GE.tEND .AND. tEND.NE.0.D0 ) l_stop_time=.true.
       !PERFOFF
       !PERFON('logics')
       !-- Checking logic for output: 
       l_graph= l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
            &                  n_graph_step,n_graphs,n_t_graph,t_graph,0) .OR. &
            &             n_graph_signal.EQ.1
       !l_graph=.FALSE.
       n_graph_signal=0   ! reset interrupt signal !
       l_spectrum=                                                  &
            &         l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
            &           n_spec_step,n_specs,n_t_spec,t_spec,0) .OR.            &
            &           n_spec_signal.EQ.1
       l_frame= l_movie .and. (                                     &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &        n_movie_step,n_movie_frames,n_t_movie,t_movie,0) .OR.     &
            &              n_time_steps_go.EQ.1 )
           IF ( l_mag .OR. l_mag_LF ) THEN
              l_dtB=( l_frame .and. l_dtBmovie ) .or.                   &
     &              ( l_log .and. l_DTrMagSpec ) 
           END IF
       l_HT  = l_frame .and. l_HTmovie

       lTOframe=l_TOmovie .AND.                                     &
            &     l_correct_step(n_time_step-1,time,timeLast,n_time_steps,     &
            &     n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)

       l_Bpot=l_storeBpot .AND. (                                   &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_Bpot_step,n_Bpots,n_t_Bpot,t_Bpot,0).OR. &
            &            n_time_steps.EQ.1 )
       l_Vpot=l_storeVpot .AND. (                                   &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_Vpot_step,n_Vpots,n_t_Vpot,t_Vpot,0).OR. &
            &            n_time_steps.EQ.1 )
       l_Tpot=l_storeTpot .AND. (                                   &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_Tpot_step,n_Tpots,n_t_Tpot,t_Tpot,0).OR. &
            &            n_time_steps.EQ.1 )

       l_cour=.TRUE.

       l_new_rst_file=                                              &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_rst_step,n_rsts,n_t_rst,t_rst,0) .OR.    &
            &        n_rst_signal.EQ.1
       n_rst_signal=0
       l_store= l_new_rst_file .OR.                                 &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       0,n_stores,0,t_rst,0)

       l_log= l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_log_step,n_logs,n_t_log,t_log,0)
       l_cmb= l_cmb_field .and.                                     &
            &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
            &                       n_cmb_step,n_cmbs,n_t_cmb,t_cmb,0)
       l_logNext=.FALSE.
       IF ( n_time_step+1.LE.n_time_steps+1 )                       &
            &        l_logNext=                                                &
            &        l_correct_step(n_time_step,time+dt,timeLast,              &
            &              n_time_steps,n_log_step,n_logs,n_t_log,t_log,0)
       l_logNext2=.FALSE.
       IF ( n_time_step+2.LE.n_time_steps+1 )                       &
            &        l_logNext2=                                               &
            &        l_correct_step(n_time_step+1,time+2*dt,timeLast,          &
            &         n_time_steps,n_log_step,n_logs,n_t_log,t_log,0)
       lTOCalc= n_time_step.GT.2 .AND. l_TO .AND.                   &
            &          l_correct_step(n_time_step-1,time,timeLast,             &
            &          n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
       lTOnext     =.FALSE.
       lTOframeNext=.FALSE.
       IF ( n_time_step+1.LE.n_time_steps+1 ) THEN
          lTONext= l_TO .AND.                                       &
               &           l_correct_step(n_time_step,time+dt,timeLast,           &
               &            n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
          lTOframeNext= l_TOmovie .AND.                             &
               &           l_correct_step(n_time_step,time+dt,timeLast,           &
               &          n_time_steps,n_TOmovie_step,n_TOmovie_frames,           &
               &                               n_t_TOmovie,t_TOmovie,0)
       END IF
       lTONext      =lTOnext.OR.lTOframeNext 
       lTONext2     =.FALSE.
       lTOframeNext2=.FALSE.
       IF ( n_time_step+2.LE.n_time_steps+1 ) THEN
          lTONext2= l_TO .AND.                                      &
               &           l_correct_step(n_time_step+1,time+2*dt,timeLast,       &
               &                                    n_time_steps,n_TO_step,       &
               &                                       n_TOs,n_t_TO,t_TO,0)
          lTOframeNext2= l_TOmovie .AND.                            &
               &           l_correct_step(n_time_step+1,time+2*dt,timeLast,       &
               &                               n_time_steps,n_TOmovie_step,       &
               &                  n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)
       END IF
       lTONext2=lTOnext2.OR.lTOframeNext2 
       lTOZhelp= n_time_step.GT.2 .AND. l_TO .AND.                  &
            &                l_correct_step(n_time_step-1,time,timeLast,       &
            &            n_time_steps,n_TOZ_step,n_TOZs,n_t_TOZ,t_TOZ,0)
       IF ( lTOZhelp ) lTOZwrite=.TRUE.

       lRmsCalc=l_RMS.AND.l_log.AND.(n_time_step.GT.1)
       IF ( l_mag .OR. l_mag_LF ) l_dtB   =l_dtB.OR.lRmsCalc
       lRmsNext=l_RMS.AND.l_logNext ! Used for storing in update routines !

       IF ( n_time_step.EQ.1 ) l_log=.TRUE.

       IF ( l_stop_time ) THEN                  ! Programm stopped by kill -30
          l_new_rst_file=.TRUE.                 ! Write rst-file and some
          IF ( n_stores.GT.0 ) l_store=.TRUE.   ! diagnostics before dying ! 
          l_log=.TRUE.
          lRmsNext=.FALSE.
       END IF
       lHelCalc=l_hel.AND.l_log

       IF ( l_graph ) THEN  ! write graphic output !
          n_graph=n_graph+1     ! increase counter for graphic file
#ifdef WITH_MPI
          IF ( l_graph_time ) THEN 
             CALL dble2str(time,string)
             IF ( ngform.NE.0 ) THEN
                graph_file='g_t='//TRIM(string)//'.'//tag_wo_rank
             ELSE
                graph_file='G_t='//trim(string)//'.'//tag_wo_rank
             END IF
          ELSE
             write(string, *) n_graph
             IF ( ngform.NE.0 ) THEN
                graph_file='g_'//trim(adjustl(string))//'.'//tag_wo_rank
             ELSE
                graph_file='G_'//trim(adjustl(string))//'.'//tag_wo_rank
             END IF
          END IF
#else
          IF ( l_graph_time ) THEN 
             CALL dble2str(time,string)
             IF ( ngform.NE.0 ) THEN
                graph_file='g_t='//trim(string)//'.'//tag
             ELSE
                graph_file='G_t='//trim(string)//'.'//tag
             END IF
          ELSE
             write(string, *) n_graph
             IF ( ngform.NE.0 ) THEN
                graph_file='g_'//trim(adjustl(string))//'.'//tag
             ELSE
                graph_file='G_'//trim(adjustl(string))//'.'//tag
             END IF
          END IF
#endif
          IF (rank.EQ.0) THEN
             WRITE(*,'(1p,/,A,/,A,D20.10,/,A,i15,/,A,A)')&
                  &" ! Storing graphic file:",&
                  &"             at time=",timeScaled,&
                  &"            step no.=",n_time_step,&
                  &"           into file=",graph_file
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(1p,/,A,/,A,d20.10,/,A,i15,/,A,A)') &
                  &" ! Storing graphic file:",&
                  &"             at time=",timeScaled,&
                  &"            step no.=",n_time_step,&
                  &"           into file=",graph_file
             CALL safeClose(nLF)
          END IF
#ifdef WITH_MPI
          CALL MPI_File_open(MPI_COMM_WORLD,graph_file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,graph_mpi_fh,ierr)
          !CALL MPI_ERROR_STRING(ierr,error_string,length_of_error,ierr)
          !PRINT*,"MPI_FILE_OPEN returned: ",TRIM(error_string)
#else
          IF ( ngform.NE.0 ) THEN
             OPEN(n_graph_file,FILE=graph_file,STATUS='NEW',FORM='FORMATTED')
          ELSE
             OPEN(n_graph_file,FILE=graph_file,STATUS='NEW',FORM='UNFORMATTED')
          END IF
#endif

       END IF

#ifdef WITH_MPI

       IF (DEBUG_OUTPUT) THEN
          DO nLMB=1+rank*nLMBs_per_rank,MIN((rank+1)*nLMBs_per_rank,nLMBs)
             lmStart=lmStartB(nLMB)
             lmStop=lmStopB(nLMB)
             lmStart_00  =MAX(2,lmStart)

             WRITE(*,"(A,I3,6ES20.12)") "start w: ",nLMB,&
                  & GET_GLOBAL_SUM( w_LMloc(lmStart:lmStop,:) ),&
                  & GET_GLOBAL_SUM( dw_LMloc(lmStart:lmStop,:) ),&
                  & GET_GLOBAL_SUM( ddw_LMloc(lmStart:lmStop,:) )
             WRITE(*,"(A,I3,4ES20.12)") "start z: ",nLMB,&
                  &GET_GLOBAL_SUM( z_LMloc(lmStart:lmStop,:) ),&
                  &GET_GLOBAL_SUM( dz_LMloc(lmStart:lmStop,:) )
             WRITE(*,"(A,I3,4ES20.12)") "start s: ",nLMB,&
                  &GET_GLOBAL_SUM( s_LMloc(lmStart:lmStop,:) ),&
                  &GET_GLOBAL_SUM( ds_LMloc(lmStart:lmStop,:) )
             WRITE(*,"(A,I3,4ES20.12)") "start p: ",nLMB,&
                  &GET_GLOBAL_SUM( p_LMloc(lmStart:lmStop,:) ),&
                  &GET_GLOBAL_SUM( dp_LMloc(lmStart_00:lmStop,:) )
             IF (l_mag) THEN
                WRITE(*,"(A,I3,8ES20.12)") "start b: ",nLMB,&
                     &GET_GLOBAL_SUM( b_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( db_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( ddb_LMloc(lmStart:lmStop,:) )
                WRITE(*,"(A,I3,8ES20.12)") "start aj: ",nLMB,&
                     &GET_GLOBAL_SUM( aj_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( dj_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( ddj_LMloc(lmStart:lmStop,:) )
                WRITE(*,"(A,I3,8ES20.12)") "start b_ic: ",nLMB,&
                     &GET_GLOBAL_SUM( b_ic_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( db_ic_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( ddb_ic_LMloc(lmStart:lmStop,:) )
                WRITE(*,"(A,I3,8ES20.12)") "start aj_ic: ",nLMB,&
                     &GET_GLOBAL_SUM( aj_ic_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( dj_ic_LMloc(lmStart:lmStop,:) ),&
                     &GET_GLOBAL_SUM( ddj_ic_LMloc(lmStart:lmStop,:) )
             END IF
          END DO
       END IF
       
       ! =================================== BARRIER ======================
       PERFON('barr_1')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       ! ==================================================================
       ! Here now comes the block where the LM distributed fields
       ! are redistributed to Rloc distribution which is needed for the radialLoop.
       ! s,ds
       ! z,dz
       ! w,dw,ddw,p,dp
       ! b,db,ddb,aj,dj,ddj
       ! b_ic,db_ic, ddb_ic,aj_ic,dj_ic,ddj_ic

       ! Waiting for the completion before we continue to the radialLoop
       IF (l_heat) THEN
          CALL lo2r_redist_wait(lo2r_s)
          CALL lo2r_redist_wait(lo2r_ds)
       END IF
       IF (l_conv) THEN
          call lo2r_redist_wait(lo2r_z)
          call lo2r_redist_wait(lo2r_dz)
          CALL lo2r_redist_wait(lo2r_w)
          CALL lo2r_redist_wait(lo2r_dw)
          CALL lo2r_redist_wait(lo2r_ddw)
          CALL lo2r_redist_wait(lo2r_p)
          CALL lo2r_redist_wait(lo2r_dp)
       END IF

       if (l_mag) then
          CALL lo2r_redist_wait(lo2r_b)
          CALL lo2r_redist_wait(lo2r_db)
          CALL lo2r_redist_wait(lo2r_ddb)

          CALL lo2r_redist_wait(lo2r_aj)
          CALL lo2r_redist_wait(lo2r_dj)
       end if


       ! Broadcast omega_ic and omega_ma
       call MPI_Bcast(omega_ic,1,MPI_DOUBLE_PRECISION,rank_with_l1m0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(omega_ma,1,MPI_DOUBLE_PRECISION,rank_with_l1m0,MPI_COMM_WORLD,ierr)

#endif

       !--- Now the real work starts with the radial loop that calculates
       !    the nonlinear terms:
       IF ( lVerbose ) THEN 
          WRITE(*,*)
          WRITE(*,*) '! Starting radial loop!'
       END IF

       !PERFOFF
       ! =============================== BARRIER ===========================
       PERFON('barr_2')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       ! ===================================================================

       CALL wallTime(runTimeRstart)
       CALL radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,      &
            &           lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,      &
            &           dsdt_Rloc,dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dbdt_Rloc,djdt_Rloc,&
            &           dVxBhLM_Rloc,dVSrLM_Rloc, &
            &           lorentz_torque_ic,lorentz_torque_ma,      &
            &           br_vt_lm_cmb,br_vp_lm_cmb,      &
            &           br_vt_lm_icb,br_vp_lm_icb,      &
            &           HelLMr_Rloc,Hel2LMr_Rloc,HelnaLMr_Rloc,Helna2LMr_Rloc,&
            &           uhLMr_Rloc,duhLMr_Rloc,dtrkc_Rloc,dthkc_Rloc)


       IF ( lVerbose ) WRITE(*,*) '! r-loop finished!'
       IF ( .NOT.l_log ) THEN
          CALL wallTime(runTimeRstop)
          IF ( .NOT.lNegTime(runTimeRstart,runTimeRstop) ) THEN
             nTimeR=nTimeR+1
             CALL subTime(runTimeRstart,runTimeRstop,runTimePassed)
             CALL addTime(runTimeR,runTimePassed)
          END IF
       END IF

       !---------------------------------------
       !--- Gather all r-distributed arrays ---
       !---------------------------------------
       ! we have here the following arrays from radialLoopG:
       !  dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM
       !  TstrRLM,TadvRLM,TomeRLM, 
       !  HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,dtrkc,dthkc
       !
       ! gather in-place, we need allgatherV because auf the unequal
       ! number of points on the processes (last block is one larger)


#ifdef WITH_MPI
       IF (DEBUG_OUTPUT) THEN
          nR_i1=MAX(2,nRstart)
          nR_i2=MIN(n_r_max-1,nRstop)
          WRITE(*,"(A,10ES20.12)") "middl: dwdt,dsdt,dzdt,dpdt = ",&
               & get_global_sum( dwdt_Rloc(:,nR_i1:nR_i2) ),&
               & get_global_sum( dsdt_Rloc(:,nR_i1:nR_i2) ),&
               & get_global_sum( dzdt_Rloc(:,nR_i1:nR_i2) ),&
               & get_global_sum( dpdt_Rloc(:,nR_i1:nR_i2) ),get_global_sum( dVSrLM_Rloc )
          IF (l_mag) THEN
             WRITE(*,"(A,6ES20.12)") "middl: dbdt,djdt,dVxBhLM = ",&
                  & get_global_SUM( dbdt_Rloc(:,nR_i1:nR_i2) ),&
                  & get_global_SUM( djdt_Rloc(:,nR_i1:nR_i2) ),get_global_sum( dVxBhLM_Rloc )
          END IF
       END IF

       ! ===================================== BARRIER =======================
       PERFON('barr_rad')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       ! =====================================================================
       !CALL r2lm_redist(dsdt_Rloc,dsdt_LMloc)
       CALL r2lo_redist(dsdt_Rloc,dsdt_LMloc)
       !CALL r2lm_redist(dwdt_Rloc,dwdt_LMloc)
       CALL r2lo_redist(dwdt_Rloc,dwdt_LMloc)
       !CALL r2lm_redist(dzdt_Rloc,dzdt_LMloc)
       CALL r2lo_redist(dzdt_Rloc,dzdt_lo)
       !CALL r2lm_redist(dpdt_Rloc,dpdt_LMloc)
       call r2lo_redist(dpdt_Rloc,dpdt_LMloc)
       !CALL r2lm_redist(dVSrLM_Rloc,dVSrLM_LMloc)
       CALL r2lo_redist(dVSrLM_Rloc,dVSrLM_LMloc)
       if (l_mag) then
          CALL r2lo_redist(dbdt_Rloc,dbdt_LMloc)
          CALL r2lo_redist(djdt_Rloc,djdt_LMloc)
          CALL r2lo_redist(dVxBhLM_Rloc,dVxBhLM_LMloc)
       end if


       PERFON('gather')

       !IF (lHelCalc) THEN
       !   sendcount  = (nRstop-nRstart+1)*(l_max+1)
       !   recvcounts = nr_per_rank*(l_max+1)
       !   recvcounts(n_procs-1) = (nr_per_rank+1)*(l_max+1)
       !   DO i=0,n_procs-1
       !      displs(i) = i*nr_per_rank*(l_max+1)
       !   END DO
       !   CALL MPI_AllGatherV(HelLMr_Rloc,sendcount,MPI_DOUBLE_PRECISION,&
       !        & HelLMr,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       !   CALL MPI_AllGatherV(Hel2LMr_Rloc,sendcount,MPI_DOUBLE_PRECISION,&
       !        & Hel2LMr,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       !   CALL MPI_AllGatherV(HelnaLMr_Rloc,sendcount,MPI_DOUBLE_PRECISION,&
       !        & HelnaLMr,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       !   CALL MPI_AllGatherV(Helna2LMr_Rloc,sendcount,MPI_DOUBLE_PRECISION,&
       !        & Helna2LMr,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       !END IF

       ! ------------------
       ! also exchange the lorentz_torques which are only set at the boundary points
       ! but are needed on all processes.
       CALL MPI_Bcast(lorentz_torque_ic,1,MPI_DOUBLE_PRECISION,n_procs-1,MPI_COMM_WORLD,ierr)
       CALL MPI_Bcast(lorentz_torque_ma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


       IF (DEBUG_OUTPUT) THEN
          WRITE(*,"(A,8ES20.12)") "lo_arr middl: dzdt_lo,z_LMloc,dz_LMloc,dzdtLast_lo = ",&
               & GET_GLOBAL_SUM( dzdt_lo(:,2:n_r_max-1) ),&
               & GET_GLOBAL_SUM( z_LMloc ),GET_GLOBAL_SUM( dz_LMloc ),get_global_sum( dzdtLast_lo )
          WRITE(*,"(A,8ES20.12)") "lo_arr middl: dsdt,s,ds,dsdtLast = ",&
               & GET_GLOBAL_SUM( dsdt_LMloc(:,2:n_r_max-1) ),&
               & GET_GLOBAL_SUM( s_LMloc ),GET_GLOBAL_SUM( ds_LMloc ),get_global_sum( dsdtLast_LMloc )
          WRITE(*,"(A,10ES20.12)") "lo_arr middl: dwdt,w,dw,ddw,dwdtLast = ",&
               & GET_GLOBAL_SUM( dwdt_LMloc(:,2:n_r_max-1) ),&
               & GET_GLOBAL_SUM( w_LMloc ),GET_GLOBAL_SUM( dw_LMloc ),&
               & GET_GLOBAL_SUM( ddw_LMloc ),get_global_sum( dwdtLast_LMloc )
          IF (l_mag) THEN
             WRITE(*,"(A,10ES20.12)") "lo_arr middl: dbdt,b,db,ddb,dbdtLast = ",&
                  & GET_GLOBAL_SUM( dbdt_LMloc(:,2:n_r_max-1) ),&
                  & GET_GLOBAL_SUM( b_LMloc ),GET_GLOBAL_SUM( db_LMloc ),&
                  & GET_GLOBAL_SUM( ddb_LMloc ),get_global_sum( dbdtLast_LMloc )
             WRITE(*,"(A,12ES20.12)") "lo_arr middl: djdt,aj,dj,ddj,djdtLast,dVxBhLM = ",&
                  & GET_GLOBAL_SUM( djdt_LMloc(:,2:n_r_max-1) ),&
                  & GET_GLOBAL_SUM( aj_LMloc ),GET_GLOBAL_SUM( dj_LMloc ),&
                  & GET_GLOBAL_SUM( ddj_LMloc ),get_global_sum( djdtLast_LMloc ), &
                  & get_global_sum( dVxBhLM_LMloc )
          END IF
       END IF

       PERFOFF

#endif
       !--- Output before update of fields in LMLoop:
       ! =================================== BARRIER ======================
       !PERFON('barr_4')
       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       !PERFOFF
       ! ==================================================================
       PERFON('output')
       IF (nRstart.LE.n_r_cmb) THEN
          ptr_dbdt_CMB => dbdt_Rloc(:,n_r_cmb)
       ELSE
          NULLIFY(ptr_dbdt_CMB)
       END IF
       CALL output(time,dt,dtNew,n_time_step,l_stop_time,                &
            &      l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,          &
            &      l_store,l_new_rst_file,                               &
            &      l_spectrum,lTOCalc,lTOframe,lTOZwrite,                &
            &      l_frame,n_frame,l_cmb,n_cmb_sets,                     &
            &      lorentz_torque_ic,lorentz_torque_ma,ptr_dbdt_CMB,     &
            &      HelLMr_Rloc,Hel2LMr_Rloc,HelnaLMr_Rloc,Helna2LMr_Rloc,&
            &      uhLMr_Rloc,duhLMr_Rloc)
       PERFOFF

       IF ( l_graph ) THEN
#ifdef WITH_MPI
          CALL MPI_File_close(graph_mpi_fh,ierr)
          !CALL MPI_ERROR_STRING(ierr,error_string,length_of_error,ierr)
          !PRINT*,"MPI_FILE_CLOSE returned: ",TRIM(error_string)

#else
          CLOSE(n_graph_file)  ! close graphic output file !
#endif
       END IF
       ! =================================== BARRIER ======================
       !PERFON('barr_5')
       !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       !PERFOFF
       ! ==================================================================

       !----- Finish time stepping, the last step is only for output!
       IF ( l_stop_time ) GOTO 6000  ! END OF TIME INTEGRATION


       !------ Nonlinear magnetic boundary conditions:
       !       For stressfree conducting boundaries
       PERFON('nl_m_bnd')
       IF ( l_b_nl_cmb ) THEN
          b_nl_cmb(1) =CMPLX(1.D0,1.D0,kind=KIND(b_nl_cmb))
          aj_nl_cmb(1)=CMPLX(1.D0,1.D0,kind=KIND(aj_nl_cmb))
          CALL get_b_nl_bcs('CMB', br_vt_lm_cmb,br_vp_lm_cmb,              &
               &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_cmb(2:lm_max))
       END IF
       IF ( l_b_nl_icb ) THEN
          aj_nl_icb(1)=CMPLX(1.D0,1.D0,kind=KIND(aj_nl_icb))
          CALL get_b_nl_bcs('ICB', br_vt_lm_icb,br_vp_lm_icb,              &
               &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_icb(2:lm_max))
       END IF
       PERFOFF

       PERFON('t_check')
       !---- Time-step check and change if needed (l_new_dtNext=.TRUE.)
       !     I anticipate the dt change here that is only used in 
       !     the next time step cause its needed for coex=(alpha-1)/w2 !
       dtLast=dt
       dt=dtNew        ! Update to new time step
       w2=w2New        ! Weight for time-derivatives of last time step
       w1=1.D0-w2      ! Weight for currect time step
       l_new_dt    =l_new_dtNext
       l_new_dtNext=.FALSE.

       !------ Checking Courant criteria, l_new_dt and dtNew are output
       IF ( l_cour ) THEN
          !PRINT*,"dtrkc_Rloc = ",dtrkc_Rloc
          CALL dt_courant(dtr,dth,l_new_dtNext,dt,dtNew,dtMax,dtrkc_Rloc,dthkc_Rloc)
       END IF

       !------ Checking whether we have to hit a specific time for output,
       !       dtNew is changed accordingly and l_new_dt is set .TRUE.
       IF ( l_true_time .AND. l_time_hits ) THEN
          CALL check_time_hits(l_new_dtHit,time,dt,dtNew)
          l_new_dtNext=l_new_dtNext .OR. l_new_dtHit
       END IF

       !----- Stop if time step has become too small:
       IF ( dtNew.LT.dtMin ) THEN
          IF (rank.EQ.0) THEN
             WRITE(*,'(1p,/,A,d14.4,/,A)') &
                  &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                  &" ! I THUS STOP THE RUN !"
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(1p,/,A,d14.4,/,A)') &
                  &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                  &" ! I THUS STOP THE RUN !"
             CALL safeClose(nLF)
          END IF
          STOP
       END IF
       IF ( l_new_dtNext ) THEN
          !------ Writing info and getting new weights:
          w2New=-0.5D0*dtNew/dt ! Weight I will be using for next update !
          n_dt_changed=0
          lCourChecking=.TRUE.
          IF (rank.EQ.0) THEN
             WRITE(*,'(1p,/,A,d18.10,/,A,i9,/,A,d15.8,/,A,d15.8)') &
                  &" ! Changing time step at time=",(time+dt)*tScale,&
                  &"                 time step no=",n_time_step+1,&
                  &"                      last dt=",dt*tScale,&
                  &"                       new dt=",dtNew*tScale
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(1p,/,A,d18.10,/,A,i9,/,A,d15.8,/,A,d15.8)') &
                  &" ! Changing time step at time=",(time+dt)*tScale,&
                  &"                 time step no=",n_time_step+1,&
                  &"                      last dt=",dt*tScale,&
                  &"                       new dt=",dtNew*tScale
             CALL safeClose(nLF)
          END IF
       ELSE
          w2New=-0.5D0 ! Normal weight if timestep is not changed !
          n_dt_changed=n_dt_changed+1
          IF ( n_dt_changed.LE.n_dt_check  ) THEN
             lCourChecking=.TRUE.
          ELSE
             lCourChecking=.FALSE.
          END IF
       END IF
       PERFOFF


       !----- UPDATING THE FIELDS:
       !      This is the second parallel part. Here we parallize over lm.

       !----- Advancing time:
       coex  =(alpha-1.D0)/w2New
       timeLast        =time               ! Time of last time step
       time            =time+dt            ! Update time
       timeScaled      =time*tScale

       lMat=.FALSE.
       IF ( l_new_dt ) THEN
          !----- Calculate matricies for new time step if dt.NE.dtLast
          lMat=.TRUE.
          IF (rank.EQ.0) THEN
             WRITE(*,'(1p,/,'' ! BUILDING MATRICIES AT STEP/TIME:'',   &
                  &              i8,d16.6)') n_time_step,timeScaled
          END IF
       END IF

       ! =================================== BARRIER ======================
       PERFON('barr_6')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       ! ==================================================================

       CALL wallTime(runTimeRstart)
       IF ( lVerbose ) WRITE(*,*) '! starting lm-loop!'
       IF ( lVerbose ) WRITE(*,*) '! No of lm-blocks:',nLMBs

       CALL LMLoop(w1,coex,time,dt,lMat,lRmsNext, &
            &      dVxBhLM_LMloc,dVSrLM_LMloc,&
            &      dsdt_LMloc,dwdt_LMloc,dzdt_lo,dpdt_LMloc,dbdt_LMloc,djdt_LMloc,      &
            &      lorentz_torque_ma,lorentz_torque_ic,      &
            &      b_nl_cmb,aj_nl_cmb,aj_nl_icb,n_time_step)

       IF ( lVerbose ) WRITE(*,*) '! lm-loop finished!'
       CALL wallTime(runTimeRstop)

       IF ( .NOT.lNegTime(runTimeRstart,runTimeRstop) ) THEN
          nTimeLM=nTimeLM+1
          CALL subTime(runTimeRstart,runTimeRstop,runTimePassed)
          CALL addTime(runTimeLM,runTimePassed)
       END IF

       IF (DEBUG_OUTPUT) THEN
          WRITE(*,"(A,6ES20.12)") "lo_arr end: z_LMloc,dz_LMloc,dzdtLast_lo = ",&
               & GET_GLOBAL_SUM( z_LMloc ),GET_GLOBAL_SUM( dz_LMloc ),get_global_sum( dzdtLast_lo )
          WRITE(*,"(A,6ES20.12)") "lo_arr end: s,ds,dsdtLast = ",&
               & GET_GLOBAL_SUM( s_LMloc ),GET_GLOBAL_SUM( ds_LMloc ),get_global_sum( dsdtLast_LMloc )
          WRITE(*,"(A,6ES20.12)") "lo_arr end: w,dw,dwdtLast = ",&
               & GET_GLOBAL_SUM( w_LMloc ),GET_GLOBAL_SUM( dw_LMloc ),get_global_sum( dwdtLast_LMloc )
          WRITE(*,"(A,6ES20.12)") "lo_arr end: p,dpdtLast = ",&
               & GET_GLOBAL_SUM( p_LMloc ),get_global_sum( dpdtLast_LMloc )
          WRITE(*,"(A,6ES20.12)") "lo_arr end: b,db,dbdtLast = ",&
               & GET_GLOBAL_SUM( b_LMloc ),GET_GLOBAL_SUM( db_LMloc ),get_global_sum( dbdtLast_LMloc )
          WRITE(*,"(A,6ES20.12)") "lo_arr end: aj,dj,djdtLast = ",&
               & GET_GLOBAL_SUM( aj_LMloc ),GET_GLOBAL_SUM( dj_LMloc ),get_global_sum( djdtLast_LMloc )
       END IF

       !----- Timing and info of advancement:
       ! =================================== BARRIER ======================
       !start_time=MPI_Wtime()
       PERFON('barr_lm')
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
       PERFOFF
       !end_time=MPI_Wtime()
       !WRITE(*,"(A,I4,A,F10.6,A)") " lm barrier on rank ",rank,&
       !     & " takes ",end_time-start_time," s."
       ! ==================================================================
       CALL wallTime(runTimeTstop)
       IF ( .NOT.lNegTime(runTimeTstart,runTimeTstop) ) THEN
          CALL subTime(runTimeTstart,runTimeTstop,runTimePassed)
          IF ( lMat .AND. .NOT.l_log ) THEN
             nTimeTM=nTimeTM+1
             CALL addTime(runTimeTM,runTimePassed)
          ELSE IF ( .NOT.lMat .AND. l_log ) THEN
             nTimeTL=nTimeTL+1
             CALL addTime(runTimeTL,runTimePassed)
          ELSE IF ( .NOT.lMat .AND. .NOT.l_log ) THEN
             nTimeT=nTimeT+1
             CALL addTime(runTimeT,runTimePassed)
          END IF
          CALL addTime(runTime,runTimePassed)
       END IF
       IF ( DBLE(n_time_step)+tenth_n_time_steps*DBLE(nPercent).GE.DBLE(n_time_steps) &
            & .OR. n_time_steps.LT.31 ) THEN
          WRITE(message,'(" ! Time step finished:",i6)') n_time_step
          call logWrite(message)
          IF ( DBLE(n_time_step)+tenth_n_time_steps*DBLE(nPercent).GE.DBLE(n_time_steps) &
               & .AND. n_time_steps.GE.10 ) THEN
             WRITE(message,'(" ! This is           :",i3,"%")') (10-nPercent)*10
             call logWrite(message)
             nPercent=nPercent-1
          END IF
          DO n=1,4
             runTimePassed(n)=runTimeT(n)
          END DO
          IF ( nTimeT.GT.0 ) THEN
             CALL meanTime(runTimePassed,nTimeT)
             IF (rank.EQ.0) THEN
                CALL writeTime(6,'! Mean wall time for time step:', runTimePassed)
                CALL safeOpen(nLF,log_file)
                CALL writeTime(nLF,'! Mean wall time for time step:', runTimePassed)
                CALL safeClose(nLF)
             END IF
          END IF
       END IF

    END DO ! end of time stepping !
    PERFOFF

6000 CONTINUE  ! jump point for termination upon "kill -30" signal

    IF ( l_movie ) THEN
       IF (rank.EQ.0) THEN
          IF (n_frame.GT.0) THEN
             WRITE(*,'(1p,/,/,A,i10,3(/,A,d16.6))') &
                  &" !  No of stored movie frames: ",n_frame,&
                  &" !     starting at time: ",t_movieS(1)*tScale,&
                  &" !       ending at time: ",t_movieS(n_frame)*tScale,&
                  &" !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(1p,/,/,A,i10,3(/,A,d16.6))') &
                  &" !  No of stored movie frames: ",n_frame,&
                  &" !     starting at time: ",t_movieS(1)*tScale,&
                  &" !       ending at time: ",t_movieS(n_frame)*tScale,&
                  &" !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
             CALL safeClose(nLF)
          ELSE
             WRITE(*,'(1p,/,/,A,i10,3(/,A,d16.6))') &
                  &" !  No of stored movie frames: ",n_frame,&
                  &" !     starting at time: ",0.0D0,&
                  &" !       ending at time: ",0.0D0,&
                  &" !      with step width: ",0.0D0
             CALL safeOpen(nLF,log_file)
             WRITE(nLF,'(1p,/,/,A,i10,3(/,A,d16.6))') &
                  &" !  No of stored movie frames: ",n_frame,&
                  &" !     starting at time: ",0.0D0,&
                  &" !       ending at time: ",0.0D0,&
                  &" !      with step width: ",0.0D0
             CALL safeClose(nLF)
          END IF
       END IF
    END IF

    IF ( l_cmb_field ) THEN
       WRITE(message,'(A,i9)') " !  No of stored sets of b at CMB: ",n_cmb_sets
       CALL logWrite(message)
    END IF

    CALL meanTime(runTimeR, nTimeR)
    CALL meanTime(runTimeLM,nTimeLM)
    CALL meanTime(runTimeTM,nTimeTM)
    CALL meanTime(runTimeTL,nTimeTL)
    CALL meanTime(runTimeT,nTimeT)
    IF (rank.EQ.0) THEN
       CALL writeTime(6,                                               &
            &     '! Mean wall time for r Loop                :',runTimeR)
       CALL writeTime(6,                                               &
            &     '! Mean wall time for LM Loop               :',runTimeLM)
       CALL writeTime(6,                                               &
            &     '! Mean wall time for t-step with matix calc:',runTimeTM)
       CALL writeTime(6,                                               &
            &     '! Mean wall time for t-step with log output:',runTimeTL)
       CALL writeTime(6,                                               &
            &     '! Mean wall time for pure t-step           :',runTimeT)
       CALL safeOpen(nLF,log_file)
       CALL writeTime(nLF,                                              &
            &     '! Mean wall time for r Loop                :',runTimeR)
       CALL writeTime(nLF,                                             &
            &     '! Mean wall time for LM Loop               :',runTimeLM)
       CALL writeTime(nLF,                                             &
            &     '! Mean wall time for t-step with matix calc:',runTimeTM)
       CALL writeTime(nLF,                                             &
            &     '! Mean wall time for t-step with log output:',runTimeTL)
       CALL writeTime(nLF,                                             &
            &     '! Mean wall time for pure t-step           :',runTimeT)
       CALL safeClose(nLF)
    END IF
    !-- Write output for variable conductivity test:
    !       IF ( imagcon.EQ.-10 ) THEN
    !          message='testVarCond.'//tag
    !          OPEN(99,FILE=message,STATUS='UNKNOWN')
    !           DO nR=n_r_max,1,-1             ! Diffusive toroidal field
    !             WRITE(99,*) r(nR),REAL(aj(4,nR)),jVarCon(nR)
    !           END DO
    !          CLOSE(99)
    !       END IF


    !-- WORK IS DONE !
#ifdef WITH_MPI
    CALL MPI_INFO_FREE(info,ierr)
    DEALLOCATE(recvcounts,displs)
#endif

    RETURN
  END SUBROUTINE step_time

  !----------------------------------------------------------------------------
  !-- End of subroutine step_time
  !----------------------------------------------------------------------------

END MODULE step_time_mod
