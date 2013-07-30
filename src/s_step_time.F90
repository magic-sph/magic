!$Id$
#include "perflib_preproc.cpp"
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
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


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
        USE output_data
        use output_mod
        USE const
        USE timing
        USE charmanip, only: capitalize,dble2str
        USE usefull, only: l_correct_step

        IMPLICIT NONE 

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
        LOGICAL :: lHelCalc         ! Calculate helocoty for output
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

        CHARACTER(len=76) :: message
        CHARACTER(len=76) :: SIG    

!--- Courant criteria/diagnosis:
        REAL(kind=8) :: dtr,dth
        REAL(kind=8) :: dtrkc(n_r_max),dthkc(n_r_max)   ! Saves values for time step

!--- Explicit part of time stepping, calculated in s_radialLoopG.f and
!    passed to s_LMLoop.f where the time step is preformed.
!    Note that the respective arrays for the changes in inner-core
!    magnetic field are calculated in s_updateB.f and are only
!    needed there.
        COMPLEX(kind=8) :: dwdt(lm_max,n_r_max)
        COMPLEX(kind=8) :: dzdt(lm_max,n_r_max)
        COMPLEX(kind=8) :: dpdt(lm_max,n_r_max)
        COMPLEX(kind=8) :: dsdt(lm_max,n_r_max)
        COMPLEX(kind=8) :: dbdt(lm_maxMag,n_r_maxMag)
        COMPLEX(kind=8) :: djdt(lm_maxMag,n_r_maxMag)

!--- Help arrays for calculating dsdt and djdt (see s_updateS.f,s_updateB.f):
        COMPLEX(kind=8) :: dVxBhLM(lm_maxMag,n_r_maxMag)
        COMPLEX(kind=8) :: dVSrLM(lm_max,n_r_max)

!--- Lorentz torques:
        REAL(kind=8) :: lorentz_torque_ma,lorentz_torque_ic

!--- Help arrays for magnetic field stretching and advection (dtB):
        COMPLEX(kind=8) :: TstrRLM(lm_max_dtB,n_r_max_dtB)
        COMPLEX(kind=8) :: TadvRLM(lm_max_dtB,n_r_max_dtB)
        COMPLEX(kind=8) :: TomeRLM(lm_max_dtB,n_r_max_dtB)

!--- Help arrays for calculating axsymmetric helicity in s_outMisc.f:
        REAL(kind=8) :: HelLMr(l_max+1,n_r_max)
        REAL(kind=8) :: Hel2LMr(l_max+1,n_r_max)
        REAL(kind=8) :: HelnaLMr(l_max+1,n_r_max)
        REAL(kind=8) :: Helna2LMr(l_max+1,n_r_max)
        REAL(kind=8) :: uhLMr(l_max+1,n_r_max)
        REAL(kind=8) :: duhLMr(l_max+1,n_r_max)

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


!-- end of declaration
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+


        IF ( lVerbose ) WRITE(*,'(/,'' ! STARTING STEP_TIME !'')')

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
        message='signal'//'.'//tag
        OPEN(19,FILE=message,STATUS='unknown')
        WRITE(19,'(A3)') 'NOT'
        CLOSE(19)

!-- STARTING THE TIME STEPPING LOOP:
        WRITE(*,*)
        WRITE(*,*) '! Starting time integration!'
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
        PERFON('tloop')
        DO n_time_step=1,n_time_steps_go 
           n_time_cour=n_time_cour+1

           IF ( lVerbose) THEN 
              WRITE(*,*)
              WRITE(*,*) '! Starting time step ',n_time_step
           END IF
        
!----- Signalling via file signal:
           message='signal'//'.'//tag
           OPEN(19,FILE=message,STATUS='old')
           READ(19,*) SIG
           CLOSE(19)
           IF ( len(trim(SIG)).GT.0 ) THEN ! Non blank string ?
              CALL capitalize(SIG)
              IF ( index(SIG,'END')/=0 ) n_stop_signal=1
              n_graph_signal=0
              IF ( index(SIG,'GRA')/=0 ) THEN 
                 n_graph_signal=1
                 OPEN(19,FILE=message,STATUS='unknown')
                 WRITE(19,'(A3)') 'NOT'
                 CLOSE(19)
              END IF
              n_rst_signal=0
              IF ( index(SIG,'RST')/=0 ) THEN
                 n_rst_signal=1
                 OPEN(19,FILE=message,STATUS='unknown')
                 WRITE(19,'(A3)') 'NOT'
                 CLOSE(19)
              END IF
              n_spec_signal=0
              IF ( index(SIG,'SPE')/=0 ) THEN
                 n_spec_signal=1
                 OPEN(19,FILE=message,STATUS='unknown')
                 WRITE(19,'(A3)') 'NOT'
                 CLOSE(19)
              END IF
           END IF
           
           CALL wallTime(runTimeTstart)

!--- Various reasons to stop the time integration:
           IF ( l_runTimeLimit ) THEN
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

!-- Checking logic for output: 
           l_graph=                                                     &
     &         l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
     &           n_graph_step,n_graphs,n_t_graph,t_graph,0) .OR.        &
     &             n_graph_signal.EQ.1
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
           l_HT=l_frame .and. l_HTmovie

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


!----- Checking whether Courant condition should be checked:
!           IF ( mode.NE.5 .AND. n_cour_step.GT.0 .AND. (
!     &          MOD(n_time_cour,n_cour_step).EQ.0 .OR.
!     &          lCourChecking ) ) THEN
!              l_cour=.TRUE.
!           ELSE
!              l_cour=.FALSE.
!           END IF
            l_cour=.TRUE.

           l_new_rst_file=                                              &
     &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
     &                       n_rst_step,n_rsts,n_t_rst,t_rst,0) .OR.    &
     &        n_rst_signal.EQ.1
           n_rst_signal=0
           l_store= l_new_rst_file .OR.                                 &
     &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
     &                       0,n_stores,0,t_rst,0)
!           IF ( l_store .AND. lCourChecking ) THEN
!              l_store     =.FALSE.
!              l_storeSched=.TRUE.
!           ELSE IF ( l_storeSched .AND. .NOT.lCourChecking ) THEN
!              l_store     =.TRUE.
!              l_storeSched=.FALSE.
!           END IF
!           IF ( l_store ) n_time_cour=-1 ! Courant check after following update

           l_log=                                                       &
     &        l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
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
              IF ( ngform.NE.0 ) THEN
                 OPEN(n_graph_file,                                     &
     &                FILE=graph_file,STATUS='NEW',FORM='FORMATTED')
              ELSE
                 OPEN(n_graph_file,                                     &
     &                FILE=graph_file,STATUS='NEW',FORM='UNFORMATTED')
              END IF
           END IF


!--- Now the real work starts with the radial loop that calucates
!    the nonlinear terms:
           IF ( lVerbose ) THEN 
              WRITE(*,*)
              WRITE(*,*) '! Starting radial loop!'
           END IF
           CALL wallTime(runTimeRstart)

           CALL radialLoopG(l_graph,l_cour,l_frame,time,dt,dtLast,      &
     &                 lTOCalc,lTONext,lTONext2,lHelCalc,lRmsCalc,      &
     &               dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM,      &
     &                        lorentz_torque_ic,lorentz_torque_ma,      &
     &                                  br_vt_lm_cmb,br_vp_lm_cmb,      &
     &                                  br_vt_lm_icb,br_vp_lm_icb,      &
     &                                    TstrRLM,TadvRLM,TomeRLM,      &
     &             HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,      &
     &                                                dtrkc,dthkc)

           IF ( lVerbose ) WRITE(*,*) '! r-loop finished!'
           IF ( .NOT.l_log ) THEN
              CALL wallTime(runTimeRstop)
              IF ( .NOT.lNegTime(runTimeRstart,runTimeRstop) ) THEN
                 nTimeR=nTimeR+1
                 CALL subTime(runTimeRstart,runTimeRstop,runTimePassed)
                 CALL addTime(runTimeR,runTimePassed)
              END IF
           END IF

!--- Output before update of fields in LMLoop:
           CALL output(time,dt,dtNew,n_time_step,l_stop_time,           &
     &           l_Bpot,l_Vpot,l_Tpot,l_log,l_graph,lRmsCalc,           &
     &                                l_store,l_new_rst_file,           &
     &                 l_spectrum,lTOCalc,lTOframe,lTOZwrite,           &
     &                      l_frame,n_frame,l_cmb,n_cmb_sets,           &
     &              lorentz_torque_ic,lorentz_torque_ma,dbdt,           &
     &                               TstrRLM,TadvRLM,TomeRLM,           &
     &        HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr)

           IF ( l_graph ) CLOSE(n_graph_file)  ! close graphic output file !

!----- Finish time stepping, the last step is only for output!
           IF ( l_stop_time ) GOTO 6000  ! END OF TIME INTEGRATION


!------ Nonlinear magnetic boundary conditions:
!       For stressfree conducting boundaries
           IF ( l_b_nl_cmb ) THEN
              b_nl_cmb(1) =CMPLX(1.D0,1.D0,KIND=KIND(0d0))
              aj_nl_cmb(1)=CMPLX(1.D0,1.D0,KIND=KIND(0d0))
              CALL get_b_nl_bcs('CMB',                                  &
     &                          br_vt_lm_cmb,br_vp_lm_cmb,              &
     &                          2,lm_max,b_nl_cmb,aj_nl_cmb)
           END IF
           IF ( l_b_nl_icb ) THEN
              aj_nl_icb(1)=CMPLX(1.D0,1.D0,KIND=KIND(0d0))
              CALL get_b_nl_bcs('ICB',                                  &
     &                          br_vt_lm_icb,br_vp_lm_icb,              &
     &                          2,lm_max,b_nl_cmb,aj_nl_icb)
           END IF


!---- Time-step check and change if needed (l_new_dtNext=.TRUE.)
!     I anticipate the dt change here that is only used in 
!     the next time step cause its needed for coex=(alpha-1)/w2 !
           dtLast=dt
           dt=dtNew        ! Update to new time step
           w2=w2New        ! Weight for time-derivatives of last time step
           w1=1.D0-w2      ! Weight for currect time step
           l_new_dt    =l_new_dtNext
           l_new_dtNext=.FALSE.

!------ Checking Courant criteria, l_new_dt and dtNew areoutput
           IF ( l_cour )                                                &
     &        CALL dt_courant(dtr,dth,l_new_dtNext,dt,dtNew,            &
     &                          dtMax,dtrkc,dthkc)

!------ Checking whether we have to hit a specific time for output,
!       dtNew is changed accordingly and l_new_dt is set .TRUE.
           IF ( l_true_time .AND. l_time_hits ) THEN
              CALL check_time_hits(l_new_dtHit,time,dt,dtNew)
              l_new_dtNext=l_new_dtNext .OR. l_new_dtHit
           END IF

!----- Stop if time step has become too small:
           IF ( dtNew.LT.dtMin ) THEN
              WRITE(*,'(1p,/,A,d14.4,/,A)') &
                   &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                   &" ! I THUS STOP THE RUN !"
              CALL safeOpen(nLF,log_file)
              WRITE(nLF,'(1p,/,A,d14.4,/,A)') &
                   &" ! TIME STEP TOO SMALL, dt=",dtNew,&
                   &" ! I THUS STOP THE RUN !"
              CALL safeClose(nLF)
              STOP
           END IF
           IF ( l_new_dtNext ) THEN
!------ Writing info and getting new weights:
              w2New=-0.5D0*dtNew/dt ! Weight I will be using for next update !
              n_dt_changed=0
              lCourChecking=.TRUE.
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
           ELSE
              w2New=-0.5D0 ! Normal weight if timestep is not changed !
              n_dt_changed=n_dt_changed+1
              IF ( n_dt_changed.LE.n_dt_check  ) THEN
                 lCourChecking=.TRUE.
              ELSE
                 lCourChecking=.FALSE.
              END IF
           END IF


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
              WRITE(*,'(1p,/,'' ! BUILDING MATRICIES AT STEP/TIME:'',   &
     &              i8,d16.6)') n_time_step,timeScaled
           END IF

           CALL wallTime(runTimeRstart)
           IF ( lVerbose ) WRITE(*,*) '! starting lm-loop!'
           IF ( lVerbose ) WRITE(*,*) '! No of lm-blocks:',nLMBs

           CALL LMLoop(w1,coex,time,dt,lMat,lRmsNext,                   &
     &               dVxBhLM,dVSrLM,dsdt,dwdt,dzdt,dpdt,dbdt,djdt,      &
     &                        lorentz_torque_ma,lorentz_torque_ic,      &
     &                               b_nl_cmb,aj_nl_cmb,aj_nl_icb)

           IF ( lVerbose ) WRITE(*,*) '! lm-loop finished!'
           CALL wallTime(runTimeRstop)
           IF ( .NOT.lNegTime(runTimeRstart,runTimeRstop) ) THEN
              nTimeLM=nTimeLM+1
              CALL subTime(runTimeRstart,runTimeRstop,runTimePassed)
              CALL addTime(runTimeLM,runTimePassed)
           END IF

!----- Timing and info of advancement:
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
           IF ( DBLE(n_time_step)+tenth_n_time_steps*DBLE(nPercent)     &
     &          .GE.DBLE(n_time_steps) .OR. n_time_steps.LT.31          &
     &        ) THEN
              WRITE(*,'(/," ! Time step finished:",i6)') n_time_step
              CALL safeOpen(nLF,log_file)
              WRITE(nLF,'(/," ! Time step finished:",i6)') n_time_step
              CALL safeClose(nLF)
              IF ( DBLE(n_time_step)+tenth_n_time_steps*DBLE(nPercent)  &
     &           .GE.DBLE(n_time_steps) .AND. n_time_steps.GE.10 ) THEN
                 WRITE(*,'(" ! This is           :",i3,"%")')           &
     &                 (10-nPercent)*10
                 CALL safeOpen(nLF,log_file)
                 WRITE(nLF,'(" ! This is           :",i3,"%")')         &
     &                 (10-nPercent)*10
                 CALL safeClose(nLF)
                 nPercent=nPercent-1
              END IF
              DO n=1,4
                 runTimePassed(n)=runTimeT(n)
              END DO
              IF ( nTimeT.GT.0 ) THEN
                 CALL meanTime(runTimePassed,nTimeT)
                 CALL writeTime(6,'! Mean wall time for time step:',    &
     &                          runTimePassed)
                 CALL safeOpen(nLF,log_file)
                 CALL writeTime(nLF,'! Mean wall time for time step:',  &
     &                          runTimePassed)
                 CALL safeClose(nLF)
              END IF
           END IF


        END DO ! end of time stepping !
        PERFOFF

6000    CONTINUE  ! jump point for termination upon "kill -30" signal

        IF ( l_movie ) THEN
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

        IF ( l_cmb_field ) THEN
           WRITE(*,'(//," !  No of stored sets of b at CMB: ",i9)')     &
     &        n_cmb_sets
           CALL safeOpen(nLF,log_file)
           WRITE(nLF,'(//," !  No of stored sets of b at CMB: ",i9)')   &
     &        n_cmb_sets
           CALL safeClose(nLF)
           CALL logWrite(message)
        END IF

        CALL meanTime(runTimeR, nTimeR)
        CALL meanTime(runTimeLM,nTimeLM)
        CALL meanTime(runTimeTM,nTimeTM)
        CALL meanTime(runTimeTL,nTimeTL)
        CALL meanTime(runTimeT,nTimeT)
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

        RETURN
        END

!----------------------------------------------------------------------------
!-- End of subroutine step_time
!----------------------------------------------------------------------------
