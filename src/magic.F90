!$Id$
#include "perflib_preproc.cpp"
!***********************************************************************
    PROGRAM magic3_44
!***********************************************************************

!--+-------------+----------------+------------------------------------+

!     A dynamic dynamo model driven by thermal convection
!     in a rotating spherical fluid shell.
!     This version is restricted to Boussinesq fluids and
!     non-dimensional variables are used throughout.

!     The following set of equations is solved:

!     E {dv/dt + v.grad(v)} = -grad(p) - 2e_z x v
!         +1/Pm rot(B) x B + RaE/Pr g/g_o T

!     dB/dt = rot(v x B) + 1/Pm Lapl(B)

!     dT/dt + v.grad(T) = 1/Pr Lapl(T) + epsc0

!     div(v)=0          div(B)=0

!       subject to the following boundary conditions
!       at the inner and outer radii:

!       v_r=0, and either no slip or stress free
!       T=0 / T=1  or fixed heat flux (the latter not tested!)
!       B fitted to exterior potential fields, or parts of B
!       specified on the boundaries

!     List of symbols:

!     v: velocity          p: pressure        B: magnetic field
!     g: gravity           g_o: reference value at outer radius
!     T: temperature       epsc0: rate of internal heating
!     e_z: unit vector parallel to the rotation axis
!     d/dt: partial time derivative  Lapl: Laplace operator

!     Scaling properties:

!     nu: kinematic viscosity         d: shell width
!     omega: angular frequency        alpha: thermal expansion coeff
!     delta_T: temperature contrast   kappa: thermal diffusivity
!     eta: magnetic diffusivity       rho: density
!     mu_o: magnetic permeability

!     Scaling:

!     Length:   d              time:      d^2/nu
!     Velocity: nu/d           pressure:  rho*nu*omega
!     Temperature: delta_T     mag.field: sqrt(rho*mu_o*eta*omega)


!     Non-dimensional numbers:

!     E: Ekman number     E= nu*d^2/omega
!     Ra: Rayleigh number Ra = alpha*g_o*delta_T*d^3/(kappa*nu)
!     Pr: Prandtl number  Pr = nu/kappa
!     Pm: Magnetic Prandtl number    Pm=nu/eta


!     Numerical simulations via a nonlinear, multimode,
!     initial-boundary value problem.

! *** entropy boundary conditions (tops and bots on input)
!     if ktops = 1, entropy specified on outer boundary
!     if ktops = 2, radial heat flux specified on outer boundary
!     if kbots = 1, entropy specified on inner boundary
!     if kbots = 2, radial heat flux specified on inner boundary
!     for example: ktops=1,
!           the spherically-symmetric temperature
!           on the outer boundary relative to the reference state

! *** velocity boundary condtions
!     if ktopv = 1, stress-free outer boundary
!     if ktopv = 2, non-slip outer boundary
!     if kbotv = 1, stress-free inner boundary
!     if kbotv = 2, non-slip inner boundary

! *** magnetic boundary condtions
!     if ktopb = 1, insulating outer boundary (mantle)
!     if kbotb = 1, insulating inner boundary (core)
!     if ktopb = 2, perfectly conducting outer boundary
!     if kbotb = 2, perfectly conducting inner boundary
!     if ktopb = 3, finitely conducting mantle
!     if kbotb = 3, finitely conducting inner core

! *** magneto-convection
!     amp_b1 = max amplitude of imposed magnetic field
!     if imagcon .eq. 1, imposed toroidal field via inner bc on J(l=2,m=0)
!     if imagcon .eq.10, imposed tor. field on both icb and cmb J(l=2,m=0)
!     if imagcon .eq.11, imposed tor. field on both icb and cmb J(l=2,m=0)
!                        opposite sign
!     if imagcon .eq.12, imposed tor. field on both icb and cmb J(l=1,m=0)
!     if imagcon .lt. 0, imposed poloidal field via inner bc on B(l=1,m=0)


!     if l_start_file=.true. initial fields are read from file $start_file$
!     if l_start_file=.false. start fields are initialised
!     according to init_s1,init_s2,init_b1,init_v1.
!     (see subroutines init_s, init_b and init_v.

!     Resolution is defined in truncation !

!     Subroutine step_time performes n_time_steps time steps.
!--+-------------------------------------------------------------------+

    USE truncation
    USE physical_parameters
    USE radial_functions
    USE num_param
    USE torsional_oscillations
    use init_fields
    use Grenoble
    use blocking
    use horizontal_data
    USE logic
    use matrices
    use fields
    use fieldsLast
    use movie_data
    use RMS,only: initialize_RMS
    use dtB_mod
    USE radial_data,only: initialize_radial_data
    use radialLoop
    USE lmLoop_data,ONLY: initialize_LMLoop_data
    USE LMLoop_mod,only: initialize_LMLoop
    USE kinetic_energy
    use magnetic_energy
    use fields_average_mod
    use Egeos_mod
    use spectrum_average_mod
    use spectrumC_average_mod
    USE output_data
    USE output_mod,only: initialize_output
    use outPV3, only:initialize_outPV3
    USE outTO_mod,only: initialize_outTO_mod
    USE parallel_mod
    use Namelists
    use step_time_mod, only: initialize_step_time, step_time
    USE timing, only: writeTime,wallTime
    USE communications, only:initialize_communications
    USE power, only: initialize_output_power
    use outPar_mod
    !USE rIterThetaBlocking_mod,ONLY: initialize_rIterThetaBlocking
#ifdef WITH_LIKWID
#   include "likwid_f90.h"
#endif
    IMPLICIT NONE

!-- Local variables:
    INTEGER :: n_time_step
    INTEGER :: n_time_step_start   ! storing initial time step no
    INTEGER :: n                   ! counter
    INTEGER :: nO                  ! output unit
    REAL(kind=8) :: time
    REAL(kind=8) :: dt
    REAL(kind=8) :: dtNew

    INTEGER :: n_stop_signal       ! signal returned from step_time
    COMMON/stop_signals/n_stop_signal

    CHARACTER(len=76) :: message

    ! MPI specific variables
    INTEGER :: required_level,provided_level

!-- end of declaration
!------------------------------------------------------------------------
#ifdef WITHOMP
    required_level=MPI_THREAD_MULTIPLE
    CALL mpi_init_thread(required_level,provided_level,ierr)
    If (provided_level.LT.required_level) THEN
       PRINT*,"We need at least thread level ",required_level,", but have ",provided_level
       stop
    END IF
#else
    call mpi_init(ierr)
#endif

    PERFINIT
    PERFON('main')
    LIKWID_INIT
    !LIKWID_ON('main')
    CALL parallel

!--- Read starting time
    IF ( rank.eq.0 ) THEN
       !CALL get_resetTime(resetTime)
        CALL wallTime(runTimeStart)
        WRITE(*,*)
        message='!--- PROGRAM MAGIC3.44, 31 Aug. 2010, BY JW ---!'
        WRITE(*,*) message
        CALL writeTime(6,'! Started at:',runTimeStart)
    END IF

!--- Read input parameters:
    CALL readNamelists  ! includes sent to other procs !

    call initialize_const
    !--- Check parameters and write info to SDTOUT
    CALL checkTruncation
    call openFiles

    call initialize_blocking
    call initialize_radial_data
    call initialize_radial_functions
    call initialize_radialLoop
    !CALL initialize_rIterThetaBlocking
    call initialize_LMLoop_data
    call initialize_LMLoop

    call initialize_num_param
    CALL initialize_TO
    IF (l_TO) THEN
       CALL initialize_outTO_mod
    END IF
    call initialize_init_fields
    call initialize_Grenoble
    call initialize_horizontal_data
    call initialize_matrices
    call initialize_fields
    call initialize_fieldsLast
    call initialize_movie_data
    if (l_RMS) call initialize_RMS
    call initialize_dtB_mod
    call initialize_kinetic_energy
    call initialize_magnetic_energy
    call initialize_fields_average_mod
    call initialize_Egeos_mod
    call initialize_spectrum_average_mod
    call initialize_spectrumC_average_mod
    if (l_PV) call initialize_outPV3
    call initialize_step_time
    call initialize_communications
    call initialize_outPar_mod
    IF ( l_power ) call initialize_output_power
    call initialize_output

!--- Open output files:
    !CALL openFiles

!--- Do pre-calculations:
    !CALL getBlocking
    CALL preCalc
    IF ( rank.eq.0 ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        CALL writeNamelists(6)
        CALL writeNamelists(n_log_file)
        IF ( l_save_out ) close(n_log_file)
    END IF

!--- Now read start-file or initialize fields:
    CALL getStartFields(time,dt,dtNew,n_time_step)

!--- Second pre-calculation:
    CALL preCalcTimes(time,n_time_step)

!--- Write info to STDOUT and log-file:
    IF ( rank.eq.0 ) THEN
        CALL writeInfo(6)
        CALL writeInfo(n_log_file)
    END IF

!--- Check whether movie data are required:
    CALL checkMovie

!--- AND NOW FOR THE TIME INTEGRATION:

!--- Write starting time to SDTOUT and logfile:
    IF ( rank.eq.0 ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
        END IF
        DO n=1,2
            IF ( n == 1 ) nO=6
            IF ( n == 2 ) nO=n_log_file
            WRITE(nO,'(/,'' ! STARTING TIME INTEGRATION AT:'')')
            WRITE(nO,'(''   start_time ='',1p,d18.10)') tScale*time
            WRITE(nO,'(''   step no    ='',i10)') n_time_step
            WRITE(nO,'(''   start dt   ='',1p,d16.4)') dt
            WRITE(nO,'(''   start dtNew='',1p,d16.4)') dtNew
        END DO
        IF ( l_save_out ) CLOSE(n_log_file)
    END IF
    timeStart        =time
    n_time_step_start=n_time_step

!--- Call time-integration routine:
    PERFON('steptime')
    CALL step_time(time,dt,dtNew,n_time_step)
    PERFOFF
!--- Write stop time to SDTOUR and logfile:
    IF ( rank.eq.0 ) THEN
        IF ( l_save_out ) THEN
            OPEN(n_log_file,FILE=log_file,STATUS='UNKNOWN', &
                 POSITION='APPEND')
        END IF

        DO n=1,2
            IF ( n == 1 ) nO=6
            IF ( n == 2 ) nO=n_log_file
            WRITE(nO,'(/,'' ! STOPPING TIME INTEGRATION AT:'')')
            WRITE(nO,'(''   stop time ='',1p,d18.10)') tScale*time
            WRITE(nO,'(''   stop step ='',i10)') &
                  n_time_step+n_time_step_start
            WRITE(nO,'(''   steps gone='',i10)') (n_time_step-1)
            WRITE(nO,*)
            IF ( n_stop_signal > 0 ) THEN
                WRITE(nO,*) '!!! MAGIC TERMINATED BY STOP SIGNAL !!!'
            ELSE
                WRITE(nO,*) '!!! REGULAR END OF PROGRAM MAGIC !!!'
            END IF
            WRITE(nO,*)
            !WRITE(nO,'('' max. thread number in  R-loop='',i3)') &
            !      nThreadsRmax
            !WRITE(nO,'('' max. thread number in LM-loop='',i3)') &
            !      nThreadsLMmax
            WRITE(nO,*)
            CALL writeTime(nO,'! Total run time:',runTime)
            WRITE(nO,*)
            WRITE(nO,*) ' !***********************************!'
            WRITE(nO,*) ' !---- THANK YOU FOR USING MAGIC ----!'
            WRITE(nO,*) ' !---- ALWAYS HAPPY TO PLEASE YOU ---!'
            WRITE(nO,*) ' !--------  CALL BACK AGAIN ---------!'
            WRITE(nO,*) ' !- GET YOUR NEXT DYNAMO WITH MAGIC -!'
            WRITE(nO,*) ' !***********************************!'
            WRITE(nO,*) '                                  JW  '
            WRITE(nO,*)
        END DO
        IF ( l_save_out ) CLOSE(n_log_file)
    END IF

!--- Closing the output files:
    CALL closeFiles


    PERFOFF
    PERFOUT('main')
    !LIKWID_OFF('main')
    LIKWID_CLOSE
!-- EVERYTHING DONE ! THE END !
#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
    END PROGRAM

!--------------------------------------------------------------------
!--- FINISHED !
!--------------------------------------------------------------------
