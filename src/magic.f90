#include "perflib_preproc.cpp"
program magic
!  A dynamic dynamo model driven by thermal convection
!  in a rotating spherical fluid shell.
!  This version can solve for both Boussinesq and anelastic fluids and
!  non-dimensional variables are used throughout the whole code.
!

!  In the Boussinesq limit, the following set of equations is solved:
!
!  ``
!  E {dv/dt + v.grad(v)} = -grad(p) - 2e_z x v
!      +1/Pm rot(B) x B + RaE/Pr g/g_o T
!
!  dB/dt = rot(v x B) + 1/Pm Lapl(B)
!
!  dT/dt + v.grad(T) = 1/Pr Lapl(T) + epsc0
!
!  div(v)=0          div(B)=0
!  ``
!
!  List of symbols:
!
!  v: velocity          p: pressure        B: magnetic field
!  g: gravity           g_o: reference value at outer radius
!  T: temperature       epsc0: rate of internal heating
!  e_z: unit vector parallel to the rotation axis
!  d/dt: partial time derivative  Lapl: Laplace operator
!
!  Scaling properties:
!
!  nu: kinematic viscosity         d: shell width
!  omega: angular frequency        alpha: thermal expansion coeff
!  delta_T: temperature contrast   kappa: thermal diffusivity
!  eta: magnetic diffusivity       rho: density
!  mu_o: magnetic permeability
!
!  Scaling:
!
!  Length:   d              time:      d^2/nu
!  Velocity: nu/d           pressure:  rho*nu*omega
!  Temperature: delta_T     mag.field: sqrt(rho*mu_o*eta*omega)
!
!  Non-dimensional numbers:
!
!  E: Ekman number     E= nu*d^2/omega
!  Ra: Rayleigh number Ra = alpha*g_o*delta_T*d^3/(kappa*nu)
!  Pr: Prandtl number  Pr = nu/kappa
!  Pm: Magnetic Prandtl number    Pm=nu/eta
!
!  Numerical simulations via a nonlinear, multimode,
!  initial-boundary value problem.
!
!    * entropy boundary conditions (tops and bots on input)
!      if ktops = 1, entropy specified on outer boundary
!      if ktops = 2, radial heat flux specified on outer boundary
!      if kbots = 1, entropy specified on inner boundary
!      if kbots = 2, radial heat flux specified on inner boundary
!      for example: ktops=1,
!            the spherically-symmetric temperature
!            on the outer boundary relative to the reference state
!
!    * velocity boundary condtions
!      if ktopv = 1, stress-free outer boundary
!      if ktopv = 2, non-slip outer boundary
!      if kbotv = 1, stress-free inner boundary
!      if kbotv = 2, non-slip inner boundary
!
!    * magnetic boundary condtions
!      if ktopb = 1, insulating outer boundary (mantle)
!      if kbotb = 1, insulating inner boundary (core)
!      if ktopb = 2, not implemented
!      if kbotb = 2, perfectly conducting inner boundary
!      if ktopb = 3, finitely conducting mantle
!      if kbotb = 3, finitely conducting inner core
!      if ktopb = 4, pseudo vacuum outer boundary (B=Br)
!      if kbotb = 4, pseudo vacuum inner boundary (B=Br)
!
!    * magneto-convection
!      amp_b1 = max amplitude of imposed magnetic field
!      if imagcon  ==  1, imposed toroidal field via inner bc on J(l=2,m=0)
!      if imagcon  == 10, imposed tor. field on both icb and cmb J(l=2,m=0)
!      if imagcon  == 11, imposed tor. field on both icb and cmb J(l=2,m=0)
!                         opposite sign
!      if imagcon  == 12, imposed tor. field on both icb and cmb J(l=1,m=0)
!      if imagcon  <  0, imposed poloidal field via inner bc on B(l=1,m=0)
!

   use truncation
   use precision_mod
   use physical_parameters
   use radial_functions
   use num_param
   use torsional_oscillations
   use init_fields
   use Grenoble
   use blocking
   use horizontal_data
   use logic
   use matrices
   use fields
   use fieldsLast
   use constants, only: codeVersion
   use movie_data, only: initialize_movie_data, finalize_movie_data
   use RMS, only: initialize_RMS
   use dtB_mod, only: initialize_dtB_mod
   use radial_data, only: initialize_radial_data
   use radialLoop
   use lmLoop_data, only: initialize_LMLoop_data
   use LMLoop_mod,only: initialize_LMLoop
   use preCalculations
   use start_fields, only: getStartFields
   use kinetic_energy
   use magnetic_energy
   use fields_average_mod
   use Egeos_mod
   use spectra, only: initialize_spectra
   use output_data
   use output_mod, only: initialize_output
   use outPV3, only:initialize_outPV3
   use outTO_mod,only: initialize_outTO_mod
   use parallel_mod
   use Namelists
   use step_time_mod, only: initialize_step_time, step_time
   use timing, only: writeTime,wallTime
   use communications, only:initialize_communications
   use power, only: initialize_output_power
   use outPar_mod, only: initialize_outPar_mod
   !use rIterThetaBlocking_mod,ONLY: initialize_rIterThetaBlocking
#ifdef WITH_LIKWID
#  include "likwid_f90.h"
#endif

#ifdef WITH_SHTNS
   use shtns
#endif
   implicit none

   !-- Local variables:
   integer :: n_time_step
   integer :: n_time_step_start   ! storing initial time step no
   integer :: n                   ! counter
   integer :: nO                  ! output unit
   real(cp) :: time
   real(cp) :: dt
   real(cp) :: dtNew

   integer :: n_stop_signal=0     ! signal returned from step_time

   ! MPI specific variables
#ifdef WITHOMP
   integer :: required_level,provided_level
#endif

#ifdef WITH_MPI
#ifdef WITHOMP
   required_level=MPI_THREAD_MULTIPLE
   call mpi_init_thread(required_level,provided_level,ierr)
   if (provided_level < required_level) then
      print*,"We need at least thread level ",required_level, &
      &      ", but have ",provided_level
      stop
   end if
#else
   call mpi_init(ierr)
#endif
#endif

   PERFINIT
   PERFON('main')
   LIKWID_INIT
   !LIKWID_ON('main')
   call parallel

   !--- Read starting time
   if ( rank == 0 ) then
      !call get_resetTime(resetTime)
      call wallTime(runTimeStart)
      write(*,*)
      write(*,*) '!--- Program MAGIC ', trim(codeVersion), ' ---!'
      call writeTime(6,'! Started at:',runTimeStart)
   end if

   !--- Read input parameters:
   call readNamelists  ! includes sent to other procs !

   call initialize_output

   !--- Check parameters and write info to SDTOUT
   call checkTruncation

   !--- Open output files:
   call openFiles

   if ( rank == 0 ) then
      if ( l_save_out ) then
         open(n_log_file, file=log_file, status='unknown', position='append')
      end if
      write(n_log_file,*) '!------------------------------------------------------!'
      write(n_log_file,*) '!                                                      !'
      write(n_log_file,*) '!       Program MAGIC ', trim(codeVersion),  &
           &              '                              !'
      write(n_log_file,*) '!                                                      !'
      write(n_log_file,*) '!                                                      !'
      write(n_log_file,*) '!------------------------------------------------------!'
      write(n_log_file,*)
      if ( l_save_out ) close(n_log_file)
   end if

   call initialize_blocking
   call initialize_radial_data
   call initialize_radial_functions
   call initialize_radialLoop
   !call initialize_rIterThetaBlocking
   call initialize_LMLoop_data
   call initialize_LMLoop

   call initialize_num_param
   if ( l_TO ) call initialize_TO
   if ( l_TO ) call initialize_outTO_mod
   call initialize_init_fields
   call initialize_Grenoble
   call initialize_horizontal_data
   call initialize_matrices
   call initialize_fields
   call initialize_fieldsLast
   if ( ldtBmem == 1 ) call initialize_dtB_mod
   call initialize_kinetic_energy
   call initialize_magnetic_energy
   call initialize_fields_average_mod
   if ( l_par ) call initialize_Egeos_mod
   call initialize_spectra
   if ( l_PV ) call initialize_outPV3
   call initialize_step_time
   call initialize_communications
   call initialize_outPar_mod
   if ( l_power ) call initialize_output_power

   !--- Do pre-calculations:
   call preCalc

   if ( l_movie ) call initialize_movie_data !Needs to be called after preCalc to get correct coordinate values

   if ( l_RMS ) call initialize_RMS

   if ( rank == 0 ) then
      if ( l_save_out ) then
         open(n_log_file, file=log_file, status='unknown', position='append')
      end if
      call writeNamelists(6)
      call writeNamelists(n_log_file)
      if ( l_save_out ) close(n_log_file)
   end if

   !--- Now read start-file or initialize fields:
   call getStartFields(time,dt,dtNew,n_time_step)

   !--- Second pre-calculation:
   call preCalcTimes(time,n_time_step)

#ifdef WITH_SHTNS
   call init_shtns()
#endif
   !--- Write info to STDOUT and log-file:
   if ( rank == 0 ) then
      call writeInfo(6)
      call writeInfo(n_log_file)
   end if

   !--- AND NOW FOR THE TIME INTEGRATION:

   !--- Write starting time to SDTOUT and logfile:
   if ( rank == 0 ) then
      if ( l_save_out ) then
         open(n_log_file, file=log_file, status='unknown',  position='append')
      end if
      do n=1,2
         if ( n == 1 ) nO=6
         if ( n == 2 ) nO=n_log_file
         write(nO,'(/,'' ! STARTING TIME INTEGRATION AT:'')')
         write(nO,'(''   start_time ='',1p,ES18.10)') tScale*time
         write(nO,'(''   step no    ='',i10)') n_time_step
         write(nO,'(''   start dt   ='',1p,ES16.4)') dt
         write(nO,'(''   start dtNew='',1p,ES16.4)') dtNew
      end do
      if ( l_save_out ) close(n_log_file)
   end if
   timeStart        =time
   n_time_step_start=n_time_step

   !--- Call time-integration routine:
   PERFON('steptime')
   call step_time(time,dt,dtNew,n_time_step)
   PERFOFF
   !--- Write stop time to SDTOUR and logfile:
   if ( rank == 0 ) then
      if ( l_save_out ) then
         open(n_log_file, file=log_file, status='unknown', position='append')
      end if

      do n=1,2
         if ( n == 1 ) nO=6
         if ( n == 2 ) nO=n_log_file
         write(nO,'(/,'' ! STOPPING TIME INTEGRATION AT:'')')
         write(nO,'(''   stop time ='',1p,ES18.10)') tScale*time
         write(nO,'(''   stop step ='',i10)') n_time_step+n_time_step_start
         write(nO,'(''   steps gone='',i10)') (n_time_step-1)
         write(nO,*)
         if ( n_stop_signal > 0 ) then
            write(nO,*) '!!! MAGIC terminated by STOP signal !!!'
         else
            write(nO,*) '!!! regular end of program MAGIC !!!'
         end if
         write(nO,*)
         !write(nO,'('' max. thread number in  R-loop='',i3)') &
         !      nThreadsRmax
         !write(nO,'('' max. thread number in LM-loop='',i3)') &
         !      nThreadsLMmax
         write(nO,*)
         call writeTime(nO,'! Total run time:',runTime)
         write(nO,*)
         write(nO,*) ' !***********************************!'
         write(nO,*) ' !---- THANK YOU FOR USING MAGIC ----!'
         write(nO,*) ' !---- ALWAYS HAPPY TO PLEASE YOU ---!'
         write(nO,*) ' !--------  call BACK AGAIN ---------!'
         write(nO,*) ' !- GET YOUR NEXT DYNAMO WITH MAGIC -!'
         write(nO,*) ' !***********************************!'
         write(nO,*) '                                  JW  '
         write(nO,*)
      end do
      if ( l_save_out ) close(n_log_file)
   end if

   !--- Closing the movie files (if any)
   call finalize_movie_data
   !--- Closing the output files:
   call closeFiles
   


   PERFOFF
   PERFOUT('main')
   !LIKWID_OFF('main')
   LIKWID_CLOSE
!-- EVERYTHING DONE ! THE END !
#ifdef WITH_MPI
   call mpi_finalize(ierr)
#endif
end program magic
