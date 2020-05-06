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
   use iso_fortran_env, only: output_unit
   use courant_mod, only: initialize_courant, finalize_courant
   use radial_der, only: initialize_der_arrays, finalize_der_arrays
   use radial_functions, only: initialize_radial_functions, &
       &                       finalize_radial_functions
   use num_param
   use torsional_oscillations
   use init_fields
   use special, only: initialize_Grenoble, finalize_Grenoble
   use blocking, only: initialize_blocking, finalize_blocking
   use timing, only: timer_type
   use horizontal_data
   use logic
   use fields
   use fieldsLast
   use constants, only: codeVersion
   use movie_data, only: initialize_movie_data, finalize_movie_data
   use RMS, only: initialize_RMS, finalize_RMS
   use dtB_mod, only: initialize_dtB_mod, finalize_dtB_mod
   use radialLoop, only: initialize_radialLoop, finalize_radialLoop
   use LMLoop_mod,only: initialize_LMLoop, finalize_LMLoop
   use preCalculations
   use start_fields, only: getStartFields
   use kinetic_energy
   use magnetic_energy
   use fields_average_mod
   use geos_mod, only: initialize_geos_mod, finalize_geos_mod
   use spectra, only: initialize_spectra, finalize_spectra
   use output_data, only: tag, log_file, n_log_file
   use output_mod, only: initialize_output, finalize_output
   use outTO_mod,only: initialize_outTO_mod, finalize_outTO_mod
   use parallel_mod
   use Namelists
   use step_time_mod, only: initialize_step_time, step_time
   use communications, only:initialize_communications, finalize_communications
   use power, only: initialize_output_power, finalize_output_power
   use outPar_mod, only: initialize_outPar_mod, finalize_outPar_mod
   use outMisc_mod, only: initialize_outMisc_mod, finalize_outMisc_mod
   use outRot, only: initialize_outRot, finalize_outRot
   use mem_alloc
   use useful, only: abortRun
   use probe_mod, only: initialize_probes, finalize_probes
   use time_schemes, only: type_tscheme
   use LMmapping, only: initialize_mapping, finalize_mapping
   use fft, only: initialize_fft_phi, finalize_fft_phi

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
   integer(lip) :: local_bytes_used
   integer :: values(8)
   character(len=72) :: date
   real(cp) :: time
   class(type_tscheme), pointer :: tscheme
   type(timer_type) :: run_time, run_time_start

   integer :: n_stop_signal=0     ! signal returned from step_time


   ! MPI specific variables
#ifdef WITHOMP
   integer :: required_level,provided_level
   character(len=100) :: message
   character(len=14) :: str, str_1
#endif

#ifdef WITH_MPI
#ifdef WITHOMP
   required_level=MPI_THREAD_FUNNELED
   call MPI_Init_Thread(required_level,provided_level,ierr)
   if (provided_level < required_level) then
      print*,"We need at least thread level ",required_level, &
      &      ", but have ",provided_level
      write(str,*) required_level
      write(str_1,*) provided_level
      message = 'We need at least thread level'//trim(adjustl(str))//&
      &         ', but have'//trim(adjustl(str_1))
      call abortRun(message)
   end if
#else
   call MPI_Init(ierr)
#endif
#endif

   PERFINIT
   PERFON('main')
   LIKWID_INIT
   !LIKWID_ON('main')
   call initialize_parallel()

   call run_time%initialize()
   call run_time%start_count()
   call run_time_start%initialize()
   call run_time_start%start_count()

   !--- Read starting time
   if ( l_master_rank ) then
      write(output_unit,*)
      write(output_unit,*) '!--- Program MagIC ', trim(codeVersion), ' ---!'
      call date_and_time(values=values)
#if defined(GIT_VERSION)
      write(output_unit, '(A,A)') ' !  Git version:  ', GIT_VERSION
#else
      write(output_unit, '(A)') ' !  Git version: unknown'
#endif
#if defined(BUILD_DATE)
      write(output_unit, '(A,A)') ' !  Build date:  ', BUILD_DATE
#else
      write(output_unit, '(A)') ' !  Build date: unknown'
#endif
      write(date, '(i4,''/'',i0.2,''/'',i0.2,'' '', i0.2,'':'',i0.2,'':'',i0.2)') &
      &     values(1), values(2), values(3), values(5), values(6), values(7)
      write(output_unit, '(A,A)') ' !  Start date:  ', date
   end if

   !--- Read input parameters:
   call readNamelists(tscheme)  ! includes sent to other procs !

   call initialize_output()

   !--- Check parameters and write info to SDTOUT
   call checkTruncation()

   log_file='log.'//tag

   if ( l_master_rank ) then
      open(newunit=n_log_file, file=log_file, status='new')

      write(n_log_file,*) '!      __  __             _____ _____   _____ ___        '
      write(n_log_file,*) '!     |  \/  |           |_   _/ ____| | ____/ _ \       '
      write(n_log_file,*) '!     | \  / | __ _  __ _  | || |      | |__| (_) |      '
      write(n_log_file,*) '!     | |\/| |/ _` |/ _` | | || |      |___ \\__, |      '
      write(n_log_file,*) '!     | |  | | (_| | (_| |_| || |____   ___) | / /       '
      write(n_log_file,*) '!     |_|  |_|\__,_|\__, |_____\_____| |____(_)_/        '
      write(n_log_file,*) '!                    __/ |                               '
      write(n_log_file,*) '!                   |___/                                '
      write(n_log_file,*) '!                                                        '
      write(n_log_file,*) '!                                                        '
      write(n_log_file,*) '!                          /^\     .                     '
      write(n_log_file,*) '!                     /\   "V"                           '
      write(n_log_file,*) '!                    /__\   I      O  o                  '
      write(n_log_file,*) '!                   //..\\  I     .                      '
      write(n_log_file,*) '!                   \].`[/  I                            '
      write(n_log_file,*) '!                   /l\/j\  (]    .  O                   '
      write(n_log_file,*) '!                  /. ~~ ,\/I          .                 '
      write(n_log_file,*) '!                  \\L__j^\/I       o                    '
      write(n_log_file,*) '!                   \/--v}  I     o   .                  '
      write(n_log_file,*) '!                   |    |  I   _________                '
      write(n_log_file,*) '!                   |    |  I c(`       ")o              '
      write(n_log_file,*) '!                   |    l  I   \.     ,/                '
      write(n_log_file,*) '!                 _/j  L l\_!  _//^---^\\_               '
      write(n_log_file,*) '!              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            '
      write(n_log_file,*) '!                                                        '
      write(n_log_file,*) '!                                                        '

#if defined(GIT_VERSION)
      write(n_log_file, '(A,A)') ' ! Git version:  ', GIT_VERSION
#else
      write(n_log_file, '(A)') ' ! Git version: unknown'
#endif
#if defined(BUILD_DATE)
      write(n_log_file, '(A,A)') ' ! Build date:  ', BUILD_DATE
#else
      write(n_log_file, '(A)') ' ! Build date: unknown'
#endif
      write(n_log_file, '(A,A)') ' ! Start date:  ', date

      if ( l_save_out ) close(n_log_file)
   end if

   call initialize_mpi_map
   call initialize_distributed_geometry
   
   call initialize_memory_counter()

   !-- Blocking/radial/horizontal
   call initialize_mapping()
   call initialize_blocking()
   local_bytes_used=bytes_allocated
   call initialize_radial_data(n_r_max)
   call initialize_radial_functions()
   call initialize_horizontal_data()
   local_bytes_used=bytes_allocated-local_bytes_used
   call memWrite('radial/horizontal', local_bytes_used)
   call initialize_fft_phi()

   !-- Initialize time scheme
   call tscheme%initialize(time_scheme, courfac, intfac, alffac)

   !-- Radial/LM Loop
   call initialize_radialLoop()
   call initialize_LMLoop(tscheme)

   call initialize_num_param()
   call initialize_init_fields()
   call initialize_Grenoble()

   local_bytes_used=bytes_allocated
   call initialize_fields()
   call initialize_fieldsLast(tscheme%nold, tscheme%nexp, tscheme%nimp)
   local_bytes_used=bytes_allocated-local_bytes_used
   call memWrite('fields/fieldsLast', local_bytes_used)

   call initialize_step_time()
   call initialize_communications()

   call initialize_der_arrays(n_r_max,n_mlo_loc)

   !-- Array allocation for I/O
   local_bytes_used=bytes_allocated
   call initialize_kinetic_energy()
   if ( l_mag ) call initialize_magnetic_energy()
   call initialize_spectra()
   call initialize_outPar_mod()
   call initialize_outMisc_mod()
   call initialize_outRot()
   if ( l_power ) call initialize_output_power()
   call initialize_fields_average_mod()
   if ( l_TO ) call initialize_TO()



#ifdef WITH_SHTNS
   call init_shtns()
! ! ! ! ! ! ! ! ! ! ! ! !    call test_shtns()
#endif

   if ( l_master_rank ) then
      call tscheme%print_info(n_log_file)
   end if

   !--- Do pre-calculations:
   call preCalc(tscheme)

   if ( l_par ) call initialize_geos_mod(l_par) ! Needs to be called after preCalc, r_icb needed
   if ( l_TO ) call initialize_outTO_mod() ! Needs to be called after preCalc, r_icb needed
   if ( l_movie ) call initialize_movie_data() !Needs to be called after preCalc to get correct coordinate values
   if ( ldtBmem == 1 ) call initialize_dtB_mod() ! Needs to be called after movie to make sure l_dtBmovie has been set
   if (l_probe) call initialize_probes()       !Needs to be called after preCalc to get correct coordinate values
   if ( l_RMS ) call initialize_RMS()
   local_bytes_used=bytes_allocated-local_bytes_used

   !local_bytes_used=bytes_allocated-local_bytes_used
   call memWrite('Total I/O', local_bytes_used)

   if (coord_r == 0) print*, '-----> coord_r 0 has', bytes_allocated, ' B allocated'


   call finalize_memory_counter()

   if ( l_master_rank ) then
      if ( l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if
      call writeNamelists(output_unit)
      call writeNamelists(n_log_file)
      if ( l_save_out ) close(n_log_file)
   end if

   !--- Now read start-file or initialize fields:
   call getStartFields(time, tscheme, n_time_step)

   !-- Open time step file
   call initialize_courant(time, tscheme%dt(1), tag)

   !--- Second pre-calculation:
   call preCalcTimes(time,n_time_step)

   !--- Write info to STDOUT and log-file:
   if ( l_master_rank ) then
      if ( l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if
      call writeInfo(output_unit)
      call writeInfo(n_log_file)
      if ( l_save_out ) close(n_log_file)
   end if

   !--- AND NOW FOR THE TIME INTEGRATION:

   !--- Write starting time to SDTOUT and logfile:
   if ( l_master_rank ) then
      if ( l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown',  &
         &    position='append')
      end if
      do n=1,2
         if ( n == 1 ) nO=6
         if ( n == 2 ) nO=n_log_file
         write(nO,'(/,'' ! STARTING TIME INTEGRATION AT:'')')
         write(nO,'(''   start_time ='',1p,ES18.10)') tScale*time
         write(nO,'(''   step no    ='',i10)') n_time_step
         write(nO,'(''   start dt   ='',1p,ES16.4)') tscheme%dt(1)
      end do
      if ( l_save_out ) close(n_log_file)
   end if
   timeStart        =time
   n_time_step_start=n_time_step

   !-- Stop measuring initiatisation of MagIC
   call run_time_start%stop_count()

   !--- Call time-integration routine:
   call step_time(time, tscheme, n_time_step, run_time_start)

   !-- Stop counting time and print
   call run_time%stop_count()
   if ( l_save_out ) then
      open(newunit=n_log_file, file=log_file, status='unknown', &
      &    position='append')
   end if
   call run_time%finalize('! Total run time:', n_log_file)
   if ( l_save_out ) close(n_log_file)

   !--- Write stop time to SDTOUR and logfile:
   if ( l_master_rank ) then
      if ( l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if

      do n=1,2
         if ( n == 1 ) nO=output_unit
         if ( n == 2 ) nO=n_log_file
         write(nO,'(/,'' ! STOPPING TIME INTEGRATION AT:'')')
         write(nO,'(''   stop time ='',1p,ES18.10)') tScale*time
         write(nO,'(''   stop step ='',i10)') n_time_step+n_time_step_start
         write(nO,'(''   steps gone='',i10)') (n_time_step-1)
         write(nO,*)
         if ( n_stop_signal > 0 ) then
            write(nO,*) '!!! MagIC terminated by STOP signal !!!'
         else
            write(nO,*) '!!! regular end of program MagIC !!!'
         end if
         write(nO,*)

         write(nO,*)
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
   if ( l_RMS ) call finalize_RMS
   if ( l_TO ) call finalize_outTO_mod
   if ( l_TO ) call finalize_TO
   if ( l_par ) call finalize_geos_mod(l_par)
   if ( ldtBmem == 1 ) call finalize_dtB_mod
   call finalize_fields_average_mod()
   if ( l_power ) call finalize_output_power()
   call finalize_outRot()
   call finalize_outMisc_mod()
   call finalize_outPar_mod()
   call finalize_spectra()
   if ( l_mag ) call finalize_magnetic_energy()
   call finalize_kinetic_energy()
   if ( l_probe ) call finalize_probes()

   call finalize_courant()

   call finalize_communications()
   call finalize_fieldsLast()
   call finalize_fields()

   call finalize_Grenoble()
   call finalize_init_fields()
   call finalize_num_param()
   call finalize_LMLoop(tscheme)
   call finalize_radialLoop()

   call finalize_der_arrays()

   call finalize_horizontal_data()
   call finalize_radial_functions()
   call finalize_blocking()
   call finalize_mapping()
   call finalize_fft_phi()
   call finalize_radial_data()

   call tscheme%finalize()
   call finalize_output()

   if ( l_master_rank .and. .not. l_save_out )  close(n_log_file)

   PERFOFF
   PERFOUT('main')
   !LIKWID_OFF('main')
   LIKWID_CLOSE
!-- EVERYTHING DONE ! THE END !
#ifdef WITH_MPI

#ifdef WITHPERF
   do n=0,n_ranks-1
      if (rank == n) then
         print *, "PerfOut for rank: ", rank
         print *, "======================================================="
         PERFOUT('main')
      end if
      flush(6)
      call mpi_barrier(mpi_comm_world,ierr)
   end do
#endif
   call mpi_finalize(ierr)
#endif

end program magic
