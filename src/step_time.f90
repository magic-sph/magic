#include "perflib_preproc.cpp"
module step_time_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use iso_fortran_env, only: output_unit
   use fields
   use fieldsLast
   use parallel_mod
   use precision_mod
   use constants, only: zero, one, half
   use truncation, only: n_r_max, l_max, l_maxMag, lm_max, lmP_max, fd_order, &
       &                 fd_order_bound
   use num_param, only: n_time_steps, run_time_limit, tEnd, dtMax, &
       &                dtMin, tScale, dct_counter, nl_counter,    &
       &                solve_counter, lm2phy_counter, td_counter, &
       &                phy2lm_counter, nl_counter, f_exp_counter
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag, &
       &                  n_r_icb, n_r_cmb
   use radial_der, only: get_dr_Rloc, get_ddr_Rloc, exch_ghosts, bulk_to_ghost
   use radial_functions, only: rscheme_oc
   use logic, only: l_mag, l_mag_LF, l_dtB, l_RMS, l_hel, l_TO,        &
       &            l_TOmovie, l_r_field, l_cmb_field, l_HTmovie,      &
       &            l_DTrMagSpec, lVerbose, l_b_nl_icb, l_par,         &
       &            l_b_nl_cmb, l_FluxProfs, l_ViscBcCalc, l_perpPar,  &
       &            l_HT, l_dtBmovie, l_heat, l_conv, l_movie,         &
       &            l_runTimeLimit, l_save_out, l_bridge_step,         &
       &            l_dt_cmb_field, l_chemical_conv, l_mag_kin,        &
       &            l_power, l_double_curl, l_PressGraph, l_probe,     &
       &            l_AB1, l_finite_diff, l_cond_ic, l_single_matrix,  &
       &            l_packed_transp, l_parallel_solve, l_mag_par_solve
   use init_fields, only: omega_ic1, omega_ma1
   use movie_data, only: t_movieS
   use radialLoop, only: radialLoopG
   use LMLoop_mod, only: LMLoop, finish_explicit_assembly, assemble_stage, &
       &                 finish_explicit_assembly_Rdist, LMLoop_Rdist,     &
       &                 assemble_stage_Rdist
   use signals_mod, only: initialize_signals, check_signals
   use graphOut_mod, only: open_graph_file, close_graph_file
   use output_data, only: tag, n_graph_step, n_graphs, n_t_graph, t_graph, &
       &                  n_spec_step, n_specs, n_t_spec, t_spec,          &
       &                  n_movie_step, n_movie_frames, n_t_movie, t_movie,&
       &                  n_TOmovie_step, n_TOmovie_frames, n_t_TOmovie,   &
       &                  t_TOmovie, n_pot_step, n_pots, n_t_pot, t_pot,   &
       &                  n_rst_step, n_rsts, n_t_rst, t_rst, n_stores,    &
       &                  n_log_step, n_logs, n_t_log, t_log, n_cmb_step,  &
       &                  n_cmbs, n_t_cmb, t_cmb, n_r_field_step,          &
       &                  n_r_fields, n_t_r_field, t_r_field, n_TO_step,   &
       &                  n_TOs, n_t_TO, t_TO, n_probe_step, n_probe_out,  &
       &                  n_t_probe, t_probe, log_file, n_log_file,        &
       &                  n_time_hits
   use updateB_mod, only: get_mag_rhs_imp, get_mag_ic_rhs_imp, b_ghost, aj_ghost, &
       &                  get_mag_rhs_imp_ghost, fill_ghosts_B
   use updateWP_mod, only: get_pol_rhs_imp, get_pol_rhs_imp_ghost, w_ghost, &
       &                   fill_ghosts_W
   use updateWPS_mod, only: get_single_rhs_imp
   use updateS_mod, only: get_entropy_rhs_imp, get_entropy_rhs_imp_ghost, s_ghost, &
       &                  fill_ghosts_S
   use updateXI_mod, only: get_comp_rhs_imp, get_comp_rhs_imp_ghost, xi_ghost, &
       &                   fill_ghosts_Xi
   use updateZ_mod, only: get_tor_rhs_imp, get_rot_rates, get_tor_rhs_imp_ghost, &
       &                  z_ghost, fill_ghosts_Z
   use output_mod, only: output
   use time_schemes, only: type_tscheme
   use useful, only: l_correct_step, logWrite
   use communications, only: lo2r_field, lo2r_flow, scatter_from_rank0_to_lo, &
       &                     lo2r_xi,  r2lo_flow, r2lo_s, r2lo_xi,r2lo_field, &
       &                     lo2r_s, lo2r_press, lo2r_one, r2lo_one
   use courant_mod, only: dt_courant
   use nonlinear_bcs, only: get_b_nl_bcs
   use timing ! Everything is needed
   use probe_mod

   implicit none

   private

   public :: initialize_step_time, step_time

contains

   subroutine initialize_step_time()

      call initialize_signals()

   end subroutine initialize_step_time
!-------------------------------------------------------------------------------
   subroutine step_time(time, tscheme, n_time_step, run_time_start)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !-- Input from initialization:
      !   time and n_time_step updated and returned to magic.f
      real(cp),            intent(inout) :: time
      class(type_tscheme), intent(inout) :: tscheme
      integer,             intent(inout) :: n_time_step
      type(timer_type),    intent(in) :: run_time_start

      !--- Local variables:

      !--- Logicals controlling output/calculation:
      logical :: l_graph          !
      logical :: l_spectrum
      logical :: l_store          ! Store output in restart file
      logical :: l_new_rst_file   ! Use new rst file
      logical :: l_log            ! Log output
      logical :: l_stop_time      ! Stop time stepping
      logical :: l_frame          ! Movie frame output
      logical :: lTOframe         ! TO movie frame output
      logical :: l_cmb            ! Store set of b at CMB
      logical :: l_r              ! Store coeff at various depths
      logical :: lHelCalc         ! Calculate helicity for output
      logical :: lPowerCalc       ! Calculate viscous heating in the physical space
      logical :: lviscBcCalc      ! Calculate horizontal velocity and (grad T)**2
      logical :: lFluxProfCalc    ! Calculate radial flux components
      logical :: lPerpParCalc     ! Calculate perpendicular and parallel Ekin
      logical :: lGeosCalc        ! Calculate geos.TAG outputs
      logical :: lTOCalc          ! Calculate TO stuff
      logical :: lTONext,lTONext2 ! TO stuff for next steps
      logical :: lTOframeNext,lTOframeNext2
      logical :: l_logNext, l_pot
      logical :: lRmsCalc,lRmsNext, l_pure, l_mat_time
      logical :: lPressCalc,lPressNext
      logical :: lMat, lMatNext   ! update matrices
      logical :: l_probe_out      ! Sensor output

      !-- Timers:
      type(timer_type) :: rLoop_counter, lmLoop_counter, comm_counter
      type(timer_type) :: mat_counter, tot_counter, io_counter, pure_counter
      real(cp) :: run_time_passed, dt_new

      !--- Counter:
      integer :: n_frame          ! No. of movie frames
      integer :: n_cmb_sets       ! No. of stored sets of b at CMB
      integer :: n_stage

      !--- Stuff needed to construct output files:
      character(len=255) :: message

      !--- Courant criteria/diagnosis:
      real(cp) :: dtr, dth, tTot
      !-- Saves values for time step
      real(cp) :: dtrkc_Rloc(nRstart:nRstop), dthkc_Rloc(nRstart:nRstop)

      !--- Explicit part of time stepping partly calculated in radialLoopG and
      !    passed to LMLoop where the time step is preformed.
      !    Note that the respective arrays for the changes in inner-core
      !    magnetic field are calculated in updateB and are only
      !    needed there.

      !--- Lorentz torques:
      real(cp) :: lorentz_torque_ma,lorentz_torque_ic

      !-- Arrays for outMisc.f90 and outPar.f90
      real(cp) :: HelASr_Rloc(2,nRstart:nRstop),Hel2ASr_Rloc(2,nRstart:nRstop)
      real(cp) :: HelnaASr_Rloc(2,nRstart:nRstop),Helna2ASr_Rloc(2,nRstart:nRstop)
      real(cp) :: HelEAASr_Rloc(nRstart:nRstop),viscASr_Rloc(nRstart:nRstop)
      real(cp) :: uhASr_Rloc(nRstart:nRstop),duhASr_Rloc(nRstart:nRstop)
      real(cp) :: gradsASr_Rloc(nRstart:nRstop),fconvASr_Rloc(nRstart:nRstop)
      real(cp) :: fkinASr_Rloc(nRstart:nRstop),fviscASr_Rloc(nRstart:nRstop)
      real(cp) :: fpoynASr_Rloc(nRstartMag:nRstopMag),fresASr_Rloc(nRstartMag:nRstopMag)
      real(cp) :: EperpASr_Rloc(nRstart:nRstop),EparASr_Rloc(nRstart:nRstop)
      real(cp) :: EperpaxiASr_Rloc(nRstart:nRstop),EparaxiASr_Rloc(nRstart:nRstop)

      !--- Nonlinear magnetic boundary conditions needed in s_updateB.f :
      complex(cp) :: br_vt_lm_cmb(lmP_max)    ! product br*vt at CMB
      complex(cp) :: br_vp_lm_cmb(lmP_max)    ! product br*vp at CMB
      complex(cp) :: br_vt_lm_icb(lmP_max)    ! product br*vt at ICB
      complex(cp) :: br_vp_lm_icb(lmP_max)    ! product br*vp at ICB
      complex(cp) :: b_nl_cmb(lm_max)         ! nonlinear bc for b at CMB
      complex(cp) :: aj_nl_cmb(lm_max)        ! nonlinear bc for aj at CMB
      complex(cp) :: aj_nl_icb(lm_max)        ! nonlinear bc for dr aj at ICB

      !--- Various stuff for time control:
      real(cp) :: timeLast, timeStage, dtLast
      integer :: n_time_steps_go
      logical :: l_finish_exp_early, l_last_RMS
      logical :: l_new_dt         ! causes call of matbuild !
      integer :: nPercent         ! percentage of finished time stepping
      real(cp) :: tenth_n_time_steps

      !-- Interupt procedure:
      integer :: signals(5)
      integer :: n_stop_signal     ! =1 causes run to stop
      integer :: n_graph_signal    ! =1 causes output of graphic file
      integer :: n_rst_signal      ! =1 causes output of rst file
      integer :: n_spec_signal     ! =1 causes output of a spec file
      integer :: n_pot_signal      ! =1 causes output for pot files

      if ( lVerbose ) write(output_unit,'(/,'' ! STARTING STEP_TIME !'')')

      run_time_passed=0.0_cp
      l_log       =.false.
      l_last_RMS  = l_RMS
      l_stop_time =.false.
      l_new_dt    =.true.   ! Invokes calculation of t-step matrices
      lMatNext    =.true.
      timeLast    =time
      timeStage   =time

      l_finish_exp_early = ( l_finite_diff .and. rscheme_oc%order==2 .and. &
      &                      rscheme_oc%order_boundary==2 )

      tenth_n_time_steps=real(n_time_steps,kind=cp)/10.0_cp
      nPercent=9

      !---- Set Lorentz torques to zero:
      lorentz_torque_ic=0.0_cp
      lorentz_torque_ma=0.0_cp

      !---- Counter for output files/sets:
      n_frame   =0    ! No. of movie frames
      n_cmb_sets=0    ! No. of store dt_b sets at CMB

      !---- Prepare signalling via file signal
      signals=0
      n_stop_signal =0     ! Stop signal returned to calling program
      n_graph_signal=0     ! Graph signal returned to calling program
      n_spec_signal=0      ! Spec signal
      n_rst_signal=0       ! Rst signal
      n_pot_signal=0       ! Potential file signal

      !-- STARTING THE TIME STEPPING LOOP:
      if ( rank == 0 ) then
         write(output_unit,*)
         write(output_unit,*) '! Starting time integration!'
      end if
      call comm_counter%initialize()
      call rLoop_counter%initialize()
      call lmLoop_counter%initialize()
      call mat_counter%initialize()
      call tot_counter%initialize()
      call pure_counter%initialize()
      call io_counter%initialize()

      !!!!! Time loop starts !!!!!!
      if ( n_time_steps == 1 ) then
         n_time_steps_go=1 ! Output only, for example G-file/movie etc.
      else if ( n_time_steps == 2 ) then
         n_time_steps_go=2 !
      else
         n_time_steps_go=n_time_steps+1  ! Last time step for output only !
      end if

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

      PERFON('tloop')
      !LIKWID_ON('tloop')
      outer: do n_time_step=1,n_time_steps_go

         if ( lVerbose ) then
            write(output_unit,*)
            write(output_unit,*) '! Starting time step ',n_time_step
         end if

         !-- Start time counters
         call mat_counter%start_count()
         call tot_counter%start_count()
         call pure_counter%start_count()
         l_pure=.false.
         l_mat_time=.false.

#ifdef WITH_MPI
         ! Broadcast omega_ic and omega_ma
         if ( l_parallel_solve ) then
            call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,n_procs-1, MPI_COMM_WORLD,ierr)
            call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,0, MPI_COMM_WORLD,ierr)
         else
            call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,rank_with_l1m0, &
                 &         MPI_COMM_WORLD,ierr)
            call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,rank_with_l1m0, &
                 &         MPI_COMM_WORLD,ierr)
         end if
#endif

         !----------------
         !-- This handling of the signal files is quite expensive
         !-- as the file can be read only on one rank and the result
         !-- must be distributed to all other ranks.
         !----------------
         call check_signals(run_time_passed, signals)
         n_stop_signal =signals(1)
         n_graph_signal=signals(2)
         n_rst_signal  =signals(3)
         n_spec_signal =signals(4)
         n_pot_signal  =signals(5)

         !--- Various reasons to stop the time integration:
         if ( l_runTimeLimit ) then
            tTot = tot_counter%tTot+run_time_start%tTot
#ifdef WITH_MPI
            call MPI_Allreduce(MPI_IN_PLACE, tTot, 1, MPI_DEF_REAL, MPI_MAX, &
                 &             MPI_COMM_WORLD, ierr)
#endif
            if ( tTot > run_time_limit ) then
               if ( .not. l_last_RMS ) then
                  write(message,'("! Run time limit exeeded !")')
                  call logWrite(message)
               end if
               l_stop_time=.true.
            end if

         end if
         !-- Handle an extra iteration in case RMS outputs are requested
         if ( (n_stop_signal > 0) .or. (l_RMS .and. (.not. l_last_RMS)) ) then
            l_stop_time=.true.   ! last time step !
         end if
         if ( n_time_step == n_time_steps_go ) then
            l_stop_time=.true.   ! last time step !
            l_last_RMS =.false.
         end if

         !--- Another reasons to stop the time integration:
         if ( time >= tEND .and. tEND /= 0.0_cp ) l_stop_time=.true.

         !-- Checking logic for output:
         l_graph= l_correct_step(n_time_step-1,time,timeLast,n_time_steps,       &
         &                       n_graph_step,n_graphs,n_t_graph,t_graph,0) .or. &
         &                  n_graph_signal == 1
         n_graph_signal=0   ! reset interrupt signal !
         l_spectrum=                                                             &
         &              l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
         &                n_spec_step,n_specs,n_t_spec,t_spec,0) .or.            &
         &                n_spec_signal == 1
         l_frame= l_movie .and. (                                                &
         &             l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &             n_movie_step,n_movie_frames,n_t_movie,t_movie,0) .or.     &
         &                   n_time_steps_go == 1 )
         if ( l_mag .or. l_mag_LF ) then
            l_dtB=( l_frame .and. l_dtBmovie ) .or. ( l_log .and. l_DTrMagSpec )
         end if

         lTOframe=l_TOmovie .and.                                                &
         &          l_correct_step(n_time_step-1,time,timeLast,n_time_steps,     &
         &          n_TOmovie_step,n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)

         l_probe_out=l_probe .and.                                               &
         &          l_correct_step(n_time_step-1,time,timeLast,n_time_steps,     &
         &          n_probe_step,n_probe_out,n_t_probe,t_probe,0)

         !-- Potential files
         l_pot= l_correct_step(n_time_step-1,time,timeLast,n_time_steps, &
         &                       n_pot_step,n_pots,n_t_pot,t_pot,0) .or. &
         &                  n_pot_signal == 1
         n_pot_signal=0   ! reset interrupt signal !

         l_new_rst_file=                                                         &
         &             l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &                            n_rst_step,n_rsts,n_t_rst,t_rst,0) .or.    &
         &             n_rst_signal == 1
         n_rst_signal=0
         l_store= l_new_rst_file .or.                                            &
         &             l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &                            0,n_stores,0,t_rst,0)

         l_log= l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &                            n_log_step,n_logs,n_t_log,t_log,0)
         l_cmb= l_cmb_field .and.                                                &
         &             l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &                            n_cmb_step,n_cmbs,n_t_cmb,t_cmb,0)
         l_r= l_r_field .and.                                                    &
         &             l_correct_step(n_time_step-1,time,timeLast,n_time_steps,  &
         &                            n_r_field_step,n_r_fields,n_t_r_field,     &
         &                            t_r_field,0)
         l_logNext=.false.
         if ( n_time_step+1 <= n_time_steps+1 )                                  &
         &             l_logNext=                                                &
         &             l_correct_step(n_time_step,time+tscheme%dt(1),timeLast,   &
         &                   n_time_steps,n_log_step,n_logs,n_t_log,t_log,0)
         lTOCalc= n_time_step > 2 .and. l_TO .and.                   &
         &               l_correct_step(n_time_step-1,time,timeLast, &
         &               n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
         lTOnext     =.false.
         lTOframeNext=.false.
         if ( n_time_step+1 <= n_time_steps+1 ) then
            lTONext= l_TO .and.                                            &
            &                l_correct_step(n_time_step,time+tscheme%dt(1),&
            &                timeLast,n_time_steps,n_TO_step,n_TOs,n_t_TO,t_TO,0)
            lTOframeNext= l_TOmovie .and.                                   &
            &                l_correct_step(n_time_step,time+tscheme%dt(1), &
            &                timeLast,n_time_steps,n_TOmovie_step,          &
            &                n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)
         end if
         lTONext      =lTOnext.or.lTOframeNext
         lTONext2     =.false.
         lTOframeNext2=.false.
         if ( n_time_step+2 <= n_time_steps+1 ) then
            lTONext2= l_TO .and.                                                 &
            &                l_correct_step(n_time_step+1,time+2*tscheme%dt(1),  &
            &                                timeLast,n_time_steps,n_TO_step,    &
            &                                            n_TOs,n_t_TO,t_TO,0)
            lTOframeNext2= l_TOmovie .and.                                      &
            &                l_correct_step(n_time_step+1,time+2*tscheme%dt(1), &
            &                             timeLast,n_time_steps,n_TOmovie_step, &
            &                       n_TOmovie_frames,n_t_TOmovie,t_TOmovie,0)
         end if
         lTONext2=lTOnext2.or.lTOframeNext2

         lRmsCalc=(l_RMS .and. l_log .and. (n_time_step > 1)) .or. &
         &        (l_RMS .and. l_stop_time)
         if ( l_mag .or. l_mag_LF ) l_dtB = l_dtB .or. lRmsCalc
         lRmsNext=l_RMS .and. l_logNext ! Used for storing in update routines !

         if ( n_time_step == 1 ) l_log=.true.

         !-- Compute one more iteration to properly terminate computations of
         !-- viscosity and pressure in the FD setup
         if ( l_last_RMS .and. l_stop_time ) then
            lRmsNext   =.true.
            l_logNext  =.true.
            l_last_RMS =.false.
            lRmsCalc   =.false.
            l_dtB      =.false.
            l_stop_time=.false.
         end if

         if ( l_stop_time ) then                  ! Programm stopped by kill -30
            l_new_rst_file=.true.                 ! Write rst-file and some
            if ( n_stores > 0 ) l_store=.true.    ! diagnostics before dying !
            l_log=.true.
            lRmsNext=.false.
            if ( n_specs > 0 ) l_spectrum=.true.
         end if

         lHelCalc     =l_hel        .and. l_log
         lPowerCalc   =l_power      .and. l_log
         lPerpParCalc =l_perpPar    .and. l_log
         lGeosCalc    =l_par        .and. l_log
         lFluxProfCalc=l_FluxProfs  .and. l_log
         lViscBcCalc  =l_ViscBcCalc .and. l_log

         l_HT  = (l_frame .and. l_movie) .or. lViscBcCalc
         lPressCalc=lRmsCalc .or. ( l_PressGraph .and. l_graph )  &
         &            .or. lFluxProfCalc
         lPressNext=( l_RMS .or. l_FluxProfs ) .and. l_logNext

         if ( l_graph ) call open_graph_file(n_time_step, time)

         tscheme%istage = 1

         do n_stage=1,tscheme%nstages

            !--- Now the real work starts with the radial loop that calculates
            !    the nonlinear terms:
            if ( lVerbose ) then
               write(output_unit,*)
               write(output_unit,*) '! Starting radial loop!'
            end if

            !------------------------
            !-- Storage or special calculatons computed in radial loop need to be
            !-- only done on the first sub-stage
            !------------------------
            l_graph       = l_graph       .and. (tscheme%istage==1)
            l_frame       = l_frame       .and. (tscheme%istage==1)
            lTOCalc       = lTOCalc       .and. (tscheme%istage==1)
            lTONext       = lTONext       .and. (tscheme%istage==1)
            lTONext2      = lTONext2      .and. (tscheme%istage==1)
            lHelCalc      = lHelCalc      .and. (tscheme%istage==1)
            lPowerCalc    = lPowerCalc    .and. (tscheme%istage==1)
            lRmsCalc      = lRmsCalc      .and. (tscheme%istage==1)
            lPressCalc    = lPressCalc    .and. (tscheme%istage==1)
            lViscBcCalc   = lViscBcCalc   .and. (tscheme%istage==1)
            lFluxProfCalc = lFluxProfCalc .and. (tscheme%istage==1)
            lPerpParCalc  = lPerpParCalc  .and. (tscheme%istage==1)
            lGeosCalc     = lGeosCalc     .and. (tscheme%istage==1)
            l_probe_out   = l_probe_out   .and. (tscheme%istage==1)

            if ( tscheme%l_exp_calc(n_stage) ) then

               !----------------
               !- Mloc -> Rloc transposes
               !----------------
               call transp_LMloc_to_Rloc(comm_counter, l_finish_exp_early, &
                    &                    lPressCalc,l_HT)

               !---------------
               !- Radial loop
               !---------------
               call rLoop_counter%start_count()
               if ( l_parallel_solve ) then
                  if ( l_mag_par_solve ) then
                  call radialLoopG(l_graph, l_frame,time,timeStage,tscheme,            &
                       &           dtLast,lTOCalc,lTONext,lTONext2,lHelCalc,           &
                       &           lPowerCalc,lRmsCalc,lPressCalc,lPressNext,          &
                       &           lViscBcCalc,lFluxProfCalc,lPerpParCalc,lGeosCalc,   &
                       &           l_probe_out,dsdt%expl(:,:,tscheme%istage),          &
                       &           dwdt%expl(:,:,tscheme%istage),                      &
                       &           dzdt%expl(:,:,tscheme%istage),                      &
                       &           dpdt%expl(:,:,tscheme%istage),                      &
                       &           dxidt%expl(:,:,tscheme%istage),                     &
                       &           dbdt%expl(:,:,tscheme%istage),                      &
                       &           djdt%expl(:,:,tscheme%istage),dVxVhLM_Rloc,         &
                       &           dVxBhLM_Rloc,dVSrLM_Rloc,dVXirLM_Rloc,              &
                       &           lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb,   &
                       &           br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,HelASr_Rloc, &
                       &           Hel2ASr_Rloc,HelnaASr_Rloc,Helna2ASr_Rloc,          &
                       &           HelEAASr_Rloc,viscASr_Rloc,uhASr_Rloc,duhASr_Rloc,  &
                       &           gradsASr_Rloc,fconvASr_Rloc,fkinASr_Rloc,           &
                       &           fviscASr_Rloc,fpoynASr_Rloc,fresASr_Rloc,           &
                       &           EperpASr_Rloc,EparASr_Rloc,EperpaxiASr_Rloc,        &
                       &           EparaxiASr_Rloc,dtrkc_Rloc,dthkc_Rloc)
                  else
                  call radialLoopG(l_graph, l_frame,time,timeStage,tscheme,            &
                       &           dtLast,lTOCalc,lTONext,lTONext2,lHelCalc,           &
                       &           lPowerCalc,lRmsCalc,lPressCalc,lPressNext,          &
                       &           lViscBcCalc,lFluxProfCalc,lPerpParCalc,lGeosCalc,   &
                       &           l_probe_out,dsdt%expl(:,:,tscheme%istage),          &
                       &           dwdt%expl(:,:,tscheme%istage),                      &
                       &           dzdt%expl(:,:,tscheme%istage),                      &
                       &           dpdt%expl(:,:,tscheme%istage),                      &
                       &           dxidt%expl(:,:,tscheme%istage),                     &
                       &           dbdt_Rloc,djdt_Rloc,dVxVhLM_Rloc,                   &
                       &           dVxBhLM_Rloc,dVSrLM_Rloc,dVXirLM_Rloc,              &
                       &           lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb,   &
                       &           br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,HelASr_Rloc, &
                       &           Hel2ASr_Rloc,HelnaASr_Rloc,Helna2ASr_Rloc,          &
                       &           HelEAASr_Rloc,viscASr_Rloc,uhASr_Rloc,duhASr_Rloc,  &
                       &           gradsASr_Rloc,fconvASr_Rloc,fkinASr_Rloc,           &
                       &           fviscASr_Rloc,fpoynASr_Rloc,fresASr_Rloc,           &
                       &           EperpASr_Rloc,EparASr_Rloc,EperpaxiASr_Rloc,        &
                       &           EparaxiASr_Rloc,dtrkc_Rloc,dthkc_Rloc)
                  end if
               else
                  call radialLoopG(l_graph, l_frame,time,timeStage,tscheme,            &
                       &           dtLast,lTOCalc,lTONext,lTONext2,lHelCalc,           &
                       &           lPowerCalc,lRmsCalc,lPressCalc,lPressNext,          &
                       &           lViscBcCalc,lFluxProfCalc,lPerpParCalc,lGeosCalc,   &
                       &           l_probe_out,dsdt_Rloc,dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,&
                       &           dxidt_Rloc,dbdt_Rloc,djdt_Rloc,dVxVhLM_Rloc,        &
                       &           dVxBhLM_Rloc,dVSrLM_Rloc,dVXirLM_Rloc,              &
                       &           lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb,   &
                       &           br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,HelASr_Rloc, &
                       &           Hel2ASr_Rloc,HelnaASr_Rloc,Helna2ASr_Rloc,          &
                       &           HelEAASr_Rloc,viscASr_Rloc,uhASr_Rloc,duhASr_Rloc,  &
                       &           gradsASr_Rloc,fconvASr_Rloc,fkinASr_Rloc,           &
                       &           fviscASr_Rloc,fpoynASr_Rloc,fresASr_Rloc,           &
                       &           EperpASr_Rloc,EparASr_Rloc,EperpaxiASr_Rloc,        &
                       &           EparaxiASr_Rloc,dtrkc_Rloc,dthkc_Rloc)
               end if
               call rLoop_counter%stop_count()

               if ( lVerbose ) write(output_unit,*) '! r-loop finished!'

#ifdef WITH_MPI
               ! ------------------
               ! also exchange the lorentz_torques which are only
               ! set at the boundary points  but are needed on all processes.
               ! ------------------
               call MPI_Bcast(lorentz_torque_ic,1,MPI_DEF_REAL, &
                    &         n_procs-1,MPI_COMM_WORLD,ierr)
               call MPI_Bcast(lorentz_torque_ma,1,MPI_DEF_REAL, &
                    &         0,MPI_COMM_WORLD,ierr)
#endif

               !---------------
               ! Finish assembing the explicit terms
               !---------------
               if ( l_finish_exp_early ) then
                  call f_exp_counter%start_count()
                  if ( l_parallel_solve ) then
                     if ( l_mag_par_solve ) then
                     call finish_explicit_assembly_Rdist(omega_ic,w_Rloc,b_ic_LMloc,   &
                          &                      aj_ic_LMloc,dVSrLM_RLoc,dVXirLM_RLoc, &
                          &                      dVxVhLM_Rloc,dVxBhLM_Rloc,            &
                          &                      lorentz_torque_ma,lorentz_torque_ic,  &
                          &                      dsdt%expl(:,:,tscheme%istage),        &
                          &                      dxidt%expl(:,:,tscheme%istage),       &
                          &                      dwdt%expl(:,:,tscheme%istage),        &
                          &                      djdt%expl(:,:,tscheme%istage),        &
                          &                      dbdt_ic, djdt_ic,                     &
                          &                      domega_ma_dt, domega_ic_dt,           &
                          &                      lorentz_torque_ma_dt,                 &
                          &                      lorentz_torque_ic_dt, tscheme)
                     else
                     call finish_explicit_assembly_Rdist(omega_ic,w_Rloc,b_ic_LMloc,   &
                          &                      aj_ic_LMloc,dVSrLM_RLoc,dVXirLM_RLoc, &
                          &                      dVxVhLM_Rloc,dVxBhLM_Rloc,            &
                          &                      lorentz_torque_ma,lorentz_torque_ic,  &
                          &                      dsdt%expl(:,:,tscheme%istage),        &
                          &                      dxidt%expl(:,:,tscheme%istage),       &
                          &                      dwdt%expl(:,:,tscheme%istage),        &
                          &                      djdt_Rloc, dbdt_ic, djdt_ic,          &
                          &                      domega_ma_dt, domega_ic_dt,           &
                          &                      lorentz_torque_ma_dt,                 &
                          &                      lorentz_torque_ic_dt, tscheme)
                     end if
                  else
                     call finish_explicit_assembly_Rdist(omega_ic,w_Rloc,b_ic_LMloc,   &
                          &                      aj_ic_LMloc,dVSrLM_RLoc,dVXirLM_RLoc, &
                          &                      dVxVhLM_Rloc,dVxBhLM_Rloc,            &
                          &                      lorentz_torque_ma,lorentz_torque_ic,  &
                          &                      dsdt_Rloc, dxidt_Rloc, dwdt_Rloc,     &
                          &                      djdt_Rloc, dbdt_ic, djdt_ic,          &
                          &                      domega_ma_dt, domega_ic_dt,           &
                          &                      lorentz_torque_ma_dt,                 &
                          &                      lorentz_torque_ic_dt, tscheme)
                  end if
                  call f_exp_counter%stop_count()
               end if

               !----------------
               !-- Rloc to Mloc transposes
               !----------------
               call transp_Rloc_to_LMloc(comm_counter,tscheme%istage, &
                    &                    l_finish_exp_early, lPressNext)

               !------ Nonlinear magnetic boundary conditions:
               !       For stress-free conducting boundaries
               PERFON('nl_m_bnd')
               if ( l_b_nl_cmb .and. (nRStart <= n_r_cmb) ) then
                  b_nl_cmb(1) =zero
                  aj_nl_cmb(1)=zero
                  call get_b_nl_bcs('CMB', br_vt_lm_cmb,br_vp_lm_cmb,              &
                       &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_cmb(2:lm_max))
               end if
               !-- Replace by scatter from rank to lo (and in updateB accordingly)
               if ( l_b_nl_cmb ) then
#ifdef WITH_MPI
                  call MPI_Bcast(b_nl_cmb,lm_max,MPI_DEF_COMPLEX,0, &
                       &         MPI_COMM_WORLD,ierr)
                  call MPI_Bcast(aj_nl_cmb,lm_max,MPI_DEF_COMPLEX,0, &
                       &         MPI_COMM_WORLD,ierr)
#endif
               end if
               if ( l_b_nl_icb .and. (nRstop >= n_r_icb) ) then
                  aj_nl_icb(1)=zero
                  call get_b_nl_bcs('ICB', br_vt_lm_icb,br_vp_lm_icb,              &
                       &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_icb(2:lm_max))
               end if
               if ( l_b_nl_icb ) then
#ifdef WITH_MPI
                  call MPI_Bcast(aj_nl_icb,lm_max,MPI_DEF_COMPLEX,n_procs-1, &
                       &         MPI_COMM_WORLD,ierr)
#endif
               end if
               PERFOFF

               !---------------
               ! Finish assembing the explicit terms
               !---------------
               call lmLoop_counter%start_count()
               if ( .not. l_finish_exp_early ) then
                  call f_exp_counter%start_count()
                  call finish_explicit_assembly(omega_ic,w_LMloc,b_ic_LMloc,         &
                       &                        aj_ic_LMloc,                         &
                       &                        dVSrLM_LMLoc(:,:,tscheme%istage),    &
                       &                        dVXirLM_LMLoc(:,:,tscheme%istage),   &
                       &                        dVxVhLM_LMloc(:,:,tscheme%istage),   &
                       &                        dVxBhLM_LMloc(:,:,tscheme%istage),   &
                       &                        lorentz_torque_ma,lorentz_torque_ic, &
                       &                        dsdt, dxidt, dwdt, djdt, dbdt_ic,    &
                       &                        djdt_ic, domega_ma_dt, domega_ic_dt, &
                       &                        lorentz_torque_ma_dt,                &
                       &                        lorentz_torque_ic_dt, tscheme)
                  call f_exp_counter%stop_count()
               end if
               call lmLoop_counter%stop_count(l_increment=.false.)
            end if

            !------------
            !--- Output before update of fields in LMLoop:
            !------------
            if ( tscheme%istage == 1 ) then
               if ( lVerbose ) write(output_unit,*) "! start output"

               if ( l_cmb .and. l_dt_cmb_field ) then
                  call scatter_from_rank0_to_lo(dbdt_Rloc(:,n_r_cmb), dbdt_CMB_LMloc)
               end if

               if ( lVerbose ) write(output_unit,*) "! start real output"
               call io_counter%start_count()
               if ( l_parallel_solve .and. (l_log .or. l_spectrum .or. lTOCalc .or. &
               &    l_dtB .or. l_cmb .or. l_r .or. l_pot .or. l_store .or. l_frame) ) then
                  call transp_Rloc_to_LMloc_IO()
               end if
               call output(time,tscheme,n_time_step,l_stop_time,l_pot,l_log,      &
                    &      l_graph,lRmsCalc,l_store,l_new_rst_file,               &
                    &      l_spectrum,lTOCalc,lTOframe,                           &
                    &      l_frame,n_frame,l_cmb,n_cmb_sets,l_r,                  &
                    &      lorentz_torque_ic,lorentz_torque_ma,dbdt_CMB_LMloc,    &
                    &      HelASr_Rloc,Hel2ASr_Rloc,HelnaASr_Rloc,Helna2ASr_Rloc, &
                    &      HelEAASr_Rloc,viscASr_Rloc,uhASr_Rloc,duhASr_Rloc,     &
                    &      gradsASr_Rloc,fconvASr_Rloc,fkinASr_Rloc,fviscASr_Rloc,&
                    &      fpoynASr_Rloc,fresASr_Rloc,EperpASr_Rloc,EparASr_Rloc, &
                    &      EperpaxiASr_Rloc,EparaxiASr_Rloc)
               call io_counter%stop_count()
               if ( lVerbose ) write(output_unit,*) "! output finished"

               if ( l_graph ) call close_graph_file()

               !----- Finish time stepping, the last step is only for output!
               if ( l_stop_time ) exit outer  ! END OF TIME INTEGRATION

               dtLast = tscheme%dt(1) ! Old time step (needed for some TO outputs)

               !---------------------
               !-- Checking Courant criteria, l_new_dt and dt_new are output
               !---------------------
               call dt_courant(dtr,dth,l_new_dt,tscheme%dt(1),dt_new,dtMax, &
                    &          dtrkc_Rloc,dthkc_Rloc,time)

               !--------------------
               !-- Set weight arrays
               !--------------------
               call tscheme%set_dt_array(dt_new,dtMin,time,n_log_file,n_time_step,&
                    &                    l_new_dt)

               !-- Store the old weight factor of matrices
               !-- if it changes because of dt factors moving
               !-- matrix also needs to be rebuilt
               call tscheme%set_weights(lMatNext)

               !----- Advancing time:
               timeLast=time               ! Time of the previous time step
               time    =time+tscheme%dt(1) ! Update time

            end if

            call tscheme%get_time_stage(timeLast, timeStage)

            lMat = lMatNext
            if ( (l_new_dt .or. lMat) .and. (tscheme%istage==1) ) then
               !----- Calculate matrices for new time step if dt /= dtLast
               lMat=.true.
               if ( rank == 0 ) then
                  write(output_unit,'(1p,'' ! Building matrices at time step:'', &
                  &                   i8,ES16.6)') n_time_step,time
               end if
            end if
            lMatNext = .false.

            !-- If the scheme is a multi-step scheme that is not Crank-Nicolson
            !-- we have to use a different starting scheme
            call start_from_another_scheme(timeLast, l_bridge_step, n_time_step, tscheme)

            !---------------
            !-- LM Loop (update routines)
            !---------------
            if ( (.not. tscheme%l_assembly) .or. (tscheme%istage/=tscheme%nstages) ) then
               if ( lVerbose ) write(output_unit,*) '! starting lm-loop!'
               call lmLoop_counter%start_count()
               if ( l_parallel_solve ) then
                  call LMLoop_Rdist(timeStage,time,tscheme,lMat,lRmsNext,lPressNext,dsdt,&
                       &            dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dbdt_ic,djdt_ic,      &
                       &            domega_ma_dt,domega_ic_dt,lorentz_torque_ma_dt,      &
                       &            lorentz_torque_ic_dt,b_nl_cmb,aj_nl_cmb,aj_nl_icb)
               else
                  call LMLoop(timeStage,time,tscheme,lMat,lRmsNext,lPressNext,dsdt,  &
                       &      dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dbdt_ic,djdt_ic,        &
                       &      domega_ma_dt,domega_ic_dt,lorentz_torque_ma_dt,        &
                       &      lorentz_torque_ic_dt,b_nl_cmb,aj_nl_cmb,aj_nl_icb)
               end if

               if ( lVerbose ) write(output_unit,*) '! lm-loop finished!'

               !-- Timer counters
               call lmLoop_counter%stop_count()
               if ( tscheme%istage == 1 .and. lMat ) l_mat_time=.true.
               if ( tscheme%istage == 1 .and. .not. lMat .and. &
               &    .not. l_log ) l_pure=.true.

               ! Increment current stage
               tscheme%istage = tscheme%istage+1
            end if
         end do

         !----------------------------
         !-- Assembly stage of IMEX-RK (if needed)
         !----------------------------
         if ( tscheme%l_assembly ) then
            if ( l_parallel_solve ) then
               call assemble_stage_Rdist(time, omega_ic, omega_ic1, omega_ma, omega_ma1,&
                    &                    dwdt, dzdt, dpdt, dsdt, dxidt, dbdt, djdt,     &
                    &                    dbdt_ic, djdt_ic, domega_ic_dt, domega_ma_dt,  &
                    &                    lorentz_torque_ic_dt, lorentz_torque_ma_dt,    &
                    &                    lPressNext, lRmsNext, tscheme)
            else
               call assemble_stage(time, omega_ic, omega_ic1, omega_ma, omega_ma1,     &
                    &              dwdt, dzdt, dpdt, dsdt, dxidt, dbdt, djdt, dbdt_ic, &
                    &              djdt_ic, domega_ic_dt, domega_ma_dt,                &
                    &              lorentz_torque_ic_dt, lorentz_torque_ma_dt,         &
                    &              lPressNext, lRmsNext, tscheme)
            end if
         end if

         !-- Update counters
         if ( l_mat_time ) call mat_counter%stop_count()
         if ( l_pure ) call pure_counter%stop_count()
         call tot_counter%stop_count()

         !-----------------------
         !----- Timing and info of advance:
         !-----------------------
         run_time_passed = tot_counter%tTot/real(tot_counter%n_counts,cp)
         if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >=  &
         &    real(n_time_steps,cp)  .or. n_time_steps < 31 ) then
            write(message,'(" ! Time step finished:",i6)') n_time_step
            call logWrite(message)
            if ( real(n_time_step,cp)+tenth_n_time_steps*real(nPercent,cp) >= &
            &    real(n_time_steps,cp) .and. n_time_steps >= 10 ) then
               write(message,'(" ! This is           :",i3,"%")') (10-nPercent)*10
               call logWrite(message)
               nPercent=nPercent-1
            end if
            !tot_counter%tTtop%
            if ( rank == 0 ) then
               call formatTime(output_unit,' ! Mean wall time for time step:',  &
               &               run_time_passed)
               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if
               call formatTime(n_log_file,' ! Mean wall time for time step:', &
               &               run_time_passed)
               if ( l_save_out ) close(n_log_file)
            end if
         end if

      end do outer ! end of time stepping !

      !LIKWID_OFF('tloop')
      PERFOFF

      if ( l_movie ) then
         if ( rank == 0 ) then
            if (n_frame > 0) then
               write(output_unit,'(1p,/,/,A,i10,3(/,A,ES16.6))')          &
               &     " !  No of stored movie frames: ",n_frame,           &
               &     " !     starting at time: ",t_movieS(1)*tScale,      &
               &     " !       ending at time: ",t_movieS(n_frame)*tScale,&
               &     " !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if
               write(n_log_file,'(1p,/,/,A,i10,3(/,A,ES16.6))')           &
               &     " !  No of stored movie frames: ",n_frame,           &
               &     " !     starting at time: ",t_movieS(1)*tScale,      &
               &     " !       ending at time: ",t_movieS(n_frame)*tScale,&
               &     " !      with step width: ",(t_movieS(2)-t_movieS(1))*tScale
               if ( l_save_out ) close(n_log_file)
            else
               write(output_unit,'(1p,/,/,A,i10,3(/,A,ES16.6))')  &
               &     " !  No of stored movie frames: ",n_frame,   &
               &     " !     starting at time: ",0.0_cp,          &
               &     " !       ending at time: ",0.0_cp,          &
               &     " !      with step width: ",0.0_cp
               if ( l_save_out ) then
                  open(newunit=n_log_file, file=log_file, status='unknown', &
                  &    position='append')
               end if
               write(n_log_file,'(1p,/,/,A,i10,3(/,A,ES16.6))') &
               &     " !  No of stored movie frames: ",n_frame, &
               &     " !     starting at time: ",0.0_cp,        &
               &     " !       ending at time: ",0.0_cp,        &
               &     " !      with step width: ",0.0_cp
               if ( l_save_out ) close(n_log_file)
            end if
         end if
      end if

      if ( l_cmb_field ) then
         write(message,'(A,i9)') " !  No of stored sets of b at CMB: ",n_cmb_sets
         call logWrite(message)
      end if

      if ( rank == 0 ) then
         write(output_unit,*)
         call logWrite('')
      end if

      if ( rank == 0 .and. l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if
      call rLoop_counter%finalize('! Mean wall time for r Loop                 :', &
           &                      n_log_file)
      call phy2lm_counter%finalize('!    - Time taken for Spat->Spec            :',&
           &                       n_log_file)
      call lm2phy_counter%finalize('!    - Time taken for Spec->Spat            :',&
           &                       n_log_file)
      call nl_counter%finalize('!    - Time taken for nonlinear terms       :',&
           &                       n_log_file)
      call td_counter%finalize('!    - Time taken for time derivative terms :',&
           &                       n_log_file)
      call lmLoop_counter%finalize('! Mean wall time for LM Loop                :',&
           &                       n_log_file)
      call f_exp_counter%finalize('!     - Time taken to compute r-der of adv. :', &
           &                      n_log_file)
      call dct_counter%finalize('!     - Time taken for DCTs and r-der       :',   &
           &                    n_log_file)
      call solve_counter%finalize('!     - Time taken for linear solves        :', &
           &                      n_log_file)
      call comm_counter%finalize('! Mean wall time for MPI communications     :',  &
           &                     n_log_file)
      call mat_counter%finalize('! Mean wall time for t-step with matrix calc:',   &
           &                    n_log_file)
      call io_counter%finalize('! Mean wall time for output routine         :',  &
           &                   n_log_file)
      call pure_counter%finalize('! Mean wall time for one pure time step     :', &
           &                     n_log_file)
      call tot_counter%finalize('! Mean wall time for one time step          :', &
           &                    n_log_file)
      if ( rank==0 .and. l_save_out ) close(n_log_file)

      !-- WORK IS DONE !

   end subroutine step_time
!------------------------------------------------------------------------------
   subroutine start_from_another_scheme(time, l_bridge_step, n_time_step, tscheme)
      !
      ! This subroutine is used to initialize multisteps schemes whenever previous
      ! steps are not known. In that case a CN/AB2 scheme is used to bridge the
      ! missing steps.
      !

      !-- Input variables
      real(cp),            intent(in) :: time
      logical,             intent(in) :: l_bridge_step
      integer,             intent(in) :: n_time_step

      !-- Output variable
      class(type_tscheme), intent(inout) :: tscheme

      !-- If the scheme is a multi-step scheme that is not Crank-Nicolson
      !-- we have to use a different starting scheme
      if ( l_bridge_step .and. tscheme%time_scheme /= 'CNAB2' .and.  &
           n_time_step <= tscheme%nold-1 .and.                       &
           tscheme%family=='MULTISTEP' ) then

         if ( l_single_matrix ) then
            call get_single_rhs_imp(s_LMloc, ds_LMloc, w_LMloc, dw_LMloc,     &
                 &                  ddw_LMloc, p_LMloc, dp_LMloc, dsdt, dwdt, &
                 &                  dpdt, tscheme, 1, .true., .false.)
         else
            if ( l_parallel_solve ) then
               call bulk_to_ghost(w_Rloc, w_ghost, 2, nRstart, nRstop, lm_max, 1, lm_max)
               call exch_ghosts(w_ghost, lm_max, nRstart, nRstop, 2)
               call fill_ghosts_W(w_ghost)
               call get_pol_rhs_imp_ghost(w_ghost, dw_Rloc, ddw_Rloc, p_Rloc, dp_Rloc, &
                    &                     dwdt, tscheme, 1, .true., .false., .false.,  &
                    &                     dwdt%expl(:,:,1)) ! Work array
            else
               call get_pol_rhs_imp(s_LMloc, xi_LMloc, w_LMloc, dw_LMloc, ddw_LMloc,  &
                    &               p_LMloc, dp_LMloc, dwdt, dpdt, tscheme, 1,        &
                    &               .true., .false., .false., work_LMloc)
            end if
            if ( l_heat ) then
               if ( l_parallel_solve ) then
                  call bulk_to_ghost(s_Rloc, s_ghost, 1, nRstart, nRstop, lm_max, 1, &
                       &             lm_max)
                  call exch_ghosts(s_ghost, lm_max, nRstart, nRstop, 1)
                  call fill_ghosts_S(s_ghost)
                  call get_entropy_rhs_imp_ghost(s_ghost, ds_Rloc, dsdt, 1, .true.)
               else
                  call get_entropy_rhs_imp(s_LMloc, ds_LMloc, dsdt, 1, .true.)
               end if
            end if
         end if

         call get_rot_rates(omega_ma, lorentz_torque_ma_dt%old(1))
         call get_rot_rates(omega_ic, lorentz_torque_ic_dt%old(1))

         if ( l_parallel_solve ) then
            call bulk_to_ghost(z_Rloc, z_ghost, 1, nRstart, nRstop, lm_max, 1, lm_max)
            call exch_ghosts(z_ghost, lm_max, nRstart, nRstop, 1)
            call fill_ghosts_Z(z_ghost)
            call get_tor_rhs_imp_ghost(time, z_ghost, dz_Rloc, dzdt, domega_ma_dt,  &
                 &                     domega_ic_dt, omega_ic, omega_ma, omega_ic1, &
                 &                     omega_ma1, tscheme, 1, .true., .false.)
         else
            call get_tor_rhs_imp(time, z_LMloc, dz_LMloc, dzdt, domega_ma_dt,  &
                 &               domega_ic_dt, omega_ic, omega_ma, omega_ic1,  &
                 &               omega_ma1, tscheme, 1, .true., .false.)
         end if

         if ( l_chemical_conv ) then
            if ( l_parallel_solve ) then
                  call bulk_to_ghost(xi_Rloc, xi_ghost, 1, nRstart, nRstop, lm_max, &
                       &             1, lm_max)
                  call exch_ghosts(xi_ghost, lm_max, nRstart, nRstop, 1)
                  call fill_ghosts_Xi(xi_ghost)
                  call get_comp_rhs_imp_ghost(xi_ghost, dxidt, 1, .true.)
            else
               call get_comp_rhs_imp(xi_LMloc, dxi_LMloc, dxidt, 1, .true.)
            end if
         end if

         if ( l_mag ) then
            if ( l_mag_par_solve ) then
               call bulk_to_ghost(b_Rloc, b_ghost, 1, nRstart, nRstop, lm_max, 1, &
                    &             lm_max)
               call bulk_to_ghost(aj_Rloc, aj_ghost, 1, nRstart, nRstop, lm_max, 1, &
                    &             lm_max)
               call exch_ghosts(aj_ghost, lm_max, nRstart, nRstop, 1)
               call exch_ghosts(b_ghost, lm_max, nRstart, nRstop, 1)
               call fill_ghosts_B(b_ghost, aj_ghost)
               call get_mag_rhs_imp_ghost(b_ghost, db_Rloc, ddb_RLoc, aj_ghost,    &
                    &                     dj_Rloc, ddj_Rloc,  dbdt, djdt, tscheme, &
                    &                     1, .true., .false.)
            else
               call get_mag_rhs_imp(b_LMloc, db_LMloc, ddb_LMLoc, aj_LMLoc, dj_LMloc, &
                    &               ddj_LMloc, dbdt, djdt, tscheme, 1, .true., .false.)
            end if
         end if

         if ( l_cond_ic ) call get_mag_ic_rhs_imp(b_ic_LMloc, db_ic_LMloc,     &
                               &                  ddb_ic_LMLoc, aj_ic_LMLoc,   &
                               &                  dj_ic_LMloc, ddj_ic_LMloc,   &
                               &                  dbdt_ic, djdt_ic, 1, .true.)

         call tscheme%bridge_with_cnab2()

      end if

      if ( l_AB1 .and. n_time_step == 1 ) then
         call tscheme%start_with_ab1()
         l_AB1 = .false.
      end if

   end subroutine start_from_another_scheme
!--------------------------------------------------------------------------------
   subroutine transp_LMloc_to_Rloc(comm_counter, l_Rloc, lPressCalc, lHTCalc)
      ! Here now comes the block where the LM distributed fields
      ! are redistributed to Rloc distribution which is needed for
      ! the radialLoop.

      !-- Input variables
      logical, intent(in) :: l_Rloc, lPressCalc, lHTCalc

      !-- Output variable
      type(timer_type), intent(inout) :: comm_counter

      call comm_counter%start_count()
      if ( l_packed_transp ) then
         if ( l_Rloc ) then
            if ( (.not. l_parallel_solve) .or. (l_mag .and. .not. l_mag_par_solve) ) then
               call lo2r_flow%transp_lm2r(flow_LMloc_container, flow_Rloc_container)
            end if
            if ( l_heat .and. lHTCalc .and. (.not. l_parallel_solve) ) then
               call get_dr_Rloc(s_Rloc, ds_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
            if ( l_chemical_conv .and. (.not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(xi_LMloc,xi_Rloc)
            end if
            if ( (l_conv .or. l_mag_kin) .and. (.not. l_parallel_solve) ) then
               call get_ddr_Rloc(w_Rloc, dw_Rloc, ddw_Rloc, lm_max, nRstart, nRstop, &
                    &            n_r_max, rscheme_oc)
               call get_dr_Rloc(z_Rloc, dz_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
            if ( lPressCalc .and. ( .not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(p_LMloc, p_Rloc)
               call get_dr_Rloc(p_Rloc, dp_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
            if ( l_mag .and. ( .not. l_mag_par_solve ) ) then
               call get_ddr_Rloc(b_Rloc, db_Rloc, ddb_Rloc, lm_max, nRstart, nRstop, &
                    &            n_r_max, rscheme_oc)
               call get_dr_Rloc(aj_Rloc, dj_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
         else
            if ( l_heat ) then
               !if ( .not. l_parallel_solve ) then
               call lo2r_one%transp_lm2r(s_LMloc, s_Rloc)
               if ( lHTCalc ) call lo2r_one%transp_lm2r(ds_LMloc, ds_Rloc)
            end if
            if ( l_chemical_conv ) call lo2r_one%transp_lm2r(xi_LMloc,xi_Rloc)
            if ( l_conv .or. l_mag_kin ) then
               call lo2r_flow%transp_lm2r(flow_LMloc_container,flow_Rloc_container)
            end if
            if ( lPressCalc ) then
               call lo2r_press%transp_lm2r(press_LMloc_container,press_Rloc_container)
            end if
            if ( l_mag ) then
               call lo2r_field%transp_lm2r(field_LMloc_container,field_Rloc_container)
            end if
         end if
      else
         if ( l_Rloc ) then
            if ( l_heat .and. (.not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(s_LMloc, s_Rloc)
               if ( lHTCalc ) then
                  call get_dr_Rloc(s_Rloc, ds_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                       &           rscheme_oc)
               end if
            end if
            if ( l_chemical_conv .and. (.not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(xi_LMloc,xi_Rloc)
            end if
            if ( (l_conv .or. l_mag_kin) .and. (.not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(w_LMloc, w_Rloc)
               call get_ddr_Rloc(w_Rloc, dw_Rloc, ddw_Rloc, lm_max, nRstart, nRstop, &
                    &            n_r_max, rscheme_oc)
               call lo2r_one%transp_lm2r(z_LMloc, z_Rloc)
               call get_dr_Rloc(z_Rloc, dz_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
            if ( lPressCalc .and. (.not. l_parallel_solve) ) then
               call lo2r_one%transp_lm2r(p_LMloc, p_Rloc)
               call get_dr_Rloc(p_Rloc, dp_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
            if ( l_mag .and. ( .not. l_mag_par_solve ) ) then
               call lo2r_one%transp_lm2r(b_LMloc, b_Rloc)
               call get_ddr_Rloc(b_Rloc, db_Rloc, ddb_Rloc, lm_max, nRstart, nRstop, &
                    &            n_r_max, rscheme_oc)
               call lo2r_one%transp_lm2r(aj_LMloc, aj_Rloc)
               call get_dr_Rloc(aj_Rloc, dj_Rloc, lm_max, nRstart, nRstop, n_r_max, &
                    &           rscheme_oc)
            end if
         else
            if ( l_heat ) then
               call lo2r_one%transp_lm2r(s_LMloc, s_Rloc)
               if ( lHTCalc ) call lo2r_one%transp_lm2r(ds_LMloc, ds_Rloc)
            end if
            if ( l_chemical_conv ) call lo2r_one%transp_lm2r(xi_LMloc,xi_Rloc)
            if ( l_conv .or. l_mag_kin ) then
               call lo2r_one%transp_lm2r(w_LMloc, w_Rloc)
               call lo2r_one%transp_lm2r(dw_LMloc, dw_Rloc)
               call lo2r_one%transp_lm2r(ddw_LMloc, ddw_Rloc)
               call lo2r_one%transp_lm2r(z_LMloc, z_Rloc)
               call lo2r_one%transp_lm2r(dz_LMloc, dz_Rloc)
            end if
            if ( lPressCalc ) then
               call lo2r_one%transp_lm2r(p_LMloc, p_Rloc)
               call lo2r_one%transp_lm2r(dp_LMloc, dp_Rloc)
            end if
            if ( l_mag ) then
               call lo2r_one%transp_lm2r(b_LMloc, b_Rloc)
               call lo2r_one%transp_lm2r(db_LMloc, db_Rloc)
               call lo2r_one%transp_lm2r(ddb_LMloc, ddb_Rloc)
               call lo2r_one%transp_lm2r(aj_LMloc, aj_Rloc)
               call lo2r_one%transp_lm2r(dj_LMloc, dj_Rloc)
            end if
         end if
      end if
      call comm_counter%stop_count(l_increment=.false.)

   end subroutine transp_LMloc_to_Rloc
!--------------------------------------------------------------------------------
   subroutine transp_Rloc_to_LMloc(comm_counter, istage, lRloc, lPressNext)
      !
      !- MPI transposition from r-distributed to LM-distributed
      !

      !-- Input variable
      logical, intent(in) :: lRloc
      logical, intent(in) :: lPressNext
      integer, intent(in) :: istage

      !-- Output variable
      type(timer_type), intent(inout) :: comm_counter

      if ( lVerbose ) write(output_unit,*) "! start r2lo redistribution"

      call comm_counter%start_count()
      if ( l_packed_transp ) then
         if ( lRloc ) then
            if ( (.not. l_parallel_solve) .or. ( l_mag .and. .not. l_mag_par_solve) ) then
               call r2lo_flow%transp_r2lm(dflowdt_Rloc_container, &
                    &                     dflowdt_LMloc_container(:,:,:,istage))
            end if
            if ( (l_conv .or. l_mag_kin) .and. (.not. l_parallel_solve) ) then
               if ( .not. l_double_curl .or. lPressNext ) then
                  call r2lo_one%transp_r2lm(dpdt_Rloc,dpdt%expl(:,:,istage))
               end if
            end if
            if ( l_chemical_conv .and. ( .not. l_parallel_solve ) ) then
               call r2lo_one%transp_r2lm(dxidt_Rloc,dxidt%expl(:,:,istage))
            end if
         else
            if ( l_conv .or. l_mag_kin ) then
               call r2lo_flow%transp_r2lm(dflowdt_Rloc_container,  &
                    &                     dflowdt_LMloc_container(:,:,:,istage))
            end if
            !if ( l_heat .and. (.not. l_parallel_solve) ) then
            if ( l_heat  ) then
               call r2lo_s%transp_r2lm(dsdt_Rloc_container,&
                    &                  dsdt_LMloc_container(:,:,:,istage))
            end if
            if ( l_chemical_conv ) then
               call r2lo_xi%transp_r2lm(dxidt_Rloc_container, &
                    &                   dxidt_LMloc_container(:,:,:,istage))
            end if
            if ( l_mag ) then
               call r2lo_field%transp_r2lm(dbdt_Rloc_container, &
                    &                      dbdt_LMloc_container(:,:,:,istage))
            end if
         end if
      else
         if ( lRloc ) then
            if ( (l_conv .or. l_mag_kin) .and. (.not. l_parallel_solve) ) then
               call r2lo_one%transp_r2lm(dwdt_Rloc,dwdt%expl(:,:,istage))
               if ( .not. l_parallel_solve ) then
                  call r2lo_one%transp_r2lm(dzdt_Rloc,dzdt%expl(:,:,istage))
               end if
               if ( (.not. l_double_curl .or. lPressNext) .and. &
               &    (.not.  l_parallel_solve) ) then
                  call r2lo_one%transp_r2lm(dpdt_Rloc,dpdt%expl(:,:,istage))
               end if
            end if
            if ( l_heat .and. (.not. l_parallel_solve) ) then
               call r2lo_one%transp_r2lm(dsdt_Rloc,dsdt%expl(:,:,istage))
            end if
            if ( l_chemical_conv .and. (.not. l_parallel_solve) ) then
               call r2lo_one%transp_r2lm(dxidt_Rloc,dxidt%expl(:,:,istage))
            end if
            if ( l_mag .and. ( .not. l_mag_par_solve ) ) then
               call r2lo_one%transp_r2lm(dbdt_Rloc,dbdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(djdt_Rloc,djdt%expl(:,:,istage))
            end if
         else
            if ( l_conv .or. l_mag_kin ) then
               call r2lo_one%transp_r2lm(dwdt_Rloc,dwdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(dzdt_Rloc,dzdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(dpdt_Rloc,dpdt%expl(:,:,istage))
               if ( l_double_curl ) then
                  call r2lo_one%transp_r2lm(dVxVhLM_Rloc,dVxVhLM_LMloc(:,:,istage))
               end if
            end if
            if ( l_heat .and. (.not. l_parallel_solve) ) then
               call r2lo_one%transp_r2lm(dsdt_Rloc,dsdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(dVSrLM_Rloc,dVSrLM_LMloc(:,:,istage))
            end if
            if ( l_chemical_conv ) then
               call r2lo_one%transp_r2lm(dxidt_Rloc,dxidt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(dVXirLM_Rloc,dVXirLM_LMloc(:,:,istage))
            end if
            if ( l_mag ) then
               call r2lo_one%transp_r2lm(dbdt_Rloc,dbdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(djdt_Rloc,djdt%expl(:,:,istage))
               call r2lo_one%transp_r2lm(dVxBhLM_Rloc,dVxBhLM_LMloc(:,:,istage))
            end if
         end if
      end if
      call comm_counter%stop_count()

      if ( lVerbose ) write(output_unit,*) "! r2lo redistribution finished"

   end subroutine transp_Rloc_to_LMloc
!--------------------------------------------------------------------------------
   subroutine transp_Rloc_to_LMloc_IO()
      !
      ! For now, most of the outputs use LM-distributed arrays as input. To handle
      ! that one has to transpose the missing fields.
      !

      if ( l_heat ) then
         call r2lo_one%transp_r2lm(s_Rloc,s_LMloc)
         call r2lo_one%transp_r2lm(ds_Rloc,ds_LMloc)
      end if
      if ( l_chemical_conv ) call r2lo_one%transp_r2lm(xi_Rloc,xi_LMloc)

      call r2lo_one%transp_r2lm(z_Rloc,z_LMloc)
      call r2lo_one%transp_r2lm(dz_Rloc,dz_LMloc)
      call r2lo_one%transp_r2lm(w_Rloc,w_LMloc)
      call r2lo_one%transp_r2lm(dw_Rloc,dw_LMloc)
      call r2lo_one%transp_r2lm(ddw_Rloc,ddw_LMloc)

      if ( l_mag .and. l_mag_par_solve ) then
         call r2lo_one%transp_r2lm(b_Rloc,b_LMloc)
         call r2lo_one%transp_r2lm(db_Rloc,db_LMloc)
         call r2lo_one%transp_r2lm(ddb_Rloc,ddb_LMloc)
         call r2lo_one%transp_r2lm(aj_Rloc,aj_LMloc)
         call r2lo_one%transp_r2lm(dj_Rloc,dj_LMloc)
         call r2lo_one%transp_r2lm(ddj_Rloc,ddj_LMloc)
      end if

   end subroutine transp_Rloc_to_LMloc_IO
!--------------------------------------------------------------------------------
end module step_time_mod
