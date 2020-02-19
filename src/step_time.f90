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
   use truncation, only: n_r_max, l_max, l_maxMag, lm_max, lmP_max
   use num_param, only: n_time_steps, run_time_limit, tEnd, dtMax, &
       &                dtMin, tScale, dct_counter, nl_counter,    &
       &                solve_counter, lm2phy_counter, td_counter, &
       &                phy2lm_counter, nl_counter
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag, &
       &                  n_r_icb, n_r_cmb
   use logic, only: l_mag, l_mag_LF, l_dtB, l_RMS, l_hel, l_TO,        &
       &            l_TOmovie, l_r_field, l_cmb_field, l_HTmovie,      &
       &            l_DTrMagSpec, lVerbose, l_b_nl_icb,                &
       &            l_b_nl_cmb, l_FluxProfs, l_ViscBcCalc, l_perpPar,  &
       &            l_HT, l_dtB, l_dtBmovie, l_heat, l_conv, l_movie,  &
       &            l_runTimeLimit, l_save_out, l_bridge_step,         &
       &            l_dt_cmb_field, l_chemical_conv, l_mag_kin,        &
       &            l_power, l_double_curl, l_PressGraph, l_probe,     &
       &            l_AB1, l_finite_diff, l_cond_ic, l_single_matrix
   use init_fields, only: omega_ic1, omega_ma1
   use movie_data, only: t_movieS
   use radialLoop, only: radialLoopG
   use LMLoop_mod, only: LMLoop, finish_explicit_assembly
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
       &                  n_TOs, n_t_TO, t_TO, n_TOZ_step, n_TOZs,         &
       &                  n_t_TOZ, t_TOZ, n_probe_step, n_probe_out,       &
       &                  n_t_probe, t_probe, log_file, n_log_file,        &
       &                  n_time_hits
   use updateB_mod, only: get_mag_rhs_imp, get_mag_ic_rhs_imp
   use updateWP_mod, only: get_pol_rhs_imp
   use updateWPS_mod, only: get_single_rhs_imp
   use updateS_mod, only: get_entropy_rhs_imp
   use updateXI_mod, only: get_comp_rhs_imp
   use updateZ_mod, only: get_tor_rhs_imp
   use output_mod, only: output
   use time_schemes, only: type_tscheme
   use useful, only: l_correct_step, logWrite
   use communications, only: lo2r_field, lo2r_flow, scatter_from_rank0_to_lo, &
       &                     lo2r_xi,  r2lo_flow, r2lo_s, r2lo_xi,r2lo_field, &
       &                     lo2r_s, lo2r_press
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
      logical :: lperpParCalc     ! Calculate perpendicular and parallel Ekin
      logical :: lTOCalc          ! Calculate TO stuff
      logical :: lTONext,lTONext2 ! TO stuff for next steps
      logical :: lTOframeNext,lTOframeNext2
      logical :: lTOZhelp,lTOZwrite
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

      !--- Explicit part of time stepping, calculated in s_radialLoopG.f and
      !    passed to LMLoop.f where the time step is preformed.
      !    Note that the respective arrays for the changes in inner-core
      !    magnetic field are calculated in s_updateB.f and are only
      !    needed there.

      !--- Lorentz torques:
      real(cp) :: lorentz_torque_ma,lorentz_torque_ic

      !-- Arrays for m_outMisc.F90 and m_outPar.F90
      real(cp) :: HelLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: Hel2LMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: HelnaLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: Helna2LMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: viscLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: uhLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: duhLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: gradsLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fconvLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fkinLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fviscLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: fpoynLMr_Rloc(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp) :: fresLMr_Rloc(l_maxMag+1,nRstartMag:nRstopMag)
      real(cp) :: EperpLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EparLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EperpaxiLMr_Rloc(l_max+1,nRstart:nRstop)
      real(cp) :: EparaxiLMr_Rloc(l_max+1,nRstart:nRstop)

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

      if ( lVerbose ) write(*,'(/,'' ! STARTING STEP_TIME !'')')

      run_time_passed=0.0_cp
      l_log       =.false.
      l_stop_time =.false.
      l_new_dt    =.true.   ! Invokes calculation of t-step matrices
      lMatNext    =.true.
      timeLast    =time
      timeStage   =time

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
         write(*,*)
         write(*,*) '! Starting time integration!'
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
            write(*,*)
            write(*,*) '! Starting time step ',n_time_step
         end if

         !-- Start time counters
         call mat_counter%start_count()
         call tot_counter%start_count()
         call pure_counter%start_count()
         l_pure=.false.
         l_mat_time=.false.

#ifdef WITH_MPI
         ! Broadcast omega_ic and omega_ma
         call MPI_Bcast(omega_ic,1,MPI_DEF_REAL,rank_with_l1m0, &
              &         MPI_COMM_WORLD,ierr)
         call MPI_Bcast(omega_ma,1,MPI_DEF_REAL,rank_with_l1m0, &
              &         MPI_COMM_WORLD,ierr)
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
               write(message,'("! Run time limit exeeded !")')
               call logWrite(message)
               l_stop_time=.true.
            end if
         end if
         if ( (n_stop_signal > 0) .or. (n_time_step == n_time_steps_go) ) then
            l_stop_time=.true.   ! last time step !
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
            l_dtB=( l_frame .and. l_dtBmovie ) .or.         &
            &                   ( l_log .and. l_DTrMagSpec )
         end if
         l_HT  = l_frame .and. l_HTmovie

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
         lTOZhelp= n_time_step > 2 .and. l_TO .and.                         &
         &                     l_correct_step(n_time_step-1,time,timeLast,  &
         &                 n_time_steps,n_TOZ_step,n_TOZs,n_t_TOZ,t_TOZ,0)
         if ( lTOZhelp ) then
            lTOZwrite=.true.
         else
            lTOZwrite=.false.
         end if

         lRmsCalc=(l_RMS .and. l_log .and. (n_time_step > 1)) .or. &
         &        (l_RMS .and. l_stop_time)
         if ( l_mag .or. l_mag_LF ) l_dtB = l_dtB .or. lRmsCalc
         lRmsNext=l_RMS .and. l_logNext ! Used for storing in update routines !

         if ( n_time_step == 1 ) l_log=.true.

         if ( l_stop_time ) then                  ! Programm stopped by kill -30
            l_new_rst_file=.true.                 ! Write rst-file and some
            if ( n_stores > 0 ) l_store=.true.    ! diagnostics before dying !
            l_log=.true.
            lRmsNext=.false.
         end if

         lHelCalc=l_hel.and.l_log

         lPowerCalc=l_power.and.l_log

         lperpParCalc=l_perpPar.and.l_log

         lFluxProfCalc =l_FluxProfs.and.l_log

         lViscBcCalc =l_ViscBcCalc.and.l_log

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
            lperpParCalc  = lperpParCalc  .and. (tscheme%istage==1)
            l_probe_out   = l_probe_out   .and. (tscheme%istage==1)

            if ( tscheme%l_exp_calc(n_stage) ) then

               !----------------------
               ! Here now comes the block where the LM distributed fields
               ! are redistributed to Rloc distribution which is needed for 
               ! the radialLoop.
               !----------------------
               call comm_counter%start_count()
               if ( l_heat ) then
                  call lo2r_s%transp_lm2r(s_LMloc_container, s_Rloc_container)
               end if
               if ( l_chemical_conv ) then
                  call lo2r_xi%transp_lm2r(xi_LMloc_container,xi_Rloc_container)
               end if
               if ( l_conv .or. l_mag_kin ) then
                  call lo2r_flow%transp_lm2r(flow_LMloc_container, &
                       &                     flow_Rloc_container)
               end if
               if ( lPressCalc ) then
                  call lo2r_press%transp_lm2r(press_LMloc_container, &
                       &                      press_Rloc_container)
               end if
               if ( l_mag ) then
                  call lo2r_field%transp_lm2r(field_LMloc_container, &
                       &                      field_Rloc_container)
               end if
               call comm_counter%stop_count(l_increment=.false.)

               call rLoop_counter%start_count()
               call radialLoopG(l_graph, l_frame,time,timeStage,tscheme,           &
                    &           dtLast,lTOCalc,lTONext,lTONext2,lHelCalc,          &
                    &           lPowerCalc,lRmsCalc,lPressCalc,                    &
                    &           lViscBcCalc,lFluxProfCalc,lperpParCalc,l_probe_out,&
                    &           dsdt_Rloc,dwdt_Rloc,dzdt_Rloc,dpdt_Rloc,dxidt_Rloc,&
                    &           dbdt_Rloc,djdt_Rloc,dVxVhLM_Rloc,dVxBhLM_Rloc,     &
                    &           dVSrLM_Rloc,dVXirLM_Rloc,                          &
                    &           lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb,  &
                    &           br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,HelLMr_Rloc,&
                    &           Hel2LMr_Rloc,HelnaLMr_Rloc,Helna2LMr_Rloc,         &
                    &           viscLMr_Rloc,uhLMr_Rloc,duhLMr_Rloc,gradsLMr_Rloc, &
                    &           fconvLMr_Rloc,fkinLMr_Rloc,fviscLMr_Rloc,          &
                    &           fpoynLMr_Rloc,fresLMr_Rloc,EperpLMr_Rloc,          &
                    &           EparLMr_Rloc,EperpaxiLMr_Rloc,EparaxiLMr_Rloc,     &
                    &           dtrkc_Rloc,dthkc_Rloc)
               call rLoop_counter%stop_count()
               phy2lm_counter%n_counts=phy2lm_counter%n_counts+1
               lm2phy_counter%n_counts=lm2phy_counter%n_counts+1
               nl_counter%n_counts=nl_counter%n_counts+1
               td_counter%n_counts=td_counter%n_counts+1

               if ( lVerbose ) write(output_unit,*) '! r-loop finished!'

               !---------------------------------------
               !- MPI transposition from r-distributed to LM-distributed
               !---------------------------------------
               if ( lVerbose ) write(output_unit,*) "! start r2lo redistribution"

               call comm_counter%start_count()
               PERFON('r2lo_dst')
               if ( l_conv .or. l_mag_kin ) then
                  call r2lo_flow%transp_r2lm(dflowdt_Rloc_container,  &
                       &             dflowdt_LMloc_container(:,:,:,tscheme%istage))
               end if

               if ( l_heat ) then
                  call r2lo_s%transp_r2lm(dsdt_Rloc_container,&
                       &             dsdt_LMloc_container(:,:,:,tscheme%istage))
               end if

               if ( l_chemical_conv ) then
                  call r2lo_xi%transp_r2lm(dxidt_Rloc_container, &
                       &             dxidt_LMloc_container(:,:,:,tscheme%istage))
               end if

               if ( l_mag ) then
                  call r2lo_field%transp_r2lm(dbdt_Rloc_container, &
                       &             dbdt_LMloc_container(:,:,:,tscheme%istage))
               end if
               call comm_counter%stop_count()

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
               if ( lVerbose ) write(output_unit,*) "! r2lo redistribution finished"


               !------ Nonlinear magnetic boundary conditions:
               !       For stressfree conducting boundaries
               PERFON('nl_m_bnd')
               if ( l_b_nl_cmb ) then
                  b_nl_cmb(1) =(1.0_cp,1.0_cp)
                  aj_nl_cmb(1)=(1.0_cp,1.0_cp)
                  call get_b_nl_bcs('CMB', br_vt_lm_cmb,br_vp_lm_cmb,              &
                       &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_cmb(2:lm_max))
               end if
               if ( l_b_nl_icb ) then
                  aj_nl_icb(1)=(1.0_cp,1.0_cp)
                  call get_b_nl_bcs('ICB', br_vt_lm_icb,br_vp_lm_icb,              &
                       &            2,lm_max,b_nl_cmb(2:lm_max),aj_nl_icb(2:lm_max))
               end if
               PERFOFF

               !---------------
               ! Finish assembing the explicit terms
               !---------------
               call lmLoop_counter%start_count()
               call finish_explicit_assembly(omega_ic,w_LMloc,b_ic_LMloc,       &
                    &                        aj_ic_LMloc,                       &
                    &                        dVSrLM_LMLoc(:,:,tscheme%istage),  &
                    &                        dVXirLM_LMLoc(:,:,tscheme%istage), &
                    &                        dVxVhLM_LMloc(:,:,tscheme%istage), &
                    &                        dVxBhLM_LMloc(:,:,tscheme%istage), &
                    &                        dsdt, dxidt, dwdt, djdt, dbdt_ic,  &
                    &                        djdt_ic, tscheme)
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
               call output(time,tscheme,n_time_step,l_stop_time,l_pot,l_log,     &
                    &      l_graph,lRmsCalc,l_store,l_new_rst_file,              &
                    &      l_spectrum,lTOCalc,lTOframe,lTOZwrite,                &
                    &      l_frame,n_frame,l_cmb,n_cmb_sets,l_r,                 &
                    &      lorentz_torque_ic,lorentz_torque_ma,dbdt_CMB_LMloc,   &
                    &      HelLMr_Rloc,Hel2LMr_Rloc,HelnaLMr_Rloc,Helna2LMr_Rloc,&
                    &      viscLMr_Rloc,uhLMr_Rloc,duhLMr_Rloc,gradsLMr_Rloc,    &
                    &      fconvLMr_Rloc,fkinLMr_Rloc,fviscLMr_Rloc,             &
                    &      fpoynLMr_Rloc,fresLMr_Rloc,EperpLMr_Rloc,EparLMr_Rloc,&
                    &      EperpaxiLMr_Rloc,EparaxiLMr_Rloc)
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
            call start_from_another_scheme(l_bridge_step, n_time_step, tscheme)

            !---------------
            !-- LM Loop (update routines)
            !---------------
            if ( lVerbose ) write(output_unit,*) '! starting lm-loop!'
            call lmLoop_counter%start_count()
            call LMLoop(timeStage,tscheme,lMat,lRmsNext,lPressNext,dsdt,  &
                 &      dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dbdt_ic,djdt_ic,   &
                 &      lorentz_torque_ma,lorentz_torque_ic,b_nl_cmb,     &
                 &      aj_nl_cmb,aj_nl_icb)

            if ( lVerbose ) write(output_unit,*) '! lm-loop finished!'

            !-- Timer counters
            call lmLoop_counter%stop_count()
            if ( tscheme%istage == 1 .and. lMat ) l_mat_time=.true.
            if (  tscheme%istage == 1 .and. .not. lMat .and. &
            &     .not. l_log ) l_pure=.true.

            ! Increment current stage
            tscheme%istage = tscheme%istage+1
         end do

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
               write(*,'(1p,/,/,A,i10,3(/,A,ES16.6))')                    &
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
               write(*,'(1p,/,/,A,i10,3(/,A,ES16.6))')          &
               &     " !  No of stored movie frames: ",n_frame, &
               &     " !     starting at time: ",0.0_cp,        &
               &     " !       ending at time: ",0.0_cp,        &
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
   subroutine start_from_another_scheme(l_bridge_step, n_time_step, tscheme)

      logical,             intent(in) :: l_bridge_step
      integer,             intent(in) :: n_time_step
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
            call get_pol_rhs_imp(s_LMloc, xi_LMloc, w_LMloc, dw_LMloc, ddw_LMloc,  &
                 &               p_LMloc, dp_LMloc, dwdt, dpdt, tscheme, 1,        &
                 &               .true., .false., .false., work_LMloc)
            if ( l_heat ) call get_entropy_rhs_imp(s_LMloc, ds_LMloc, dsdt, 1, .true.)
         end if

         call get_tor_rhs_imp(z_LMloc, dz_LMloc, dzdt, domega_ma_dt, domega_ic_dt, &
              &               omega_ic, omega_ma, omega_ic1, omega_ma1, tscheme,   &
              &               1, .true., .false.)

         if ( l_chemical_conv ) call get_comp_rhs_imp(xi_LMloc,dxi_LMloc,    &
                                     &                dxidt, 1, .true.)

         if ( l_mag ) call get_mag_rhs_imp(b_LMloc, db_LMloc, ddb_LMLoc,       &
                           &               aj_LMLoc, dj_LMloc, ddj_LMloc,      &
                           &               dbdt, djdt, tscheme, 1, .true.,     &
                           &               .false.)

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
end module step_time_mod
