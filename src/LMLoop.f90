#include "perflib_preproc.cpp"
module LMLoop_mod

   use iso_fortran_env, only: output_unit
   use fields
   use precision_mod
   use parallel_mod
   use constants, only: one
   use useful, only: abortRun, logWrite
   use num_param, only: solve_counter
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag, n_r_ic_max, lm_max, &
       &                 lm_maxMag
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
   use blocking, only: lo_map, llm, ulm, llmMag, ulmMag, st_map
   use logic, only: l_mag, l_conv, l_heat, l_single_matrix, l_double_curl, &
       &            l_chemical_conv, l_cond_ic, l_update_s, l_z10mat,      &
       &            l_parallel_solve, l_mag_par_solve, l_phase_field, l_onset
   use time_array, only: type_tarray, type_tscalar
   use time_schemes, only: type_tscheme
   use timing, only: timer_type
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod
   use updatePhi_mod

   implicit none

   private

   integer :: n_tri, n_penta ! Number of tridiagonal and pentadiagonal solvers
   integer :: block_sze, n_requests, nblocks
   integer, allocatable :: array_of_requests(:)

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop, finish_explicit_assembly, &
   &         assemble_stage, finish_explicit_assembly_Rdist, LMLoop_Rdist,         &
   &         test_LMLoop, assemble_stage_Rdist

contains

   subroutine initialize_LMLoop(tscheme)
      !
      ! This subroutine handles the memory allocation of the matrices needed
      ! in the time advance of MagIC
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variable
      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated
      if ( l_single_matrix ) then
         call initialize_updateWPS()
      else
         if ( l_heat ) call initialize_updateS()
         call initialize_updateWP(tscheme)
      end if

      if ( l_chemical_conv ) call initialize_updateXi()

      if ( l_phase_field ) then
         call initialize_updatePhi()
      else
         allocate( phi_ghost(1,1) ) ! For debugging only
      end if

      call initialize_updateZ()
      if ( l_mag ) call initialize_updateB()
      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('LMLoop.f90',local_bytes_used)

      if ( l_parallel_solve ) then
         if ( l_conv ) then
            n_tri  =1 ! z equation
            n_penta=1 ! w equation
         end if
         if ( l_heat ) n_tri = n_tri+1
         if ( l_chemical_conv ) n_tri = n_tri+1
         if ( l_mag_par_solve ) n_tri = n_tri+2

         block_sze=50
         n_requests=10
         nblocks = lm_max
         nblocks = set_block_number(nblocks)
         allocate( array_of_requests(n_requests))

      end if

   end subroutine initialize_LMLoop
!----------------------------------------------------------------------------
   subroutine test_LMLoop(tscheme)
      !
      ! This subroutine is used to solve dummy linear problem to estimate the best
      ! blocking size. This is done once at the initialisation stage of MagIC.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      !-- Local variable
      type(type_tarray) :: dummy
      type(type_tscalar) :: dum_scal

      lWPmat(:)=.false.
      if ( l_heat ) lSmat(:) =.false.
      lZmat(:) =.false.
      if ( l_mag ) lBmat(:) =.false.
      if ( l_chemical_conv ) lXimat(:)=.false.

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      call dum_scal%initialize(tscheme%nold, tscheme%nexp, tscheme%nimp)
      call dummy%initialize(1, lm_max, nRstart, nRstop, tscheme%nold, tscheme%nexp,&
           &                tscheme%nimp, l_allocate_exp=.true.)

      if ( l_heat ) call prepareS_FD(tscheme, dummy, phi_Rloc)
      if ( l_chemical_conv ) call prepareXi_FD(tscheme, dummy)
      if ( l_conv ) then
         call prepareZ_FD(0.0_cp, tscheme, dummy, omega_ma, omega_ic, dum_scal, &
              &           dum_scal)
         call prepareW_FD(tscheme, dummy, .false.)
      end if
      if ( l_mag_par_solve ) call prepareB_FD(0.0_cp, tscheme, dummy, dummy)

      call find_faster_block() ! Find the fastest blocking

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      call dum_scal%finalize()
      call dummy%finalize()

   end subroutine test_LMLoop
!----------------------------------------------------------------------------
   subroutine finalize_LMLoop(tscheme)
      !
      ! This subroutine deallocates the matrices involved in the time advance
      ! of MagIC.
      !

      !-- Input variables:
      class(type_tscheme), intent(in) :: tscheme ! time scheme

      if ( l_single_matrix ) then
         call finalize_updateWPS()
      else
         if ( l_heat ) call finalize_updateS()
         call finalize_updateWP(tscheme)
      end if
      call finalize_updateZ()

      if ( l_chemical_conv ) call finalize_updateXi()
      if ( l_phase_field ) then
         call finalize_updatePhi()
      else
         deallocate( phi_ghost )
      end if
      if ( l_mag ) call finalize_updateB()
      if ( l_parallel_solve ) deallocate(array_of_requests)

   end subroutine finalize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(time,timeNext,tscheme,lMat,lRmsNext,lPressNext,     &
              &      dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,djdt,dbdt_ic, &
              &      djdt_ic,domega_ma_dt,domega_ic_dt,                  &
              &      lorentz_torque_ma_dt,lorentz_torque_ic_dt,          &
              &      b_nl_cmb,aj_nl_cmb,aj_nl_icb)
      !
      !  This subroutine performs the actual time-stepping. It calls succesively
      !  the update routines of the various fields.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      real(cp),            intent(in) :: timeNext
      logical,             intent(in) :: lMat
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext
      complex(cp),         intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp),         intent(in)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp),         intent(in)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

      !--- Input from radialLoop:
      type(type_tarray),  intent(inout) :: dsdt, dxidt, dwdt, dpdt, dzdt, dphidt
      type(type_tarray),  intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tscalar), intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt
      !integer,     intent(in) :: n_time_step

      !--- Inner core rotation from last time step
      real(cp) :: z10(n_r_max)

      PERFON('LMloop')
      if ( lMat ) then ! update matrices:
         !---- The following logicals tell whether the respective inversion
         !     matrices have been updated. lMat=.true. when a general
         !     update is necessary.
         lZ10mat=.false.
         if ( l_single_matrix ) then
            lWPSmat(:)=.false.
         else
            lWPmat(:)=.false.
            if ( l_heat ) lSmat(:) =.false.
         end if
         lZmat(:) =.false.
         if ( l_mag ) lBmat(:) =.false.
         if ( l_chemical_conv ) lXimat(:)=.false.
         if ( l_phase_field ) lPhimat(:)=.false.
      end if

      if ( l_phase_field ) call updatePhi(phi_LMloc, dphidt, tscheme)

      if ( l_heat .and. .not. l_single_matrix ) then
         PERFON('up_S')
         call updateS( s_LMloc, ds_LMloc, dsdt, phi_LMloc, tscheme )
         PERFOFF
      end if

      if ( l_chemical_conv ) call updateXi(xi_LMloc, dxi_LMloc, dxidt, tscheme)

      if ( l_conv ) then
         PERFON('up_Z')
         call updateZ( time, timeNext, z_LMloc, dz_LMloc, dzdt, omega_ma,  &
              &        omega_ic, domega_ma_dt,domega_ic_dt,                &
              &        lorentz_torque_ma_dt,lorentz_torque_ic_dt, tscheme, &
              &        lRmsNext)
         PERFOFF

         if ( l_single_matrix ) then
            if ( rank == rank_with_l1m0 ) then
               z10(:)=real(z_LMloc(lo_map%lm2(1,0),:))
            end if
#ifdef WITH_MPI
            if ( rank_with_l1m0 >= 0 ) then ! This is -1 if m_min > 0
               call MPI_Bcast(z10,n_r_max,MPI_DEF_REAL,rank_with_l1m0, &
                    &         MPI_COMM_WORLD,ierr)
            end if
#endif
            call updateWPS( w_LMloc, dw_LMloc, ddw_LMloc, z10, dwdt,    &
                 &          p_LMloc, dp_LMloc, dpdt, s_LMloc, ds_LMloc, &
                 &          dsdt, tscheme, lRmsNext )
         else
            PERFON('up_WP')
            call updateWP( s_LMloc, xi_LMLoc, w_LMloc, dw_LMloc, ddw_LMloc, &
                 &         dwdt, p_LMloc, dp_LMloc, dpdt, tscheme,          &
                 &         lRmsNext, lPressNext )
            PERFOFF
         end if
      end if
      if ( l_mag ) then ! dwdt,dpdt used as work arrays
         PERFON('up_B')
         call updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc, &
              &        dbdt, djdt, b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc,      &
              &        aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, dbdt_ic,        &
              &        djdt_ic, b_nl_cmb, aj_nl_cmb, aj_nl_icb, time, tscheme, &
              &        lRmsNext )
         PERFOFF
      end if

      PERFOFF
   end subroutine LMLoop
!--------------------------------------------------------------------------------
   subroutine LMLoop_Rdist(time,timeNext,tscheme,lMat,lRmsNext,lPressNext,    &
              &            lP00Next,dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,    &
              &            djdt,dbdt_ic,djdt_ic,domega_ma_dt,domega_ic_dt,    &
              &            lorentz_torque_ma_dt,lorentz_torque_ic_dt,         &
              &            b_nl_cmb,aj_nl_cmb,aj_nl_icb)
      !
      !  This subroutine performs the actual time-stepping. It calls succesively
      !  the update routines of the various fields. This is used with the parallel
      !  finite difference solver.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      real(cp),            intent(in) :: timeNext
      logical,             intent(in) :: lMat
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lP00Next ! Do wee need p00 pressure on next log
      complex(cp),         intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp),         intent(in)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp),         intent(in)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB
      !--- Input from radialLoop:
      type(type_tarray),  intent(inout) :: dsdt, dxidt, dwdt, dpdt, dzdt, dphidt
      type(type_tarray),  intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tscalar), intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      !-- Local variable
      logical :: lPress

      lPress = lPressNext .or. lP00Next

      if ( lMat ) then ! update matrices:
         lZ10mat=.false.
         lWPmat(:)=.false.
         if ( l_heat ) lSmat(:) =.false.
         lZmat(:) =.false.
         if ( l_mag ) lBmat(:) =.false.
         if ( l_chemical_conv ) lXimat(:)=.false.
         if ( l_phase_field ) lPhimat(:)=.false.
      end if

      !-- Phase field needs to be computed first on its own to allow a proper
      !-- advance of temperature afterwards
      if ( l_phase_field ) then
         call preparePhase_FD(tscheme, dphidt)
         call parallel_solve_phase(block_sze)
         call fill_ghosts_Phi(phi_ghost)
         call updatePhase_FD(phi_Rloc, dphidt, tscheme)
      end if

      !-- Mainly assemble the r.h.s. and rebuild the matrices if required
      if ( l_heat ) call prepareS_FD(tscheme, dsdt, phi_Rloc)
      if ( l_chemical_conv ) call prepareXi_FD(tscheme, dxidt)
      if ( l_conv ) then
         call prepareZ_FD(time, tscheme, dzdt, omega_ma, omega_ic, domega_ma_dt, &
              &           domega_ic_dt)
         if ( l_z10mat ) call z10Mat_FD%solver_single(z10_ghost, nRstart, nRstop)
         call prepareW_FD(tscheme, dwdt, lPress)
         if ( lPress ) call p0Mat_FD%solver_single(p0_ghost, nRstart, nRstop)
      end if
      if ( l_mag_par_solve ) call prepareB_FD(time, tscheme, dbdt, djdt)

      !-----------------------------------------------------------
      !--- This is where the matrices are solved
      !-- Here comes the real deal:
      call solve_counter%start_count()
      call parallel_solve(block_sze)
      call solve_counter%stop_count()
      !-----------------------------------------------------------

      !-- Copy z10 into z_ghost after solving when needed
      if ( l_z10Mat ) z_ghost(st_map%lm2(1,0),:)=cmplx(real(z10_ghost(:)),0.0_cp,cp)

      !-- Now simply fill the ghost zones to ensure the boundary conditions
      if ( l_heat ) call fill_ghosts_S(s_ghost)
      if ( l_chemical_conv ) call fill_ghosts_Xi(xi_ghost)
      if ( l_conv ) then
         call fill_ghosts_Z(z_ghost)
         call fill_ghosts_W(w_ghost, p0_ghost, lPress)
      end if
      if ( l_mag_par_solve ) call fill_ghosts_B(b_ghost, aj_ghost)

      !-- Finally build the radial derivatives and the arrays for next iteration
      if ( l_heat ) call updateS_FD(s_Rloc, ds_Rloc, dsdt, phi_Rloc, tscheme)
      if ( l_chemical_conv ) call updateXi_FD(xi_Rloc, dxidt, tscheme)

      call updateZ_FD(time, timeNext, z_Rloc, dz_Rloc, dzdt, omega_ma, omega_ic, &
           &          domega_ma_dt, domega_ic_dt, lorentz_torque_ma_dt,          &
           &          lorentz_torque_ic_dt, tscheme, lRmsNext)
      call updateW_FD(w_Rloc, dw_Rloc, ddw_Rloc, dwdt, p_Rloc, dp_Rloc, dpdt, tscheme, &
           &          lRmsNext, lPressNext, lP00Next)

      if ( l_mag_par_solve ) then
         call updateB_FD(b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc, dj_Rloc, ddj_Rloc, dbdt, &
              &          djdt, tscheme, lRmsNext)
      end if

      if ( l_mag .and. (.not. l_mag_par_solve) ) then ! dwdt,dpdt used as work arrays
         call updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc, &
              &        dbdt, djdt, b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc,      &
              &        aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, dbdt_ic,        &
              &        djdt_ic, b_nl_cmb, aj_nl_cmb, aj_nl_icb, time, tscheme, &
              &        lRmsNext )
      end if

   end subroutine LMLoop_Rdist
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly(omega_ic, w, b_ic, aj_ic, dVSr_LMloc,      &
              &                        dVXir_LMloc, dVxVh_LMloc, dVxBh_LMloc,     &
              &                        lorentz_torque_ma, lorentz_torque_ic,      &
              &                        dsdt, dxidt, dwdt, djdt, dbdt_ic,          &
              &                        djdt_ic, domega_ma_dt, domega_ic_dt,       &
              &                        lorentz_torque_ma_dt, lorentz_torque_ic_dt,&
              &                        tscheme)
      !
      ! This subroutine is used to finish the computation of the explicit terms.
      ! This is only possible in a LM-distributed space since it mainly involves
      ! computation of radial derivatives.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: omega_ic
      real(cp),            intent(in) :: lorentz_torque_ic
      real(cp),            intent(in) :: lorentz_torque_ma
      complex(cp),         intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(inout) :: dVSr_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVXir_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVxVh_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVxBh_LMloc(llmMag:ulmMag,n_r_maxMag)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dsdt, dxidt, djdt, dwdt
      type(type_tarray),   intent(inout) :: dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      if ( l_chemical_conv ) then
         call finish_exp_comp(w, dVXir_LMloc, dxidt%expl(:,:,tscheme%istage))
      end if

      if ( l_single_matrix ) then
         call finish_exp_smat(dVSr_LMloc, dsdt%expl(:,:,tscheme%istage))
      else
         if ( l_heat ) then
            call finish_exp_entropy(w, dVSr_LMloc, dsdt%expl(:,:,tscheme%istage))
         end if
         if ( l_double_curl ) then
            call finish_exp_pol(dVxVh_LMloc, dwdt%expl(:,:,tscheme%istage))
         end if
      end if

      if ( .not. l_onset ) then
         call finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic,     &
              &              domega_ma_dt%expl(tscheme%istage),        &
              &              domega_ic_dt%expl(tscheme%istage),        &
              &              lorentz_torque_ma_dt%expl(tscheme%istage),&
              &              lorentz_torque_ic_dt%expl(tscheme%istage))
      end if

      if ( l_mag ) then
         call finish_exp_mag(dVxBh_LMloc, djdt%expl(:,:,tscheme%istage))
      end if

      if ( l_cond_ic ) then
         call finish_exp_mag_ic(b_ic, aj_ic, omega_ic,            &
              &                 dbdt_ic%expl(:,:,tscheme%istage), &
              &                 djdt_ic%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly_Rdist(omega_ic, w, b_ic, aj_ic, dVSr_Rloc, &
              &                        dVXir_Rloc, dVxVh_Rloc, dVxBh_Rloc,        &
              &                        lorentz_torque_ma, lorentz_torque_ic,      &
              &                        dsdt_Rloc, dxidt_Rloc, dwdt_Rloc,          &
              &                        djdt_Rloc, dbdt_ic, djdt_ic, domega_ma_dt, &
              &                        domega_ic_dt, lorentz_torque_ma_dt,        &
              &                        lorentz_torque_ic_dt,tscheme)
      !
      ! This subroutine is used to finish the computation of the explicit terms.
      ! This is the version that handles R-distributed arrays used when FD are
      ! employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: omega_ic
      real(cp),            intent(in) :: lorentz_torque_ic
      real(cp),            intent(in) :: lorentz_torque_ma
      complex(cp),         intent(in) :: w(lm_max,nRstart:nRstop)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(inout) :: dVSr_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVXir_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVxVh_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVxBh_Rloc(lm_maxMag,nRstartMag:nRstopMag)

      !-- Output variables
      complex(cp),         intent(inout) :: dxidt_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dsdt_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: dwdt_Rloc(lm_max,nRstart:nRstop)
      complex(cp),         intent(inout) :: djdt_Rloc(lm_max,nRstart:nRstop)
      type(type_tarray),   intent(inout) :: dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      if ( l_chemical_conv ) call finish_exp_comp_Rdist(w, dVXir_Rloc, dxidt_Rloc)

      if ( l_single_matrix ) then
         call finish_exp_smat_Rdist(dVSr_Rloc, dsdt_Rloc)
      else
         if ( l_heat ) call finish_exp_entropy_Rdist(w, dVSr_Rloc, dsdt_Rloc)
         if ( l_double_curl ) call finish_exp_pol_Rdist(dVxVh_Rloc, dwdt_Rloc)
      end if

      if ( .not. l_onset ) then
         call finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic,     &
              &              domega_ma_dt%expl(tscheme%istage),        &
              &              domega_ic_dt%expl(tscheme%istage),        &
              &              lorentz_torque_ma_dt%expl(tscheme%istage),&
              &              lorentz_torque_ic_dt%expl(tscheme%istage))
      end if

      if ( l_mag ) call finish_exp_mag_Rdist(dVxBh_Rloc, djdt_Rloc)

      if ( l_cond_ic ) then
         call finish_exp_mag_ic(b_ic, aj_ic, omega_ic,            &
              &                 dbdt_ic%expl(:,:,tscheme%istage), &
              &                 djdt_ic%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly_Rdist
!--------------------------------------------------------------------------------
   subroutine assemble_stage(time, omega_ic, omega_ic1, omega_ma, omega_ma1,        &
              &              dwdt, dzdt, dpdt, dsdt, dxidt, dphidt, dbdt, djdt,     &
              &              dbdt_ic, djdt_ic, domega_ic_dt, domega_ma_dt,          &
              &              lorentz_torque_ic_dt, lorentz_torque_ma_dt, lPressNext,&
              &              lRmsNext, tscheme)
      !
      ! This routine is used to call the different assembly stage of the different
      ! equations. This is only used for a special subset of IMEX-RK schemes that
      ! have ``tscheme%l_assembly=.true.``
      !

      !-- Input variables
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time

      !-- Output variables
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt
      real(cp),            intent(inout) :: omega_ic, omega_ma, omega_ic1, omega_ma1
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dpdt, dsdt, dxidt, dphidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic

      if ( l_chemical_conv )  call assemble_comp(xi_LMloc, dxi_LMloc, dxidt, tscheme)

      if ( l_phase_field )  call assemble_phase(phi_LMloc, dphidt, tscheme)

      if ( l_single_matrix ) then
         call assemble_single(s_LMloc, ds_LMloc, w_LMloc, dw_LMloc, ddw_LMloc, &
              &               dsdt, dwdt, dpdt, tscheme,lRmsNext)
      else
         if ( l_heat )  call assemble_entropy(s_LMloc, ds_LMloc, dsdt, phi_LMloc, &
                             &                tscheme)
         call assemble_pol(s_LMloc, xi_LMloc, w_LMloc, dw_LMloc, ddw_LMloc, p_LMloc, &
              &            dp_LMloc, dwdt, dpdt, dpdt%expl(:,:,1), tscheme,          &
              &            lPressNext, lRmsNext)
      end if

      call assemble_tor(time, z_LMloc, dz_LMloc, dzdt, domega_ic_dt, domega_ma_dt, &
           &            lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic,      &
           &            omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)

      if ( l_mag ) call assemble_mag(b_LMloc, db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc,  &
                        &            ddj_LMloc, b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc,  &
                        &            aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, dbdt, djdt,&
                        &            dbdt_ic, djdt_ic, lRmsNext, tscheme)

   end subroutine assemble_stage
!--------------------------------------------------------------------------------
   subroutine assemble_stage_Rdist(time, omega_ic, omega_ic1, omega_ma, omega_ma1,    &
              &                    dwdt, dzdt, dpdt, dsdt, dxidt, dphidt, dbdt,       &
              &                    djdt, dbdt_ic, djdt_ic, domega_ic_dt, domega_ma_dt,&
              &                    lorentz_torque_ic_dt, lorentz_torque_ma_dt,        &
              &                    lPressNext, lRmsNext, tscheme)
      !
      ! This routine is used to call the different assembly stage of the different
      ! equations. This is only used for a special subset of IMEX-RK schemes that
      ! have ``tscheme%l_assembly=.true.``
      !

      !-- Input variables
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time

      !-- Output variables
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt
      real(cp),            intent(inout) :: omega_ic, omega_ma, omega_ic1, omega_ma1
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dsdt, dxidt, dpdt, dphidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic

      if ( l_phase_field )  call assemble_phase_Rloc(phi_Rloc, dphidt, tscheme)
      if ( l_chemical_conv )  call assemble_comp_Rloc(xi_Rloc, dxidt, tscheme)
      if ( l_heat )  call assemble_entropy_Rloc(s_Rloc, ds_Rloc, dsdt, phi_Rloc, tscheme)

      call assemble_pol_Rloc(block_sze, nblocks, w_Rloc, dw_Rloc, ddw_Rloc, p_Rloc, &
           &                 dp_Rloc, dwdt, dpdt%expl(:,:,1), tscheme, lPressNext,  &
           &                 lRmsNext)

      call assemble_tor_Rloc(time, z_Rloc, dz_Rloc, dzdt, domega_ic_dt, domega_ma_dt, &
           &                 lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic, &
           &                 omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)

      if ( l_mag ) then
         if ( l_mag_par_solve ) then
            call assemble_mag_Rloc(b_Rloc, db_Rloc, ddb_Rloc, aj_Rloc, dj_Rloc,   &
                 &                 ddj_Rloc, dbdt, djdt, lRmsNext, tscheme)
         else
            call assemble_mag(b_LMloc, db_LMloc, ddb_LMloc, aj_LMloc, dj_LMloc,   &
                 &            ddj_LMloc, b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc,   &
                 &            aj_ic_LMloc, dj_ic_LMloc, ddj_ic_LMloc, dbdt, djdt, &
                 &            dbdt_ic, djdt_ic, lRmsNext, tscheme)
         end if
      end if

   end subroutine assemble_stage_Rdist
!--------------------------------------------------------------------------------
   subroutine parallel_solve_phase(block_sze)
      !
      ! This subroutine handles the parallel solve of the phase field matrices.
      ! This needs to be updated before the temperature.
      !
      integer, intent(in) :: block_sze ! Size ot the LM blocks

      !-- Local variables
      integer :: req
      integer :: start_lm, stop_lm, tag, nlm_block, lms_block

#ifdef WITH_MPI
      array_of_requests(:)=MPI_REQUEST_NULL
#endif
      !$omp parallel default(shared) private(tag, req, start_lm, stop_lm, nlm_block, lms_block)
      tag = 0
      req=1
      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier

         call phiMat_FD%solver_up(phi_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
              &                   array_of_requests, req, lms_block, nlm_block)
         tag = tag+1
      end do

      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier

         call phiMat_FD%solver_dn(phi_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
           &                      array_of_requests, req, lms_block, nlm_block)
         tag = tag+1
      end do

      !$omp master
      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze

         call phiMat_FD%solver_finish(phi_ghost, lms_block, nlm_block, nRstart, &
              &                       nRstop, tag, array_of_requests, req)
         tag = tag+1
      end do

#ifdef WITH_MPI
      call MPI_Waitall(req-1, array_of_requests(1:req-1), MPI_STATUSES_IGNORE, ierr)
      if ( ierr /= MPI_SUCCESS ) call abortRun('MPI_Waitall failed in LMLoop')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      !$omp end master
      !$omp barrier

      !$omp end parallel

   end subroutine parallel_solve_phase
!--------------------------------------------------------------------------------
   subroutine parallel_solve(block_sze)
      !
      ! This subroutine handles the parallel solve of the time-advance matrices.
      ! This works with R-distributed arrays (finite differences).
      !
      integer, intent(in) :: block_sze ! Size ot the LM blocks

      !-- Local variables
      integer :: req
      integer :: start_lm, stop_lm, tag, nlm_block, lms_block

#ifdef WITH_MPI
      array_of_requests(:)=MPI_REQUEST_NULL
#endif

      !$omp parallel default(shared) private(tag, req, start_lm, stop_lm, nlm_block, lms_block)
      tag = 0
      req=1

      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier

         if ( l_heat ) then
            call sMat_FD%solver_up(s_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if

         if ( l_chemical_conv ) then
            call xiMat_FD%solver_up(xi_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                  array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if

         if ( l_conv ) then
            call zMat_FD%solver_up(z_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
            call wMat_FD%solver_up(w_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+2
         end if

         if ( l_mag_par_solve ) then
            call bMat_FD%solver_up(b_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
            call jMat_FD%solver_up(aj_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if
      end do

      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze
         start_lm=lms_block; stop_lm=lms_block+nlm_block-1
         call get_openmp_blocks(start_lm,stop_lm)
         !$omp barrier

         if ( l_heat ) then
            call sMat_FD%solver_dn(s_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
              &                    array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if

         if ( l_chemical_conv ) then
            call xiMat_FD%solver_dn(xi_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                  array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if

         if ( l_conv ) then
            call zMat_FD%solver_dn(z_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
            call wMat_FD%solver_dn(w_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+2
         end if

         if ( l_mag_par_solve ) then
            call bMat_FD%solver_dn(b_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
            call jMat_FD%solver_dn(aj_ghost, start_lm, stop_lm, nRstart, nRstop, tag, &
                 &                 array_of_requests, req, lms_block, nlm_block)
            tag = tag+1
         end if
      end do

      !$omp master
      do lms_block=1,lm_max,block_sze
         nlm_block = lm_max-lms_block+1
         if ( nlm_block > block_sze ) nlm_block=block_sze

         if ( l_heat ) then
            call sMat_FD%solver_finish(s_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+1
         end if

         if ( l_chemical_conv ) then
            call xiMat_FD%solver_finish(xi_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                      tag, array_of_requests, req)
            tag = tag+1
         end if

         if ( l_conv ) then
            call zMat_FD%solver_finish(z_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+1

            call wMat_FD%solver_finish(w_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+2
         end if

         if ( l_mag_par_solve ) then
            call bMat_FD%solver_finish(b_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+1
            call jMat_FD%solver_finish(aj_ghost, lms_block, nlm_block, nRstart, nRstop, &
                 &                     tag, array_of_requests, req)
            tag = tag+1
         end if
      end do

#ifdef WITH_MPI
      call MPI_Waitall(req-1, array_of_requests(1:req-1), MPI_STATUSES_IGNORE, ierr)
      if ( ierr /= MPI_SUCCESS ) call abortRun('MPI_Waitall failed in LMLoop')
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      !$omp end master
      !$omp barrier

      !$omp end parallel

   end subroutine parallel_solve
!--------------------------------------------------------------------------------
   integer function set_block_number(nb) result(nbo)
      !
      ! This routine returns the number of lm-blocks for solving the LM loop.
      ! This is adapted from xshells.
      !
      integer, intent(inout) :: nb ! Number of blocks

      !-- Local variable
      integer :: nlm_4

      nlm_4=(lm_max+3)/4 ! Don't split cache lines
      if ( nb < 1 ) nb=1
      if ( (nlm_4+nb-1)/nb < 16) nb=(nlm_4+15)/16 ! Size not less than 1KB (64 cplx)
      if ( nb > 2048 ) nb=2048 ! No more than 2048 blocks
      block_sze=((nlm_4+nb-1)/nb)*4 ! Block size (rounded)

      nblocks = (lm_max+block_sze-1)/block_sze
      n_requests = ( n_tri*4+n_penta*8)*nblocks ! Max number of MPI request per block
      nbo = nb

      if ( allocated(array_of_requests) ) then
         deallocate(array_of_requests)
         allocate(array_of_requests(n_requests))
      end if

   end function set_block_number
!--------------------------------------------------------------------------------
   subroutine find_faster_block
      !
      ! This routine is used to find the best lm-block size for MPI communication.
      ! This is adapted from xshells.
      !
      real(cp), parameter :: alpha=0.61803_cp
      real(cp), parameter :: beta=one-alpha
      integer :: bmin, bmax, nb, nblk, b(0:2), nloops, bb
      real(cp) :: t(0:2), mul, tt
      character(len=72) :: str

      call logWrite('')
      call logWrite('! Time the LM Loop to get the optimal blocking:')

      t(:)=0.0_cp; b(:)=1
      nloops=1
      bmin = 1 ! Minimum number of blocks
      nb = lm_max
      bmax =  set_block_number(nb) ! Maximum number of blocks
      nb = (n_procs+1)/2 ! Start with a default block number
      nblk = set_block_number(nb)

      if ( n_procs > 1 ) then
         nb=int(bmin**beta * bmax**alpha)
         nblk = set_block_number(nb)
      end if
      b(1) = set_block_number(nblk)

      !-- First estimate the number of loops to get accurate timing
      !-- To do so measure the walltime until this is stable to a 2% accuracy
      t(1) = time_LM_loop(b(1), nloops) ! Starting value
      tt=t(1)+one
      do while( abs(t(1)-tt) > 0.02*abs(t(1)) ) ! 2% difference
         nloops = 2*nloops
         tt = t(1)
         t(1) = time_LM_loop(b(1), nloops)
         if ( t(1)*(nloops/2) > 10 ) exit ! Longer than 10sec
      end do
      nloops=nloops/2
      nloops=max(nloops,2) ! At least 2 loops
      write(str,'(A,I5,A)') '! I am computing', nloops, ' calls'
      call logWrite(str)

      if ( n_procs > 1) then
         if ( bmax*bmin > b(1)*b(1) ) then
            mul = real(bmax,cp)
         else
            mul = real(bmin,cp)
         end if
         mul = (mul/b(1))**beta
         nb = int(b(1)*mul)
         b(0)=set_block_number(nb)
         t(0)=time_LM_loop(b(0), nloops)
         if ( t(0) < t(1) ) then ! Exchange indices
            bb=b(0); b(0)=b(1); b(1)=bb
            tt=t(0); t(0)=t(1); t(1)=tt
         end if

         if ( b(1) > b(0) ) then
            mul = real(bmax,cp)
         else
            mul = real(bmin,cp)
         end if
         mul = (mul/b(1))**beta
         b(2)=b(1)

         do while ( (b(1) > bmin) .and. ( b(1) < bmax) )
            nb=int(b(2)*mul)
            b(2)=set_block_number(nb)
            t(2)=time_LM_loop(b(2), nloops)
            if ( t(2) < t(1) ) then
               if ( t(1) > t(2)*1.02) then ! If we are within 2% of the minimum time
                  !-- Change the interval boundaries
                  t(0)=t(1); b(0)=b(1)
               end if
               t(1)=t(2); b(1)=b(2)
            end if
            if ( (b(2)==bmin) .or. (b(2)==bmax) .or. (t(2)>t(1)*1.02 ) ) exit
         end do
         write(str,'(A,1x,I0,A,1x,I0)') '! I am braketing the block number between', &
         &                              b(0), ' and', b(2)
         call logWrite(str)

         !-- At this stage, the minimum wall time has been stored in t(1) and is bracketed
         !-- between in t(0) and t(2)
         !if ( rank == 0 ) print*, 'b', 't', b, t

         !-- Determine the largest interval
         if ( abs(log(real(b(2)))-log(real(b(0)))) > abs(log(real(b(0)))-log(real(b(1)))) ) then
            nb=int(b(1)**alpha*b(2)**beta)
            bb=set_block_number(nb)
            tt=time_LM_loop(bb, nloops)
         else
            bb=b(1); tt=t(1)
            nb=int(b(1)**alpha*b(0)**beta)
            t(1)=time_LM_loop(b(1), nloops)
         end if

         !-- Refined block determination
         if ( rank == 0 ) write(output_unit,*) ''
         do while( (t(2)>min(t(1),tt)*1.02) .and. (t(0)>min(t(1),tt)*1.02) &
         &          .and. (abs(b(2)-b(0))>1) .and. (maxval(b) < 4) )
            if ( tt < t(1) ) then ! This is better than before
               t(0)=t(1); b(0)=b(1);
               t(1)=tt; b(1)=bb
               nb=int(b(1)**alpha*b(2)**beta)
               bb=set_block_number(nb)
               tt=time_LM_loop(bb, nloops)
            else
               t(2)=tt; b(2)=bb
               tt=t(1); bb=b(1)
               nb=int(b(1)**alpha*b(0)**beta)
               b(1)=set_block_number(nb)
               t(1)=time_LM_loop(b(1), nloops)
            end if
            if ( rank==0 ) write(output_unit, '(A,I0,1x,I0,1x,I0,1x,I0,A,4ES9.2,A)')  &
            &              'Searching for number of blocks b=(',b(0), b(1), bb, b(2), &
            &              '), t=(',t(0), t(1), tt, t(2),')'
         end do
         if ( rank == 0 ) write(output_unit,*) ''

         if ( tt < t(1) ) t(1)=tt; b(1)=bb
#ifdef WITH_MPI
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

         !- At this stage this should be done: b(1) should be the optimal block_number
         nb=b(1)
         nblk=set_block_number(nb)
      end if

      !print*, 'b,t=', b, t
      if ( nblk == 1 ) then
         write(str,'(A,1x,I0,A,ES9.2)') '! Found 1 block of size',  &
         &                                 block_sze, ', timing was:', t(1)
      else
         write(str,'(A,1x,I0,A,1x,I0,A,ES9.2)') '! Found', nblk, ' blocks of size',  &
         &                                      block_sze, ', timing was:', t(1)
      end if
      call logWrite(str)

   end subroutine find_faster_block
!--------------------------------------------------------------------------------
   real(cp) function time_LM_loop(nblk, nloops) result(time)
      !
      ! This function is defined to measure the time of a call to the parallel solvers
      ! when a a number of block nblk is used. This takes the average over nloop calls.
      !

      integer, intent(inout) :: nblk ! Number of lm blocks
      integer, intent(in) :: nloops ! Number of calls

      !-- Local variables
      type(timer_type) :: tcount
      integer :: n_l

      time = 0.0_cp
      nblk = set_block_number(nblk)
      call tcount%initialize()

#ifdef WITH_MPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      do n_l=1,nloops
         call tcount%start_count()
         call parallel_solve(block_sze)
         call tcount%stop_count()
      end do
      call tcount%finalize()

      time=tcount%tTot

   end function time_LM_loop
!--------------------------------------------------------------------------------
end module LMLoop_mod
