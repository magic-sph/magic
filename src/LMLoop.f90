#include "perflib_preproc.cpp"
module LMLoop_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use precision_mod
   use parallel_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag, n_r_ic_max
   use blocking, only: lo_map, llm, ulm, llmMag, ulmMag
   use logic, only: l_mag, l_conv, l_heat, l_single_matrix, l_double_curl, &
       &            l_chemical_conv, l_cond_ic
   use time_array, only: type_tarray, type_tscalar
   use time_schemes, only: type_tscheme
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod
   
   ! DELETEMEEEE
   use communications
   use mpi_thetap_mod

   implicit none

   private

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop, finish_explicit_assembly, &
   &         assemble_stage

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
         if ( l_heat ) call initialize_updateS_dist()
         call initialize_updateWP(tscheme)
      end if

      if ( l_chemical_conv ) call initialize_updateXi()
      if ( l_chemical_conv ) call initialize_updateXi_dist()

      call initialize_updateZ()
      if ( l_mag ) call initialize_updateB()
      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('LMLoop.f90',local_bytes_used)

   end subroutine initialize_LMLoop
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
      if ( l_mag ) call finalize_updateB()

   end subroutine finalize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(time,timeNext,tscheme,lMat,lRmsNext,lPressNext,   &
              &      dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dbdt_ic,      &
              &      djdt_ic,dsdt_dist,dwdt_dist,dzdt_dist,dpdt_dist,  &
              &      dxidt_dist,dbdt_dist,djdt_dist,dbdt_ic_dist,      &
              &      djdt_ic_dist,domega_ma_dt,domega_ic_dt,           &
              &      lorentz_torque_ma_dt,lorentz_torque_ic_dt,        &
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
      type(type_tarray),  intent(inout) :: dsdt, dxidt, dwdt, dpdt, dzdt
      type(type_tarray),  intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic
      type(type_tarray),  intent(inout) :: dsdt_dist, dxidt_dist, dwdt_dist
      type(type_tarray),  intent(inout) :: dpdt_dist, dzdt_dist
      type(type_tarray),  intent(inout) :: dbdt_dist, djdt_dist, dbdt_ic_dist, djdt_ic_dist
      type(type_tscalar), intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt
      !integer,     intent(in) :: n_time_step

      !--- Inner core rotation from last time step
      real(cp) :: z10(n_r_max)

      PERFON('LMloop')
      !LIKWID_ON('LMloop')
      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin of porting point
      call transform_old2new(s_LMloc, s_LMdist)
      call dsdt%slice_all(dsdt_dist) ! <-- needed otherwise the expl array is not sliced
      !call test_field(s_LMdist, s_LMloc, 'entropy_')

      if ( l_chemical_conv ) then
         call transform_old2new(xi_LMloc, xi_LMdist)
         call dxidt%slice_all(dxidt_dist)
         !call test_field(xi_LMdist, xi_LMloc, 'comp_')
      end if
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin of porting point

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
            if ( l_heat ) lSmat_dist(:) =.false.
         end if
         lZmat(:) =.false.
         if ( l_mag ) lBmat(:) =.false.
         if ( l_chemical_conv ) lXimat(:)=.false.
      end if

      if ( l_heat .and. .not. l_single_matrix ) then
         PERFON('up_S')
         call updateS( s_LMloc, ds_LMloc, dsdt, tscheme )
         call updateS_dist( s_LMdist, ds_LMdist, dsdt_dist, tscheme )
         
         call test_field(s_LMdist, s_LMloc, 'entropy_')
         call test_field(ds_LMdist, ds_LMloc, 'dentropydr_')
         PERFOFF
      end if
      
      if ( l_chemical_conv ) then
         call updateXi(xi_LMloc, dxi_LMloc, dxidt, tscheme)
         call updateXi_dist(xi_LMdist, dxi_LMdist, dxidt_dist, tscheme)

         call test_field(xi_LMdist, xi_LMloc, 'comp_')
         call test_field(dxi_LMdist, dxi_LMloc, 'dcompdr_')
      end if

      if ( l_conv ) then
         PERFON('up_Z')
         call updateZ( time, timeNext, z_LMloc, dz_LMloc, dzdt, omega_ma,  &
              &        omega_ic, domega_ma_dt,domega_ic_dt,                &
              &        lorentz_torque_ma_dt,lorentz_torque_ic_dt, tscheme, &
              &        lRmsNext)
         PERFOFF

         if ( l_single_matrix ) then
            if ( coord_r == rank_with_l1m0 ) then
               z10(:)=real(z_LMloc(lo_map%lm2(1,0),:))
            end if
#ifdef WITH_MPI
            call MPI_Bcast(z10,n_r_max,MPI_DEF_REAL,rank_with_l1m0, &
                 &         comm_r,ierr)
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
         !LIKWID_OFF('up_B')
      end if

      !LIKWID_OFF('LMloop')
      PERFOFF
   end subroutine LMLoop
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
         call finish_exp_comp(dVXir_LMloc, dxidt%expl(:,:,tscheme%istage))
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

      call finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic,     &
           &              domega_ma_dt%expl(tscheme%istage),        &
           &              domega_ic_dt%expl(tscheme%istage),        &
           &              lorentz_torque_ma_dt%expl(tscheme%istage),&
           &              lorentz_torque_ic_dt%expl(tscheme%istage))

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
   subroutine assemble_stage(time, w, dw, ddw, p, dp, z, dz, s, ds, xi, dxi, b, db, &
              &              ddb, aj, dj, ddj, b_ic, db_ic, ddb_ic, aj_ic, dj_ic,   &
              &              ddj_ic, omega_ic, omega_ic1, omega_ma, omega_ma1,      &
              &              dwdt, dzdt, dpdt, dsdt, dxidt, dbdt, djdt, dbdt_ic,    &
              &              djdt_ic, domega_ic_dt, domega_ma_dt,                   &
              &              lorentz_torque_ic_dt, lorentz_torque_ma_dt, lPressNext,&
              &              lRmsNext, tscheme)
      !
      ! This routine is used to call the different assembly stage of the different
      ! equations. This is only used for a special subset of IMEX-RK schemes that
      ! have tscheme%l_assembly=.true.
      !

      !-- Input variables
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time

      !-- Output variables
      complex(cp),         intent(inout) :: w(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: dw(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: ddw(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: z(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: dz(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: p(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dp(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: ds(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: xi(llm:ulm,n_r_max)
      complex(cp),         intent(out) :: dxi(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: b(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: db(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: ddb(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(inout) :: aj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: dj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(out) :: ddj(llmMag:ulmMag,n_r_maxMag)
      complex(cp),         intent(inout) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(out) :: db_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(out) :: ddb_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(inout) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(out) :: dj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(out) :: ddj_ic(llmMag:ulmMag,n_r_ic_max)

      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ma_dt
      real(cp),            intent(inout) :: omega_ic, omega_ma, omega_ic1, omega_ma1
      type(type_tarray),   intent(inout) :: dwdt, dzdt, dpdt, dsdt, dxidt
      type(type_tarray),   intent(inout) :: dbdt, djdt, dbdt_ic, djdt_ic

      if ( l_chemical_conv )  call assemble_comp(xi, dxi, dxidt, tscheme)

      if ( l_single_matrix ) then
         call assemble_single(s, ds, w, dw, ddw, dsdt, dwdt, dpdt, tscheme,lRmsNext)
      else
         if ( l_heat )  call assemble_entropy(s, ds, dsdt, tscheme)
         call assemble_pol(s, xi, w, dw, ddw, p, dp, dwdt, dpdt, dpdt%expl(:,:,1), &
              &            tscheme, lPressNext, lRmsNext)
      end if

      call assemble_tor(time, z, dz, dzdt, domega_ic_dt, domega_ma_dt,        &
           &            lorentz_torque_ic_dt, lorentz_torque_ma_dt, omega_ic, &
           &            omega_ma, omega_ic1, omega_ma1, lRmsNext, tscheme)

      if ( l_mag ) call assemble_mag(b, db, ddb, aj, dj, ddj, b_ic, db_ic,     &
                        &            ddb_ic, aj_ic, dj_ic, ddj_ic, dbdt, djdt, &
                        &            dbdt_ic, djdt_ic, lRmsNext, tscheme)

   end subroutine assemble_stage
!--------------------------------------------------------------------------------
end module LMLoop_mod
