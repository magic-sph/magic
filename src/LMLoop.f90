#include "perflib_preproc.cpp"
module LMLoop_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use precision_mod
   use parallel_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: lm_max, n_r_max, n_r_maxMag, n_r_ic_max, nRstartMag, &
       &                 n_mlo_loc, n_mloMag_loc, nRstart, nRstop, nRstopMag, &
       &                 n_lm_loc, n_lmMag_loc
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag, n_r_ic_max, lm_max, &
       &                 lm_maxMag
   use logic, only: l_mag, l_conv, l_heat, l_single_matrix, l_double_curl, &
       &            l_chemical_conv, l_cond_ic
   use LMmapping, only: map_mlo
   use time_array, only: type_tarray, type_tscalar
   use time_schemes, only: type_tscheme
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod
   
   implicit none

   private

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop, finish_explicit_assembly, &
   &         assemble_stage, finish_explicit_assembly_Rdist

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
              &      dsdt_dist,dwdt_dist,dzdt_dist,dpdt_dist,          &
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
      type(type_tarray),  intent(inout) :: dsdt_dist, dxidt_dist, dwdt_dist
      type(type_tarray),  intent(inout) :: dpdt_dist, dzdt_dist
      type(type_tarray),  intent(inout) :: dbdt_dist, djdt_dist, dbdt_ic_dist, djdt_ic_dist
      type(type_tscalar), intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar), intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt
      !integer,     intent(in) :: n_time_step

      !--- Inner core rotation from last time step
      real(cp) :: z10(n_r_max)

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
      end if

      if ( l_heat .and. .not. l_single_matrix ) then
         call updateS( s_LMdist, ds_LMdist, dsdt_dist, tscheme )
      end if
      
      if ( l_chemical_conv ) then
         call updateXi(xi_LMdist, dxi_LMdist, dxidt_dist, tscheme)
      end if

      if ( l_conv ) then
         call updateZ( time, timeNext, z_LMdist, dz_LMdist, dzdt_dist, omega_ma,  &
              &        omega_ic, domega_ma_dt,domega_ic_dt,                &
              &        lorentz_torque_ma_dt,lorentz_torque_ic_dt, tscheme, &
              &        lRmsNext)

         if ( l_single_matrix ) then
            if ( map_mlo%ml2i(0,1) > 0 ) then
               z10(:)=real(z_LMdist(map_mlo%ml2i(0,1),:))
            end if
#ifdef WITH_MPI
            !@> TODO: probably overkill here: ask Rafael whether he has an idea
            call MPI_Bcast(z10,n_r_max,MPI_DEF_REAL,map_mlo%ml2rnk(0,1),MPI_COMM_WORLD,ierr)
#endif
            call updateWPS( w_LMdist, dw_LMdist, ddw_LMdist, z10, dwdt_dist,     &
                 &          p_LMdist, dp_LMdist, dpdt_dist, s_LMdist, ds_LMdist, &
                 &          dsdt_dist, tscheme, lRmsNext )
         else
            call updateWP( s_LMdist, xi_LMdist, w_LMdist, dw_LMdist, ddw_LMdist,&
                 &         dwdt_dist, p_LMdist, dp_LMdist, dpdt_dist, tscheme,  &
                 &         lRmsNext, lPressNext )
         end if
      end if
      if ( l_mag ) then ! dwdt,dpdt used as work arrays
         call updateB( b_LMdist,db_LMdist,ddb_LMdist,aj_LMdist,dj_LMdist,  &
              &        ddj_LMdist, dbdt_dist, djdt_dist, b_ic_LMdist,      &
              &        db_ic_LMdist, ddb_ic_LMdist, aj_ic_LMdist,          &
              &        dj_ic_LMdist, ddj_ic_LMdist, dbdt_ic_dist,          &
              &        djdt_ic_dist, b_nl_cmb, aj_nl_cmb, aj_nl_icb, time, &
              &        tscheme, lRmsNext )
      end if

   end subroutine LMLoop
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly(omega_ic, w, b_ic, aj_ic, dVSr_LMdist,     &
              &                        dVXir_LMdist, dVxVh_LMdist, dVxBh_LMdist,  &
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
      complex(cp),         intent(in) :: w(n_mlo_loc,n_r_max)
      complex(cp),         intent(in) :: b_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(inout) :: dVSr_LMdist(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: dVXir_LMdist(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: dVxVh_LMdist(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: dVxBh_LMdist(n_mloMag_loc,n_r_maxMag)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dsdt, dxidt, djdt, dwdt
      type(type_tarray),   intent(inout) :: dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      if ( l_chemical_conv ) then
         call finish_exp_comp(dVXir_LMdist, dxidt%expl(:,:,tscheme%istage))
      end if

      if ( l_single_matrix ) then
         call finish_exp_smat(dVSr_LMdist, dsdt%expl(:,:,tscheme%istage))
      else
         if ( l_heat ) then
            call finish_exp_entropy(w, dVSr_LMdist, dsdt%expl(:,:,tscheme%istage))
         end if
         if ( l_double_curl ) then
            call finish_exp_pol(dVxVh_LMdist, dwdt%expl(:,:,tscheme%istage))
         end if
      end if

      call finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic,     &
           &              domega_ma_dt%expl(tscheme%istage),        &
           &              domega_ic_dt%expl(tscheme%istage),        &
           &              lorentz_torque_ma_dt%expl(tscheme%istage),&
           &              lorentz_torque_ic_dt%expl(tscheme%istage))

      if ( l_mag ) then
         call finish_exp_mag(dVxBh_LMdist, djdt%expl(:,:,tscheme%istage))
      end if

      if ( l_cond_ic ) then
         call finish_exp_mag_ic(b_ic, aj_ic, omega_ic,            &
              &                 dbdt_ic%expl(:,:,tscheme%istage), &
              &                 djdt_ic%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly_Rdist(omega_ic, w, b_ic, aj_ic, dVSr_Rdist,  &
              &                        dVXir_Rdist, dVxVh_Rdist, dVxBh_Rdist,       &
              &                        lorentz_torque_ma, lorentz_torque_ic,        &
              &                        dsdt_Rdist, dxidt_Rdist, dwdt_Rdist,         &
              &                        djdt_Rdist, dbdt_ic, djdt_ic, domega_ma_dt,  &
              &                        domega_ic_dt, lorentz_torque_ma_dt,          &
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
      complex(cp),         intent(in) :: w(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(in) :: b_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(in) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(inout) :: dVSr_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVXir_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVxVh_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: dVxBh_Rdist(n_lmMag_loc,nRstartMag:nRstopMag)

      !-- Output variables
      complex(cp),         intent(inout) :: dxidt_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: dsdt_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: dwdt_Rdist(n_lm_loc,nRstart:nRstop)
      complex(cp),         intent(inout) :: djdt_Rdist(n_lm_loc,nRstart:nRstop)
      type(type_tarray),   intent(inout) :: dbdt_ic, djdt_ic
      type(type_tscalar),  intent(inout) :: domega_ic_dt, domega_ma_dt
      type(type_tscalar),  intent(inout) :: lorentz_torque_ic_dt, lorentz_torque_ma_dt

      if ( l_chemical_conv ) call finish_exp_comp_Rdist(dVXir_Rdist, dxidt_Rdist)

      if ( l_single_matrix ) then
         call finish_exp_smat_Rdist(dVSr_Rdist, dsdt_Rdist)
      else
         if ( l_heat ) call finish_exp_entropy_Rdist(w, dVSr_Rdist, dsdt_Rdist)
         if ( l_double_curl ) call finish_exp_pol_Rdist(dVxVh_Rdist, dwdt_Rdist)
      end if

      call finish_exp_tor(lorentz_torque_ma, lorentz_torque_ic,     &
           &              domega_ma_dt%expl(tscheme%istage),        &
           &              domega_ic_dt%expl(tscheme%istage),        &
           &              lorentz_torque_ma_dt%expl(tscheme%istage),&
           &              lorentz_torque_ic_dt%expl(tscheme%istage))

      if ( l_mag ) call finish_exp_mag_Rdist(dVxBh_Rdist, djdt_Rdist)

      if ( l_cond_ic ) then
         call finish_exp_mag_ic(b_ic, aj_ic, omega_ic,            &
              &                 dbdt_ic%expl(:,:,tscheme%istage), &
              &                 djdt_ic%expl(:,:,tscheme%istage))
      end if

   end subroutine finish_explicit_assembly_Rdist
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
      ! have ``tscheme%l_assembly=.true.``
      !

      !-- Input variables
      logical,             intent(in) :: lPressNext
      logical,             intent(in) :: lRmsNext
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time

      !-- Output variables
      complex(cp),         intent(inout) :: w(n_mlo_loc,n_r_max)
      complex(cp),         intent(out) :: dw(n_mlo_loc,n_r_max)
      complex(cp),         intent(out) :: ddw(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: z(n_mlo_loc,n_r_max)
      complex(cp),         intent(out) :: dz(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: p(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: dp(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: s(n_mlo_loc,n_r_max)
      complex(cp),         intent(out) :: ds(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: xi(n_mlo_loc,n_r_max)
      complex(cp),         intent(out) :: dxi(n_mlo_loc,n_r_max)
      complex(cp),         intent(inout) :: b(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(out) :: db(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(out) :: ddb(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(inout) :: aj(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(out) :: dj(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(out) :: ddj(n_mloMag_loc,n_r_maxMag)
      complex(cp),         intent(inout) :: b_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(out) :: db_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(out) :: ddb_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(inout) :: aj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(out) :: dj_ic(n_mloMag_loc,n_r_ic_max)
      complex(cp),         intent(out) :: ddj_ic(n_mloMag_loc,n_r_ic_max)

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
