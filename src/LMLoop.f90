#include "perflib_preproc.cpp"
module LMLoop_mod

#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

   use fields
   use fieldsLast
   use omp_lib
   use precision_mod
   use parallel_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag,          &
       &            n_r_ic_max, n_r_icb, n_r_cmb
   use blocking, only: lo_map, llm, ulm, llmMag, ulmMag
   use logic, only: l_mag, l_conv, lVerbose, l_heat, l_single_matrix, &
       &            l_double_curl, l_chemical_conv, l_save_out,       &
       &            l_cond_ic
   use output_data, only: n_log_file, log_file
   use debugging,  only: debug_write
   use time_array, only: type_tarray
   use time_schemes, only: type_tscheme
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod

   implicit none

   private

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop, finish_explicit_assembly

contains

   subroutine initialize_LMLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

      if ( l_single_matrix ) then
         call initialize_updateWPS()
      else
         if ( l_heat ) call initialize_updateS()
         call initialize_updateWP()
      end if

      if ( l_chemical_conv ) call initialize_updateXi()

      call initialize_updateZ()
      if ( l_mag ) call initialize_updateB()

      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('LMLoop.f90',local_bytes_used)

   end subroutine initialize_LMLoop
!----------------------------------------------------------------------------
   subroutine finalize_LMLoop

      if ( l_single_matrix ) then
         call finalize_updateWPS()
      else
         if ( l_heat ) call finalize_updateS()
         call finalize_updateWP()
      end if
      call finalize_updateZ()

      if ( l_chemical_conv ) call finalize_updateXi()
      if ( l_mag ) call finalize_updateB()

   end subroutine finalize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(time,tscheme,lMat,lRmsNext,lPressNext,            &
              &      dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dbdt_ic,      &
              &      djdt_ic,lorentz_torque_ma,lorentz_torque_ic,      &
              &      b_nl_cmb,aj_nl_cmb,aj_nl_icb)
      !
      !  This subroutine performs the actual time-stepping.
      !
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: time
      logical,             intent(in) :: lMat
      logical,             intent(in) :: lRmsNext
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      complex(cp),         intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp),         intent(in)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp),         intent(in)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

      !--- Input from radialLoop:
      type(type_tarray), intent(inout) :: dsdt
      type(type_tarray), intent(inout) :: dxidt
      type(type_tarray), intent(inout) :: dwdt
      type(type_tarray), intent(inout) :: dpdt
      type(type_tarray), intent(inout) :: dzdt
      type(type_tarray), intent(inout) :: dbdt
      type(type_tarray), intent(inout) :: djdt
      type(type_tarray), intent(inout) :: dbdt_ic
      type(type_tarray), intent(inout) :: djdt_ic
      !integer,     intent(in) :: n_time_step

      !--- Local counter
      integer :: l,ierr

      !--- Inner core rotation from last time step
      real(cp) :: z10(n_r_max)


      PERFON('LMloop')
      !LIKWID_ON('LMloop')
      if ( lVerbose .and. l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if

      if ( lMat ) then ! update matrices:
      !---- The following logicals tell whether the respective inversion
      !     matrices have been updated. lMat=.true. when a general
      !     update is necessary. These logicals are THREADPRIVATE and
      !     stored in the module matrices in m_mat.F90:
         lZ10mat=.false.
         do l=0,l_max
            if ( l_single_matrix ) then
               lWPSmat(l)=.false.
            else
               lWPmat(l)=.false.
               if ( l_heat ) lSmat(l) =.false.
            end if
            lZmat(l) =.false.
            if ( l_mag ) lBmat(l) =.false.
            if ( l_chemical_conv ) lXimat(l)=.false.
         end do
      end if

      if ( l_heat .and. .not. l_single_matrix ) then
         PERFON('up_S')
         call updateS( s_LMloc, ds_LMloc, dsdt, tscheme )
         PERFOFF
      end if

      if ( l_chemical_conv ) call updateXi(xi_LMloc, dxi_LMloc, dxidt, tscheme)

      if ( l_conv ) then
         PERFON('up_Z')
         call updateZ( z_LMloc, dz_LMloc, dzdt, time, omega_ma, omega_ic, &
              &        lorentz_torque_ma, lorentz_torque_ic, tscheme,lRmsNext)
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

      if ( lVerbose .and. l_save_out ) close(n_log_file)

      !LIKWID_OFF('LMloop')
      PERFOFF
   end subroutine LMLoop
!--------------------------------------------------------------------------------
   subroutine finish_explicit_assembly(omega_ic, w, b_ic, aj_ic, dVSr_LMloc,  &
              &                        dVXir_LMloc, dVxVh_LMloc, dVxBh_LMloc, &
              &                        dsdt, dxidt, dwdt, djdt, dbdt_ic,      &
              &                        djdt_ic, tscheme)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      real(cp),            intent(in) :: omega_ic
      complex(cp),         intent(in) :: w(llm:ulm,n_r_max)
      complex(cp),         intent(in) :: b_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(in) :: aj_ic(llmMag:ulmMag,n_r_ic_max)
      complex(cp),         intent(inout) :: dVSr_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVXir_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVxVh_LMloc(llm:ulm,n_r_max)
      complex(cp),         intent(inout) :: dVxBh_LMloc(llmMag:ulmMag,n_r_maxMag)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dsdt
      type(type_tarray),   intent(inout) :: dxidt
      type(type_tarray),   intent(inout) :: djdt
      type(type_tarray),   intent(inout) :: dwdt
      type(type_tarray),   intent(inout) :: dbdt_ic
      type(type_tarray),   intent(inout) :: djdt_ic

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
end module LMLoop_mod
