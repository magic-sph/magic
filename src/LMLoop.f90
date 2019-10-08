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
   use truncation, only: l_max, lm_max, n_r_max, n_r_maxMag
   use radial_data, only: n_r_icb, n_r_cmb
   use blocking, only: lo_map, llm, ulm, llmMag, ulmMag
   use logic, only: l_mag, l_conv, l_anelastic_liquid, lVerbose, l_heat, &
       &            l_single_matrix, l_chemical_conv, l_save_out
   use output_data, only: n_log_file, log_file
   use debugging,  only: debug_write
   use updateS_mod
   use updateZ_mod
   use updateWP_mod
   use updateWPS_mod
   use updateB_mod
   use updateXi_mod

   implicit none

   private

   public :: LMLoop, initialize_LMLoop, finalize_LMLoop

contains

   subroutine initialize_LMLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

      if ( l_single_matrix ) then
         call initialize_updateWPS
      else
         call initialize_updateS
         call initialize_updateWP
      end if

      if ( l_chemical_conv ) call initialize_updateXi

      call initialize_updateZ
      if ( l_mag ) call initialize_updateB

      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('LMLoop.f90',local_bytes_used)

   end subroutine initialize_LMLoop
!----------------------------------------------------------------------------
   subroutine finalize_LMLoop

      if ( l_single_matrix ) then
         call finalize_updateWPS
      else
         call finalize_updateS
         call finalize_updateWP
      end if

      if ( l_chemical_conv ) call finalize_updateXi

      call finalize_updateZ
      if ( l_mag ) call finalize_updateB

   end subroutine finalize_LMLoop
!----------------------------------------------------------------------------
   subroutine LMLoop(w1,coex,time,dt,lMat,lRmsNext,lPressNext,dVxVhLM, &
              &      dVxBhLM,dVSrLM,dVXirLM,dsdt,dwdt,                 &
              &      dzdt,dpdt,dxidt,dbdt,djdt,lorentz_torque_ma,      &
              &      lorentz_torque_ic,b_nl_cmb,aj_nl_cmb,             &
              &      aj_nl_icb)
      !
      !  This subroutine performs the actual time-stepping.
      !
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1,coex
      real(cp),    intent(in) :: dt,time
      logical,     intent(in) :: lMat
      logical,     intent(in) :: lRmsNext
      logical,     intent(in) :: lPressNext

      !--- Input from radialLoop:
      !    These fields are provided in the R-distributed space!
      ! for djdt in update_b
      complex(cp), intent(inout) :: dVxBhLM(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: dVxVhLM(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)   ! for dsdt in update_s
      complex(cp), intent(inout) :: dVXirLM(llm:ulm,n_r_max)  ! for dxidt in update_xi
      !integer,     intent(in) :: n_time_step

      !--- Input from radialLoop and then redistributed:
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dxidt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dwdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dzdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dpdt(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dbdt(llmMag:ulmMag,n_r_maxMag)
      complex(cp), intent(inout) :: djdt(llmMag:ulmMag,n_r_maxMag)
      real(cp),    intent(in) :: lorentz_torque_ma,lorentz_torque_ic
      complex(cp), intent(in) :: b_nl_cmb(lm_max)   ! nonlinear bc for b at CMB
      complex(cp), intent(in)  :: aj_nl_cmb(lm_max)  ! nonlinear bc for aj at CMB
      complex(cp), intent(in)  :: aj_nl_icb(lm_max)  ! nonlinear bc for dr aj at ICB

      !--- Local counter
      integer :: l,nR,ierr

      !--- Inner core rotation from last time step
      real(cp), save :: omega_icLast
      real(cp) :: z10(n_r_max)


      PERFON('LMloop')
      !LIKWID_ON('LMloop')
      if ( lVerbose .and. l_save_out ) then
         open(newunit=n_log_file, file=log_file, status='unknown', &
         &    position='append')
      end if

      omega_icLast=omega_ic

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
               lSmat(l) =.false.
            end if
            lZmat(l) =.false.
            if ( l_mag ) lBmat(l) =.false.
            if ( l_chemical_conv ) lXimat(l)=.false.
         end do
      end if

      if ( l_heat ) then ! dp,workA usead as work arrays
         if ( .not. l_single_matrix ) then
            PERFON('up_S')
            if ( l_anelastic_liquid ) then
               call updateS_ala(s_LMloc, ds_LMloc, w_LMloc, dVSrLM,dsdt,  &
                    &           dsdtLast_LMloc, w1, coex, dt)
            else
               call updateS( s_LMloc, ds_LMloc, w_LMloc, dVSrLM,dsdt, &
                    &        dsdtLast_LMloc, w1, coex, dt )
            end if
            PERFOFF
         end if

      end if

      if ( l_chemical_conv ) then ! dp,workA usead as work arrays
         call updateXi(xi_LMloc,dxi_LMloc,dVXirLM,dxidt,dxidtLast_LMloc, &
              &        w1,coex,dt)
      end if

      if ( l_conv ) then
         PERFON('up_Z')
         call updateZ( z_LMloc, dz_LMloc, dzdt, dzdtLast_lo, time, &
              &        omega_ma,d_omega_ma_dtLast,                 &
              &        omega_ic,d_omega_ic_dtLast,                 &
              &        lorentz_torque_ma,lorentz_torque_maLast,    &
              &        lorentz_torque_ic,lorentz_torque_icLast,    &
              &        w1,coex,dt,lRmsNext)
         PERFOFF

         if ( l_single_matrix ) then
            if ( rank == rank_with_l1m0 ) then
               do nR=1,n_r_max
                  z10(nR)=real(z_LMloc(lo_map%lm2(1,0),nR))
               end do
            end if
#ifdef WITH_MPI
            call MPI_Bcast(z10,n_r_max,MPI_DEF_REAL,rank_with_l1m0, &
                 &         MPI_COMM_WORLD,ierr)
#endif
            call updateWPS( w_LMloc, dw_LMloc, ddw_LMloc, z10, dwdt,    &
                 &          dwdtLast_LMloc, p_LMloc, dp_LMloc, dpdt,    &
                 &          dpdtLast_LMloc, s_LMloc, ds_LMloc, dVSrLM,  &
                 &          dsdt, dsdtLast_LMloc, w1, coex, dt,         &
                 &          lRmsNext )
         else
            PERFON('up_WP')
            call updateWP( w_LMloc, dw_LMloc, ddw_LMloc, dVxVhLM, dwdt,     &
                 &         dwdtLast_LMloc, p_LMloc, dp_LMloc, dpdt,         &
                 &         dpdtLast_LMloc, s_LMloc, xi_LMloc, w1, coex, dt, &
                 &         lRmsNext, lPressNext)
            PERFOFF
         end if
      end if
      if ( l_mag ) then ! dwdt,dpdt used as work arrays
         PERFON('up_B')
         call updateB( b_LMloc,db_LMloc,ddb_LMloc,aj_LMloc,dj_LMloc,ddj_LMloc, &
              &        dVxBhLM, dbdt, dbdtLast_LMloc, djdt, djdtLast_LMloc,    &
              &        b_ic_LMloc, db_ic_LMloc, ddb_ic_LMloc, aj_ic_LMloc,     &
              &        dj_ic_LMloc, ddj_ic_LMloc, dbdt_icLast_LMloc,           &
              &        djdt_icLast_LMloc, b_nl_cmb, aj_nl_cmb, aj_nl_icb,      &
              &        omega_icLast, w1, coex, dt, time, lRmsNext )
         PERFOFF
         !LIKWID_OFF('up_B')
      end if

      lorentz_torque_maLast=lorentz_torque_ma
      lorentz_torque_icLast=lorentz_torque_ic

      if ( lVerbose .and. l_save_out ) close(n_log_file)

      !LIKWID_OFF('LMloop')
      PERFOFF
   end subroutine LMLoop
!--------------------------------------------------------------------------------
end module LMLoop_mod
