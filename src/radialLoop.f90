module radialLoop

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: lm_max, lm_maxMag, l_max, l_maxMag
   use radial_data,only: nRstart, nRstop, nRstartMag, nRstopMag
   use time_schemes, only: type_tscheme
   use rIteration, only: rIter_t
   use rIter_mod, only: rIter_single_t

   implicit none

   private

   public :: initialize_radialLoop, finalize_radialLoop, radialLoopG

   class(rIter_t), pointer :: rIter

contains

   subroutine initialize_radialLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated
      allocate( rIter_single_t :: rIter )
      call rIter%initialize()
      local_bytes_used = bytes_allocated-local_bytes_used

      call memWrite('radialLoop.f90', local_bytes_used)

   end subroutine initialize_radialLoop
!----------------------------------------------------------------------------
   subroutine finalize_radialLoop

      call rIter%finalize()
      deallocate(rIter)

   end subroutine finalize_radialLoop
!----------------------------------------------------------------------------
   subroutine radialLoopG(l_graph,l_frame,time,timeStage,tscheme,dtLast,   &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,     &
              &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,       &
              &          lFluxProfCalc,lPerpParCalc,lGeosCalc,lHemiCalc,   &
              &          l_probe_out,dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,&
              &          djdt,dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,              &
              &          lorentz_torque_ic,lorentz_torque_ma,br_vt_lm_cmb, &
              &          br_vp_lm_cmb,br_vt_lm_icb,br_vp_lm_icb,dtrkc,dthkc)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !--- Input of variables:
      logical,             intent(in) :: l_graph,l_frame
      logical,             intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,             intent(in) :: lPowerCalc,lGeosCalc
      logical,             intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,             intent(in) :: lRmsCalc,lHemiCalc
      logical,             intent(in) :: l_probe_out
      logical,             intent(in) :: lPressCalc
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: time,timeStage,dtLast
      class(type_tscheme), intent(in) :: tscheme

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dphidt(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)
      complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(lm_max) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(lm_max) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(lm_max) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(lm_max) ! product br*vp at ICB
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for Courant criteria:
      real(cp),    intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)

      call rIter%radialLoop(l_graph,l_frame,time,timeStage,tscheme,dtLast,      &
           &             lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,          &
           &             lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,            &
           &             lFluxProfCalc,lPerpParCalc,lGeosCalc,lHemiCalc,        &
           &             l_probe_out,dsdt,dwdt,dzdt,dpdt,dxidt,dphidt,dbdt,djdt,&
           &             dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,lorentz_torque_ic,      &
           &             lorentz_torque_ma, br_vt_lm_cmb,br_vp_lm_cmb,          &
           &             br_vt_lm_icb,br_vp_lm_icb,dtrkc,dthkc)

   end subroutine radialLoopG
!----------------------------------------------------------------------------
end module radialLoop
