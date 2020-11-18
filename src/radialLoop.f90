module radialLoop

   use precision_mod
   use mem_alloc, only: memWrite, bytes_allocated
   use truncation, only: n_lm_loc, n_lmMag_loc, nRstart, nRstop, nRstartMag, &
       &                 nRstopMag,  n_lmP_loc, rIter_type
   use time_schemes, only: type_tscheme
   use rIteration, only: rIter_t
   use rIter_split, only: rIter_split_t
   use rIter_single, only: rIter_single_t

   implicit none

   private

   public :: initialize_radialLoop, finalize_radialLoop, radialLoopG

   class(rIter_t), pointer :: rIter

contains

   subroutine initialize_radialLoop

      integer(lip) :: local_bytes_used

      local_bytes_used = bytes_allocated

      if ( index(rIter_type, 'SINGLE') /= 0 ) then
         allocate( rIter_single_t :: rIter )
      else
         allocate( rIter_split_t :: rIter )
      end if
      
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
   subroutine radialLoopG(l_graph,l_frame,time,timeStage,tscheme,dtLast, &
              &          lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,   &
              &          lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,     &
              &          lFluxProfCalc,lPerpParCalc,l_probe_out,dsdt,    &
              &          dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxVhLM,dVxBhLM, &
              &          dVSrLM,dVXirLM,lorentz_torque_ic,               &
              &          lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,    &
              &          br_vt_lm_icb,br_vp_lm_icb,HelASr,Hel2ASr,       &
              &          HelnaASr,Helna2ASr,HelEAASr,viscAS,uhASr,       &
              &          duhASr,gradsASr,fconvASr,fkinASr,fviscASr,      &
              &          fpoynASr,fresASr,EperpASr,EparASr,              &
              &          EperpaxiASr,EparaxiASr,dtrkc,dthkc)
      !
      !  This subroutine performs the actual time-stepping.
      !

      !--- Input of variables:
      logical,             intent(in) :: l_graph,l_frame
      logical,             intent(in) :: lTOcalc,lTONext,lTONext2,lHelCalc
      logical,             intent(in) :: lPowerCalc
      logical,             intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc
      logical,             intent(in) :: lRmsCalc
      logical,             intent(in) :: l_probe_out
      logical,             intent(in) :: lPressCalc
      logical,             intent(in) :: lPressNext
      real(cp),            intent(in) :: time,timeStage,dtLast
      class(type_tscheme), intent(in) :: tscheme

      !---- Output of explicit time step:
      !---- dVSrLM and dVxBhLM are output of contributions to explicit time step that
      !     need a further treatment (radial derivatives required):
      complex(cp), intent(out) :: dwdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dzdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dpdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dsdt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dxidt(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVSrLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVXirLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dbdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: djdt(n_lmMag_loc,nRstartMag:nRstopMag)
      complex(cp), intent(out) :: dVxVhLM(n_lm_loc,nRstart:nRstop)
      complex(cp), intent(out) :: dVxBhLM(n_lmMag_loc,nRstartMag:nRstopMag)
      real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

      !---- Output for axisymmetric helicity:
      real(cp),    intent(out) :: HelASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: Hel2ASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: HelnaASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: Helna2ASr(2,nRstart:nRstop)
      real(cp),    intent(out) :: HelEAASr(nRstart:nRstop)
      real(cp),    intent(out) :: uhASr(nRstart:nRstop)
      real(cp),    intent(out) :: duhASr(nRstart:nRstop)
      real(cp),    intent(out) :: viscAS(nRstart:nRstop)
      real(cp),    intent(out) :: gradsASr(nRstart:nRstop)
      real(cp),    intent(out) :: fkinASr(nRstart:nRstop)
      real(cp),    intent(out) :: fconvASr(nRstart:nRstop)
      real(cp),    intent(out) :: fviscASr(nRstart:nRstop)
      real(cp),    intent(out) :: fresASr(nRstartMag:nRstopMag)
      real(cp),    intent(out) :: fpoynASr(nRstartMag:nRstopMag)
      real(cp),    intent(out) :: EperpASr(nRstart:nRstop)
      real(cp),    intent(out) :: EparASr(nRstart:nRstop)
      real(cp),    intent(out) :: EperpaxiASr(nRstart:nRstop)
      real(cp),    intent(out) :: EparaxiASr(nRstart:nRstop)

      !---- Output of nonlinear products for nonlinear
      !     magnetic boundary conditions (needed in s_updateB.f):
      complex(cp), intent(out) :: br_vt_lm_cmb(n_lmP_loc) ! product br*vt at CMB
      complex(cp), intent(out) :: br_vp_lm_cmb(n_lmP_loc) ! product br*vp at CMB
      complex(cp), intent(out) :: br_vt_lm_icb(n_lmP_loc) ! product br*vt at ICB
      complex(cp), intent(out) :: br_vp_lm_icb(n_lmP_loc) ! product br*vp at ICB

      !---- Output for Courant criteria:
      real(cp),    intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)

      call rIter%radialLoop(l_graph,l_frame,time,timeStage,tscheme,dtLast,  &
              &             lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,   &
              &             lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,     &
              &             lFluxProfCalc,lPerpParCalc,l_probe_out,dsdt,    &
              &             dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxVhLM,dVxBhLM, &
              &             dVSrLM,dVXirLM,lorentz_torque_ic,               &
              &             lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,    &
              &             br_vt_lm_icb,br_vp_lm_icb,HelASr,Hel2ASr,       &
              &             HelnaASr,Helna2ASr,HelEAASr,viscAS,uhASr,       &
              &             duhASr,gradsASr,fconvASr,fkinASr,fviscASr,      &
              &             fpoynASr,fresASr,EperpASr,EparASr,              &
              &             EperpaxiASr,EparaxiASr,dtrkc,dthkc)

   end subroutine radialLoopG
!----------------------------------------------------------------------------
end module radialLoop
