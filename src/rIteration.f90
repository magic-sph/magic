module rIteration
   !
   ! This module is used to define an abstract class for the radial loop
   !

   use precision_mod
   use radial_data, only: nRstart, nRstop, nRstartMag, nRstopMag
   use truncation, only: lm_max, lm_maxMag
   use time_schemes, only: type_tscheme

   implicit none

   private

   type, abstract, public :: rIter_t
   contains
      procedure(initialize_if), deferred :: initialize
      procedure(finalize_if), deferred :: finalize
      procedure(radialLoop_if), deferred :: radialLoop
   end type rIter_t

   interface

      subroutine initialize_if(this)
         import
         class(rIter_t) :: this
      end subroutine initialize_if

      subroutine finalize_if(this)
         import
         class(rIter_t) :: this
      end subroutine finalize_if

      subroutine radialLoop_if(this,l_graph,l_frame,time,timeStage,tscheme,dtLast, &
              &                lTOCalc,lTONext,lTONext2,lHelCalc,lPowerCalc,       &
              &                lRmsCalc,lPressCalc,lPressNext,lViscBcCalc,         &
              &                lFluxProfCalc,lPerpParCalc,l_probe_out,dsdt,        &
              &                dwdt,dzdt,dpdt,dxidt,dbdt,djdt,dVxVhLM,dVxBhLM,     &
              &                dVSrLM,dVXirLM,lorentz_torque_ic,                   &
              &                lorentz_torque_ma,br_vt_lm_cmb,br_vp_lm_cmb,        &
              &                br_vt_lm_icb,br_vp_lm_icb,HelAS,Hel2AS,             &
              &                HelnaAS,Helna2AS,HelEAAS,viscAS,uhAS,               &
              &                duhAS,gradsAS,fconvAS,fkinAS,fviscAS,               &
              &                fpoynAS,fresAS,EperpAS,EparAS,                      &
              &                EperpaxiAS,EparaxiAS,dtrkc,dthkc)
         import
         class(rIter_t) :: this
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

         complex(cp), intent(out) :: dwdt(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dzdt(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dpdt(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dsdt(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dxidt(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dVSrLM(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dVXirLM(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dbdt(lm_maxMag,nRstartMag:nRstopMag)
         complex(cp), intent(out) :: djdt(lm_maxMag,nRstartMag:nRstopMag)
         complex(cp), intent(out) :: dVxVhLM(lm_max,nRstart:nRstop)
         complex(cp), intent(out) :: dVxBhLM(lm_maxMag,nRstartMag:nRstopMag)
         real(cp),    intent(out) :: lorentz_torque_ma,lorentz_torque_ic

         !---- inoutput for axisymmetric helicity:
         real(cp),    intent(inout) :: HelAS(2,nRstart:nRstop)
         real(cp),    intent(inout) :: Hel2AS(2,nRstart:nRstop)
         real(cp),    intent(inout) :: HelnaAS(2,nRstart:nRstop)
         real(cp),    intent(inout) :: Helna2AS(2,nRstart:nRstop)
         real(cp),    intent(inout) :: HelEAAS(nRstart:nRstop)
         real(cp),    intent(inout) :: uhAS(nRstart:nRstop)
         real(cp),    intent(inout) :: duhAS(nRstart:nRstop)
         real(cp),    intent(inout) :: viscAS(nRstart:nRstop)
         real(cp),    intent(inout) :: gradsAS(nRstart:nRstop)
         real(cp),    intent(inout) :: fkinAS(nRstart:nRstop)
         real(cp),    intent(inout) :: fconvAS(nRstart:nRstop)
         real(cp),    intent(inout) :: fviscAS(nRstart:nRstop)
         real(cp),    intent(inout) :: fresAS(nRstartMag:nRstopMag)
         real(cp),    intent(inout) :: fpoynAS(nRstartMag:nRstopMag)
         real(cp),    intent(inout) :: EperpAS(nRstart:nRstop)
         real(cp),    intent(inout) :: EparAS(nRstart:nRstop)
         real(cp),    intent(inout) :: EperpaxiAS(nRstart:nRstop)
         real(cp),    intent(inout) :: EparaxiAS(nRstart:nRstop)

         !---- inoutput of nonlinear products for nonlinear
         !     magnetic boundary conditions (needed in s_updateB.f):
         complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
         complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
         complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
         complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB

         !---- inoutput for Courant criteria:
         real(cp),    intent(out) :: dtrkc(nRstart:nRstop),dthkc(nRstart:nRstop)
      end subroutine radialLoop_if

   end interface

end module rIteration
