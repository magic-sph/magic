module rIteration_mod

   use precision_mod
   use time_schemes, only: type_tscheme

   implicit none

   private
 
   type, abstract, public :: rIteration_t
      integer :: nR,nBc
      logical :: lTOCalc,lTOnext,lTOnext2
      logical :: lDeriv,lRmsCalc,lHelCalc,l_frame, lMagNlBc
      logical :: lPowerCalc, l_probe_out
      logical :: l_graph,lPerpParCalc,lViscBcCalc,lFluxProfCalc,lPressCalc
      logical :: isRadialBoundaryPoint, lPressNext
      real(cp) :: dtrkc,dthkc
 
   contains
 
      procedure(empty_if), deferred :: initialize
      procedure(empty_if), deferred :: finalize
      procedure :: set_steering_variables
      procedure(do_iteration_if), deferred :: do_iteration
      procedure(getType_if), deferred :: getType
 
   end type rIteration_t
 
   interface 
 
      subroutine empty_if(this)
         import
         class(rIteration_t) :: this
      end subroutine empty_if
   !-----------------------------------------------------------------------------
      subroutine do_iteration_if(this,nR,nBc,time,timeStage,tscheme,dtLast,&
                 &               dsdt,dwdt,dzdt,dpdt,dxidt,dbdt,djdt,      &
                 &               dVxVhLM,dVxBhLM,dVSrLM,dVXirLM,           &
                 &               br_vt_lm_cmb,br_vp_lm_cmb,br_vt_lm_icb,   &
                 &               br_vp_lm_icb,lorentz_torque_ic,           &
                 &               lorentz_torque_ma,HelLMr,Hel2LMr,         &
                 &               HelnaLMr,Helna2LMr,viscLMr,uhLMr,duhLMr,  &
                 &               gradsLMr,fconvLMr,fkinLMr,fviscLMr,       &
                 &               fpoynLMr,fresLMr,EperpLMr,EparLMr,        &
                 &               EperpaxiLMr,EparaxiLMr)
         import
         class(rIteration_t) :: this
 
         !-- Input variables
         integer,             intent(in) :: nR,nBc
         class(type_tscheme), intent(in) :: tscheme
         real(cp),            intent(in) :: time,timeStage,dtLast
     
         !-- Output variables
         complex(cp), intent(out) :: dwdt(:), dzdt(:), dpdt(:), dsdt(:), dVSrLM(:)
         complex(cp), intent(out) :: dbdt(:), djdt(:), dVxVhLM(:), dVxBhLM(:)
         complex(cp), intent(out) :: dxidt(:), dVXirLM(:)
         !---- Output of nonlinear products for nonlinear
         !     magnetic boundary conditions (needed in s_updateB.f):
         complex(cp), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
         complex(cp), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
         complex(cp), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
         complex(cp), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
         real(cp),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
         real(cp),    intent(out) :: HelLMr(:), Hel2LMr(:)
         real(cp),    intent(out) :: HelnaLMr(:), Helna2LMr(:)
         real(cp),    intent(out) :: viscLMr(:)
         real(cp),    intent(out) :: uhLMr(:), duhLMr(:), gradsLMr(:)
         real(cp),    intent(out) :: fconvLMr(:), fkinLMr(:), fviscLMr(:)
         real(cp),    intent(out) :: fpoynLMr(:),fresLMr(:)
         real(cp),    intent(out) :: EperpLMr(:), EparLMr(:)
         real(cp),    intent(out) :: EperpaxiLMr(:), EparaxiLMr(:)
 
      end subroutine do_iteration_if
   !-----------------------------------------------------------------------------
      function getType_if(this)
         import
         class(rIteration_t) :: this
         character(len=100) :: getType_if
      end function getType_if
   !-----------------------------------------------------------------------------
   end interface

contains

   subroutine set_steering_variables(this,lTOCalc,lTOnext,lTOnext2,        &
              &                      lDeriv,lRmsCalc,lHelCalc,lPowerCalc,  &
              &                      l_frame,lMagNlBc,l_graph,lViscBcCalc, &
              &                      lFluxProfCalc,lPerpParCalc,lPressCalc,&
              &                      lPressNext,l_probe_out)

      class(rIteration_t) :: this
      logical, intent(in) :: lDeriv, lRmsCalc
      logical, intent(in) :: lHelCalc, lPowerCalc, l_frame, l_probe_out
      logical, intent(in) :: lTOCalc, lTOnext, lTOnext2, lMagNlBc, l_graph
      logical, intent(in) :: lViscBcCalc, lFluxProfCalc, lPerpParCalc
      logical, intent(in) :: lPressCalc, lPressNext

      this%lTOCalc = lTOCalc
      this%lTOnext = lTOnext
      this%lTOnext2 = lTOnext2
      this%lDeriv = lDeriv
      this%lRmsCalc = lRmsCalc
      this%lHelCalc = lHelCalc
      this%lPowerCalc = lPowerCalc
      this%l_frame = l_frame
      this%lMagNlBc = lMagNlBc
      this%l_graph = l_graph
      this%lPerpParCalc = lPerpParCalc
      this%lFluxProfCalc = lFluxProfCalc
      this%lViscBcCalc = lViscBcCalc
      this%lPressCalc = lPressCalc
      this%lPressNext = lPressNext
      this%l_probe_out = l_probe_out

   end subroutine set_steering_variables
!------------------------------------------------------------------------------
end module rIteration_mod
!------------------------------------------------------------------------------
module rIteration_boundary_mod

   use rIteration_mod

   implicit none

   type, abstract, public, extends(rIteration_t) :: rIteration_boundary_t
   end type rIteration_boundary_t

end module rIteration_boundary_mod
!------------------------------------------------------------------------------
module rIteration_inner_mod

   use rIteration_mod

   implicit none

   type, abstract, public, extends(rIteration_t) :: rIteration_inner_t
   end type rIteration_inner_t

end module rIteration_inner_mod
