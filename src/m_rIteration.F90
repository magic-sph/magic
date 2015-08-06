!$Id$
module rIteration_mod

   implicit none

   private
 
   type, abstract, public :: rIteration_t
      integer :: nR,nBc
      logical :: l_cour
      logical :: lTOCalc,lTOnext,lTOnext2
      logical :: lDeriv,lRmsCalc,lHelCalc,l_frame, lMagNlBc
      logical :: l_graph,lPerpParCalc,lViscBcCalc,lFluxProfCalc
      logical :: isRadialBoundaryPoint
      real(kind=8) :: dtrkc,dthkc
 
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
      subroutine do_iteration_if(this,nR,nBc,time,dt,dtLast,                      &
        &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM,           &
        &                 br_vt_lm_cmb,br_vp_lm_cmb,                              &
        &                 br_vt_lm_icb,br_vp_lm_icb,                              &
        &                 lorentz_torque_ic,lorentz_torque_ma,                    &
        &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr,&
        &                 fconvLMr,fkinLMr,fviscLMr,fpoynLMr,fresLMr,             &
        &                 EperpLMr,EparLMr,EperpaxiLMr,EparaxiLMr)
        import
        class(rIteration_t) :: this
 
        !-- Input variables
        integer,      intent(in) :: nR,nBc
        real(kind=8), intent(in) :: time,dt,dtLast
     
        !-- Output variables
        complex(kind=8), intent(out) :: dwdt(:), dzdt(:), dpdt(:), dsdt(:), dVSrLM(:)
        complex(kind=8), intent(out) :: dbdt(:), djdt(:), dVxBhLM(:)
        !---- Output of nonlinear products for nonlinear
        !     magnetic boundary conditions (needed in s_updateB.f):
        complex(kind=8), intent(out) :: br_vt_lm_cmb(:) ! product br*vt at CMB
        complex(kind=8), intent(out) :: br_vp_lm_cmb(:) ! product br*vp at CMB
        complex(kind=8), intent(out) :: br_vt_lm_icb(:) ! product br*vt at ICB
        complex(kind=8), intent(out) :: br_vp_lm_icb(:) ! product br*vp at ICB
        real(kind=8),    intent(out) :: lorentz_torque_ic,lorentz_torque_ma
        real(kind=8),    intent(out) :: HelLMr(:), Hel2LMr(:)
        real(kind=8),    intent(out) :: HelnaLMr(:), Helna2LMr(:)
        real(kind=8),    intent(out) :: uhLMr(:), duhLMr(:), gradsLMr(:)
        real(kind=8),    intent(out) :: fconvLMr(:), fkinLMr(:), fviscLMr(:)
        real(kind=8),    intent(out) :: fpoynLMr(:),fresLMr(:)
        real(kind=8),    intent(out) :: EperpLMr(:), EparLMr(:)
        real(kind=8),    intent(out) :: EperpaxiLMr(:), EparaxiLMr(:)
 
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

   subroutine set_steering_variables(this,l_cour,lTOCalc,lTOnext,lTOnext2,&
        & lDeriv,lRmsCalc,lHelCalc,l_frame,lMagNlBc,l_graph,lViscBcCalc,  &
        & lFluxProfCalc,lPerpParCalc)

      class(rIteration_t) :: this
      logical, intent(in) :: l_cour,lDeriv,lRmsCalc,lHelCalc,l_frame
      logical, intent(in) :: lTOCalc,lTOnext,lTOnext2, lMagNlBc,l_graph
      logical, intent(in) :: lViscBcCalc,lFluxProfCalc,lPerpParCalc

      this%l_cour = l_cour
      this%lTOCalc = lTOCalc
      this%lTOnext = lTOnext
      this%lTOnext2 = lTOnext2
      this%lDeriv = lDeriv
      this%lRmsCalc = lRmsCalc
      this%lHelCalc = lHelCalc
      this%l_frame = l_frame
      this%lMagNlBc = lMagNlBc
      this%l_graph = l_graph
      this%lPerpParCalc = lPerpParCalc
      this%lFluxProfCalc = lFluxProfCalc
      this%lViscBcCalc = lViscBcCalc
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
