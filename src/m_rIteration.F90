MODULE rIteration_mod
  IMPLICIT NONE

  TYPE,ABSTRACT,PUBLIC :: rIteration_t
     INTEGER :: nR,nBc
     LOGICAL :: l_cour
     LOGICAL :: lTOCalc,lTOnext,lTOnext2
     LOGICAL :: lDeriv,lRmsCalc,lHelCalc,l_frame, lMagNlBc
     logical :: l_graph
     logical :: isRadialBoundaryPoint
     REAL(kind=8) :: dtrkc,dthkc

   CONTAINS
     PROCEDURE(empty_if), DEFERRED :: initialize
     PROCEDURE(empty_if), DEFERRED :: finalize
     PROCEDURE :: set_steering_variables
     PROCEDURE(do_iteration_if), DEFERRED :: do_iteration
     PROCEDURE(getType_if), deferred :: getType
  END TYPE rIteration_t
  INTERFACE 
     SUBROUTINE empty_if(this)
       IMPORT
       class(rIteration_t) :: this
     END SUBROUTINE empty_if
     SUBROUTINE do_iteration_if(this,nR,nBc,time,dt,dtLast,&
       &                 dsdt,dwdt,dzdt,dpdt,dbdt,djdt,dVxBhLM,dVSrLM, &
       &                 br_vt_lm_cmb,br_vp_lm_cmb,   &
       &                 br_vt_lm_icb,br_vp_lm_icb,&
       &                 lorentz_torque_ic,lorentz_torque_ma,&
       &                 HelLMr,Hel2LMr,HelnaLMr,Helna2LMr,uhLMr,duhLMr,gradsLMr)
       IMPORT
       CLASS(rIteration_t) :: this
       INTEGER,INTENT(IN) :: nR,nBc
       REAL(kind=8),INTENT(IN) :: time,dt,dtLast
    
       COMPLEX(kind=8),INTENT(OUT),DIMENSION(:) :: dwdt,dzdt,dpdt,dsdt,dVSrLM
       COMPLEX(kind=8),INTENT(OUT),DIMENSION(:) :: dbdt,djdt,dVxBhLM
       !---- Output of nonlinear products for nonlinear
       !     magnetic boundary conditions (needed in s_updateB.f):
       COMPLEX(kind=8),INTENT(OUT) :: br_vt_lm_cmb(:) ! product br*vt at CMB
       COMPLEX(kind=8),INTENT(OUT) :: br_vp_lm_cmb(:) ! product br*vp at CMB
       COMPLEX(kind=8),INTENT(OUT) :: br_vt_lm_icb(:) ! product br*vt at ICB
       COMPLEX(kind=8),INTENT(OUT) :: br_vp_lm_icb(:) ! product br*vp at ICB
       REAL(kind=8),INTENT(OUT) :: lorentz_torque_ic,lorentz_torque_ma
       REAL(kind=8),INTENT(OUT),DIMENSION(:) :: HelLMr,Hel2LMr,HelnaLMr,Helna2LMr
       REAL(kind=8),INTENT(OUT),DIMENSION(:) :: uhLMr,duhLMr,gradsLMr

     END SUBROUTINE do_iteration_if

     FUNCTION getType_if(this)
       import
       class(rIteration_t) :: this
       character(len=100) :: getType_if
     END FUNCTION getType_if
  END INTERFACE
CONTAINS
  SUBROUTINE set_steering_variables(this,l_cour,lTOCalc,lTOnext,lTOnext2,&
       & lDeriv,lRmsCalc,lHelCalc,l_frame,lMagNlBc,l_graph)
    class(rIteration_t) :: this
    LOGICAL,INTENT(IN) :: l_cour,lDeriv,lRmsCalc,lHelCalc,l_frame
    LOGICAL,INTENT(IN) :: lTOCalc,lTOnext,lTOnext2, lMagNlBc,l_graph
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
  END SUBROUTINE set_steering_variables
END MODULE rIteration_mod

MODULE rIteration_boundary_mod
  use rIteration_mod
  IMPLICIT NONE

  TYPE,ABSTRACT,PUBLIC,EXTENDS(rIteration_t) :: rIteration_boundary_t
  END TYPE rIteration_boundary_t
END MODULE rIteration_boundary_mod

MODULE rIteration_inner_mod
  use rIteration_mod
  implicit none

  TYPE,ABSTRACT,PUBLIC,EXTENDS(rIteration_t) :: rIteration_inner_t
  END TYPE rIteration_inner_t
END MODULE rIteration_inner_mod

