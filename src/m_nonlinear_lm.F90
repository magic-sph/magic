#include "intrinsic_sizes.h"
MODULE nonlinear_lm_mod
  !USE cutils
  implicit none

  TYPE :: nonlinear_lm_t
     !----- Nonlinear terms in lm-space: 
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: AdvrLM, AdvtLM, AdvpLM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: LFrLM,  LFtLM,  LFpLM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: VxBrLM, VxBtLM, VxBpLM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: VSrLM,  VStLM,  VSpLM
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: ViscHeatLM, OhmLossLM
   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: output
     procedure :: set_zero
  END TYPE nonlinear_lm_t

contains
  SUBROUTINE initialize(this,lmP_max)
    CLASS(nonlinear_lm_t) :: this
    INTEGER,intent(IN) :: lmP_max
    integer :: size_in_bytes

    ALLOCATE( this%AdvrLM(lmP_max) )   
    ALLOCATE( this%AdvtLM(lmP_max) )   
    ALLOCATE( this%AdvpLM(lmP_max) )   
    ALLOCATE( this%LFrLM(lmP_max) )    
    ALLOCATE( this%LFtLM(lmP_max) )    
    ALLOCATE( this%LFpLM(lmP_max) )    
    ALLOCATE( this%VxBrLM(lmP_max) )   
    ALLOCATE( this%VxBtLM(lmP_max) )   
    ALLOCATE( this%VxBpLM(lmP_max) )   
    ALLOCATE( this%VSrLM(lmP_max) )    
    ALLOCATE( this%VStLM(lmP_max) )    
    ALLOCATE( this%VSpLM(lmP_max) )    
    ALLOCATE( this%ViscHeatLM(lmP_max) )
    ALLOCATE( this%OhmLossLM(lmP_max) )
    !size_in_bytes=14*lmP_max*SIZEOF_DOUBLE_COMPLEX
    !WRITE(*,"(A,I15,A)") "nonlinear_lm: allocated ",size_in_bytes,"B."
    !CALL this%set_zero()

    !WRITE(*,"(A,I5,A,I10,A)") "cache info for first element, size is lmP_max=",lmP_max,&
    !     &" = ",lmP_max*16,"B"
    !CALL print_cache_info_dcmplx("nl_lm%AdvrLM"//C_NULL_CHAR,this%AdvrLM(1))
    !CALL print_cache_info_dcmplx("nl_lm%AdvtLM"//C_NULL_CHAR,this%AdvtLM(1))
    !CALL print_cache_info_dcmplx("nl_lm%AdvpLM"//C_NULL_CHAR,this%AdvpLM(1))

    !CALL print_address("nl_lm%VStLM"//C_NULL_CHAR,this%VStLM(1))
    !CALL print_address("nl_lm%VSpLM"//C_NULL_CHAR,this%VSpLM(1))
  END SUBROUTINE initialize

  SUBROUTINE finalize(this)
    CLASS(nonlinear_lm_t) :: this
    integer :: size_in_bytes

    DEALLOCATE( this%AdvrLM )   
    DEALLOCATE( this%AdvtLM )   
    DEALLOCATE( this%AdvpLM )   
    DEALLOCATE( this%LFrLM )    
    DEALLOCATE( this%LFtLM )    
    DEALLOCATE( this%LFpLM )    
    DEALLOCATE( this%VxBrLM )   
    DEALLOCATE( this%VxBtLM )   
    DEALLOCATE( this%VxBpLM )   
    DEALLOCATE( this%VSrLM )    
    DEALLOCATE( this%VStLM )    
    DEALLOCATE( this%VSpLM )    
    DEALLOCATE( this%ViscHeatLM )
    DEALLOCATE( this%OhmLossLM )
  END SUBROUTINE finalize

  SUBROUTINE set_zero(this)
    CLASS(nonlinear_lm_t) :: this
    
    this%AdvrLM=cmplx(0.0,0.0,kind=8)   
    this%AdvtLM=cmplx(0.0,0.0,kind=8)   
    this%AdvpLM=cmplx(0.0,0.0,kind=8)   
    this%LFrLM=cmplx(0.0,0.0,kind=8)    
    this%LFtLM=cmplx(0.0,0.0,kind=8)    
    this%LFpLM=cmplx(0.0,0.0,kind=8)    
    this%VxBrLM=cmplx(0.0,0.0,kind=8)   
    this%VxBtLM=cmplx(0.0,0.0,kind=8)   
    this%VxBpLM=cmplx(0.0,0.0,kind=8)   
    this%VSrLM=cmplx(0.0,0.0,kind=8)    
    this%VStLM=cmplx(0.0,0.0,kind=8)    
    this%VSpLM=cmplx(0.0,0.0,kind=8)    
    this%ViscHeatLM=CMPLX(0.0,0.0,kind=8)
    this%OhmLossLM=cmplx(0.0,0.0,kind=8)
  END SUBROUTINE set_zero

  SUBROUTINE output(this)
    CLASS(nonlinear_lm_t) :: this
    
    WRITE(*,"(A,6ES20.12)") "AdvrLM,AdvtLM,AdvpLM = ",&
         & SUM(this%AdvrLM),&
         & SUM(this%AdvtLM),&
         & SUM(this%AdvpLM)
  END SUBROUTINE output
END MODULE nonlinear_lm_mod
