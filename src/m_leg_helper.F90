MODULE leg_helper_mod
  IMPLICIT NONE
  
  TYPE :: leg_helper_t
     !----- Help arrays for Legendre transform calculated in legPrepG:
     !      Parallelizatio note: these are the R-distributed versions
     !      of the field scalars.
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: dLhw, dLhdw, dLhz, dLhb, dLhj
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: vhG, vhC, dvhdrG, dvhdrC
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: bhG, bhC, cbhG, cbhC
     !----- R-distributed versions of scalar fields (see c_fields.f):
     COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE :: sR, dsR, preR, dpR
     REAL(kind=8),DIMENSION(:), ALLOCATABLE :: zAS, dzAS, ddzAS ! used in TO
     REAL(kind=8) :: omegaIC,omegaMA
     COMPLEX(kind=8),ALLOCATABLE :: bCMB(:)
   CONTAINS
     procedure :: initialize
  END TYPE leg_helper_t
CONTAINS
  SUBROUTINE initialize(this,lm_max,lm_maxMag,l_max)
    CLASS(leg_helper_t) :: this
    INTEGER,INTENT(IN) :: lm_max,lm_maxMag,l_max

    ALLOCATE( this%dLhw(lm_max) )
    ALLOCATE( this%dLhdw(lm_max) )
    ALLOCATE( this%dLhz(lm_max) )
    ALLOCATE( this%dLhb(lm_max) )
    ALLOCATE( this%dLhj(lm_max) )
    ALLOCATE( this%vhG(lm_max) )
    ALLOCATE( this%vhC(lm_max) )
    ALLOCATE( this%dvhdrG(lm_max) )
    ALLOCATE( this%dvhdrC(lm_max) )
    ALLOCATE( this%bhG(lm_maxMag) )
    ALLOCATE( this%bhC(lm_maxMag) )
    ALLOCATE( this%cbhG(lm_maxMag) )
    ALLOCATE( this%cbhC(lm_maxMag) )
    !----- R-distributed versions of scalar fields (see c_fields.f):
    ALLOCATE( this%sR(lm_max),this%dsR(lm_max) )
    ALLOCATE( this%preR(lm_max),this%dpR(lm_max) )
    ALLOCATE( this%zAS(l_max+1),this%dzAS(l_max+1),&
         & this%ddzAS(l_max+1) ) ! used in TO

    ALLOCATE( this%bCMB(lm_maxMag) )
  END SUBROUTINE initialize
END MODULE leg_helper_mod
