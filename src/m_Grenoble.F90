!$Id$
!--------------------------------------------------------------------
!  Contains all variables for the case of an imposed IC dipole
!--------------------------------------------------------------------

MODULE Grenoble
  use truncation
  IMPLICIT NONE

  LOGICAL :: lGrenoble
  REAL(kind=8) :: BIC

  !-- Magnetic field potentials of imposed inner core field:
  !   This is calculated only once and is constant in time.
  REAL(kind=8),allocatable :: b0(:)
  REAL(kind=8),allocatable :: db0(:)
  REAL(kind=8),allocatable :: ddb0(:)

  !COMMON/Grenoble/b0,db0,ddb0,BIC,lGrenoble
contains
  SUBROUTINE initialize_Grenoble
    ALLOCATE( b0(n_r_maxMag) )
    ALLOCATE( db0(n_r_maxMag) )
    ALLOCATE( ddb0(n_r_maxMag) )
  END SUBROUTINE initialize_Grenoble
END MODULE Grenoble
