!$Id$
module Grenoble
!--------------------------------------------------------------------
!  Contains all variables for the case of an imposed IC dipole
!--------------------------------------------------------------------

  implicit none

  logical :: lGrenoble
  real(kind=8) :: BIC

  !-- Magnetic field potentials of imposed inner core field:
  !   This is calculated only once and is constant in time.
  real(kind=8),allocatable :: b0(:)
  real(kind=8),allocatable :: db0(:)
  real(kind=8),allocatable :: ddb0(:)

contains

   subroutine initialize_Grenoble

      use truncation, only: n_r_maxMag

      allocate( b0(n_r_maxMag) )
      allocate( db0(n_r_maxMag) )
      allocate( ddb0(n_r_maxMag) )

   end subroutine initialize_Grenoble

end module Grenoble
