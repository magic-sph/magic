!$Id$
module Grenoble
!--------------------------------------------------------------------
!  Contains all variables for the case of an imposed IC dipole
!--------------------------------------------------------------------

  use truncation, only: n_r_maxMag
  use precision_mod, only: cp

  implicit none

  private

  logical, public :: lGrenoble
  real(cp), public :: BIC

  !-- Magnetic field potentials of imposed inner core field:
  !   This is calculated only once and is constant in time.
  real(cp), public, allocatable :: b0(:)
  real(cp), public, allocatable :: db0(:)
  real(cp), public, allocatable :: ddb0(:)

  public :: initialize_Grenoble

contains

   subroutine initialize_Grenoble

      allocate( b0(n_r_maxMag) )
      allocate( db0(n_r_maxMag) )
      allocate( ddb0(n_r_maxMag) )

   end subroutine initialize_Grenoble

end module Grenoble
