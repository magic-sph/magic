module Grenoble
   !
   ! This module contains all variables for the case of an imposed IC dipole
   !

   use truncation, only: n_r_maxMag
   use precision_mod
   use mem_alloc, only: bytes_allocated

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
      bytes_allocated = bytes_allocated + 3*n_r_maxMag*SIZEOF_DEF_REAL

   end subroutine initialize_Grenoble

end module Grenoble
