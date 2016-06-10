module special
   !
   ! This module contains all variables for the case of an imposed IC dipole
   !

   use truncation, only: n_r_maxMag
   use precision_mod

   implicit none

   private

   !-- Magnetic field potentials of imposed inner core field:
   !   This is calculated only once and is constant in time.
   logical, public :: lGrenoble
   real(cp), public :: BIC
   real(cp), public, allocatable :: b0(:)
   real(cp), public, allocatable :: db0(:)
   real(cp), public, allocatable :: ddb0(:)

   !-- Control parameters for the external magnetic field:
   integer, public :: n_imp       ! Controls external field model
   integer, public :: l_imp       ! Mode of external field (dipole,quadrupole etc.)
   real(cp), public :: rrMP       ! Magnetopause radius
   real(cp), public :: amp_imp    ! Amplitude of the time varying osc
   real(cp), public :: expo_imp   ! Exponent for decay
   real(cp), public :: bmax_imp   ! Location of maximum in g_ext/g_int

   logical, public :: l_curr      ! Switch for current loop at the equator
   real(cp), public :: amp_curr   ! Amplitude of magnetic field of current loop
   real(cp), public, allocatable :: fac_loop(:)  ! Array of factors for computing magnetic field for loop

   !-- A toroidal boundary force of symmetry (m,m) or (m+1,m) is used.
   !-- Can be applied to inner or outer boundary.
   logical, public :: l_Ri                          !Decide whether to use forcing at all
   integer, public :: m_RiIcAsym,m_RiMaAsym         !Order of forcing at boundaries
   integer, public :: m_RiIcSym,m_RiMaSym           !Order of forcing at boundaries
   real(cp), public :: amp_RiIcAsym,omega_RiIcAsym  !Inner boundary (eq anti symm)
   real(cp), public :: amp_RiIcSym,omega_RiIcSym    !Inner boundary (eq symm)
   real(cp), public :: amp_RiMaAsym,omega_RiMaAsym  !Outer boundary (eq anti symm)
   real(cp), public :: amp_RiMaSym,omega_RiMaSym    !Outer boundary (eq symm)

   public :: initialize_Grenoble

contains

   subroutine initialize_Grenoble

      allocate( b0(n_r_maxMag) )
      allocate( db0(n_r_maxMag) )
      allocate( ddb0(n_r_maxMag) )

   end subroutine initialize_Grenoble

end module special
