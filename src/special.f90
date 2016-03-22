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

   !-- A toroidal boundary force of symmetry (m,m) is used.
   !-- Can be applied to inner or outer boundary.
   logical, public :: l_Ri                  !Decide whether to use forcing at all
   logical, public :: l_RiIc,l_RiMa         !Switches to decide which boundary
                                            !should be forced           
   integer, public :: m_RiIc,m_RiMa         !Order of forcing at boundaries
   real(cp), public :: amp_RiIc,omega_RiIc  !Inner boundary
   real(cp), public :: amp_RiMa,omega_RiMa  !Outer boundary

   public :: initialize_Grenoble

contains

   subroutine initialize_Grenoble

      allocate( b0(n_r_maxMag) )
      allocate( db0(n_r_maxMag) )
      allocate( ddb0(n_r_maxMag) )

   end subroutine initialize_Grenoble

end module special
