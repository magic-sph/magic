module special
   !
   ! This module contains all variables for the case of an imposed IC dipole,
   ! an imposed external magnetic field and a special boundary forcing to excite
   ! inertial modes

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
   real(cp), public :: loopRadRatio   ! Radius ratio of outer boundary/current loop
   real(cp), public, allocatable :: fac_loop(:)  ! Array of factors for computing magnetic field for loop

   !-- Parameters for a toroidal boundary forcing of symmetry (m,m) or (m+1,m)
   !-- to excite inertial modes of desired symmetry
   !-- Can be applied to inner or outer boundary
   integer, public :: m_RiIc, m_RiMa !Order of forcing 
   real(cp), public :: amp_RiIc, omega_RiIc  !Amplitude and frequency of forcing at the inner boundary
   real(cp), public :: amp_RiMa, omega_RiMa  !Amplitude and frequency of forcing at the outer boundary
   integer, public :: RiSymmMa, RiSymmIc     !Symmetry of forcing: 1 (0) for eq symm (antisymm)

   public :: initialize_Grenoble, finalize_Grenoble

contains

   subroutine initialize_Grenoble

      allocate( b0(n_r_maxMag) )
      allocate( db0(n_r_maxMag) )
      allocate( ddb0(n_r_maxMag) )

   end subroutine initialize_Grenoble
!------------------------------------------------------------------------------
   subroutine finalize_Grenoble

      deallocate( b0, db0, ddb0 )

   end subroutine finalize_Grenoble
!------------------------------------------------------------------------------
end module special
