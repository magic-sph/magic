module special
   !
   ! This module contains all variables for the case of an imposed IC dipole,
   ! an imposed external magnetic field and a special boundary forcing to excite
   ! inertial modes

   use precision_mod

   implicit none

   private

   !-- Magnetic field potentials of imposed inner core field:
   !   This is calculated only once and is constant in time.
   logical, public  :: lGrenoble
   real(cp), public :: BIC

   !-- Control parameters for the external magnetic field:
   integer, public  :: n_imp       ! Controls external field model
   integer, public  :: l_imp       ! Mode of external field (dipole,quadrupole etc.)
   real(cp), public :: rrMP       ! Magnetopause radius
   real(cp), public :: amp_imp    ! Amplitude of the time varying osc
   real(cp), public :: expo_imp   ! Exponent for decay
   real(cp), public :: bmax_imp   ! Location of maximum in g_ext/g_int

   logical, public  :: l_curr      ! Switch for current loop at the equator
   real(cp), public :: Le         ! Lehnert number defined by the magnetic field at the centre
   real(cp), public :: amp_curr   ! Amplitude of magnetic field of current loop to be scaled by Lehnert
   real(cp), public :: loopRadRatio   ! Radius ratio of outer boundary/current loop
   real(cp), public, allocatable :: fac_loop(:)  ! Array of factors for computing magnetic field for loop

   !-- Parameters for a toroidal boundary forcing of symmetry (m,m) or (m+1,m)
   !-- to excite inertial modes of desired symmetry
   !-- Can be applied to inner or outer boundary
   integer, public  :: m_mode_ic, m_mode_ma         !Order of forcing
   real(cp), public :: amp_mode_ic, amp_mode_ma     !Amplitude of forcing
   real(cp), public :: omega_mode_ic, omega_mode_ma !Frequency of forcing
   integer, public  :: mode_symm_ic, mode_symm_ma   !Symmetry of forcing: 1 : eq symm, 0 : eq antisymm

end module special
