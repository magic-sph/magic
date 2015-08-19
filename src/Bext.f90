!$Id$
module Bext
   !----------------------------------------------------------
   ! Common block containing the external field parameters
   !----------------------------------------------------------

   use precision_mod, only: cp

   implicit none
 
   !-- Control parameters for the external magnetic field:
   integer :: n_imp       ! Controls external field model
   integer :: l_imp       ! Mode of external field (dipole,quadrupole etc.)
   real(cp) :: rrMP         ! Magnetopause radius
   real(cp) :: amp_imp      ! Amplitude of the time varying osc
   real(cp) :: expo_imp     ! Exponent for decay
   real(cp) :: bmax_imp     ! Location of maximum in g_ext/g_int

   logical :: l_curr        ! Switch for current loop at the equator
   real(cp) :: amp_curr     ! Amplitude of magnetic field of current loop
   real(cp),allocatable :: fac_loop(:)  ! Array of factors for computing magnetic field for loop

end module Bext
