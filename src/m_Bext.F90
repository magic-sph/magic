!$Id$
module Bext
   !----------------------------------------------------------
   ! Common block containing the external field parameters
   !----------------------------------------------------------

   implicit none
 
   !-- Control parameters for the external magnetic field:
   integer :: n_imp       ! Controls external field model
   integer :: l_imp       ! Mode of external field (dipole,quadrupole etc.)
   real(kind=8) :: rrMP         ! Magnetopause radius
   real(kind=8) :: amp_imp      ! Amplitude of the time varying osc
   real(kind=8) :: expo_imp     ! Exponent for decay
   real(kind=8) :: bmax_imp     ! Location of maximum in g_ext/g_int

end module Bext
