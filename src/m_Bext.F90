!$Id$
!**************************************************************
! Common block containing the external field parameters
!**************************************************************
MODULE Bext

  IMPLICIT NONE

  !-- Control parameters for the external magnetic field:
  INTEGER :: n_imp       ! Controls external field model
  REAL(kind=8) :: rrMP         ! Magnetopause radius
  REAL(kind=8) :: amp_imp      ! Amplitude of the time varying osc
  REAL(kind=8) :: expo_imp     ! Exponent for decay
  REAL(kind=8) :: bmax_imp     ! Location of maximum in g_ext/g_int
END MODULE Bext
