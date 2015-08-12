!$Id$
module precision_mod
   !--------------------------------------------------
   ! This module controls the precision used in MagIC
   !--------------------------------------------------

   implicit none

   !-- Current precision for calculations
   integer, parameter :: cp=selected_real_kind(15)
   !-- Precision for long integers
   integer, parameter :: lip=selected_int_kind(12)
   !-- Precision for outputs in unformatted files (G files, movie files)
   integer, parameter :: outp=selected_real_kind(4)

end module precision_mod
