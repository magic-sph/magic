module precision_mod
   !
   ! This module controls the precision used in MagIC
   !

#ifdef WITH_MPI
   use mpi
#endif

   implicit none

   private

   !-- Current precision for calculations

#if (DEFAULT_PRECISION==sngl)
   integer, public, parameter :: cp=selected_real_kind(6)
#elif (DEFAULT_PRECISION==dble)
   integer, public, parameter :: cp=selected_real_kind(15)
!#elif (DEFAULT_PRECISION==quad)
!   integer, public, parameter :: cp=selected_real_kind(20)
#endif
   !-- Precision for outputs in unformatted files (G files, movie files)
#if (DEFAULT_OUTPUT_PRECISION==sngl)
   integer, public, parameter :: outp=selected_real_kind(6)
#elif (DEFAULT_OUTPUT_PRECISION==dble)
   integer, public, parameter :: outp=selected_real_kind(15)
#endif

   !-- SIZEOFs
   integer, public, parameter :: SIZEOF_DEF_COMPLEX=2*cp
   integer, public, parameter :: SIZEOF_DEF_REAL=cp
   integer, public, parameter :: SIZEOF_OUT_REAL=outp

   !-- Precision for long integers
   integer, public, parameter :: lip=selected_int_kind(12)
   integer, public, parameter :: SIZEOF_INTEGER=4
   integer, public, parameter :: SIZEOF_CHARACTER=1

#ifdef WITH_MPI
   !-- Precision for MPI communications

   !-- Real numbers
#if (DEFAULT_PRECISION==sngl)
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL4
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX8
#elif (DEFAULT_PRECISION==dble)
   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL8
   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX16
!#elif (DEFAULT_PRECISION==quad)
!   integer, public, parameter :: MPI_DEF_REAL=MPI_REAL16
!   integer, public, parameter :: MPI_DEF_COMPLEX=MPI_COMPLEX32
#endif

   !-- Output
#if (DEFAULT_OUTPUT_PRECISION==sngl)
   integer, public, parameter :: MPI_OUT_REAL=MPI_REAL4
#elif (DEFAULT_OUTPUT_PRECISION==dble)
   integer, public, parameter :: MPI_OUT_REAL=MPI_REAL8
#endif
#endif

end module precision_mod
