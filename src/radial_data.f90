module radial_data
   !
   ! This module defines the MPI decomposition in the radial direction.
   !

   use parallel_mod, only: rank, n_procs, nR_per_rank, load, getBlocks
   use logic, only: l_mag, lVerbose, l_finite_diff
 
   implicit none
 
   private

   type(load), public, allocatable :: radial_balance(:)
 
   integer, protected, public :: nRstart,nRstop,nRstartMag,nRstopMag
   integer, protected, public :: n_r_cmb,n_r_icb
 
   public :: initialize_radial_data, finalize_radial_data

contains

   subroutine initialize_radial_data(n_r_max)
      !
      ! This subroutine is used to set up the MPI decomposition in the
      ! radial direction
      !
      integer, intent(in) :: n_r_max ! Number of radial grid points

      n_r_cmb=1
      n_r_icb=n_r_max

      allocate(radial_balance(0:n_procs-1))
      call getBlocks(radial_balance, n_r_max, n_procs)   

      nRstart = radial_balance(rank)%nStart
      nRstop = radial_balance(rank)%nStop
      nR_per_rank = radial_balance(rank)%n_per_rank

      if ( l_mag ) then
         nRstartMag = nRstart
         nRstopMag  = nRstop
      else
         nRstartMag = 1
         nRstopMag  = 1
      end if

      if ( lVerbose ) then
         write(*,"(4(A,I4))") "On rank ",rank," nR is in (", &
               nRstart,",",nRstop,"), nR_per_rank is ",nR_per_rank
      end if

   end subroutine initialize_radial_data
!------------------------------------------------------------------------------
   subroutine finalize_radial_data

      deallocate( radial_balance )

   end subroutine finalize_radial_data
!------------------------------------------------------------------------------
end module radial_data

