module radial_data
   !
   ! This module defines the MPI decomposition in the radial direction.
   !

   use truncation, only: n_r_max
   use parallel_mod, only: rank, n_procs, nR_per_rank, nR_on_last_rank
   use logic, only: l_mag, lVerbose, l_finite_diff
 
   implicit none
 
   private


   type, public :: load
      integer :: nStart
      integer :: nStop
      integer :: n_per_rank
      integer :: n_points
   end type load

   type(load), public, allocatable :: radial_balance(:)
 
   integer, public :: nRstart,nRstop,nRstartMag,nRstopMag
   integer, public :: n_r_cmb,n_r_icb
 
   public :: initialize_radial_data

contains

   subroutine initialize_radial_data
      !
      ! This subroutine is used to set up the MPI decomposition in the
      ! radial direction
      !

      !-- Local variable:
      integer :: nR_remaining

      n_r_cmb=1
      n_r_icb=n_r_max

#ifdef WITH_MPI
      if ( .not. l_finite_diff ) then ! Cheb grid are restriced to odd numbers for now
         nR_per_rank = (n_r_max-1)/n_procs
         nRstart = n_r_cmb + rank*nR_per_rank
         nRstop  = n_r_cmb + (rank+1)*nR_per_rank - 1

         if ( rank == n_procs-1 ) then
            ! add the last point to the last process, which now has nR_per_rank+1
            ! radial points
            nRstop = nRstop+1
         end if
         nR_on_last_rank = nR_per_rank+1
      else ! In FD, any grid is allowed
         nR_per_rank = n_r_max/n_procs
         nRstart = n_r_cmb + rank*nR_per_rank
         nRstop  = n_r_cmb + (rank+1)*nR_per_rank - 1

         nR_remaining = n_r_max-(n_r_cmb + n_procs*nR_per_rank - 1)
         if ( rank == n_procs-1 ) then
            nRstop = nRstop+nR_remaining
         end if
         nR_on_last_rank = nR_per_rank+nR_remaining
      end if

#else
      nR_per_rank = n_r_max
      nR_on_last_rank = n_r_max
      nRstart = n_r_cmb
      nRstop  = n_r_icb
#endif
      allocate(radial_balance(0:n_procs-1))
      call getBlocks(radial_balance, n_r_max, n_procs)   

      !print*, rank, nRstart, nRstop, nR_per_rank
      !nRstart = radial_balance(rank)%nStart
      !nRstop = radial_balance(rank)%nStop
      !nR_per_rank = radial_balance(rank)%n_per_rank
      !print*, rank, nRstart, nRstop, nR_per_rank

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
   subroutine getBlocks(bal, n_points, n_procs)

      type(load), intent(inout) :: bal(0:)
      integer, intent(in) :: n_procs
      integer, intent(in) :: n_points

      integer :: n_points_loc, check, p

      n_points_loc = n_points/n_procs

      check = mod(n_points,n_procs)!-1

      bal(0)%nStart = 1

      do p =0, n_procs-1
         if ( p /= 0 ) bal(p)%nStart=bal(p-1)%nStop+1
         bal(p)%n_per_rank=n_points_loc
         if ( p == n_procs-1 ) then
            bal(p)%n_per_rank=n_points_loc+check
         end if
         bal(p)%nStop=bal(p)%nStart+bal(p)%n_per_rank-1
         bal(p)%n_points=n_points
      end do

   end subroutine getBlocks
!------------------------------------------------------------------------------
end module radial_data

