!$Id$
module radial_data

   use truncation, only: n_r_max
   use parallel_mod, only: rank, n_procs, nr_per_rank, nr_on_last_rank
   use logic, only: l_mag, lVerbose
 
   implicit none
 
   private
 
   integer, public :: nRstart,nRstop,nRstartMag,nRstopMag
   integer, public :: n_r_cmb,n_r_icb
 
   public :: initialize_radial_data

contains

   subroutine initialize_radial_data

      n_r_cmb=1
      n_r_icb=n_r_max

#ifdef WITH_MPI
      nR_per_rank = (n_r_max-1)/n_procs
      nRstart = n_r_cmb + rank*nR_per_rank
      nRstop  = n_r_cmb + (rank+1)*nR_per_rank - 1

      if ( rank == n_procs-1 ) then
         ! add the last point to the last process, which now has nR_per_rank+1
         ! radial points
         nRstop = nRstop+1
      end if
      nR_on_last_rank = nR_per_rank+1
#else
      nR_per_rank = n_r_max
      nR_on_last_rank = n_r_max
      nRstart = n_r_cmb
      nRstop  = n_r_icb
#endif
      if ( l_mag ) then
         nRstartMag = nRstart
         nRstopMag  = nRstop
      else
         nRstartMag = 1
         nRstopMag  = 1
      end if

      if ( lVerbose ) then
         write(*,"(4(A,I4))") "On rank ",rank," nR is in (", &
               nRstart,",",nRstop,"), nr_per_rank is ",nr_per_rank
      end if

   end subroutine initialize_radial_data
!------------------------------------------------------------------------------
end module radial_data

