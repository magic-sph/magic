module parallel_mod
   !
   !  This module contains the blocking information
   !

#ifdef WITH_MPI
   use MPI
#endif
   use omp_lib

   implicit none

   integer :: nThreads
   integer :: rank,n_procs
   integer :: nR_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr

   type, public :: load
      integer :: nStart
      integer :: nStop
      integer :: n_per_rank
   end type load

   public :: getBlocks, get_openmp_blocks

contains

   subroutine parallel

      !--- Get number (name) of processor
#ifdef WITH_MPI
      call MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD,n_procs, ierr)
      !write(*,"(A,I3,A,I3)") "Running MPI rank no. ",rank," out of ",n_procs
#else
      rank    = 0
      n_procs = 1
      ierr    = 0
#endif

#ifdef WITHOMP
      nThreads = omp_get_max_threads()
#else
      nThreads = 1
#endif

      chunksize=16

   end subroutine parallel
!------------------------------------------------------------------------------
   subroutine check_MPI_error(code)

      integer, intent(in) :: code
#ifdef WITH_MPI
      character(len=MPI_MAX_ERROR_STRING) :: error_str
      integer :: ierr, strlen
#endif

      if (code /= 0) then
#ifdef WITH_MPI
          call MPI_Error_string(code, error_str, strlen, ierr)
          write(*, '(A, A)') 'MPI error: ', trim(error_str)
          call MPI_Abort(MPI_COMM_WORLD, code, ierr)
#else
          write(*, '(A, I4)') 'Error code: ', code
          stop
#endif
      endif

   end subroutine check_MPI_error
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
      end do

   end subroutine getBlocks
!------------------------------------------------------------------------------
   subroutine get_openmp_blocks(nStart, nStop)

      !--Input/Outputs variables:
      integer, intent(inout) :: nStart
      integer, intent(inout) :: nStop

      !-- Local variables
      integer :: n_threads, threadid, n_points_per_thread, n_points_left
      integer :: n_points, n_glob_start
#ifdef WITHOMP
      integer :: n_max_threads
#endif

      n_points=nStop-nStart+1
      n_glob_start=nStart

#ifdef WITHOMP
      n_threads=omp_get_num_threads()
      threadid =omp_get_thread_num()
      if ( n_points < n_threads) then
         call omp_set_num_threads(n_points)
         n_points_per_thread=1
         n_points_left=0
      else
         n_points_per_thread=n_points/n_threads
         n_points_left=n_points-n_threads*n_points_per_thread
      end if
#else
      n_threads=1
      threadid =0
      n_points_per_thread=n_points
      n_points_left=0
#endif

      !-- This is a way to reshuffle the points which are not in-balance
      !-- more regularly
      if ( threadid+1 <= n_points_left ) then
         nStart = n_glob_start+threadid*n_points_per_thread+threadid
         nStop  = nStart+n_points_per_thread
      else
         nStart = n_glob_start+threadid*n_points_per_thread+n_points_left
         nStop  = nStart+n_points_per_thread-1
      end if

#ifdef WITHOMP
      if ( n_points < n_threads) then
         n_max_threads=omp_get_max_threads()
         call omp_set_num_threads(n_max_threads)
      end if
#endif

   end subroutine get_openmp_blocks
!------------------------------------------------------------------------------
end module parallel_mod
