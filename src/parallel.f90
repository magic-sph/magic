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
   integer :: nLMBs_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr


   type, public :: load
      integer :: nStart
      integer :: nStop
      integer :: n_per_rank
      integer :: n_points
   end type load

   public :: getBlocks

contains

   subroutine parallel

      !--- Get number (name) of processor
#ifdef WITH_MPI
      call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
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
         bal(p)%n_points=n_points
      end do

   end subroutine getBlocks
!------------------------------------------------------------------------------
end module parallel_mod
