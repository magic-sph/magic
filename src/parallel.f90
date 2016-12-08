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
   integer :: nr_per_rank,nr_on_last_rank
   integer :: nLMBs_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr

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
       character(len=MPI_MAX_ERROR_STRING) :: error_str
       integer :: ierr, strlen

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
end module parallel_mod
