module parallel_mod
   !
   !  This module contains the blocking information
   !

#ifdef WITH_MPI
   use mpimod
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
!------------------------------------------------------------------------------
end module parallel_mod
