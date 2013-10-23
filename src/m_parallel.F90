!$Id$
!********************************************************************
!  Module containing blocking information
!********************************************************************
MODULE parallel_mod
#ifdef WITH_MPI
  USE MPI
#endif
  USE omp_lib
  IMPLICIT NONE

  !INTEGER :: nThreadsRun,nThreads
  !INTEGER :: nThreadsRmax,nThreadsLMmax

  INTEGER :: rank,n_procs
  INTEGER :: nr_per_rank,nr_on_last_rank
  integer :: nLMBs_per_rank
  integer :: rank_with_l1m0
#ifdef WITH_MPI
  ! a common declaration of the MPI error variable
  INTEGER :: ierr
#endif

CONTAINS
  !***********************************************************************
  SUBROUTINE parallel

    !--- Get number (name) of processor
#ifdef WITH_MPI
    CALL mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
    !WRITE(*,"(A,I3,A,I3)") "Running MPI rank no. ",rank," out of ",n_procs
#else
    rank=0
    n_procs=1
#endif

    !nThreadsRmax=0
    !nThreadsLMmax=0

  END SUBROUTINE parallel
  !------------------------------------------------------------------------
END MODULE parallel_mod
