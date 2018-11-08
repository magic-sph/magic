module mpimod
#ifdef WITH_MPI
   ! include "mpif.h"
   use mpi
#endif
end module mpimod
