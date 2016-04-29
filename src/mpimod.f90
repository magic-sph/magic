module mpimod
#ifdef WITH_MPIF_H
   include "mpif.h"
#else
   use mpi
#endif
end module mpimod
