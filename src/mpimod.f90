module mpimod
#ifdef WITH_MPIF_H
   include "mpif.h"

   interface
     subroutine mpi_abort(communicator, errorcode)
       integer :: communicator, errorcode
     end subroutine
   end interface
#else
   use mpi
#endif
end module mpimod
