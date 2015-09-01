#!/bin/bash
if [ "$#" -ge "1" ]; then
    hybrid=$1
    nomp=$2
else
    hybrid="no"
fi

export F_UFMTENDIAN=big
if [ "$hybrid" == "hybrid" ]; then
    nmpi=2
    #nomp=4
    echo "Running hybrid mode with $nomp OpenMP threads per $nmpi MPI processes."
    export I_MPI_PIN_DOMAIN=socket
    # MP_BINDPROC must be set to no on hydra
    export MP_BINDPROC=no
    export MP_TASK_AFFINITY=core:$nomp
    export KMP_STACKSIZE=500m
    export KMP_AFFINITY=noverbose,granularity=core,compact
else
    nomp=1
    nmpi=8
    echo "Running pure MPI code with $nmpi MPI processes."
    export I_MPI_PIN_PROCESSOR_LIST=allcores
fi
export OMP_NUM_THREADS=$nomp

which_poe=`which poe 2>/dev/null`
if [ -x "$which_poe" ]; then
    poe ./magic.exe input.nml -procs $nmpi >stdout.out
else
    mpiexec -n $nmpi ./magic.exe input.nml >stdout.out
fi

