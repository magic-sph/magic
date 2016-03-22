#!/bin/bash
if [ $# -ge 1 ]; then
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
    export MP_BINDPROC=no
    export MP_TASK_AFFINITY=core:$nomp
    export KMP_STACKSIZE=1g
    export KMP_AFFINITY=verbose,granularity=core,compact
else
    nomp=1
    nmpi=8
    echo "Running pure MPI code with $nmpi MPI processes."
    export I_MPI_PIN_PROCESSOR_LIST=allcores
fi
export OMP_NUM_THREADS=$nomp

# Run
if grep -Fxq "USE_MPI=yes" $MAGIC_HOME/src/Makefile; then
    mpiexec -n $nmpi ./magic.exe input.nml
else
    ./magic.exe input.nml
fi

# Read coeff files
python getData.py

# Clean
rm *.start
