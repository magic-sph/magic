#!/bin/bash
if [ $# -eq 1 ]; then
    hybrid=$1
else
    hybrid="no"
fi

export F_UFMTENDIAN=big
if [ "$hybrid" == "hybrid" ]; then
    nmpi=2
    nomp=4
    echo "Running hybrid mode with $nomp OpenMP threads per $nmpi MPI processes."
    export I_MPI_PIN_DOMAIN=socket
    export KMP_STACKSIZE=500m
    export KMP_AFFINITY=noverbose,granularity=core,compact
else
    nomp=1
    nmpi=8
    echo "Running pure MPI code with $nmpi MPI processes."
    export I_MPI_PIN_PROCESSOR_LIST=allcores
fi
export OMP_NUM_THREADS=$nomp

# switch off the threading in the input file
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" inputStart.nml >inputStart_omp1.nml
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

# First run
mpiexec -n $nmpi ./magic.exe  inputStart_omp1.nml

# Restart
mpiexec -n $nmpi ./magic.exe  input_omp1.nml

# Clean
rm *.start
rm inputStart_omp1.nml
rm input_omp1.nml
