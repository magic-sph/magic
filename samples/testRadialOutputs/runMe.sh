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

sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

# Run
mpiexec -n $nmpi ./magic.exe input_omp1.nml

# Concatenate the different output file in one single e_kin.test file
cat eKinR.start eMagR.start parR.start powerR.start bLayersR.start fluxesR.start perpParR.start > e_kin.test

# Clean
rm *.start
rm input_omp1.nml
