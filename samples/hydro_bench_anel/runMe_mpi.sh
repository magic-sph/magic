#!/bin/bash
export OMP_NUM_THREADS=1
export F_UFMTENDIAN=big
export KMP_STACKSIZE=500m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# switch off the threading in the input file
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

max_procs=`cat /proc/cpuinfo | grep processor | wc -l`
if [ "$max_procs" -ge 16 ]; then
    n_procs=16
elif [ "$max_procs" -ge 8 ]; then
    n_procs=8
fi

echo "Starting magic with ${n_procs} MPI processes"
mpiexec -n ${n_procs} ./magic.exe input_omp1.nml

rm input_omp1.nml
