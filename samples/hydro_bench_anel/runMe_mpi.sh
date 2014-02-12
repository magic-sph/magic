#!/bin/bash
export OMP_NUM_THREADS=1
export F_UFMTENDIAN=big
export KMP_STACKSIZE=500m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# switch off the threading in the input file
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

max_procs=`cat /proc/cpuinfo | grep processor | wc -l`
echo "Starting magic with ${max_procs} MPI processes"
mpiexec -n ${max_procs} ./magic.exe input_omp1.nml

rm input_omp1.nml
