#!/bin/bash
export OMP_NUM_THREADS=1
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# switch off the threading in the input file
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" inputStart.nml >inputStart_omp1.nml
sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

# First run
mpiexec -n 8 ./magic.exe  inputStart_omp1.nml

# Restart
mpiexec -n 8 ./magic.exe  input_omp1.nml

# Clean
rm *.start
rm inputStart_omp1.nml
rm input_omp1.nml
