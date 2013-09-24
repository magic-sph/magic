#!/bin/bash
export OMP_NUM_THREADS=8
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

sed "s/nThreadsRun *= *[0-9]*/nThreadsRun = 1/" input.nml >input_omp1.nml

# Run
mpiexec -n 8 ./magic.exe input_omp1.nml

# Concatenate the different output file in one single e_kin.test file
cat e_kin.start e_mag_oc.start e_mag_ic.start dipole.start misc.start par.start power.start u_square.start > e_kin.test
#mv dtBrms.start myBrms.dat

# Clean
rm *.start
rm input_omp1.nml

