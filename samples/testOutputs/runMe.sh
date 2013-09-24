#!/bin/bash
export OMP_NUM_THREADS=8
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# Run
./magic.exe input.nml

# Concatenate the different output file in one single e_kin.test file
cat e_kin.start e_mag_oc.start e_mag_ic.start dipole.start misc.start par.start power.start u_square.start > e_kin.test

# Clean
rm *.start
