#!/bin/bash
export OMP_NUM_THREADS=8
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# First run
./magic.exe  inputStart.nml

# Restart
./magic.exe  input.nml

# Clean
rm *.start
