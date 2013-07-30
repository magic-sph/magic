#!/bin/bash
export OMP_NUM_THREADS=8
export F_UFMTENDIAN=big
export KMP_STACKSIZE=500m
export KMP_AFFINITY=noverbose,granularity=core,scatter
./magic.exe input.nml
