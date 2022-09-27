#!/bin/sh

module load craype/2.7.17 PrgEnv-cray/8.3.3

module load perftools-base/21.12.0 cce/14.0.2 craype-network-ofi craype-x86-rome cray-fftw/3.3.8.12 cray-libsci/22.06.1.1 libfabric/1.13.1 craype-accel-amd-gfx908 rocm/4.5.0 cray-mpich/8.1.13 gdb4hpc/4.14.2 perftools

export FC=ftn && export CRAY_ACC_DEBUG=1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/rocm-5.0.2/hip/lib/:/opt/rocm-4.5.0/lib/

if [ -e "build" ];then rm -rf "build" ; fi

mkdir -p build

cd build

cmake .. -DUSE_SHTNS=yes -DUSE_GPU=yes

make -j16

cp magic.exe ../samples/boussBenchSat

cp magic.exe ../samples/dynamo_benchmark

cd ../samples/boussBenchSat

sbatch submitScript.sh

cd ../dynamo_benchmark

sbatch submitScript.sh

cd ../boussBenchSat

squeue
