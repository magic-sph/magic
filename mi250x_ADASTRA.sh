#!/bin/sh
module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-trento libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.21
module load cpe/22.11 craype-accel-amd-gfx90a cray-fftw/3.3.10.3 cray-libsci/22.11.1.2 rocm/5.2.3
module use -a ${SHAREDHOMEDIR}/modulefiles
module load rocm/5.2.3
module load hipfort
module list

export HIPFORT_PATH=/lus/home/NAT/cpa2204/SHARED/INSTALL/hipfort/0.4-6f6ae98e/cpe-cray/15.0.0/rocm/5.2.3/cmpich/8.1.21
export LDFLAGS+=""
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.3/x86_genoa/lib/"
export HIP_LAUNCH_BLOCKING=1
export OMP_NUM_THREADS=1
export AMD_LOG_LEVEL=0 #For debug
export CRAY_ACC_DEBUG=0  #For debug
export FC=ftn

if [ -e "build" ];then rm -rf "build" ; fi
mkdir build && cd build

cmake .. -DUSE_SHTNS=yes -DUSE_GPU=yes && make -j16

source ../namelists_hpc/run_middle_ADASTRA/GPU/ clear.sh
source ../namelists_hpc/run_big_ADASTRA/GPU/ clear.sh
source ../samples/boussBenchSat/ clear.sh

cp magic.exe ../namelists_hpc/run_middle_ADASTRA/GPU
cp magic.exe ../namelists_hpc/run_big_ADASTRA/GPU
cp magic.exe ../samples/boussBenchSat

cd ../namelists_hpc/run_middle_ADASTRA/GPU
sbatch submit.sh

cd ../../run_big_ADASTRA/GPU
sbatch submit.sh

cd ../../../samples/boussBenchSat/
sbatch genoa_ADATRA_submit.sh

squeue | grep MagIC
