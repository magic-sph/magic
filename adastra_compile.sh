#!/bin/sh

### Env for ADASTRA

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-trento libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.21
module load cpe/22.11 craype-accel-amd-gfx90a cray-fftw/3.3.10.2 cray-libsci/22.11.1.2 rocm/5.2.3

export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.2/x86_trento/lib/"
export HIP_LAUNCH_BLOCKING=1
export AMD_LOG_LEVEL=1
export OMP_NUM_THREADS=1
export FC=ftn
export CRAY_ACC_DEBUG=0

module use -a ${SHAREDHOMEDIR}/modulefiles
module load rocm/5.2.3
module load hipfort
export HIPFORT_PATH=/lus/home/NAT/cpa2204/SHARED/INSTALL/hipfort/0.4-6f6ae98e/cpe-cray/15.0.0/rocm/5.2.3/cmpich/8.1.21

### Create a build directory
if [ -e "build" ];then rm -rf "build" ; fi
mkdir -p build
cd build

### Compile (assume SHTns is installed in $HOME/local
cmake .. -DUSE_SHTNS=yes -DUSE_GPU=yes
make -j16

### Submit a job for the boussBenchSat sample
cd ../samples/boussBenchSat
chmod +x clear.sh
./clear.sh
cp ../../build/magic.exe .
sbatch submit_adastra.sh

squeue

cd -
