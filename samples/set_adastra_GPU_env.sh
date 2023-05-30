#!/bin/bash

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-trento libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.24
module load cpe/23.02 craype-accel-amd-gfx90a cray-fftw/3.3.10.3 cray-libsci/23.02.1.1 rocm/5.3.0
module use -a ${SHAREDHOMEDIR}/modulefiles
module load rocm/5.2.3
module load hipfort
module list

export FC=ftn
export LDFLAGS=""
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.3/x86_trento/lib/"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/fftw/3.3.10.3/x86_trento/lib/
export KMP_STACKSIZE=1g
export ROMIO_FSTYPE_FORCE="ufs:"
export OMP_NUM_THREADS=1
export AMD_LOG_LEVEL=0
export CRAY_ACC_DEBUG=0
export HIP_LAUNCH_BLOCKING=1

# GPU direct communication can work in this context.
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_OFI_NIC_VERBOSE=1

