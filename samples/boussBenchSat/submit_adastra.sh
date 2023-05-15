#!/bin/sh

#SBATCH -J MagIC

#SBATCH --nodes=1

#SBATCH --cpus-per-task=32    # Thread OpenMP per MPI

#SBATCH --ntasks=4           # Total number of MPI

#SBATCH --time=00:10:00

#SBATCH --exclusive

#SBATCH --gpus-per-task=1

#SBATCH --constraint=MI250

#SBATCH -A cpa

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-trento libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.21
module load cpe/22.11 craype-accel-amd-gfx90a cray-fftw/3.3.10.2 cray-libsci/22.11.1.2 rocm/5.2.3
module use -a ${SHAREDHOMEDIR}/modulefiles
module load rocm/5.2.3
module load hipfort
export HIPFORT_PATH=/lus/home/NAT/cpa2204/SHARED/INSTALL/hipfort/0.4-6f6ae98e/cpe-cray/15.0.0/rocm/5.2.3/cmpich/8.1.21

export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.2/x86_trento/lib/"

export KMP_STACKSIZE=1g

ulimit -s unlimited

# Following export is needed to bypass the FS detection fail since CRAY-MPI has not been built with XFS support on the HPE dev machine
export ROMIO_FSTYPE_FORCE="ufs:"

echo "Running on: $SLURM_NODELIST"
echo "Total tasks: $TOTAL_NTASKS"

export OMP_NUM_THREADS=1
export AMD_LOG_LEVEL=2
export CRAY_ACC_DEBUG=1
export HIP_LAUNCH_BLOCKING=1
export CRAY_ACC_DEBUG=1

srun -l --mpi=cray_shasta -c 32 --cpu-bind=verbose,cores --gpu-bind=verbose,closest ./magic.exe input.nml
