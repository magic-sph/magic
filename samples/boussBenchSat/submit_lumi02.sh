#!/bin/sh

#SBATCH -J MagIC

#SBATCH --nodes=1

#SBATCH --cpus-per-task=1    # Thread OpenMP per MPI

#SBATCH --ntasks=2           # Total number of MPI

#SBATCH --time=00:10:00

#SBATCH --gpus-per-task=1

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-rome libfabric/1.13.1 craype-network-ofi cray-mrnet/5.0.4 cray-mpich/8.1.21
module load cce/15.0.0 craype-accel-amd-gfx90a cray-fftw/3.3.10.1 cray-libsci/21.08.1.2 rocm/5.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/mpich/8.1.21/gtl/lib/
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.1/x86_rome/lib/"
module use -a /common/magic/modulefiles/
module load hipfort
export HIPFORT_PATH=/common/magic/INSTALL/hipfort/0.4-6f6ae98e/cpe-cray/15.0.0/rocm/5.2.0/gfx90a/cmpich/8.1.21

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

srun -l -c 1 ./magic.exe input.nml
