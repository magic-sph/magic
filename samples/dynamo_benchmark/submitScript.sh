#!/bin/sh

#SBATCH -J MagIC_Dyna

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1  # MPI per node

#SBATCH --cpus-per-task=1    # Thread OpenMP per MPI

#SBATCH --threads-per-core=1

#SBATCH --ntasks=1           # Total number of MPI

#SBATCH --time=03:00:00

#SBATCH  --gres=gpu:1

module purge

module load craype/2.7.17 PrgEnv-cray/8.3.3

module load perftools-base/21.12.0 cce/14.0.2 craype-network-ofi craype-x86-rome cray-fftw/3.3.8.12 cray-libsci/22.06.1.1 libfabric/1.13.1 craype-accel-amd-gfx908 rocm/4.5.0 cray-mpich/8.1.13 gdb4hpc/4.14.2 perftools

export FC=ftn && export CRAY_ACC_DEBUG=1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/rocm-5.0.2/hip/lib/:/opt/rocm-4.5.0/lib/

module list

export KMP_STACKSIZE=1g

ulimit -s unlimited

# Following export is needed to bypass the FS detection fail since CRAY-MPI has not been built with XFS support on the HPE dev machine
export ROMIO_FSTYPE_FORCE="ufs:"


echo "Running on: $SLURM_NODELIST"

export TOTAL_NTASKS=$(($SLURM_NNODES*$SLURM_NTASKS_PER_NODE))

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export CRAY_ACC_DEBUG=1

export HIP_LAUNCH_BLOCKING=1

echo "Total tasks: $TOTAL_NTASKS"

srun -l ./magic.exe input.nml
