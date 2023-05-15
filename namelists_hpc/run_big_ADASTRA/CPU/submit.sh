#!/bin/sh

#SBATCH -J MagIC

#SBATCH --nodes=8

#SBATCH --cpus-per-task=8    # Thread OpenMP per MPI

#SBATCH --ntasks=192           # Total number of MPI

#SBATCH --ntasks-per-node=24

#SBATCH --threads-per-core=1

#SBATCH --time=00:50:00

#SBATCH --exclusive

#SBATCH --constraint=GENOA

#SBATCH -A cpa

#SBATCH --reservation=eolen_cpu

##SBATCH --kill-on-bad-exit

#SBATCH --output=Magic_8N192MPI8OMP_%j.out

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-genoa libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.21
module load cpe/22.11 cray-fftw/3.3.10.3 cray-libsci/22.11.1.2
module list

export KMP_STACKSIZE=2g

ulimit -s unlimited

# Following export is needed to bypass the FS detection fail since CRAY-MPI has not been built with XFS support on the HPE dev machine
export ROMIO_FSTYPE_FORCE="ufs:"

echo "Running on: $SLURM_NODELIST"
echo "Total tasks: $TOTAL_NTASKS"

export OMP_NUM_THREADS=8

srun --kill-on-bad-exit -l --mpi=cray_shasta -c ${SLURM_CPUS_PER_TASK} --cpu-bind=verbose,cores  ./magic.exe input_big.nml
