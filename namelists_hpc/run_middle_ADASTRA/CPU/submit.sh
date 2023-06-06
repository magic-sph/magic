#!/bin/sh

#SBATCH -J MagIC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3    # Thread OpenMP per MPI
#SBATCH --ntasks=64           # Total number of MPI
#SBATCH --threads-per-core=1
#SBATCH --time=00:50:00
#SBATCH --exclusive
#SBATCH --constraint=GENOA
#SBATCH -A cpa
#SBATCH --output=Magic_1N64MPI3OMP_%j.out

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-genoa libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.24
module load cpe/23.02 cray-fftw/3.3.10.3 cray-libsci/23.02.1.1
module list

export KMP_STACKSIZE=1g

ulimit -s unlimited

# Following export is needed to bypass the FS detection fail since CRAY-MPI has not been built with XFS support on the HPE dev machine
export ROMIO_FSTYPE_FORCE="ufs:"

echo "Running on: $SLURM_NODELIST"
echo "Total tasks: $TOTAL_NTASKS"

export OMP_NUM_THREADS=3

srun -l --mpi=cray_shasta -c ${SLURM_CPUS_PER_TASK} --cpu-bind=verbose,cores  ./magic.exe input_middle.nml
