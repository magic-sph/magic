#!/bin/bash -l


############################################
#  This is a submit script for the DRACO   # 
#  cluster in Garching. The job scheduler  #
#  is SLURM. Can be used for other SLURM   #
#  implementations as well. Make sure that #
#  this script is executable.              #
############################################


#SBATCH -o out.%j                       # Path to output file. Default directory is ./
#SBATCH -e err.%j                       # Path to error file. Default directory is ./
#SBATCH -D ./                           # Working directory. Default is ./ 
#SBATCH -J jobname
#SBATCH --partition=general             # Defines the queue
#SBATCH --nodes=6                       # Number of nodes
#SBATCH --ntasks-per-node=2             # Tasks per node
#SBATCH --cpus-per-task=16              # OMP threads per node, must match OMP_NUM_THREADS
#SBATCH --mail-type=BEGIN,END           # Valid values are BEGIN, END, ALL, FAIL etc.
#SBATCH --mail-user=user@domain
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=16
export OMP_PLACES=cores                 # Binds threads to cores
export KMP_AFFINITY=verbose,granularity=core,compact,1   # Sets thread topology
export KMP_STACKSIZE=1g

srun path/to/magic.exe input.nml
