#!/bin/bash
# MOAB/Torque submission script for SciNet GPC (hybrid job)
#
#PBS -l nodes=4:ppn=8,walltime=00:20:00
#PBS -N test
#PBS -m abe

# This is an example batch script for the GPC cluster (Nehalem 8CPUs/core)
# This example uses 4 OpenMP threads, 2 MPI per rank and a total number of
# 8 MPI ranks (so in total 4*8=32 CPUs will be used)

# load modules
module load intel/15.0.2 gcc/4.8.1 cmake/3.4.0 intelmpi/5.0.3.048

# Env variables
export F_UFMTENDIAN=big
export KMP_STACKSIZE=1g

# Run dir
cd $PBS_O_WORKDIR

# Threads
export OMP_NUM_THREADS=4

# MPI per node
export mpi_per_node=2

# Total of MPI ranks requested
export mpi_sum=8

# PIN THE MPI DOMAINS ACCORDING TO OMP
export I_MPI_PIN_DOMAIN=omp

# EXECUTION COMMAND; -np = nodes*ppn
mpirun -ppn ${mpi_per_node} -np ${mpi_sum} ./magic.exe input.nml
