#############################################################
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
#$ -M $USER@ipp.mpg.de 
#$ -N MAGIC_test
#$ -pe impi_hydra_cuda 16
#$ -l h_rt=01:00:00

export MACHINE=odin

source ../sourceme.sh

module load intel gcc mkl
module load perflib/2.2
export PERFLIB_OUTPUT_FORMAT=xml
#export I_MPI_PIN_PROCESSOR_LIST=allcores
export KMP_AFFINITY=verbose,granularity=core,compact,1

### set OMP_NUM_THREADS
#export OMP_NUM_THREADS=16

### set MKL to serial mode
#export MKL_SERIAL=yes

### set MPI algorithms for Intel MPI ###
#export I_MPI_ADJUST_REDUCE=4
#export F_UFMTENDIAN=big

### start program
#mpiexec -np $NSLOTS ./gene_odin

#rm *_BIS

./magic_checks.pl --all --clean


