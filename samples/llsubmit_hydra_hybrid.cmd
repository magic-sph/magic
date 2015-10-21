# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error   = $(job_name).err.$(jobid)
# @ output  = $(job_name).out.$(jobid)
# @ job_type = parallel
# @ job_name = MagIC
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 2
# @ resources = ConsumableCpus(10) ConsumableMemory(28000)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 00:15:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue

#
# run the program
#

export MACHINE=hydra

export FC=mpiifort
export CC=mpiicc

source ../sourceme.sh

module load intel gcc mkl perflib/3.0

#export PERFLIB_OUTPUT_FORMAT=xml
#export I_MPI_PIN_PROCESSOR_LIST=allcores
export KMP_AFFINITY=verbose,granularity=core,compact,1

### set OMP_NUM_THREADS
export OMP_NUM_THREADS=10

### set MKL to serial mode
#export MKL_SERIAL=yes

### set MPI algorithms for Intel MPI ###
#export I_MPI_ADJUST_REDUCE=4
#export F_UFMTENDIAN=big

### start program

./magic_checks.pl --all --clean --hybrid --use-cmake --use-mkl


