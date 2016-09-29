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
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 01:00:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ queue

module load intel gcc mkl cmake
module load mpi.ibm
module load python27/python/2.7
module load python27/scipy/2015.10

export MACHINE=hydra
export FC=mpiifort
export CC=mpiicc
export KMP_AFFINITY=verbose,granularity=core,compact,1

source ../sourceme.sh

./magic_wizard.py --use-mpi --use-mkl --nranks 16 --mpicmd mpiexec
