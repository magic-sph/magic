# @ shell=/bin/bash
#
# Sample script for LoadLeveler
#
# @ error   = $(job_name).e_$(jobid)
# @ output  = $(job_name).o_$(jobid)
# @ job_type = parallel
# @ job_name = MagIC
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 2
## @ first_node_tasks = 0
# @ task_affinity=core(10)
# @ resources = ConsumableCpus(10) ConsumableMemory(28000)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 24:00:00
# @ notification = complete
# @ notify_user = $(user)@rzg.mpg.de
# @ node_topology = island
# @ island_count = 1
# @ queue

#
# run the program
#

export MACHINE=hydra

export KMP_AFFINITY=noverbose,granularity=core,compact
export KMP_STACKSIZE=1g

export MP_CSS_INTERRUPT=yes
export MP_WAIT_MODE=nopoll
#export MP_POLLING_INTERVAL=100
#export MP_PRINTENV=yes
export MP_BULK_MIN_MSG_SIZE=8k
export MP_EAGER_LIMIT=8k

### set OMP_NUM_THREADS
n_omp=10
export OMP_NUM_THREADS=${n_omp}

### set MKL to serial mode
export MKL_SERIAL=yes

### set MPI algorithms for Intel MPI ###
#export I_MPI_ADJUST_REDUCE=4
export F_UFMTENDIAN=big

# EXECUTABLE and input_file must be set
EXECUTABLE=
input_file=

poe $EXECUTABLE ${input_file} 
 
