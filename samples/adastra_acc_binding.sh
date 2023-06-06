#!/bin/bash

set -eu

function Adastra_MI250_8TasksWith8ThreadsAnd1GPU() {
    AFFINITY_NUMACTL=('48-55' '56-63' '16-23' '24-31' '0-7' '8-15' '32-39' '40-47')
    AFFINITY_GPU=('0' '1' '2' '3' '4' '5' '6' '7')
    export MPICH_OFI_NIC_POLICY=NUMA
}

function Adastra_GENOA_24TasksWith8Threads() {
    AFFINITY_NUMACTL=('0-7' '8-15' '16-23' '24-31' '32-39' '40-47' '48-55' '56-63' '64-71' '72-79' '80-87' '88-95' '96-103' '104-111' '112-119' '120-127' '128-135' '136-143' '144-151' '152-159' '160-167' '168-175' '176-183' '184-191')
}

function GPUAffinityAMD() {
    export ROCR_VISIBLE_DEVICES="${AFFINITY_GPU[${1}]}"
}

function CPUAffinityNumactl() {
    START_COMMAND="exec numactl --localalloc --physcpubind=${AFFINITY_NUMACTL[${1}]} --"
}

Adastra_MI250_8TasksWith8ThreadsAnd1GPU
# Adastra_GENOA_24TasksWith8Threads

LOCAL_RANK=$((${SLURM_LOCALID} % ${#AFFINITY_NUMACTL[@]}))

GPUAffinityAMD ${LOCAL_RANK}
CPUAffinityNumactl ${LOCAL_RANK}

echo "Starting local rank: ${LOCAL_RANK} with: '${START_COMMAND}'"

${START_COMMAND} "${@}"
