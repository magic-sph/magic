#!/bin/bash
if [ $# -eq 1 ]; then
    hybrid=$1
else
    hybrid="no"
fi

export F_UFMTENDIAN=big
if [ "$hybrid" == "hybrid" ]; then
    nmpi=2
    nomp=4
    echo "Running hybrid mode with $nomp OpenMP threads per $nmpi MPI processes."
    export I_MPI_PIN_DOMAIN=socket
    export KMP_STACKSIZE=500m
    export KMP_AFFINITY=noverbose,granularity=core,compact
else
    nomp=1
    nmpi=8
    echo "Running pure MPI code with $nmpi MPI processes."
    export I_MPI_PIN_PROCESSOR_LIST=allcores
fi
export OMP_NUM_THREADS=$nomp

# switch off the threading in the input file
cp input.nml input_tmp.nml

if grep -Fxq "USE_HDF5=yes" $MAGIC_HOME/src/Makefile; then
    sed -i 's/start_file  ="rst_end.start"/start_file  ="h5_rst_end.start"/' input_tmp.nml
fi

# First run
mpiexec -n $nmpi ./magic.exe  inputStart.nml

# Restart
mpiexec -n $nmpi ./magic.exe  input_tmp.nml

# Clean
rm *.start
rm input_tmp.nml
