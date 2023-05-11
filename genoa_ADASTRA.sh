#!/bin/sh

#---------- MagIC compilation & execution for CPU Genoa on ADASTRA : SHTns must be installed in $HOME/local -------------

module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-genoa libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.21
module load cpe/22.11 cray-fftw/3.3.10.3 cray-libsci/22.11.1.2
module list

export LDFLAGS+=""
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.3/x86_genoa/lib/"
export FC=ftn

if [ -e "build" ];then rm -rf "build" ; fi

mkdir build && cd build
cmake .. -DUSE_SHTNS=yes -DUSE_GPU=no 
make -j16

cd ../namelists_hpc/run_middle_ADASTRA/CPU/ && source clear.sh
cd ../../run_big_ADASTRA/CPU/ && source clear.sh
cd ../../../samples/boussBenchSat/ && source clear.sh

cd ../../build
cp magic.exe ../namelists_hpc/run_middle_ADASTRA/CPU
cp magic.exe ../namelists_hpc/run_big_ADASTRA/CPU
cp magic.exe ../samples/boussBenchSat

cd ../namelists_hpc/run_middle_ADASTRA/CPU
sbatch submit.sh

cd ../../run_big_ADASTRA/CPU
sbatch submit.sh

cd ../../../samples/boussBenchSat/
sbatch genoa_ADATRA_submit.sh

squeue | grep MagIC


