#!/bin/sh
#----------SHTns compilation et & installation for CPU Genoa on ADASTRA-------------
module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-genoa libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.24
module load cpe/23.02 cray-fftw/3.3.10.3 cray-libsci/23.02.1.1
export LDFLAGS=""
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.3/x86_genoa/lib/"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/fftw/3.3.10.3/x86_genoa/lib/
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/schaeffn/shtns.git 
cd shtns
git checkout hip-rtc
#env CC=cc ./configure --enable-openmp --prefix=/lus/home/NAT/cpa2204/SHARED/INSTALL/shtns #(ex.: --prefix=PATH_TO_INSTALL_DIR)
./configure --enable-openmp --prefix=$HOME/local #(ex.: --prefix=PATH_TO_INSTALL_DIR)
make -j
make install #Install libshtns_omp.a in $PATH_TO_INSTALL_DIR/lib
make time_SHT
#Run time_SHT sample on CPU
./time_SHT 255 -nlorder=2 -quickinit -vector -batch
./time_SHT  127 -quickinit -vector -nlorder=2
./time_SHT 1022 -quickinit -vector -nlorder=2 -iter=1
./time_SHT 1022 -quickinit -vector -nlorder=2
./time_SHT 1022 -quickinit -vector -nlorder=2 -batch
./time_SHT 2300 -quickinit -vector -nlorder=2 -iter=8 -batch
./time_SHT 4095 -quickinit -vector -iter=1 -mres=5
./time_SHT 8191 -quickinit -vector -iter=1 -mmax=0
#----- SHTns compilation et & installation for GPU (with CPU Trento) on ADASTRA----------
make clean
module purge
module load craype/2.7.19 PrgEnv-cray/8.3.3
module load craype-x86-trento libfabric/1.15.2.0 craype-network-ofi cray-dsmml/0.2.2 cray-mpich/8.1.24
module load craype-accel-amd-gfx90a rocm/5.2.3
module load cray-fftw/3.3.10.3 cray-libsci/23.02.1.1
export LDFLAGS=""
export LDFLAGS+="-L/opt/cray/pe/fftw/3.3.10.3/x86_trento/lib/"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/fftw/3.3.10.3/x86_trento/lib/
env CC=cc ./configure --enable-hip=mi250 --prefix=$HOME/local #(ex.: --prefix=PATH_TO_INSTALL_DIR)
srun --account=cpa --constraint=MI250 --ntasks=1 --gres=gpu:1 --exclusive --time=00:10:00 make
make install #Install libshtns_cuda.a in $PATH_TO_INSTALL_DIR/lib
#Run time_SHT sample on GPU
srun --account=cpa --constraint=MI250 --ntasks=1 --gres=gpu:1 --exclusive --time=00:10:00 /bin/bash << END_CMD
make time_SHT
./time_SHT 255 -nlorder=2 -quickinit -vector -batch
./time_SHT  127 -quickinit -vector -nlorder=2
./time_SHT 1022 -quickinit -vector -nlorder=2 -iter=1
./time_SHT 1022 -quickinit -vector -nlorder=2
./time_SHT 1022 -quickinit -vector -nlorder=2 -batch
./time_SHT 2300 -quickinit -vector -nlorder=2 -iter=8 -batch
./time_SHT 4095 -quickinit -vector -iter=1 -mres=5
./time_SHT 8191 -quickinit -vector -iter=1 -mmax=0
END_CMD
