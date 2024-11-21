FROM ubuntu:22.04

# This is a demo container for testing magic!
# docker build -t magic .

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    fftw3-dev fftw3 \
    gfortran \
    git \
    libmpich-dev \
    python3 \
    python3-matplotlib \
    python3-f2py \
    wget && \
    ln -s /usr/bin/python3 /usr/bin/python
        
WORKDIR /code
COPY . /code

# install SHTns with fft, then build magic and run (and build) samples
RUN /bin/bash -c ". sourceme.sh && \
    cd ./bin && \
    CC=gcc ./install-shtns.sh && \
    cd ../ && mkdir -p ./build && cd ./build && \
    FC=mpif90 cmake .. -DUSE_SHTNS=yes && make && \
    cd ../samples && \
   ./magic_wizard.py --use-mpi --nranks 4 --mpicmd mpiexec"

ENTRYPOINT ["/code/entrypoint.sh"]
