#!/bin/sh

if [[ -n $MAGIC_HOME ]]
then
  cd $MAGIC_HOME/bin
else
  echo "error, please set MAGIC_HOME first (e.g. by running 'source sourceme.sh')"
  exit 1
fi

if command -v hg; then
  hg clone https://bitbucket.org/nschaeff/shtns
else
  wget https://bitbucket.org/nschaeff/shtns/downloads/shtns-3.0-r618.tar.gz
  tar -xvf shtns-3.0-r618.tar.gz
  rm shtns-3.0-r618.tar.gz
  mv shtns-3.0-r618 shtns
fi

opts="--enable-magic-layout"
if [[ -n $MKLROOT ]]
then
   echo "MKL found, installing with MKL"
   opts="$opts -enable-mkl"
else
   echo "MKL not found, will try to install with FFTW"
fi

cd shtns
## compile without OpenMP
./configure $opts
make -j
## compile with OpenMP
./configure --enable-openmp $opts
make -j

