#!/bin/sh

if test ! -d $HOME/local; then
    mkdir $HOME/local
fi

if command -v hg; then
  hg clone https://bitbucket.org/nschaeff/shtns
  cd shtns
else
  wget https://bitbucket.org/nschaeff/shtns/downloads/shtns-3.0-r618.tar.gz
  tar -xvf shtns-3.0-r618.tar.gz
  rm shtns-3.0-r618.tar.gz
  cd shtns-3.0-r618
fi

./configure --prefix=$HOME/local --enable-openmp --enable-magic-layout
make install -j
