#!/bin/sh

if test ! -d $HOME/local; then
    mkdir $HOME/local
fi

wget 'https://bitbucket.org/bputigny/shtns-magic/downloads/shtns-magic.tar.bz2'

tar xvjf shtns-magic.tar.bz2 && \
rm -f shtns-magic.tar.bz2 && \
cd shtns-magic && \
./configure --prefix=$HOME/local --enable-openmp && make install -j4

