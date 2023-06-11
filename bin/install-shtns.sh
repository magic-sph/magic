#!/bin/bash

ver="3.5.2"

if test ! -d $HOME/local; then
    mkdir $HOME/local
fi

wget https://bitbucket.org/nschaeff/shtns/downloads/shtns-$ver.tar.gz
tar -xvf shtns-$ver.tar.gz
rm shtns-$ver.tar.gz

if [ -d "shtns-$ver" ]
then
    mv shtns-$ver shtns
fi

cd shtns

shtns_version=`grep configure -e "PACKAGE_VERSION=" | sed -e "s/.*=//" -e "s/'//g"`
shtns_major="$(cut -d'.' -f1 <<< "$shtns_version")"
shtns_minor="$(cut -d'.' -f2 <<< "$shtns_version")"
shtns_minor_sub="$(cut -d'.' -f3 <<< "$shtns_version")"

if [ ${shtns_major} -lt 3 ] # Lower than 3 always "sed"
then
   sed -i '/magic_layout/{N;N;N;s/enable_many_core=yes/enable_many_core=no/;}' configure
else
   if [ ${shtns_major} -eq 3 ]
   then
     if [ ${shtns_minor} -lt 3 ] # Lower than 3.3
     then
        if [ ${shtns_minor} -eq 2 ] # 3.2....
        then
           if [ ${shtns_minor_sub} -lt 2 ]
           then
              sed -i '/magic_layout/{N;N;N;N;s/enable_many_core=yes/enable_many_core=no/;}' configure
           fi
         else # 3.1 or so
           sed -i '/magic_layout/{N;N;N;s/enable_many_core=yes/enable_many_core=no/;}' configure 
         fi
     fi  
   fi
fi

opts="--enable-magic-layout --prefix=$HOME/local --enable-ishioka"
if [[ -n $MKLROOT ]]
then
   echo "MKL found, installing with MKL"
   opts="$opts --enable-mkl"
else
   echo "MKL not found, will try to install with FFTW"
fi

# Compile without OpenMP
./configure $opts

if [ `echo $CC` ]
then
    sed -i "s/shtcc=gcc/shtcc=${CC}/" Makefile
fi

make -j
make install -j

# Compile with OpenMP
./configure --enable-openmp $opts

if [ `echo $CC` ]
then
    sed -i "s/shtcc=gcc/shtcc=${CC}/" Makefile
fi

make -j
make install -j

