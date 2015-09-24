#!/bin/bash

#
# This script is used to remove the pre-processore flags from magic and copy them in
# the sphinx directory
#

dir=$MAGIC_HOME/doc/sphinx/.fortranCleanSrc

if [ ! -d "$dir" ]; then
  mkdir $dir

  for filepath in $MAGIC_HOME/src/*.f90; do
    filename=`basename $filepath`
    targetFile="$dir/$filename"
    fpp -free -P -DWITH_PRECOND_BJ -DWITH_PRECOND_S -DWITH_PRECOND_Z -DWITH_PRECOND_S0 -DWITH_PRECOND_Z10 -DJW -DFFTLIB=JW -Ddble -DDEFAULT_PRECISION=dble $filepath $targetFile
  done

  rm $dir/fft_mkl.f90

fi

