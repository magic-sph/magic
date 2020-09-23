#!/bin/bash

#
# This script is used to remove the pre-processore flags from magic and copy them in
# the sphinx directory
#

dir=$MAGIC_HOME/doc/sphinx/.fortranCleanSrc

if [[ $OSTYPE == *"darwin"* ]]; then
    cppCommand='cpp -P -traditional-cpp -DWITH_MPI -DWITH_PRECOND_BJ -DWITH_PRECOND_S -DWITH_PRECOND_Z -DWITH_PRECOND_S0 -DWITH_PRECOND_Z10 -Ddble -DDEFAULT_PRECISION=dble'
else
    cppCommand='cpp -free -P -traditional-cpp -DWITH_MPI -DWITH_PRECOND_BJ -DWITH_PRECOND_S -DWITH_PRECOND_Z -DWITH_PRECOND_S0 -DWITH_PRECOND_Z10 -Ddble -DDEFAULT_PRECISION=dble'
fi

if [ ! -d "$dir" ]; then
  mkdir $dir

  for filepath in $MAGIC_HOME/src/*.f90; do
    filename=`basename $filepath`
    targetFile="$dir/$filename"

    $cppCommand $filepath $targetFile
  done

  rm $dir/shtns.f90
  rm $dir/*_fftw.f90
  rm $dir/*_lapack.f90

fi

