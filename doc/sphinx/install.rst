Running the code
################

Download the code
=================


Setting up the environment
==========================


Setting up compiler options and compiling
=========================================

Compiler options
----------------

Supported compilers
-------------------

Preparing a production run
==========================


Running the code without OpenMP::
  
  mpiexev -n 4 ./magic.exe input.nml

Running the code with OpenMP::
  
  export OMP_NUM_THREAD=4
