# Makefile options

### 1. Select compiler

Set a suitable compiler in the first line of the Makefile: `COMPILER = <compiler_phrase>`. The options are `intel`, `gnu` or `amd` - depending on your available compilers.

**Table** : List of default compilers
 
| Compiler Option |    Normal    |     With MPI    |
|:---------------:|:------------:|:---------------:|
|      intel      |   ifort, icc  | mpiifort, mpiicc |
|       gnu       | gfortran, gcc |   mpif90, mpicc  |
|       amd       |    openf95   |                 |

### 2. Select compiling options

* `PRODRUN` Set it to `yes` for production run, `no` for debugging.
* `OPENMP`  Set to `yes` to use openmp
* `DEBUG`   Set to `yes` to run in debugging mode. *While running in debugging mode, set* `PRODRUN` *to* `no`. The debug mode with intel compilers uses `marmot90`. 
* `USE_MPI` Set to `yes` to use MPI
* `USE_FFTLIB` This option lets you select the library you want to use for Fast Fourier Transforms. This can be set to 'JW' or 'MKL'. 'JW' refers to the inbuilt library by **J**ohannes **W**icht, while 'MKL' refers to the [Intel **M**ath **K**ernel **L**ibrary](https://software.intel.com/en-us/intel-mkl). Use 'JW' if you don't have Intel MKL installed.
* `USE_MKL` Set to `yes` if you have Intel MKL installed and want to use it for matrix operations.
* `USE_HDF5` Set to `yes` if you want the restart file to be written in HDF5 format

### 3. Architecture (Intel compilers only)

If you're using intel compilers and if your computer is capable of following specific intel instruction sets (sse3 or AVX). Then set `FFLAG_ARCH_OPT = -xsse3` or `FFLAG_ARCH_OPT = -xAVX` under intel compiler options.

### 4. MPI_INCPATH

Make sure you set the path for your mpi header file `mpif.h` in `MPI_INCPATH`. The path depends on the computer. For PCs, this is commonly `/usr/include` or `/usr/include/mpi`. Use [Open MPI](http://www.open-mpi.de/) for running MagIC on a PC. For computing clusters, please look through the documentation of the respective cluster for their MPI implementation.

### 5. Other compilers

If your available compilers are different from the options provided in the Makefile, then change them suitably using the options `COMP_FC` and `COMP_CC` for serial fortran and C compilers and `COMP_MPFC` and `COMP_MPCC` for compilers with mpi implementation.
