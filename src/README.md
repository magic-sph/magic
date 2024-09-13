# Makefile options

### 1. Select compiler

Set a suitable compiler in the first line of the Makefile: `COMPILER = <compiler_phrase>`. The options are `intel`, `gnu` or `portland` - depending on your available compilers.

**Table** : List of default compilers

| Compiler Option |    Normal     |     With MPI     |
|:---------------:|:-------------:|:----------------:|
|      intel      |   ifort, icc  | mpiifort, mpiicc |
|       gnu       | gfortran, gcc |   mpif90, mpicc  |
|     portland    |  pgf95, pgcc  |   mpif90, mpicc  |

### 2. Select compiling options

* `PRODRUN` Set it to `yes` for production run, `no` for debugging.
* `PRECISION` Set it to 'dble' for double-precision calculations or to 'sngl' for single-precision calculations
* `OUT_PREC` Set it to 'dble' for double-precision in binary outputs or to 'sngl' for single precision
* `DEBUG`   Set to `yes` to run in debugging mode. *While running in debugging mode, set* `PRODRUN` *to* `no`. The debug mode with intel compilers uses `marmot90`.
* `USE_MPI` Set to `yes` to use MPI
* `USE_OMP`  Set to `yes` to use openmp (cannot work without MPI)
* `USE_PRECOND` Set to `yes` to perform some pre-conditioning  of the matrices
* `USE_FFTLIB` This option lets you select the library you want to use for Fast Fourier Transforms. This can be set to `JW`, `FFTW` or `MKL`. `JW` refers to the inbuilt library by **J**ohannes **W**icht, while `MKL` refers to the [Intel **M**ath **K**ernel **L**ibrary](https://software.intel.com/en-us/intel-mkl). Use `JW` if you don't have Intel MKL installed.
* `USE_DCTLIB` This option lets you select the library you want to use for Discrete Cosine Transforms. This can be set to `JW`, `FFTW` or `MKL`.
* `USE_LAPACKLIB` This option allows you to select the library you want to use for LU factorisation. This can be set to `JW`, `MKL`, `LIBFLAME` or `LAPACK`.
* `USE_SHTNS` Set to `yes` to use SHTns library for spherical harmonics transforms. The helper script `install-shtns.sh` is available in the `bin` directory to help installing SHTns.

### 3. MPI_INCPATH

Make sure you set the path for your mpi header file `mpif.h` in `MPI_INCPATH`. The path depends on the computer. For PCs, this is commonly `/usr/include` or `/usr/include/mpi`. Use [Open MPI](http://www.open-mpi.de/) for running MagIC on a PC. For computing clusters, please look through the documentation of the respective cluster for their MPI implementation.

### 4. Other compilers

If your available compilers are different from the options provided in the Makefile, then change them suitably using the options `COMP_FC` and `COMP_CC` for serial fortran and C compilers and `COMP_MPFC` and `COMP_MPCC` for compilers with mpi implementation.
