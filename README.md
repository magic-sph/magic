![alt tag](https://raw.github.com/magic-sph/magic/master/doc/sphinx/.themes/magic/static/logo.png)

[![Build Status](https://travis-ci.org/magic-sph/magic.svg?branch=master)](https://travis-ci.org/magic-sph/magic)
[![DOI](https://zenodo.org/badge/22163/magic-sph/magic.svg)](https://zenodo.org/badge/latestdoi/22163/magic-sph/magic)
[![GPLv3](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl.html)

# Foreword

* **MagIC** is a numerical code that can simulate fluid dynamics in a spherical shell. MagIC solves for the Navier-Stokes equation including Coriolis force, optionally coupled with an induction equation for Magneto-Hydro Dynamics (MHD), a temperature (or entropy) equation and an equation for chemical composition under both the anelastic and the Boussinesq approximations.  

* **MagIC** uses either Chebyshev polynomials or finite differences in the radial direction and spherical harmonic decomposition in the azimuthal and latitudinal directions. The time-stepping scheme relies on a semi-implicit [Crank-Nicolson]( https://en.wikipedia.org/wiki/Crank–Nicolson_method) for the linear terms of the MHD equations and a [Adams-Bashforth](<https://en.wikipedia.org/wiki/Linear_multistep_method) scheme for the non-linear terms and the Coriolis force.  

* **MagIC** is written in Fortran and designed to be used on supercomputing clusters.  It thus relies on a hybrid parallelisation scheme using both [OpenMP](http://openmp.org/wp/) and [MPI](http://www.open-mpi.org/). Postprocessing functions written in python (requiring [matplotlib](http://matplotlib.org/) and [scipy](http://www.scipy.org/) are also provided to allow a useful data analysis.  

* **MagIC** is a free software. It can be used, modified and redistributed under the terms of the [GNU GPL v3 licence](http://www.gnu.org/licenses/gpl-3.0.en.html).


# Quickly start using MagIC

### 1) In order to check out the code, use the command

```sh
$ git clone https://github.com/magic-sph/magic.git
```
or via SSH (it requires a public key):

```sh
$ git clone ssh://git@github.com/magic-sph/magic.git
```

### 2) Go to the root directory and source the environment variables (useful for python and auto-tests)

```sh
$ cd magic
```

If you are using sh, bash or zsh as default shell (`echo $SHELL`), just use the command

```sh
$ source sourceme.sh
```

If you are using csh or tcsh, then use the following command

```sh
$ source sourceme.csh
```

### 3) Set up your compiler and compile the code


#### a) Using CMake (recommended)

Create a directory where the sources will be built

```sh
$ mkdir $MAGIC_HOME/build
$ cd $MAGIC_HOME/build
```
Set up your compilers

```sh
$ export FC=mpiifort
$ export CC=mpiicc
```

Compile and produce the executable (options can be passed to cmake using `-DOPTION=value`)

```sh
$ cmake ..
$ make -j
````
The executable `magic.exe` has been produced!

#### b) Using make (backup solution)

Go to the source directory

```sh
$ cd $MAGIC_HOME/src
```

Edit the Makefile with your favourite editor and specify your compiler 
(intel, gnu, portland) and additional 
compiler options (production run or not, debug mode, MKL library, ...)

```sh
$ make -j
```
The executable `magic.exe` has been produced!

### 4) Go to the samples directory and check that everything is fine

```sh
$ cd $MAGIC_HOME/samples
$ ./magic_wizard.py --use-mpi --nranks 4 --mpicmd mpiexec
```

If everything is correctly set, all auto-tests should pass!

### 5) You're ready for a production run

```sh
$ cd $SCRATCHDIR/run
$ cp $MAGIC_HOME/build/magic.exe .
$ cp $MAGIC_HOME/samples/hydro_bench_anel/input.nml .
```
    
Then change the input namelist to the setup you want and run the code:

```sh
$ export OMP_NUM_THREADS=2
$ export KMP_AFFINITY=verbose,granularity=core,compact,1
$ mpiexec -n 4 ./magic.exe input.nml
```

### 6) Data visualisation and postprocessing

a) Set-up your PYTHON environment ([ipython](http://ipython.org/), [scipy](http://www.scipy.org/) and [matplotlib](http://matplotlib.org/) are needed)

b) Modify `magic.cfg` according to your machine in case the auto-configuration didn't work

```sh
$ vi $MAGIC_HOME/python/magic/magic.cfg
```

c) You can now import the python classes:

```python
python> from magic import *
```

and use them to read time series, graphic files, movies, ...

```python
python> ts = MagicTs(field='e_kin', all=True)
python> s = Surf()
python> s.equat(field='vr')
python> ...
```

### 7) Modify the code and submit your modifications

a) Before commiting your modifications **always** make sure that the auto-tests
pass correctly.

b) Try to follow the same coding style rules as in the rest of the code:

1. **Never** use TABS but always SPACES instead
2. Use 3 spaces for indentation
3. Never use capital letters for variable declaration or Fortran keywords
4. Never use `dimension(len)` for declaring array but rather `real(cp) :: data(len)`
5. Always use the default precisions when introducing new variables `(cp)`

More on that topic [here](http://www.fortran90.org/src/best-practices.html)

### 8) Make sure you cite the following papers if you intend to publish scientific results using MagIC:

* Boussinesq equations: [Wicht (2002, PEPI, 132, 281-302)](http://dx.doi.org/10.1016/S0031-9201(02)00078-X)
* Anelastic equations: [Gastine & Wicht (2012, Icarus, 219, 428-442)](http://dx.doi.org/10.1016/j.icarus.2012.03.018)
* In case you use the [SHTns](https://bitbucket.org/bputigny/shtns-magic) library for the spherical harmonics transforms (MagIC 5.3 or later), please also cite: [Schaeffer (2013, GGG, 14, 751-758)](http://dx.doi.org/10.1002/ggge.20071)

MagIC has been tested and validated against several international dynamo benchmarks:
* [Christensen et al. (2001, PEPI, 128, 25-34)](http://dx.doi.org/10.1016/S0031-9201(01)00275-8)
* [Breuer et al. (2010, GJI, 183, 150-162)](http://dx.doi.org/10.1111/j.1365-246X.2010.04722.x)
* [Jones et al. (2011, Icarus, 216, 120-135)](http://dx.doi.org/10.1016/j.icarus.2011.08.014)
