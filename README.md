# Quickly start using MagIC

### 1) In order to check out the code, use the command

```sh
git clone https://github.com/tgastin/magic.git
```
or via SSH (it requires a public key):

```sh
git clone ssh://git@github.com/tgastine/magic.git
```

### 2) Go to the root directory and source the environnement variables (useful for python and auto-tests)

```sh
cd MAGIC_mpi
```

If you are using sh, bash or zsh as default shell (echo $SHELL), just use the command

```sh
source sourceme.sh
```

If you are using csh or tcsh, then use the following command

```sh
source sourceme.csh
```

### 3) Set up your compiler and compile the code

```sh
cd $MAGIC_HOME/src
```

Edit the Makefile and specify your compiler (intel or gnu) and additional 
compiler options (production run or not, debug mode, MKL library, HDF5, ...)

```sh
make -j
```
The executable magic.exe has been produced!

### 4) Go to the samples directory and check that everything is fine

```sh
-> cd $MAGIC_HOME/samples
-> ./magic_checks.pl --all --clean --no-recompile
```

If everything is correctly set, all auto-tests should pass!

### 5) You're ready for a production run

```sh
-> cd $SCRATCHDIR/run
-> cp $MAGIC_HOME/src/magic.exe .
-> cp $MAGIC_HOME/samples/hydro_bench_anel/input.nml .
```
    
Then change the input namelist to the setup you want and run the code:

```sh
-> mpiexec -n 4 ./magic.exe input.nml
```

### 6) Data visualisation and postprocessing

    a) Set-up your PYTHON environnement ([ipython](http://ipython.org/), [scipy](http://www.scipy.org/) and [matplotlib](http://matplotlib.org/) are needed)

    b) Modify `magic.cfg` according to your machine

```sh
-> vi $MAGIC_HOME/python/magic/magic.cfg
```

    c) You can now import the python class:

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

a) Before commiting your modifications ALWAYS make sure that the auto-tests
pass correctly.

b) Try to follow the same coding style rules as in the rest of the code:

1. **Never** use TABS but alwyas SPACES instead
2. Use 3 spaces for indentation
3. Never use capital letters for variable declaration
4. Never use 'dimension(len)' for declaring array but rather real(cp) :: data(len)
5. Always use the default precisions when introducing new variables (cp)


More on that topic [here](http://www.fortran90.org/src/best-practices.html)

### 8) Make sure you cite the following papers if you intend to publish scientific results using Magic:

* Boussinesq equations: Wicht (2002, PEPI, 132, 281-302)
* Anelastic equations: Gastine & Wicht (2012, Icarus, 219, 428-442)

Magic has been tested and validated against several international dynamo benchmarks:
* Christensen et al. (2001, PEPI, 128, 25-34)
* Jones et al. (2011, Icarus, 216, 120-135)
