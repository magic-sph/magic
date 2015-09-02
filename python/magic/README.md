# Setting up the python environment

To properly set up the python post-processing routines, you have to edit the file `$MAGIC_HOME/python/magic/magic.cfg`

### 1) Matplotlib setup

You can modify the default matploblib rendering backend (among the possible options: TKAgg GTKAgg, Qt4Agg). Default is

```python
backend = TkAgg
```
If LaTeX is installed on your machine, you can also use the LaTeX fonts for matplotlib labels. On clusters, however, LaTeX is however not installed most of the time, so the default is

```python
labTex = False
```
### 2) Fortran libs setup

If you want to enable all the features of the python suboutines (reading the G files with fortran, support for VTK files, ...), the fortran libraries in `$MAGIC_HOME/python/magic/fortranLib`  need to be built using `f2py`.
This can be set with the following variable

```python
buildLib = True
```
Depending on the machine and the version of python (2 or 3), the name of the executable however changes (`f2py`, `f2py2`, `f2py3`).
You can set it up with the following option:

```python
f2pyexec = f2py2
```
Finally, you have to specify your fortran compiler. The available compilers on your machine can be obtained by using

```sh
f2py2 -c --help-fcompiler
```
The option is then set via

```python
compiler = intelem
````

Once this is done you should be able to use the python class:

```python
from magic import *
````

If `buildLib=True`, `f2py` will then try to build the libraries.
