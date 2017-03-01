# Requirements

Most of the output files produced during a run of MagIC
can be treated with the python post-processing classes and functions present in this
directory. These classes depend on several python libraries
that can be usually found in most of the Linux distributions.

## Hard dependencies

* [python](https://www.python.org) 2.7/3.3 or higher.

* [matplotlib](http://matplotlib.org) 1.0 or higher.

* [scipy](http://www.scipy.org/) 0.10 or higher.

## Optional dependencies

* Although entirely optional, the installation of [ipython](http://ipython.org/)  makes the interactive use of the post-processing python functions much more pleasant. Installing it is therefore recommanded for a smoother interactive usage of the python functions.

* The installation of the [basemap toolkit](http://matplotlib.org/basemap/) is optional. If installed, additional projections will be available for visualisations ([Aitoff](https://en.wikipedia.org/wiki/Aitoff_projection), [orthographic](https://en.wikipedia.org/wiki/Orthographic_projection), [Mollweide](https://en.wikipedia.org/wiki/Mollweide_projection), etc.). Otherwise, the display of surfaces will be restricted to the [Hammer projection](https://en.wikipedia.org/wiki/Hammer_projection).


# Setting up the python environment


A file named **magic.cfg** located in `$MAGIC_HOME/python/magic/magic.cfg` should have
been created when you used the `source path/sourceme.sh` command for the first time on
your machine. At that stage, it tried to **automatically fill the best options** that
correspond to your setup. Although tested on several various machine configurations, the
auto-configuration script might however fail on your setup. The paragraph below details
the possible options that you may want to adjust in the `magic.cfg` file.


In case, the file `magic.cfg` doesn't exist in the directory `$MAGIC_HOME/python/magic`,
you can easily copy it from the default configuration `magic.cfg.default` and then adjust
the options manually:

```sh
$ cp $MAGIC_HOME/python/magic/magic.cfg.default $MAGIC_HOME/python/magic/magic.cfg
```


### 1) Matplotlib setup

You can modify the default matplotlib rendering backend (among the possible options: TkAgg GTKAgg, Qt4Agg). Default is

```python
backend = TkAgg
```
If LaTeX is installed on your machine, you can also use the LaTeX fonts for matplotlib labels. On clusters, however, LaTeX is not installed most of the time, so the default is

```python
labTex = False
```

Finally, you can change the default colormap and number of contours that will be used by default when contour plots are requestd.

```python
defaultCm = seismic
defaultLevels = 65
```

### 2) Fortran libs setup

If you want to enable all the features of the python functions (reading the G files with fortran, support for VTK files, ...), the fortran libraries in `$MAGIC_HOME/python/magic/fortranLib`  need to be built using `f2py`.
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
fcompiler = intelem
````

A similar procedure can be applied to the C compilers. On most of the machines, the default (`unix` that corresponds to gcc) is enough, but in case of problem, it can also be possibly changed to `intelem` using the option

```python
ccompiler = intelem
````

Once this is done you should be able to use the python class:

```python
from magic import *
````

If `buildLib=True`, `f2py` will then try to build the libraries.
