.. _secQuickStart:

Get MagIC and run it
####################

Download the code
=================

You can download a snapshot of the code from the `Git <https://git-scm.com/>`_ repository using

.. code-block:: bash

   $ git clone https://github.com/magic-sph/magic.git


In case you already have an account on `github.com <https://github.com/>`_ and uploaded a public SSH key on it, you could then rather use SSH:

.. code-block:: bash

   $ git clone ssh://git@github.com/magic-sph/magic.git

Setting up the environment variables
====================================

Although not mandatory, it is strongly recommended to correctly source the
environment variables of the MagIC code. It will ensure a smoother usage of the
post-processing :ref:`python classes <secPythonPostProc>` and allow to run the
:ref:`auto-test suite <secAutoTest>`.  To do that, just go to the root directory of
the MagIC code (``magic``) and source ``sourceme`` file that corresponds to your
``$SHELL`` environment variable.

In case you use `bash <http://tiswww.case.edu/php/chet/bash/bashtop.html>`_,
`ksh <http://www.kornshell.com/>`_ or `zsh <http://www.zsh.org/>`_, just use:

.. code-block:: bash
 
   $ source sourceme.sh

In case you use `csh <http://www.tcsh.org/Home>`_ or `tcsh <http://www.tcsh.org/Home>`_,
rather use

.. code-block:: csh
 
   $ source sourceme.csh

You can make sure that the environment variables have been correctly sourced by typing:

.. code-block:: bash

   $ echo $MAGIC_HOME
   $ echo $PYTHONPATH

If you don't want to ``source sourceme.[c]sh`` on each session, you can add the following
into your ``.bash_profile`` (or ``.profile`` or ``.zprofile`` or ``.cshrc``):

.. code-block:: bash

   $ source whereverYouCheckedOut/magic/sourceme.sh

To get started, you then need to compile the code

Setting up compiler options and compiling
=========================================

The **recommended way of compiling MagIC** is to use the build system `CMake
<https://cmake.org/>`_ , if available on your platform. Otherwise, a backup
solution is provided via the manual edition of a ``Makefile``.

Generic compiling options
-------------------------

For both build systems (cmake or make), several build options can be toggled using
the following available options:

* ``PRECISION`` Set it to 'dble' for double-precision calculations or to 'sngl' for single-precision calculations
* ``OUT_PREC`` Set it to 'dble' for double-precision in binary outputs or to 'sngl' for single precision
* ``USE_MPI`` Set to ``yes`` to use MPI, set it to ``no`` if you want a serial version of the code .
* ``USE_OMP``  Set it to ``yes`` to use the hybrid version of the code, or to ``no`` for a pure MPI (or serial) version.
* ``USE_PRECOND`` Set to ``yes`` to perform some pre-conditioning of the matrices.
* ``USE_FFTLIB`` This option lets you select the library you want to use for Fast Fourier Transforms. This can be set to 'JW' or 'MKL'. 'JW' refers to the inbuilt library by **J** ohannes **W** icht, while 'MKL' refers to the `Intel Math Kernel Library <https://software.intel.com/en-us/intel-mkl>`_. Use 'JW' if you don't have Intel MKL installed.
* ``USE_LAPACKLIB`` This option allows you to select the library you want to use for LU factorisation. This can be set to 'JW', 'MKL' or 'LAPACK'. 'JW' refers to the built-in library, while 'MKL' refers to the `Intel Math Kernel Library <https://software.intel.com/en-us/intel-mkl>`_ and 'LAPACK' to the `Lapack library <http://www.netlib.org/lapack>`_
* ``USE_HDF5`` Set to ``yes`` if you want the restart file to be written in the  `HDF5 <http://www.hdfgroup.org/>`_ format
* ``PRODRUN`` Set it to ``yes`` for production run, ``no`` for debugging.
* ``DEBUG``   Set to ``all`` to enable the full debug flags. *While running in debugging mode, set* ``PRODRUN`` *to* ``no``. 

.. warning:: MagIC cannot run with openMP alone, therefore a configuration of the form
          ``USE_MPI=no``, ``USE_OMP=yes`` will be overwritten to force ``USE_OMP=no``

Using ``CMake`` (recommended)
-----------------------------

`CMake <https://cmake.org/>`_  is a powerful tool that can automatically detects
and finds the best appropriate configuration for your platform. To use it, you
just need to create a directory where you want to build the sources. For instance:

.. code-block:: bash

   $ mkdir $MAGIC_HOME/build
   $ cd $MAGIC_HOME/build
   
In a second step, you might want to specify your C and Fortran compilers (in
case you skip this step, `CMake <https://cmake.org/>`_ will look for compilers
for you but it might pick up another compiler as the one you might have wanted).
For instance, in case you want to use the `Intel compilers
<https://software.intel.com/en-us/intel-compilers>`_, you can export the following
environment variables

.. code-block:: bash

   $ export FC=mpiifort
   $ export CC=mpiicc
   
for bash/ksh/zsh users and

.. code-block:: tcsh

   $ setenv FC=mpiifort
   $ setenv CC=mpiicc

for csh/tcsh users. At this stage you should be ready to build the code. If you simply use:

.. code-block:: bash

   $ cmake ..

`CMake <https://cmake.org/>`_ will try to use the best options available on your
machine (for instance it will try to locate and link the `Intel Math Kernel Library
<https://software.intel.com/en-us/intel-mkl>`_). Otherwise you can
pass the aforementioned available options to `CMake <https://cmake.org/>`_ using the
generic form `-DOPTION=value`. For instance, in case you want to make use of
the built-in libraries of MagIC and want to disable OpenMP, simply use

.. code-block:: bash

   $ cmake .. -DUSE_OMP=no -DUSE_FFTLIB=JW -DUSE_LAPACKLIB=JW

Once you're happy with your configuration, just compile the code:

.. code-block:: bash

   $ make -j
   
The executable ``magic.exe`` should have been produced in the local directory.

If you want to recompile the code from scratch do

.. code-block:: bash

   $ make clean

to remove all the files generated by the compiler.

Once the executable is built, you are now ready to run your first production run!

Using ``make`` (backup solution)
--------------------------------

In case `CMake <https://cmake.org/>`_  is not available on your platform, it is
still possible to compile the code directly.  Go to the directory where the
source files of MagIC are contained

.. code-block:: bash

   $ cd $MAGIC_HOME/src
   
**Select compiler**

Edit the file named ``Makefile`` using your favourite editor and set a suitable
compiler for your platform using the variable: ``COMPILER = value``. The possible
options are ``intel``, ``gnu`` or ``portland`` compilers.

*List of default compilers*

  +-----------------+---------------+------------------+ 
  | Compiler Option |    Normal     |     With MPI     |
  +-----------------+---------------+------------------+
  | intel           | ifort, icc    | mpiifort, mpiicc |
  +-----------------+---------------+------------------+
  | gnu             | gfortran, gcc | mpif90, mpicc    |
  +-----------------+---------------+------------------+
  | portland        | pgf95, pgcc   | mpif90, mpicc    |
  +-----------------+---------------+------------------+

.. warning::
   In case you want to use intel but ``mpiifort`` and ``mpiicc`` are not available,
   you may also need to adapt the variables ``COMP_MPFC`` and ``COMP_MPCC``.

**Select compiling options**

You can also modify the different compiling options by editing the values of
the various parameters defined in the first lines of the ``Makefile``.
For instance, in case you want to make use of
the built-in libraries and want to disable OpenMP, just define

.. code-block:: make

   USE_OMP=no
   USE_FFTLIB=JW
   USE_LAPACKLIB=JW

**MPI_INCPATH**

This variable sets the path for your MPI header file ``mpif.h``. This is in
general useless if you already use the MPI wrappers such as ``mpiifort`` or
``mpif90`` to compile the code. It might be however required to define this
path for some compiler configurations: ``MPI_INCPATH`` is usually
``/usr/include`` or ``/usr/include/mpi`` and should be found by the
``Makefile`` automatically thanks to the command ``mpif90 --showme:incdirs``.
In case this doesn't work, you may need to specify this variable manually in
the ``Makefile``. On supercomputing clusters, this variable is in general not
used.

**Other compilers**

If your available compilers are different from the options provided in the
``Makefile``, then just create a new profile for your desired compiler
by changing the options ``COMP_FC`` and
``COMP_CC`` for serial fortran and C compilers and ``COMP_MPFC`` and
``COMP_MPCC`` for the possible MPI wrappers.

Once you've set up your compiling options compile the code using

.. code-block:: bash

   $ make -j

The compiler should then produce an executable named ``magic.exe``.

If you want to recompile the code from scratch do

.. code-block:: bash

   $ make clean

to remove all the files generated by the compiler.

Once the executable is built, you are now ready to run your first production run!

Preparing a production run
==========================

After building the executable, use one of the namelists provided in the
``$MAGIC_HOME/samples`` directory (called ``input.nml``), adapt it to your
physical problem (see :ref:`here <secNamelists>` for an exhaustive
description of the possible options) and run **MagIC** as follows:

* Running a serial version of the code (``USE_MPI=no`` and ``USE_OMP=no``):

  .. code-block:: bash

     $ ./magic.exe input.nml

* Running the code without OpenMP (``USE_MPI=yes`` and ``USE_OMP=no``) with ``<n_mpi>`` 
  MPI ranks:
  
  .. code-block:: bash

     $ mpiexec -n <n_mpi> ./magic.exe input.nml

* Running the hybrid code (``USE_MPI=yes`` and ``USE_OMP=yes``) with ``<n_mpi>`` MPI ranks 
  and ``<n_omp>`` OpenMP threads:
  
  .. code-block:: bash

     $ export OMP_NUM_THREAD = <n_omp>
     $ mpiexec -n <n_mpi> ./magic.exe input.nml

Note that the :ref:`n_r_max <varn_r_max>` must be a multiple of ``<n_mpi>``,
where :ref:`n_r_max <varn_r_max>` is the number of radial grid points (see
:ref:`here <secGridNml>`). 
