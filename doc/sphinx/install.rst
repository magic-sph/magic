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

Go to the directory where the source files of MagIC are contained

.. code-block:: bash

   $ cd $MAGIC_HOME/src 
   
and edit the ``Makefile`` there to specify your specific compiler options.

Makefile options
----------------

**Select compiler**

Set a suitable compiler in the first line of the ``Makefile``: ``COMPILER =
<compiler_phrase>``. The options are ``intel``, ``gnu`` or ``amd`` - depending
on your available compilers.

*List of default compilers*

  +-----------------+---------------+------------------+ 
  | Compiler Option |    Normal     |     With MPI     |
  +-----------------+---------------+------------------+
  | intel           | ifort, icc    | mpiifort, mpiicc |
  +-----------------+---------------+------------------+
  | gnu             | gfortran, gcc | mpif90, mpicc    |
  +-----------------+---------------+------------------+
  | amd             | openf95       |                  |
  +-----------------+---------------+------------------+

**Select compiling options**

* ``PRECISION`` Set it to 'dble' for double-precision calculations or to 'sngl' for single-precision calculations
* ``OUT_PREC`` Set it to 'dble' for double-precision in binary outputs or to 'sngl' for single precision
* ``PRODRUN`` Set it to ``yes`` for production run, ``no`` for debugging.
* ``USE_MPI`` Set to ``yes`` to use MPI, set it to ``no`` if you want a serial version of the code .
* ``OPENMP``  Set it to ``yes`` to use the hybrid version of the code, or to ``no`` for a pure MPI (or serial) version.
* ``DEBUG``   Set to ``all`` to enable the full debug flags. *While running in debugging mode, set* ``PRODRUN`` *to* ``no``. 
* ``USE_FFTLIB`` This option lets you select the library you want to use for Fast Fourier Transforms. This can be set to 'JW' or 'MKL'. 'JW' refers to the inbuilt library by **J** ohannes **W** icht, while 'MKL' refers to the `Intel Math Kernel Library <https://software.intel.com/en-us/intel-mkl>`_. Use 'JW' if you don't have Intel MKL installed.
* ``USE_MKL`` Set to ``yes`` if you have Intel MKL installed and want to use it for matrix operations.
* ``USE_HDF5`` Set to ``yes`` if you want the restart file to be written in the  `HDF5 <http://www.hdfgroup.org/>`_ format

**Architecture (Intel compilers only)**

If you're using intel compilers and if your computer is capable of following
specific intel instruction sets (sse3 or AVX), then the ``Makefile``
automatically should automatically detects and sets ``FFLAG_ARCH_OPT = -xsse3``
or ``FFLAG_ARCH_OPT = -xAVX`` under intel compiler options.

**MPI_INCPATH**

This sets the path for your mpi header file ``mpif.h`` . The path depends on
the computer. For PCs, this is commonly ``/usr/include`` or
``/usr/include/mpi`` and should be found by the ``Makefile`` automatically thanks
to the command ``mpif90 --showme:incdirs``. In case this doesn't work, you may
need to specify this variable manually in the ``Makefile``. On supercomputing clusters,
this variable is in general not used since the ``mpi.mod`` file is usually find
the standard ``$PATH``.

**Other compilers**

If your available compilers are different from the options provided in the
``Makefile``, then just create a new profile for your desired compiler
by changing the options ``COMP_FC`` and
``COMP_CC`` for serial fortran and C compilers and ``COMP_MPFC`` and
``COMP_MPCC`` for compilers with mpi implementation.


Compiling the code
------------------

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

* Running a serial version of the code (``USE_MPI=no`` and ``OPENMP=no``):

  .. code-block:: bash

     $ ./magic.exe input.nml

* Running the code without OpenMP (``USE_MPI=yes`` and ``OPENMP=no``) with ``<n_mpi>`` 
  MPI ranks:
  
  .. code-block:: bash

     $ mpiexec -n <n_mpi> ./magic.exe input.nml

* Running the hybrid code (``USE_MPI=yes`` and ``OPENMP=yes``) with ``<n_mpi>`` MPI ranks 
  and ``<n_omp>`` OpenMP threads:
  
  .. code-block:: bash

     $ export OMP_NUM_THREAD = <n_omp>
     $ mpiexec -n <n_mpi> ./magic.exe input.nml

Note that the :ref:`n_r_max <varn_r_max>` must be a multiple of ``<n_mpi>``,
where :ref:`n_r_max <varn_r_max>` is the number of radial grid points (see
:ref:`here <secGridNml>`). 
