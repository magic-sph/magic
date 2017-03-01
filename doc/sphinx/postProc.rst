.. _secPythonPostProc:

Data visualisation and post-processing
######################################

Most of the :ref:`output files <secOutputFiles>` written during a run of MagIC
can be treated with the python post-processing classes and functions present in the
``$MAGIC_HOME/python/magic`` directory. These classes depend on several python libraries
that can be usually found in most of the Linux distributions.

Requirements
============

Hard dependencies
-----------------

* `python <https://www.python.org>`_ 2.7/3.3 or higher.

* `matplotlib <http://matplotlib.org>`_ 1.0 or higher.

* `scipy <http://www.scipy.org/>`_ 0.10 or higher.

Optional dependencies
---------------------

* Although entirely optional, the installation of `ipython
  <http://ipython.org/>`_  makes the interactive use of the post-processing python
  functions much more pleasant. Installing it is therefore recommanded for a smoother
  interactive usage of the python functions.

* The installation of the `basemap toolkit <http://matplotlib.org/basemap/>`_
  is optional. If installed, additional projections for the
  :py:class:`magic.Surf` (`Aitoff <https://en.wikipedia.org/wiki/Aitoff_projection>`_,
  `orthographic <https://en.wikipedia.org/wiki/Orthographic_projection>`_,
  `Mollweide <https://en.wikipedia.org/wiki/Mollweide_projection>`_, etc.) class 
  will be provided for 2-D surface plotting. Otherwise, the usage of
  :py:class:`magic.Surf` is limited to the 
  `Hammer projection <https://en.wikipedia.org/wiki/Hammer_projection>`_.


Configuration: ``magic.cfg`` file
=================================

A file named **magic.cfg** located in ``$MAGIC_HOME/python/magic/magic.cfg`` should have
been created when you used the ``source path/sourceme.sh`` command for the first time on
your machine. At that stage, it tried to **automatically fill the best options** that
correspond to your setup. Although tested on several various machine configurations, the
auto-configuration script might however fail on your setup. The paragraph below details
the possible options that you may want to adjust in the ``magic.cfg`` file.


Detailed options
----------------

In case, the file ``magic.cfg`` doesn't exist in the directory ``$MAGIC_HOME/python/magic``,
you can easily copy it from the default configuration ``magic.cfg.default`` and then adjust
the options manually:

.. code-block:: bash

   $ cp $MAGIC_HOME/python/magic/magic.cfg.default $MAGIC_HOME/python/magic/magic.cfg

In that file, you can set up the default `matplotlib  <http://matplotlib.org>`_  
rendering backend (among the possible options: ``TkAgg``, ``GTKAgg``, ``Qt4Agg``, ...).
The default configuration is

.. code-block:: python

   backend = TkAgg

.. note::

   This is usually the default configuration which is the most likely to work on
   supercomputing clusters.

If `LaTeX <http://www.latex-project.org/>`_ is installed on your work station, you
might also want to make use of the better looking LaTeX fonts for all your
displayed matplotlib figures (labels, caption, ticks, etc.). Be careful though that
most of the time LaTeX is **not installed**  on supercomputers. The default
configuration is thus:

.. code-block:: python

   labTex = False

You can change the default colormap that will be used in the plotting routines.

.. code-block:: python

   defaultCm = seismic

You cant change the default number of contours that will be used in the plotting routines.

.. code-block:: python

   defaultLevels = 65


If you want to enable all the features of the python functions (faster reading the 
:ref:`G_#.TAG <secGraphFile>`, conversion to the `VTK/VTS <http://www.vtk.org/>`_
file format, potential extrapolation of the field lines, etc.), some fortran libraries
present in the `$MAGIC_HOME/python/magic/fortranLib` directory need to be built using
the `f2py <http://www.f2py.com/>`_, which should be available on your Linux workstation
if all the required python libraries have been correctly installed. The boolean
``buildLib`` can control whether you want to try building the fortran libraries
with `f2py <http://www.f2py.com/>`_. The following configuration will try to build
the libraries:

.. code-block:: python

   buildLib = True

The exact name of the executable `f2py <http://www.f2py.com/>`_
however varies from one Linux distribution to the other. Among possible
options, one frequently finds: ``f2py``, ``f2py2``, ``f2py3``. This can be set to
your proper configuration using the ``f2pyexec`` option of the ``magic.cfg`` file.
The default configuration is:

.. code-block:: python

   f2pyexec = f2py2

You can also choose the fortran compiler you want to use on your machine. A list
of the installed compilers can be obtained by using (where ``f2py`` has to be
replaced by your own executable):

.. code-block:: bash
   
   $ f2py -c --help-fcompiler

The most frequent options are:

  * ``gnu95`` for the GNU gfortran compiler.
  * ``intelem`` for the Intel ifort compiler.
  * ``pg`` for the Portlang group pgf compiler.

Once you've decided the ideal configuration for your machine, set it up via the option
``fcompiler``:

.. code-block:: python
   
   fcompiler = intelem

Finally, he same configuration procedure can be applied to the C compiler using
the variable named ``ccompiler``. The possible options are:

  * ``unix`` for the GNU gcc compiler.
  * ``intelem`` for the Intel icc compiler.

In most of the configurations, the default configuration should do a good job:

.. code-block:: python
   
   ccompiler = unix

If you encounter any problem during the building stage, you can try playing with
this parameter though.

Ready?!
-------

Once you think you set up your ``magic.cfg`` file correctly, you can test your
configuration. If you decided to build the fortran libraries (i.e. ``buildLib=True``),
you can easily test it with any python shell by typing the following command::

   >>> from magic import *

If the build was successful, it should display::

   Please wait: building greader_single...
   Please wait: building greader_double...
   Please wait: building lmrreader_single...
   Please wait: building Legendre transforms...
   Please wait: building vtklib...

Once the libraries have been successfully installed, this message won't be displayed
again, except if you remove the ``*.so`` files that are now present in the
``$MAGIC_HOME/python/magic/`` directory.


.. _secPythonClasses:

Python functions and classes
============================

Once the python environment is correctly configured you can use the available functions and
classes to analyse and post-process your data. The following pages will give you the detailed
API of the available classes, as well as some practical examples:

.. topic:: Python classes
   
     1. To read the **log.TAG** files , see :ref:`here <secPythonLogFile>`.

     2. To read and analyse the time series, see :ref:`here <secPythonTimeSeries>`.

     3. To read and analyse the radial profiles, see :ref:`here <secPythonRadial>`.

     4. To read and analyse spectra **_spec_#.TAG**, see :ref:`here <secPythonSpectra>`.

     5. To read and analyse the **G_#.TAG** files, see :ref:`here <secPythonGraphFile>`.

     6. To read and analyse movie files **_mov.TAG**, see :ref:`here <secPythonMovieFiles>`.

     7. To read and analyse coeff files **[V|B|T]_coeff.TAG**, see :ref:`here <secPythonCoeffs>`.

     8. To read and analyse potential files **[V|B|T]_lmr_#.TAG**, see :ref:`here <secPythonPot>`.

     9. To read and analyse radial spectra **B[r|p]Spec.TAG**, see :ref:`here <secPythonRSpec>`.

     10. To read and analyse TO outputs, see :ref:`here <secPythonTO>`.

     11. To compare several runs simultaneously, see :ref:`here <secPythonComp>`.

     12. To transform the graphic files **G_#.TAG** to a file format readable by `paraview <http://www.paraview.org/>`_, `VisIt <https://visit.llnl.gov/>`_ or `mayavi <http://code.enthought.com/projects/mayavi/>`_ and do some fancy 3-D visualisation, see :ref:`here <secPythonParaview>`.

     13. For additional diagnostics (boundary layer, heat transport, interpolation on cylindrical grids, etc.), see :ref:`here <secPythonMisc>`.

     14. To take a look at the additional useful functions available (derivation, integration, interpolation, etc.), see :ref:`here <secPythonLibs>`.


.. toctree::
   :hidden:

   apiPython/timeSeries.rst
   apiPython/radial.rst
   apiPython/spectra.rst
   apiPython/graphMovCoeff.rst
   apiPython/outTO.rst
   apiPython/runComparison.rst
   apiPython/toParaview.rst
   apiPython/misc.rst
   apiPython/libs.rst
