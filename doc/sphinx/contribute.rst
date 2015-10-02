.. _secContribute:

Contributing to the code
#########################

MagIC is an open-source code, we thus value any possible contribution! There are
several ways to directly contribute to the code:


.. topic:: Contribute

   * **Do you want to contribute to the code?** Just clone the code and start modyfing it.
     Make sure that your modifications :ref:`don't alter the code <secAutoTest>`, try
     to :ref:`document your changes <secDocSphinx>` as much as you can and follow
     the recommended :ref:`Fortran coding style <secModernFortran>`.

   * **Do you want to improve the documentation?** Feel free to document some missing
     features. The documentation is stored in the directory :code:`$MAGIC_HOME/doc/sphinx`
     and relies on the documenting tool `Sphinx <http://sphinx-doc.org/>`_. Some 
     recommendations regarding documentation can be found :ref:`below <secDocSphinx>`.

   * **Did you find a bug?** Issues and feature requests should be raised in the
     `github tracker <https://github.com/magic-sph/magic/issues>`_.

.. _secAutoTest:

Checking the consistency of the code
====================================

It is frequently required to check the consistency of the code, especially **after
the implementation of new features**. For this reason, we have the
`Perl <https://www.perl.org/>`_ script ``magic_checks.pl``, located in the
directory :code:`$MAGIC_HOME/samples/`, which tests the compilation of the code
and it's results against a set of standard solutions in sample directories to
check if the code produces the correct output. It has been initially ported from the
auto-test subroutines of the `pencil-code <https://github.com/pencil-code/>`_
developed by W. Dobler and adapted to the MagIC code.

You can run it as follows:

.. code-block:: bash

 ./magic_checks.pl <options>

It supports the following options:

.. code-block:: bash
 
  -h,  --help              Show usage overview
  -c,  --clean             Clean the directories when it is finished
  -a,  --all               All auto-tests are computed
       --level=LEV         Run only tests from level LEV
       --max-level=LEV     Run all tests below with level <= LEV (default: 0)
       --no-recompile      Compile only once
       --hybrid            Run the hybrid version

The ``level=LEV`` defines the priority level of check and validation of the
code. It has the following levels of checking:

  .. tabularcolumns:: |c|p{14cm}|  

  +---------+--------------------------------------------------------+
  | Level   |  Cases to check (subdirectories)                       |
  +=========+========================================================+
  | 0       | * Boussinesq dynamo benchmark                          |
  |         |   (`Christensen et al., 2001                           |
  |         |   <http://dx.doi.org/10.1016/S0031-9201(01)00275-8>`_) |
  |         |   - start from zero (``dynamo_benchmark``)             |
  |         | * Variable transport properties (viscosity,            |
  |         |   thermal diffusivity and electrical diffusivity)      | 
  |         |   in an anelastic convective model (``varProps``)      |
  |         | * Boussinesq dynamo benchmark                          |
  |         |   (`Christensen et al., 2001                           |
  |         |   <http://dx.doi.org/10.1016/S0031-9201(01)00275-8>`_) |
  |         |   - start from a saturated state (``boussBenchSat``)   |
  +---------+--------------------------------------------------------+
  | 1       | * Test reading and writing of                          |
  |         |   restart files (``testRestart``)                      |
  |         | * Test different grid truncations (``testTruncations``)|
  |         | * Test mapping on to a new grid (``testMapping``)      |
  |         | * Test different outputs produced (``testOutputs``)    |
  |         | * Test different radial outputs -                      |
  |         |   ``*R.TAG`` (``testRadialOutputs``)                   |
  +---------+--------------------------------------------------------+
  | 2       | * Hydrodynamic anelastic benchmark                     |
  |         |   (`Jones et al., 2011                                 |
  |         |   <http://dx.doi.org/10.1016/j.icarus.2011.08.014>`_)  |
  |         |   (``hydro_bench_anel``)                               |
  +---------+--------------------------------------------------------+
  | 3       | * Heat flux perturbation (``fluxPerturbation``)        |
  |         | * Isothermal model with :math:`N_{\rho}=3`             |
  |         |   (``isothermal_nrho3``)                               |
  |         | * Boussinesq Dynamo benchmark for conducting and       |
  |         |   rotating inner core                                  |
  |         |   (``dynamo_benchmark_condICrotIC``)                   |
  |         | * Anelastic dynamo with variable conductivity          |
  |         |   (``varCond``)                                        |
  +---------+--------------------------------------------------------+



.. _secModernFortran:

Advices when contributing to the code
=====================================

* Before commiting your modifications **always** make sure that the auto-tests pass correctly.

* Try to follow the same coding style rules as in the rest of the code:

  1. **Never** use TABS but always SPACES instead

  2. Use 3 spaces for indentation
  
     .. note::
        
	These two rules can be easily set in your $HOME/.vimrc file if you use
	`vim <http://www.vim.org/>`_:

	.. code-block:: vim

	    au FileType fortran set shiftwidth=3
	    au FileType fortran set tabstop=3
	    au FileType fortran set expandtab

  3. Never use capital letters for variable declaration
  
  4. Never use :code:`dimension(len)` for declaring array but rather :code:`real(cp) :: data(len)`
 
  5. Always use the default precisions when introducing new variables :code:`(cp)`


These rules try to follow the general recommendations on modern fortran programming
that can be found on `www.fortran90.org <http://www.fortran90.org/src/best-practices.html>`_ or in the book
`Modern Fortran - style and usage <http://www.cambridge.org/us/academic/subjects/computer-science/scientific-computing-scientific-software/modern-fortran-style-and-usage>`_ by N. S. Clerman and W. Spector.


.. _secDocSphinx:

Building the documentation and contributing to it
=================================================

The documentation is generated using `Sphinx <http://sphinx-doc.org/>`_. To
build it you'll thus need to install this python module on your machine. This is in general
directly available on most of the Linux distributions under the name
``python-sphinx``. Once installed, just go to the documentation directory

.. code-block:: bash

   $ cd $MAGIC_HOME/doc/sphinx

and build the html documentation

.. code-block:: bash

   $ make html

The complete documentation will then be built in a local directory named
:code:`$MAGIC_HOME/doc/sphinx/.build/html`. 

If `LaTeX <http://www.latex-project.org/>`_ is installed on your work station, it 
is also possible to build the corresponding manual of the documentation in 
the pdf format:

.. code-block:: bash

  $ make latexpdf

The resulting pdf is then generated in a local directory named
:code:`$MAGIC_HOME/doc/sphinx/.build/latex`. 

It is pretty straightforward to contribute to the documentation by simply adding some
contents to the different :code:`rst` files. Informations about `reStructuredText <http://docutils.sourceforge.net/rst.html>`_ syntax can be found on `www.sphinx-doc.org <http://sphinx-doc.org/rest.html>`_,
while helpful CheatSheet are accessible `here <http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html>`_ or `there <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_.
