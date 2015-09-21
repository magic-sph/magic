
Contributing to the code
#########################

.. _secAutoTest:

Checking the consistency of the code
====================================

It is often necessary to check the consistency of the code, especially after one has implemented something new in the code. For this reason, we have the `Perl <https://www.perl.org/>`_ script ``magic_checks.pl``, located in the directory :code:`$MAGIC_HOME/samples/`, which tests the compilation of the code and it's results against a set of standard solutions in sample directories to check if the code produces the correct output. It has been adapted from the auto-test subroutines of the `pencil-code <https://github.com/pencil-code/>`_ developed by W. Dobler.

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

The ``level=LEV`` defines the priority level of check and validation of the code. It has the following levels of checking:

  .. tabularcolumns:: |c|p{12cm}|  

  +---------+--------------------------------------------------------+
  | Level   |  Cases to check (subdirectories)                       |
  +=========+========================================================+
  | 0       | * The dynamo benchmark (``dynamo_benchmark``)          |
  |         | *           (``varProps``)                             |
  |         | *           (``boussBenchSat``)                        |
  +---------+--------------------------------------------------------+
  | 1       | * Test reading and writing of                          |
  |         |   restart files (``testRestart``)                      |
  |         | * Test different grid truncations (``testTruncations``)|
  |         | * Test mapping on to a new grid (``testMapping``)      |
  |         | * Test different outputs produced (``testOutputs``)    |
  |         | * Test different radial outputs -                      |
  |         |   ``*R.TAG`` (``testRadialOutputs``)                   |
  +---------+--------------------------------------------------------+
  | 2       | * Hydrodynamic Anelastic Benchmark                     |
  |         |   (``hydro_bench_anel``)                               |
  +---------+--------------------------------------------------------+
  | 3       | * Heat flux perturbation (``fluxPerturbation``)        |
  |         | * Isothermal model with :math:`N_{\rho}=3`             |
  |         |   (``isothermal_nrho3``)                               |
  |         | * Dynamo benchmark for conducting and rotating         |
  |         |   inner core (``dynamo_benchmark_condICrotIC``)        |
  |         | * Anelastic dynamo with variable conductivity          |
  |         |   (``varCond``)                                        |
  +---------+--------------------------------------------------------+



Modifying/contributing to the code
==================================

* Before commiting your modifications **always** make sure that the auto-tests pass correctly.

* Try to follow the same coding style rules as in the rest of the code:

  1. **Never** use TABS but always SPACES instead
  2. Use 3 spaces for indentation
  3. Never use capital letters for variable declaration
  4. Never use :code:`dimension(len)` for declaring array but rather :code:`real(cp) :: data(len)`
  5. Always use the default precisions when introducing new variables :code:`(cp)`


  More on that topic `here <http://www.fortran90.org/src/best-practices.html>`_.
