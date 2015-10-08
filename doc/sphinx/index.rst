.. role:: mybox

:orphan:

Welcome
#######

.. container:: mybox

     **MagIC** is a numerical code that can simulate fluid dynamics in a spherical
     shell. MagIC solves for the Navier-Stokes equation including Coriolis force,
     optionally coupled with an induction equation for Magneto-Hydro Dynamics (MHD)
     and a temperature (or entropy) equation under both the anelastic and the Boussinesq 
     approximations.
  

     **MagIC** uses Chebyshev polynomials in the radial direction and spherical
     harmonic decomposition in the azimuthal and latitudinal directions. The
     time-stepping scheme relies on a semi-implicit `Crank-Nicolson
     <https://en.wikipedia.org/wiki/Crank–Nicolson_method>`_ for the linear terms of
     the MHD equations and a `Adams-Bashforth
     <https://en.wikipedia.org/wiki/Linear_multistep_method>`_ scheme for the
     non-linear terms and the Coriolis force.
     
     
     **MagIC** is written in Fortran and designed to be used on supercomputing
     clusters.  It thus relies on a hybrid parallelisation scheme using both `OpenMP
     <http://openmp.org/wp/>`_ and `MPI <http://www.open-mpi.org/>`_. Postprocessing
     functions written in python (requiring `matplotlib <http://matplotlib.org/>`_
     and `scipy <http://www.scipy.org/>`_) are also provided to allow a useful data
     analysis.
     
     
     **MagIC** is a free software. It can be used, modified and redistributed under the 
     terms of the `GNU GPL v3 licence <http://www.gnu.org/licenses/gpl-3.0.en.html>`_.

    .. only:: html

      .. raw:: html

        <div class="galleria">
           <img src=_images/benchmark.png data-title="Dynamo Benchmark" data-description="<a href='http://dx.doi.org/10.1016/S0031-9201(01)00275-8'>Christensen et al., PEPI, 2001</a>" >
           <img src=_images/3DBr.png data-title="The Jovian dynamo" data-description="<a href='http://dx.doi.org/10.1002/2014GL060814'>Gastine, T. et al., GRL, 2014</a>">
           <img src=_images/VsNS_E-4oIC-20000.jpg data-title="Inertial mode in a spherical Couette flow" data-description="<a href='http://dx.doi.org/10.1017/jfm.2013.545'>Wicht, J., JFM, 2014</a>">
        </div>

        <script>
           Galleria.loadTheme('galleria/themes/classic/galleria.classic.min.js');
           Galleria.run('.galleria', {responsive:true,height:0.4})
        </script>



Quickly starting using MagIC
============================

.. container:: mybox

   * The :ref:`quick-starting guide <secQuickStart>` will help you to download,
     set up and run your first numerical simulations using **MagIC**.
   
   * The :ref:`description of the input namelists <secNamelists>` will then help
     you to define the exact physical setup you may want to simulate.
   
   * Finally, the :ref:`python functions and classes <secPythonPostProc>` will
     allow you to do some advanced post-processing analyses on the outputs of **MagIC**.


Documentation
=============

.. container:: mybox

   * The :ref:`table of contents <contents>` gives an overview of the complete documentation.
   
   * The :ref:`formulation of the (M)HD problem <secEquations>` contains an exhaustive
     description of the equations solved by the MagIC code.

   * The :ref:`numerical methods section <secNumerics>` contains the description of the
     numerical technique.

   * The :ref:`search page <search>` allows to search the documentation.
   
   * The :ref:`fortran API <secFortranAPI>` contains a generic description of all
     Fortran variables, subroutines and modules used in **MagIC**.
   
   You can also download a :download:`PDF version <../magic_manual.pdf>` of this
   documentation generated from LaTeX Sphinx.

Contributing to the code
========================

.. container:: mybox

   If you want to contribute to **MagIC**, :ref:`the contributor
   guide<secContribute>` might be helpful for you.
