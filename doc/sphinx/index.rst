.. role:: mybox

:orphan:

Welcome
#######

.. container:: mybox

    .. only:: html

      .. raw:: html

        <div class="galleria">
           <img src=_images/benchmark.png data-title="Dynamo benchmark" data-description="<a href='http://dx.doi.org/10.1016/S0031-9201(01)00275-8'>Christensen et al., PEPI, 2001</a>" >
           <img src=_images/3DBr.png data-title="A numerical model of the Jovian dynamo" data-description="<a href='http://dx.doi.org/10.1002/2014GL060814'>Gastine, T. et al., GRL, 2014</a>">
           <img src=_images/VsNS_E-4oIC-20000.jpg data-title="Inertial mode in a spherical Couette flow" data-description="<a href='http://dx.doi.org/10.1017/jfm.2013.545'>Wicht, J., JFM, 2014</a>">
           <img src=_images/Conv.png data-title="Rayleigh-Bénard convection in non-rotating spherical shell" data-description="<a href='http://dx.doi.org/10.1017/jfm.2015.401'>Gastine, T. et al., JFM, 2015</a>">
           <img src=_images/wielandMars.png data-title="A hemispherical dynamo model to explain the Martian crustal magnetization" data-description="<a href='http://dx.doi.org/10.1016/j.pepi.2013.01.001'>Dietrich, W. et al., PEPI, 2013</a>">
           <img src=_images/MRI_sphere.png data-title="Magneto-rotational instability in spherical shell" data-description="<a href='http://dx.doi.org/10.1051/0004-6361/201425240'>Jouve, L. et al., A&amp;A, 2015</a>">
           <img src=_images/kristaEuropa.jpg data-title="Ocean-driven heating of Europa's icy shell at low latitudes" data-description="<a href='http://dx.doi.org/10.1038/ngeo2021'>Soderlund, K. et al., Nature Geoscience, 2014</a>">
           <img src=_images/AubertNature.png data-title="Numerical model of the geodynamo" data-description="<a href='http://dx.doi.org/10.1038/nature07109'>Aubert, J. et al., Nature, 2008</a>">
           <img src=_images/caoIcarus.png data-title="Spherical Couette dynamo to explain the Saturnian magnetic field" data-description="<a href='http://dx.doi.org/10.1016/j.icarus.2012.08.007'>Cao, H. et al., Icarus, 2012</a>">
           <img src=_images/yadav2015.png data-title="Formation of polar spots in a fully-convective star model" data-description="<a href='http://dx.doi.org/10.1051/0004-6361/201424589'>Yadav, R. et al., A&amp;A, 2015</a>">
           <img src=_images/heimpel05.png data-title="Numerical model of the Jovian zonal jets" data-description="<a href='http://dx.doi.org/10.1038/nature04208'>Heimpel, M. et al., Nature, 2005</a>">
        </div>

        <script>
           Galleria.loadTheme('galleria/themes/classic/galleria.classic.min.js');
           Galleria.run('.galleria', {responsive:true, height:0.4, autoplay: 3000, transition:'fadeslide', pauseOnInteraction: false});
        </script>

     **MagIC** is a numerical code that can simulate fluid dynamics in a spherical
     shell. MagIC solves for the Navier-Stokes equation including Coriolis force,
     optionally coupled with an induction equation for Magneto-Hydro Dynamics (MHD),
     a temperature (or entropy) equation and an equation for chemical composition
     under both the anelastic and the Boussinesq approximations.

     **MagIC** uses Chebyshev polynomials or finite difference in the radial 
     direction and spherical harmonic decomposition in the azimuthal and latitudinal
     directions. MagIC supports several Implicit-Explicit time schemes where the
     nonlinear terms and the Coriolis force are treated explicitly, while the
     remaining linear terms are treated implicitly.
          
     **MagIC** is written in Fortran and designed to be used on supercomputing
     clusters.  It thus relies on a hybrid parallelisation scheme using both `OpenMP
     <http://openmp.org/wp/>`_ and `MPI <http://www.open-mpi.org/>`_. Postprocessing
     functions written in python (requiring `matplotlib <http://matplotlib.org/>`_
     and `scipy <http://www.scipy.org/>`_) are also provided to allow a useful data
     analysis.

     .. figure:: figs/magic_occigen.png
        :width: 450px
        :align: center
        :alt: caption

        Walltime of MagIC on `Occigen
        <https://www.cines.fr/en/supercomputing-2/hardwares/the-supercomputer-occigen/>`_ versus number of cores
        for a Boussinesq dynamo model computed at three different numerical resolutions
        :math:`(N_r,\ell_{\text{max}})`.
     
     
     **MagIC** is a free software. It can be used, modified and redistributed under the 
     terms of the `GNU GPL v3 licence <http://www.gnu.org/licenses/gpl-3.0.en.html>`_.



Quickly starting using MagIC
============================

.. container:: mybox

   * The :ref:`quick-starting guide <secQuickStart>` will help you to download,
     set up and run your first numerical simulations using **MagIC**.
   
   * The :ref:`description of the input namelists <secNamelists>` will then help
     you to define the exact physical setup you may want to simulate.

   * The :ref:`description of the output files <secOutputFiles>` will help you to
     understand what are the diagnostic quantities computed by **MagIC**.
   
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

Giving credit
=============

.. container:: mybox

   In case you intend to publish scientific results obtained with **MagIC**
   or present them in a conference, we (the developers of MagIC) kindly
   ask to be acknowledged with a reference to the website 
   https://magic-sph.github.io/ or https://github.com/magic-sph/magic.
   
   We also suggest to give appropriate reference to one or several of the following
   papers:
   
   * Boussinesq equations: `Wicht (2002, PEPI, 132, 281-302) <http://dx.doi.org/10.1016/S0031-9201(02)00078-X>`_
   
   * Anelastic equations: `Gastine & Wicht (2012, Icarus, 219, 28-442) <http://dx.doi.org/10.1016/j.icarus.2012.03.018>`_
   
   * Boussinesq benchmark: `Christensen et al. (2001, PEPI, 128, 25-34) <http://dx.doi.org/10.1016/S0031-9201(01)00275-8>`_

   * Benchmark for double diffusive convection: `Breuer et al. (2010, GJI, 183, 150-162) <http://dx.doi.org/11.1111/j.1365-246X.2010.04722.x>`_

   * Anelastic benchmark: `Jones et al. (2011, Icarus, 216, 120-135) <http://dx.doi.org/10.1016/j.icarus.2011.08.014>`_

   * In case you use the `SHTns <https://bitbucket.org/bputigny/shtns-magic>`_ library for the spherical harmonics transforms (MagIC 5.3 or later), please also cite: `Schaeffer (2013, GGG, 14, 751-758) <http://dx.doi.org/10.1002/ggge.20071>`_

  
   .. seealso:: A (tentative) comprehensive list of the publications that have 
                been produced to date (May 2019) using **MagIC**
                is accessible `here <https://ui.adsabs.harvard.edu/public-libraries/LVt1vdaKQsC5P09In2iloA>`_.
                To date, more than **100 publications** have been-accepted in
                more than 10 different peer-reviewed journals: `PEPI
                <http://www.journals.elsevier.com/physics-of-the-earth-and-planetary-interiors/>`_
                (22), `Icarus <http://www.journals.elsevier.com/icarus/>`_ (11), `E&PSL
                <www.journals.elsevier.com/earth-and-planetary-science-letters/>`_ (7), `GJI
                <http://gji.oxfordjournals.org/>`_ (8), `A&A <http://www.aanda.org/>`_ (6), 
                `GRL <http://agupubs.onlinelibrary.wiley.com/agu/journal/10.1002/(ISSN)1944-8007/>`_ (4), 
                `JFM <http://journals.cambridge.org/action/displayJournal?jid=FLM>`_ (6), 
                `GAFD <http://www.tandfonline.com/toc/ggaf20/current>`_ (3),
                `Nature <http://www.nature.com/nature>`_ (2), etc.

   .. figure:: figs/magic_pubs.png
      :width: 600px
      :align: center
      :alt: caption

      Number of peer-reviewed publications produced using **MagIC**



