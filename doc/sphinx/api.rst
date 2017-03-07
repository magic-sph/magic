.. _secFortranAPI:

Description of the Fortran modules
##################################

The following pages contain an exhaustive description of the different
variables, subroutines and modules used in MagIC. This documentation is
automatically generated from the source code docstrings using the `Sphinx extention
for the Fortran domain <http://www.ifremer.fr/vacumm/user.desc.fdoc.dom.html>`_.

.. topic:: Fortran modules

     1. For the main program file ``magic.f90``, see :ref:`here <secFortranMagic>`.

     2. For the core modules that contain most of the global variables, see :ref:`here <secFortranBase>`.

     3. For the MPI related modules, see :ref:`here <secFortranMPI>`.

     4. For the code initialization and the pre-calculations done in the initial stage of the computation (before the time-stepping loop), see :ref:`here <secFortranInit>` and :ref:`there <secFortranPreCalc>`.

     5. For the time-stepping loop, see :ref:`here <secFortranTimeStep>`.

     6. For the calculation of the non-linear terms (in the physical space) and their time-advance, see :ref:`here <secFortranNl>`.

     7. For the calculation of the linear terms (in spectral space) and their time-advance, see :ref:`here <secFortranLin>`.

     8. For the Chebyshev, Fourier and Legendre transforms, see :ref:`here <secFortranLibs>`.

     9. For the computation of the radial derivatives (Chebyshev) and the integration, see :ref:`here <secFortranDerInt>`.

     10. For the definition of the blocking, see :ref:`here <secFortranBlocking>`.

     11. For the calculation of the standard outputs  (time-series, spectra and radial files), see :ref:`here <secFortranIOGeneral>`.

     12. For the calculation of binary outputs (graphic files, movie files, potential and coeff files), see :ref:`here <secFortranIOBinary>`.

     13. For the additional calculations of specific outputs (torsional oscillations, RMS force balance, etc.), see :ref:`here <secFortranIOAdd>`.

     14. For reading and writing the check points (restart files), see :ref:`here <secFortranRestart>`.

     15. For additional useful functions (string manipulation, etc.), see :ref:`here <secFortranMisc>`.

.. toctree::
   :hidden:

   apiFortran/magic.rst
   apiFortran/baseModules.rst
   apiFortran/parallelModules.rst
   apiFortran/initModules.rst
   apiFortran/precalcModules.rst
   apiFortran/timeAdvance
   apiFortran/timeAdvLinear
   apiFortran/timeAdvNonLinear
   apiFortran/legFourierChebAlgebra.rst
   apiFortran/derInt.rst
   apiFortran/blocking.rst
   apiFortran/IOGeneral.rst
   apiFortran/IOBinaryOutputs.rst
   apiFortran/IOAdd.rst
   apiFortran/checkPoints.rst
   apiFortran/misc.rst
