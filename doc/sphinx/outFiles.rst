Output files
############

While some information of a run is written into ``STDOUT`` to monitor its
progress, most outputs are printed into dedicated files identified by the
extension ``.TAG``. Most of the information found in ``STDOUT`` is also written
to the log-file called ``log.TAG``. In addition, this file contains all input
parameters, truncation, information on other output files, and some results
like the time averaged energies (for ``l_average=.true.``). The python
functions that can be used to process the results can extract the informations
concerning the run from ``log.TAG``.  Other output files are organised in
columns or lines (datasets). Their meaning is explained below. The number of
graphic files (``G_#.TAG``), restart files (``rst_#.tag``), and spectrum files
(``(kin|mag)_spec_#.TAG``) are determined by the input parameters in the
namelist ``&output``. A new file is produced for each output time. If time
averaging has been chosen the time-averaged graphic, potential and spectra
files will have the number 0 and a prefix of the form ``_ave``. All other files
are numbered consecutively starting with 1. Times at which these files have
been written can be found in ``log.tag``. Other files (except the log-file)
contain time series. The frequency of storage for the kinetic and magnetic
energy files, the :ref:`rot.TAG <secRotFile>` file, the :ref:`dipole.TAG
<secDipoleTag>` file, the :ref:`par.TAG <secParFile>` file and the
:ref:`misc.TAG <secMiscFile>` file are determined by the log-times. Output
times for cmb files (``B_coeff_cmb.TAG``), coeff files
(``(T|V|B)_coeff_r.TAG``), potential files, torsional oscillations files and
movie frames are chosen independently (see above). Those files are stored as
unformatted files described below.

Log file: ``log.TAG``
=====================

Transport properties for the reference state
============================================

``anel.TAG``
------------

Standard time-series outputs
============================

.. _secEkinFile:

``e_kin.TAG``
-------------

This file contains the kinetic energy of the outer core, defined by

.. math::
   \begin{aligned}
   E_k = \frac{1}{2}\int_V \tilde{\rho}u^2\,{\rm d}V & = E_{pol}+E_{tor} \\
   & = \frac{1}{2}\sum_{\ell, m} \ell(\ell+1)\int_{r_i}^{r_o}\frac{1}{\tilde{\rho}}\left[
   \frac{\ell(\ell+1)}{r^2}|W_{\ell m}|^2+\left|\frac{{\rm d} W_{\ell m}}{{\rm d} r}\right|^2
   \right]\, {\rm d}r \\ 
   & +\frac{1}{2}\sum_{\ell, m} \ell(\ell+1)
   \int_{r_i}^{r_o}\frac{1}{\tilde{\rho}}|Z_{\ell m}|^2\,{\rm d} r
   \end{aligned}
   :label: eqEkin

The detailed calculations are done in the subroutine ``get_e_kin`` in the file ``kinetic_energy.f90``.  This file contains the following informations:

   +---------------+------------------------------------------------------+
   | No. of column | Contents                                             |
   +===============+======================================================+
   | 1             | time                                                 |
   +---------------+------------------------------------------------------+
   | 2	           | poloidal energy                                      |
   +---------------+------------------------------------------------------+
   | 3             | toroidal energy                                      |
   +---------------+------------------------------------------------------+
   | 4             | axisymmetric poloidal energy                         | 
   +---------------+------------------------------------------------------+
   | 5             | axisymmetric toroidal energy                         |
   +---------------+------------------------------------------------------+
   | 6             | equatorial symmetric poloidal energy                 |
   +---------------+------------------------------------------------------+
   | 7             | equatorial symmetric toroidal energy                 |
   +---------------+------------------------------------------------------+
   | 8             | equatorial symmetric and axisymmetric poloidal energy|
   +---------------+------------------------------------------------------+
   | 9             | equatorial symmetric and axisymmetric toroidal energy|
   +---------------+------------------------------------------------------+

.. _secEmagocFile:

``e_mag_oc.TAG``
----------------

This file contains the magnetic energy of the outer core, defined by

.. math::
   \begin{aligned}
   E_m = \frac{1}{2}\int_V B^2\,{\rm d}V & = E_{pol}+E_{tor} \\
   & = \frac{1}{2}\sum_{\ell, m} \ell(\ell+1)\int_{r_i}^{r_o}\left[
   \frac{\ell(\ell+1)}{r^2}|b_{\ell m}|^2+\left|\frac{{\rm d} b_{\ell m}}{{\rm d} r}\right|^2
   \right]\, {\rm d}r \\ 
   & +\frac{1}{2}\sum_{\ell, m} \ell(\ell+1)
   \int_{r_i}^{r_o}|j_{\ell m}|^2\,{\rm d} r
   \end{aligned}
   :label: eqEmag

The detailed calculations are done in the subroutine ``get_e_mag`` in the file ``magnetic_energy.f90``.  This file contains the following informations:

   +---------------+------------------------------------------------------+
   | No. of column | Contents                                             |
   +===============+======================================================+
   | 1             | time                                                 |
   +---------------+------------------------------------------------------+
   | 2             | outer core poloidal energy                           |
   +---------------+------------------------------------------------------+
   | 3             | outer core toroidal energy                           |
   +---------------+------------------------------------------------------+
   | 4             | outer core axisymmetric poloidal energy              |
   +---------------+------------------------------------------------------+
   | 5             | outer core axisymmetric toroidal energy              |
   +---------------+------------------------------------------------------+
   | 6             | outside potential field energy                       |
   +---------------+------------------------------------------------------+
   | 7             | outside axisymmetric potential field energy          |
   +---------------+------------------------------------------------------+
   | 8             | equatorial symmetric poloidal energy                 |
   +---------------+------------------------------------------------------+
   | 9             | equatorial symmetric toroidal energy                 |
   +---------------+------------------------------------------------------+
   | 10            | equatorial symmetric and axisymmetric poloidal energy|
   +---------------+------------------------------------------------------+
   | 11            | equatorial symmetric and axisymmetric toroidal energy|
   +---------------+------------------------------------------------------+
   | 12            | outside potential field energy                       |
   +---------------+------------------------------------------------------+
   | 13            | outside potential field axisymmetric energy          |
   +---------------+------------------------------------------------------+

.. _secEmagicFile:

``e_mag_ic.TAG``
----------------

This file contains the magnetic energy of the inner core. The detailed calculations are done in the subroutine ``get_e_mag`` in the file ``magnetic_energy.f90``.  This file contains the following informations:

   +---------------+------------------------------------------+
   | No. of column | Contents                                 |
   +===============+==========================================+
   | 1             | time                                     |
   +---------------+------------------------------------------+
   | 2             | inner core poloidal energy               |
   +---------------+------------------------------------------+
   | 3             | inner core toroidal energy               |
   +---------------+------------------------------------------+
   | 4             | inner core axisymmetric poloidal energy  |
   +---------------+------------------------------------------+
   | 5             | inner core axisymmetric toroidal energy  |
   +---------------+------------------------------------------+

.. _secRotFile:

``rot.TAG``
-----------

This files contains the rotation of the inner core and the mantle. Output concerning the rotation of inner core and mantle. This file is written by the subroutine ``write_rot`` in the file ``out_Rot.f90``.

   +---------------+--------------------------------+
   | No. of column | Contents                       |
   +===============+================================+
   | 1             | time                           |
   +---------------+--------------------------------+
   | 2             | Inner core rotation rate       |
   +---------------+--------------------------------+
   | 3             | Lorentz torque on inner core   |
   +---------------+--------------------------------+
   | 4             | viscous torque on inner core   |
   +---------------+--------------------------------+
   | 5             | mantle rotation rate           |
   +---------------+--------------------------------+
   | 6             | Lorentz torque on mantle       |
   +---------------+--------------------------------+
   | 7             | viscous torque on mantle       |
   +---------------+--------------------------------+

.. _secDipoleFile:

``dipole.TAG``
--------------

This file contains several informations about the magnetic dipole. This file is written by the subroutine ``get_e_mag`` in the file ``magnetic_energy.f90``.

   +---------------+---------------------------------------------------------------------------+
   | No. of column | Contents                                                                  |
   +===============+===========================================================================+
   | 1             | time                                                                      |
   +---------------+---------------------------------------------------------------------------+
   | 2             | tilt angle (colatitude in degrees) of the dipole                          |
   +---------------+---------------------------------------------------------------------------+
   | 3             | longitude (in degress) of dipole-pole                                     |
   +---------------+---------------------------------------------------------------------------+
   | 4             | relative energy of the axisymmetric dipole                                |
   +---------------+---------------------------------------------------------------------------+
   | 5             | relative energy of the axisymmetric dipole at the CMB                     |
   +---------------+---------------------------------------------------------------------------+
   | 6             | energy of the axisymmetric dipole at the CMB normalized with the          |
   |               | total energy up to spherical harmonic degree and order 11                 |
   +---------------+---------------------------------------------------------------------------+
   | 7             | relative energy of the total (axisymmetric and equatorial) dipole         |
   +---------------+---------------------------------------------------------------------------+
   | 8             | relative energy of the total (axisymmetric and equatorial) dipole         |
   |               | in the outer core                                                         |
   +---------------+---------------------------------------------------------------------------+
   | 9             | relative energy of the total dipole (axisymmetric and equatorial)         |
   |               | at the CMB                                                                |
   +---------------+---------------------------------------------------------------------------+
   | 10            | energy of the total (axisymmetric and equatorial) dipole at the CMB       |
   +---------------+---------------------------------------------------------------------------+
   | 11            | energy of the axisymmetric dipole at the CMB                              |
   +---------------+---------------------------------------------------------------------------+
   | 12            | energy of the dipole                                                      |
   +---------------+---------------------------------------------------------------------------+
   | 13            | energy of the axisymmetric dipole                                         |
   +---------------+---------------------------------------------------------------------------+
   | 14            | magnetic energy at the CMB                                                |
   +---------------+---------------------------------------------------------------------------+
   | 15            | magnetic energy up to spherical harmonic degree and order 11              |
   +---------------+---------------------------------------------------------------------------+
   | 16            | ratio between equatorial dipole energy and equatorial poloidal energy     |
   +---------------+---------------------------------------------------------------------------+
   | 17            | difference between energy at the CMB and equatorial symmetric             |
   |               | energy at the CMB, normalized by energy at the CMB                        |
   +---------------+---------------------------------------------------------------------------+
   | 18            | difference between energy at the CMB and axisymmetric energy at           |
   |               | the CMB, normalized by energy at the CMB                                  |
   +---------------+---------------------------------------------------------------------------+
   | 19            | difference between total energy and equatorial symmetric part             |
   |               | of the total energy, normalized by the total energy                       |
   +---------------+---------------------------------------------------------------------------+
   | 20            | difference between total energy and axisymmetric part of the              |
   |               | total energy, normalized by the total energy                              |
   +---------------+---------------------------------------------------------------------------+


.. _secParFile:

``par.TAG``
-----------

.. _secMiscFile:

``misc.TAG``
------------
