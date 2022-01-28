.. _secTimeSeriesStd:

Default time-series outputs
===========================

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

The detailed calculations are done in the subroutine :f:subr:`get_e_kin <kinetic_energy/get_e_kin()>`.  This file contains the following informations:

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

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the e_kin.TAG files of the current directory
   >>> ts = MagicTs(field='e_kin', all=True)
   >>> # To only read e_kin.N0m2
   >>> ts = MagicTs(field='e_kin', tag='N0m2')

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

The detailed calculations are done in the subroutine :f:subr:`get_e_mag <magnetic_energy/get_e_mag()>`.  This file contains the following informations:

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

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the e_mag_oc.TAG files of the current directory
   >>> ts = MagicTs(field='e_mag_oc', all=True)
   >>> # To only read e_mag_oc.N0m2
   >>> ts = MagicTs(field='e_mag_oc', tag='N0m2')

.. _secEmagicFile:

``e_mag_ic.TAG``
----------------

This file contains the magnetic energy of the inner core. The detailed
calculations are done in the subroutine :f:subr:`get_e_mag <magnetic_energy/get_e_mag()>`.
This file contains the following informations:

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

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the e_mag_ic.TAG files of the current directory
   >>> ts = MagicTs(field='e_mag_ic', all=True)
   >>> # To only read e_mag_ic.N0m2
   >>> ts = MagicTs(field='e_mag_ic', tag='N0m2')


.. _secRotFile:

``rot.TAG``
-----------

This files contains the rotation of the inner core and the mantle. Output
concerning the rotation of inner core and mantle. This file is written by the
subroutine :f:subr:`write_rot <outrot/write_rot()>`.

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

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the rot.TAG files of the current directory
   >>> ts = MagicTs(field='rot', iplot=False, all=True)


.. _secDipoleFile:

``dipole.TAG``
--------------

This file contains several informations about the magnetic dipole. This file is written by the subroutine :f:subr:`get_e_mag <magnetic_energy/get_e_mag()>`. The maximum degree used to compute columns 6 and 15 is given by :ref:`l_geo <varl_earth_like>`.

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

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the dipole.TAG files of the current directory
   >>> ts = MagicTs(field='dipole', all=True)



.. _secParFile:

``par.TAG``
-----------

This files contains the outputs of several parameters that describe flow and
magnetic fields (Reynolds number, Elsasser number, flow lengthscales, etc.).
This file is written by the subroutine :f:subr:`output <output_mod/output()>`.

   +---------------+-----------------------------------------+
   | No. of column | Contents                                |
   +===============+=========================================+
   | 1             | time                                    |
   +---------------+-----------------------------------------+
   | 2             | (magnetic) Reynolds number              |
   +---------------+-----------------------------------------+
   | 3             | Elsasser number                         |
   +---------------+-----------------------------------------+
   | 4             | Local Rossby number Rol                 |
   +---------------+-----------------------------------------+
   | 5             | Realtive geostrophic kinetic energy     |
   +---------------+-----------------------------------------+
   | 6             | Total dipolarity                        |
   +---------------+-----------------------------------------+
   | 7             | CMB dipolarity                          |
   +---------------+-----------------------------------------+
   | 8             | Axial flow length scale dlV             |
   +---------------+-----------------------------------------+
   | 9             | Flow length scale dmV                   |
   +---------------+-----------------------------------------+
   | 10            | Flow length scale dpV                   |
   +---------------+-----------------------------------------+
   | 11            | Flow length scale dzV                   |
   +---------------+-----------------------------------------+
   | 12            | Dissipation length scale lvDiss         |
   +---------------+-----------------------------------------+
   | 13            | Dissipation length scale lbDiss         |
   +---------------+-----------------------------------------+
   | 14            | Magnetic length scale dlB               |
   +---------------+-----------------------------------------+
   | 15            | Magnetic length scale dmB               |
   +---------------+-----------------------------------------+
   | 16            | Elsasser number at CMB                  |
   +---------------+-----------------------------------------+
   | 17            | Local Rol based on non-ax. flow         |
   +---------------+-----------------------------------------+
   | 18            | Convective flow length scale dlVc       |
   +---------------+-----------------------------------------+
   | 19            | Peak of the poloidal kinetic energy     |
   +---------------+-----------------------------------------+
   | 20            | CMB zonal flow at the equator           |
   +---------------+-----------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the par.TAG files of the current directory
   >>> ts = MagicTs(field='par', all=True)
