
Additional optional time-series outputs
=======================================

.. _secAMFile:

``AM.TAG``
-------------

.. note:: This file is **only** written when :ref:`l_AM=.true. <varl_AM>`

This file contains the time series of the angular momentum of the inner core, the outer
core and the mantle. This file is written by the subroutine :f:subr:`write_rot <outrot/write_rot()>`.

  +---------------+-----------------------------------------------------+
  | No. of column | Contents                                            |
  +===============+=====================================================+
  | 1             | time                                                |
  +---------------+-----------------------------------------------------+
  | 2             | angular momentum of the outer core                  |
  +---------------+-----------------------------------------------------+
  | 3             | angular momentum of the inner core                  |
  +---------------+-----------------------------------------------------+
  | 4             | angular momentum of the mantle                      |
  +---------------+-----------------------------------------------------+
  | 5             | total angular momentum                              |
  +---------------+-----------------------------------------------------+
  | 6             | relative in angular momentum, per time step         |
  +---------------+-----------------------------------------------------+
  | 7             | total kinetic angular momentum                      |
  +---------------+-----------------------------------------------------+
  | 8             | relative change in kinetic energy, per time step    |
  +---------------+-----------------------------------------------------+
  | 9             | kinetic angular momentum of the inner core          |
  +---------------+-----------------------------------------------------+
  | 10            | kinetic angular momentum of the outer core          |
  +---------------+-----------------------------------------------------+
  | 11            | kinetic angular momentum of the mantle              |
  +---------------+-----------------------------------------------------+


This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the AM.TAG files of the current directory
   >>> ts = MagicTs(field='AM', all=True)


.. _secpowerFile:

``power.TAG``
-------------

.. note:: This file is **only** written when :ref:`l_power=.true. <varl_power>`

This file contains the power budget diagnostic. This file is computed by the subroutine
:f:subr:`get_power <power/get_power()>`.

   +---------------+------------------------------------------------------------------+
   | No. of column | Contents                                                         |
   +===============+==================================================================+
   | 1             | time                                                             |
   +---------------+------------------------------------------------------------------+
   | 2             | Buoyancy power: :math:`Ra\,g(r)\,\langle u_r T'\rangle_s`        |
   +---------------+------------------------------------------------------------------+
   | 3             | Viscous power at the inner boundary (ICB)                        |
   +---------------+------------------------------------------------------------------+
   | 4             | Viscous power at the outer boundary (CMB)                        |
   +---------------+------------------------------------------------------------------+
   | 5             | Viscous dissipation: :math:`\langle(\nabla \times u)^2\rangle_s` |
   +---------------+------------------------------------------------------------------+
   | 6             | Ohmic dissipation: :math:`\langle(\nabla \times B)^2\rangle_s`   |
   +---------------+------------------------------------------------------------------+
   | 7             | Total power at the CMB (viscous + Lorentz)                       |
   +---------------+------------------------------------------------------------------+
   | 8             | Total power at the ICB (viscous + Lorentz)                       |
   +---------------+------------------------------------------------------------------+
   | 9             | Total power                                                      |
   +---------------+------------------------------------------------------------------+
   | 10            | Time variation of total power                                    |
   +---------------+------------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack the files that match the pattern  ``power.N0m2*``
   >>> ts = MagicTs(field='power', tags='N0m2*')

.. _secdtEFile:

``dtE.TAG``
-----------

.. note:: This file is **only** written when :ref:`l_power=.true. <varl_power>`

This file contains the time-derivatives of the total energy. It allows to accurately
monitor how the total energy varies with time. This file is generated
by the subroutine :f:subr:`output <output_mod/output()>`.

   +---------------+------------------------------------------------------------------+
   | No. of column | Contents                                                         |
   +===============+==================================================================+
   | 1             | time                                                             |
   +---------------+------------------------------------------------------------------+
   | 2             | time-derivative of the total energy :math:`\partial E/\partial t`|
   +---------------+------------------------------------------------------------------+
   | 3             | integrated time variation of the total energy                    |
   +---------------+------------------------------------------------------------------+
   | 4             | relative time variation of the total energy                      |
   +---------------+------------------------------------------------------------------+

.. _secu_squareFile:

``u_square.TAG``
----------------

.. note:: This file is **only** written in anelastic models, i.e. either when
          :ref:`strat/=0 <varstrat>` or when :ref:`interior_model/="None" <varinterior_model>`

This file contains the square velocity of the outer core. It is actually very similar
to the :ref:`e_kin.TAG <secEkinFile>` file, except that the density background
:math:`\tilde{\rho}` is removed:

.. math::
   \begin{aligned}
   {\cal U} = \frac{1}{2}\int_V u^2\,{\rm d}V & = {\cal U}_{pol}+{\cal U}_{tor} \\
   & = \frac{1}{2}\sum_{\ell, m} \ell(\ell+1)\int_{r_i}^{r_o}\frac{1}{\tilde{\rho}^2}\left[
   \frac{\ell(\ell+1)}{r^2}|W_{\ell m}|^2+\left|\frac{{\rm d} W_{\ell m}}{{\rm d} r}\right|^2
   \right]\, {\rm d}r \\ 
   & +\frac{1}{2}\sum_{\ell, m} \ell(\ell+1)
   \int_{r_i}^{r_o}\frac{1}{\tilde{\rho}^2}|Z_{\ell m}|^2\,{\rm d} r
   \end{aligned}

The detailed calculations are done in the subroutine :f:subr:`get_u_square <kinetic_energy/get_u_square()>`.  This file contains the following informations:

  +----------------+--------------------------------------------------------------------+
  | No. of columns | Contents                                                           |
  +================+====================================================================+
  | 1	           | time                                                               |
  +----------------+--------------------------------------------------------------------+
  | 2              | poloidal part :math:`{\cal U}_{pol}`                               |
  +----------------+--------------------------------------------------------------------+
  | 3              | toroidal part :math:`{\cal U}_{pol}`                               |
  +----------------+--------------------------------------------------------------------+
  | 4              | axisymmetric contribution to the poloidal part                     |
  +----------------+--------------------------------------------------------------------+
  | 5              | axisymmetric contribution to the toroidal part                     |
  +----------------+--------------------------------------------------------------------+
  | 6              | Rossby number: :math:`Ro=E\,\sqrt{\frac{2{\cal U}}{V}}`            |
  +----------------+--------------------------------------------------------------------+
  | 7              | Magnetic Reynolds number: :math:`Rm=Pm\,\sqrt{\frac{2{\cal U}}{V}}`|
  +----------------+--------------------------------------------------------------------+
  | 8              | local Rossby number: :math:`Ro_l=Ro\frac{d}{l}`                    |
  +----------------+--------------------------------------------------------------------+
  | 9              | average flow length scale: :math:`l`                               |
  +----------------+--------------------------------------------------------------------+
  | 10             | local Rossby number based on the non-axisymmetric components       |
  |                | of the flow                                                        |
  +----------------+--------------------------------------------------------------------+
  | 11             | average flow length scale based on the non-axisymmetric            |
  |                | components of the flow                                             |
  +----------------+--------------------------------------------------------------------+


This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the u_square.TAG files of the current directory
   >>> ts = MagicTs(field='u_square', all=True)

.. _secdriftFile:

``drift[V|B][D|Q].TAG``
-----------------------

.. note:: These files are **only** written when :ref:`l_drift=.true. <varl_drift>`

These files store spherical harmonic coefficients of the toroidal (poloidal) potential of the flow (magnetic) field, only for :math:`\ell=m` or :math:`\ell=m+1` depending on the symmetry - ``D`` for **D** ipolar and ``Q`` for **Q** uadrupolar. The coefficients are stored at different three different radial levels - ``n_r1, nr_2, n_r3`` for the velocity and two different radial levels - ``n_r1`` and ``n_r2`` - for the magnetic field.



The symmetries can be summarized below:

 +---------+-----------------+-----------------+
 | Field   | Dipolar         | Quadrupolar     | 
 +=========+=================+=================+
 | Velocity| :math:`\ell=m`  | :math:`\ell=m+1`|
 +---------+-----------------+-----------------+
 | Magnetic| :math:`\ell=m+1`| :math:`\ell=m`  |
 +---------+-----------------+-----------------+

:math:`\ell+m=` even for toroidal potential refers to an equatorially antisymmetric field (*Dipolar*), while the same for a poloidal potential is associated with an equatorially symmetric field (*Quadrupolar*). The sense is opposite when :math:`\ell+m=` odd. This is the reason for the choice of selecting these specific coefficients.

The columns of the files look like follows:

For the flow field:

 * n_r1 = (1/3) * :ref:`n_r_max-1 <varn_r_max>`
 * n_r2 = (2/3) * :ref:`n_r_max-1 <varn_r_max>`
 * n_r3 = :ref:`n_r_max-1 <varn_r_max>`
 
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | Column no.|   DriftVD.TAG                                     |       DriftVQ.TAG                               |
 +===========+===================================================+=================================================+
 | 1         | Time                                              | Time                                            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 |                                                                                                                 |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 2         | :math:`z` (:ref:`minc <varminc>`, minc) at n_r1   | :math:`z` (:ref:`minc+1<varminc>`, minc) at n_r1|
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 3         | :math:`z` (2*minc, 2*minc) at n_r1                | :math:`z` (2*minc+1, 2*minc) at n_r1            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 4         | :math:`z` (3*minc, 3*minc) at n_r1                | :math:`z` (3*minc+1, 3*minc) at n_r1            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 5         | :math:`z` (4*minc, 4*minc) at n_r1                | :math:`z` (4*minc+1, 4*minc) at n_r1            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 |                                                                                                                 |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 6         | :math:`z` (minc, minc) at n_r2                    | :math:`z` (minc+1, minc) at n_r2                |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 7         | :math:`z` (2*minc, 2*minc) at n_r2                | :math:`z` (2*minc+1, 2*minc) at n_r2            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 8         | :math:`z` (3*minc, 3*minc) at n_r2                | :math:`z` (3*minc+1, 3*minc) at n_r2            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 9         | :math:`z` (4*minc, 4*minc) at n_r2                | :math:`z` (4*minc+1, 4*minc) at n_r2            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 |                                                                                                                 |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 10        | :math:`z` (minc, minc) at n_r3                    | :math:`z` (minc+1, minc) at n_r3                |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 11        | :math:`z` (2*minc, 2*minc) at n_r3                | :math:`z` (2*minc+1, 2*minc) at n_r3            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 12        | :math:`z` (3*minc, 3*minc) at n_r3                | :math:`z` (3*minc+1, 3*minc) at n_r3            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 13        | :math:`z` (4*minc, 4*minc) at n_r3                | :math:`z` (4*minc+1, 4*minc) at n_r3            |
 +-----------+---------------------------------------------------+-------------------------------------------------+

For the magnetic field:

 * n_r1 = :f:var:`n_r_ICB <n_r_cmb>`
 * n_r2 = :f:var:`n_r_CMB <n_r_icb>`
 
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | Column no.|   DriftBD.TAG                                     |       DriftBQ.TAG                               |
 +===========+===================================================+=================================================+
 | 1         | Time                                              | Time                                            |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 |                                                                                                                 |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 2         | :math:`b` (:ref:`minc+1 <varminc>`, minc) at n_r1 | :math:`b` (:ref:`minc<varminc>`, minc) at n_r1  |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 3         | :math:`b` (2*minc+1, 2*minc) at n_r1              | :math:`b` (2*minc, 2*minc) at n_r1              |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 4         | :math:`b` (3*minc+1, 3*minc) at n_r1              | :math:`b` (3*minc, 3*minc) at n_r1              |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 5         | :math:`b` (4*minc+1, 4*minc) at n_r1              | :math:`b` (4*minc, 4*minc) at n_r1              |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 |                                                                                                                 |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 6         | :math:`b` (minc+1, minc) at n_r2                  | :math:`b` (minc, minc) at n_r2                  |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 7         | :math:`b` (2*minc+1, 2*minc) at n_r2              | :math:`b` (2*minc, 2*minc) at n_r2              |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 8         | :math:`b` (3*minc+1, 3*minc) at n_r2              | :math:`b` (3*minc, 3*minc) at n_r2              |
 +-----------+---------------------------------------------------+-------------------------------------------------+
 | 9         | :math:`b` (4*minc+1, 4*minc) at n_r2              | :math:`b` (4*minc, 4*minc) at n_r2              |
 +-----------+---------------------------------------------------+-------------------------------------------------+

Analysis of these files can give you information about the drift frequency of the solution and it's symmetry.


.. _secinerFile:

``iner[P|T].TAG``
-----------------------

.. note:: These files are **only** written when :ref:`l_iner=.true. <varl_iner>` and :ref:`minc = 1 <varMinc>`.

These files contain time series of spherical harmonic coefficients upto degree,
:math:`\ell=6` at a radius :math:`r = (r_{cmb} - r_{icb})/2`. The ``inerP.TAG``
contains coefficients of the poloidal potential while the ``inerT.TAG``
contains coefficients of the toroidal potential.These files are written by 
the subroutine :f:subr:`write_rot <outrot/write_rot()>`. The oscillations of these
coefficients can be analysed to look for inertial modes. The
columns of the ``inerP.TAG`` look like follows:

  +--------------+------------------------+
  | No. of column| Coefficient            |
  +==============+========================+
  | 1            | :math:`w(\ell=1,m=1)`  |
  +--------------+------------------------+
  | 2            | :math:`w(\ell=2,m=1)`  |
  +--------------+------------------------+
  | 3            | :math:`w(\ell=2,m=2)`  |
  +--------------+------------------------+
  | 4            | :math:`w(\ell=3,m=1)`  |
  +--------------+------------------------+
  |                 ...                   |
  +--------------+------------------------+
  | 20           | :math:`w(\ell=6,m=5)`  |
  +--------------+------------------------+
  | 21           | :math:`w(\ell=6,m=6)`  |
  +--------------+------------------------+

where :math:`w(\ell,m)` is the poloidal potential with degree :math:`\ell` and order :math:`m`.

The columns of the ``inerT.TAG`` follow the following structure:

  +--------------+------------------------+
  | No. of column| Coefficient            |
  +==============+========================+
  | 1            | :math:`z(\ell=1,m=1)`  |
  +--------------+------------------------+
  | 2            | :math:`z(\ell=2,m=1)`  |
  +--------------+------------------------+
  | 3            | :math:`z(\ell=2,m=2)`  |
  +--------------+------------------------+
  | 4            | :math:`z(\ell=3,m=1)`  |
  +--------------+------------------------+
  |                 ...                   |
  +--------------+------------------------+
  | 20           | :math:`z(\ell=6,m=5)`  |
  +--------------+------------------------+
  | 21           | :math:`z(\ell=6,m=6)`  |
  +--------------+------------------------+

where :math:`z(\ell,m)` is the toroidal potential with degree :math:`\ell` and order :math:`m`.


.. _secSRFile:

``SR[IC|MA].TAG``
-------------------

.. note:: These files are **only** written for :ref:`nRotIc=-1 <varnRotIc>` (for ``SRIC.TAG``) or :ref:`nRotMa=-1 <varnRotMa>` (for ``SRMA.TAG``). In other words, these outputs are produced **only** when one of the boundaries is made to rotate at a prescribed rotation rate.

These files contain information about power due to torque from viscous and Lorentz forces at the inner core boundary (``SRIC.TAG``) or core mantle boundary (``SRMA.TAG``).The columns look like follows:

  +--------------+----------------------------------+
  | No. of column| Contents                         |
  +==============+==================================+
  | 1            | Time                             |
  +--------------+----------------------------------+
  | 2            | :math:`\Omega_{IC} | \Omega_{MA}`|
  +--------------+----------------------------------+
  | 3            | Total power = Lorentz + Viscous  |
  +--------------+----------------------------------+
  | 4            | Viscous power                    |
  +--------------+----------------------------------+
  | 5            | Lorentz force power              |
  +--------------+----------------------------------+

.. _secdtVrmsFile:

``dtVrms.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS force balance of the Navier Stokes equation. This file is
written by the subroutine :f:subr:`dtVrms <out_rms/dtvrms()>`.

   +---------------+--------------------------------------------------+
   | No. of column | Contents                                         |
   +===============+==================================================+
   | 1             | time                                             |
   +---------------+--------------------------------------------------+
   | 2             | Flow changes: inertia--advection                 |
   +---------------+--------------------------------------------------+
   | 3             | Coriolis force                                   |
   +---------------+--------------------------------------------------+
   | 4             | Lorentz force                                    |
   +---------------+--------------------------------------------------+
   | 5             | Advection term                                   |
   +---------------+--------------------------------------------------+
   | 6             | Diffusion term                                   |
   +---------------+--------------------------------------------------+
   | 7             | Buoyancy term                                    |
   +---------------+--------------------------------------------------+
   | 8             | Pressure gradient term                           |
   +---------------+--------------------------------------------------+
   | 9             | Sum of force terms: geostrophic balance          |
   +---------------+--------------------------------------------------+
   | 10            | Sum of force terms: magnetostrophic balance      |
   +---------------+--------------------------------------------------+
   | 11            | Sum of force terms: Archemidian balance          |
   +---------------+--------------------------------------------------+
   | 12            | Sum of force terms: Lorentz/Coriolis             |
   +---------------+--------------------------------------------------+
   | 13            | Sum of force terms: Pressure/Lorentz             |
   +---------------+--------------------------------------------------+
   | 14            | Sum of force terms: Coriolis/Inertia/Archimedean |
   +---------------+--------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the dtVrms.TAG files of the current directory
   >>> ts = MagicTs(field='dtVrms', all=True)


.. _secdtBrmsFile:

``dtBrms.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS terms that enter the induction equation. This file is
written by the subroutine :f:subr:`dtBrms <out_rms/dtbrms()>`.

   +---------------+-------------------------------------------------------+
   | No. of column | Contents                                              |
   +===============+=======================================================+
   | 1             | time                                                  |
   +---------------+-------------------------------------------------------+
   | 2             | Changes in magnetic field (poloidal)                  |
   +---------------+-------------------------------------------------------+
   | 3             | Changes in magnetic field (toroidal)                  |
   +---------------+-------------------------------------------------------+
   | 4             | Poloidal induction term                               |
   +---------------+-------------------------------------------------------+
   | 5             | Toroidal induction term                               |
   +---------------+-------------------------------------------------------+
   | 8             | Poloidal diffusion term                               |
   +---------------+-------------------------------------------------------+
   | 9             | Toroidal diffusion term                               |
   +---------------+-------------------------------------------------------+
   | 10            | Omega effect / toroidal induction term                |
   +---------------+-------------------------------------------------------+
   | 11            | Omega effect                                          |
   +---------------+-------------------------------------------------------+
   | 12            | Production of the dipole field                        |
   +---------------+-------------------------------------------------------+
   | 13            | Production of the axisymmetric dipole field           |
   +---------------+-------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the dtBrms.TAG files of the current directory
   >>> ts = MagicTs(field='dtBrms', all=True)


.. _secperpParFile:

``perpPar.TAG``
---------------

.. note:: This file is **only** written when :ref:`l_perpPar=.true. <varl_perpPar>`


This file contains several time series that decompose the kinetic energy into
components parallel and perpendicular to the rotation axis. This file is
calculated by the subroutine :f:subr:`outPerpPar <outpar_mod/outperppar()>`.

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | radial level                                                    |
   +---------------+-----------------------------------------------------------------+
   | 2             | Total kinetic energy perpendicular to the rotation axis:        |
   |               | :math:`\frac{1}{2}\langle u_s^2+u_\phi^2 \rangle_V`             |
   +---------------+-----------------------------------------------------------------+
   | 3             | Total kinetic energy parallel to the rotation axis:             |
   |               | :math:`\frac{1}{2}\langle u_z^2\rangle_V`                       |
   +---------------+-----------------------------------------------------------------+
   | 4             | Axisymmetric kinetic energy perpendicular to the rotation axis  |
   +---------------+-----------------------------------------------------------------+
   | 5             | Axisymmetric kinetic energy parallel to the rotation axis       |
   +---------------+-----------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the perpPar.TAG files of the current directory
   >>> ts = MagicTs(field='perpPar', all=True)
