
Additional optional time-series outputs
=======================================


.. _secHeatFile:

``heat.TAG``
------------

This files contains informations about the heat transfer (Nusselt number, entropy and
temperature at both boundaries). This file is written by the
subroutine :f:subr:`outHeat <outmisc_mod/outheat()>`.  

   +---------------+-------------------------------------------------------------+
   | No. of column | Contents                                                    |
   +===============+=============================================================+
   | 1             | time                                                        |
   +---------------+-------------------------------------------------------------+
   | 2             | Nusselt number at the inner boundary                        |
   +---------------+-------------------------------------------------------------+
   | 3             | Nusselt number at the outer boundary                        |
   +---------------+-------------------------------------------------------------+
   | 4             | Nusselt number based on :math:`\Delta T` ratio              |
   +---------------+-------------------------------------------------------------+
   | 5             | Temperature at the inner boundary                           |
   +---------------+-------------------------------------------------------------+
   | 6             | Temperature at the outer boundary                           |
   +---------------+-------------------------------------------------------------+
   | 7             | Entropy at the inner boundary                               |
   +---------------+-------------------------------------------------------------+
   | 8             | Entropy at the outer boundary                               |
   +---------------+-------------------------------------------------------------+
   | 9             | Heat flux at the inner boundary                             |
   +---------------+-------------------------------------------------------------+
   | 10            | Heat flux at the outer boundary                             |
   +---------------+-------------------------------------------------------------+
   | 11            | Pressure perturbation at the outer boundary                 |
   +---------------+-------------------------------------------------------------+
   | 12            | volume integrated mass perturbation                         |
   +---------------+-------------------------------------------------------------+
   | 13            | Sherwood number at the inner boundary                       |
   +---------------+-------------------------------------------------------------+
   | 14            | Sherwood number at the outer boundary                       |
   +---------------+-------------------------------------------------------------+
   | 15            | Sherwood number based on :math:`\Delta \xi` ratio           |
   +---------------+-------------------------------------------------------------+
   | 16            | Chemical composition at the inner boundary                  |
   +---------------+-------------------------------------------------------------+
   | 17            | Chemical composition at the outer boundary                  |
   +---------------+-------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the heat.TAG files of the current directory
   >>> ts = MagicTs(field='heat', all=True)

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
   | 3             | Chemical power: :math:`Ra_\xi\,g(r)\,\langle u_r \xi'\rangle_s`  |
   +---------------+------------------------------------------------------------------+
   | 4             | Viscous power at the inner boundary (ICB)                        |
   +---------------+------------------------------------------------------------------+
   | 5             | Viscous power at the outer boundary (CMB)                        |
   +---------------+------------------------------------------------------------------+
   | 6             | Viscous dissipation: :math:`\langle(\nabla \times u)^2\rangle_s` |
   +---------------+------------------------------------------------------------------+
   | 7             | Ohmic dissipation: :math:`\langle(\nabla \times B)^2\rangle_s`   |
   +---------------+------------------------------------------------------------------+
   | 8             | Total power at the CMB (viscous + Lorentz)                       |
   +---------------+------------------------------------------------------------------+
   | 9             | Total power at the ICB (viscous + Lorentz)                       |
   +---------------+------------------------------------------------------------------+
   | 10            | Total power                                                      |
   +---------------+------------------------------------------------------------------+
   | 11            | Time variation of total power                                    |
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


.. _secEarthLikeFile:

``earth_like.TAG``
------------------

This  contains informations about the Earth-likeness of the CMB radial magnetic
field. This file is written by the subroutine :f:subr:`get_e_mag <magnetic_energy/get_e_mag()>`.

.. note:: This file is **only** calculated when
          :ref:`l_earth_like=.true. <varl_earth_like>`.

..

   +---------------+--------------------------------------------------------------+
   | No. of column | Contents                                                     |
   +===============+==============================================================+
   | 1             | time                                                         |
   +---------------+--------------------------------------------------------------+
   | 2             | Ratio of axial dipole to non-dipole component at the CMB     |
   +---------------+--------------------------------------------------------------+
   | 3             | Equatorial symmetry of the CMB field (odd/even ratio)        |
   +---------------+--------------------------------------------------------------+
   | 4             | Zonality: zonal to non-zonal ratio of the CMB field          |
   +---------------+--------------------------------------------------------------+
   | 5             | Magnetic flux concentration at the CMB                       |
   +---------------+--------------------------------------------------------------+

The details of the calculations are given in (`Christensen et al., 2010 <http://dx.doi.org/10.1016/j.epsl.2010.06.009>`_).

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the earth_like.TAG files of the current directory
   >>> ts = MagicTs(field='earth_like', all=True)



.. _secGeosFile:

``geos.TAG``
------------

This file contains informations about the geostrophy of the flow.
This file is written by the subroutine :f:subr:`getEgeos <egeos_mod/getegeos()>`.  

.. note:: This file is **only** calculated when 
          :ref:`l_par=.true. <varl_par>`.

..

   +---------------+--------------------------------------------------------------+
   | No. of column | Contents                                                     |
   +===============+==============================================================+
   | 1             | time                                                         |
   +---------------+--------------------------------------------------------------+
   | 2             | Relative geostrophic kinetic energy                          |
   +---------------+--------------------------------------------------------------+
   | 3             | Relative kinetic energy in the northern part of the TC       |
   +---------------+--------------------------------------------------------------+
   | 4             | Relative kinetic energy in the southern part of the TC       |
   +---------------+--------------------------------------------------------------+
   | 5             | Kinetic energy (calculated on the cylindrical grid)          |
   +---------------+--------------------------------------------------------------+
   | 6             | North/South correlation of Vz, outside the TC                |
   +---------------+--------------------------------------------------------------+
   | 7             | North/South correlation of vorticity outside the TC          |
   +---------------+--------------------------------------------------------------+
   | 8             | North/South correlation of helicity outside the TC           |
   +---------------+--------------------------------------------------------------+
   | 9             | Geostrophy of axisymmetic flow                               |
   +---------------+--------------------------------------------------------------+
   | 10            | Geostrophy of zonal flow                                     |
   +---------------+--------------------------------------------------------------+
   | 11            | Geostrophy of meridional flow                                |
   +---------------+--------------------------------------------------------------+
   | 12            | Geostrophy of non-axisymmetric flow                          |
   +---------------+--------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the geos.TAG files of the current directory
   >>> ts = MagicTs(field='geos', all=True)

.. _secHelicityFile:

``helicity.TAG``
----------------

This files contains informations about the kinetic helicity in both the 
Northern and the Southern hemispheres.  This file is written by the
subroutine :f:subr:`outHelicity <outmisc_mod/outhelicity()>`.  

.. note:: This file is **only** calculated when :ref:`l_hel=.true. <varl_hel>`. 

..

   +---------------+-------------------------------------------------------------+
   | No. of column | Contents                                                    |
   +===============+=============================================================+
   | 1             | time                                                        |
   +---------------+-------------------------------------------------------------+
   | 2             | Helicity (northern hemisphere)                              |
   +---------------+-------------------------------------------------------------+
   | 3             | Helicity (southern hemisphere)                              |
   +---------------+-------------------------------------------------------------+
   | 4             | RMS helicity (northern hemisphere)                          |
   +---------------+-------------------------------------------------------------+
   | 5             | RMS helicity (southern hemisphere)                          |
   +---------------+-------------------------------------------------------------+
   | 6             | Helicity (northern hemisphere, only non-axisym. flow)       |
   +---------------+-------------------------------------------------------------+
   | 6             | Helicity (southern hemisphere, only non-axisym. flow)       |
   +---------------+-------------------------------------------------------------+
   | 8             | RMS helicity (northern hemisphere, only non-axisym. flow)   |
   +---------------+-------------------------------------------------------------+
   | 9             | RMS helicity (southern hemisphere, only non-axisym. flow)   |
   +---------------+-------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the helicity.TAG files of the current directory
   >>> ts = MagicTs(field='helicity', all=True)

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

   +---------------+--------------------------------------------------------------+
   | No. of column | Contents                                                     |
   +===============+==============================================================+
   | 1             | Time                                                         |
   +---------------+--------------------------------------------------------------+
   | 2             | Total inertia: dU/dt and advection                           |
   +---------------+--------------------------------------------------------------+
   | 3             | Coriolis force                                               |
   +---------------+--------------------------------------------------------------+
   | 4             | Lorentz force                                                |
   +---------------+--------------------------------------------------------------+
   | 5             | Advection term                                               |
   +---------------+--------------------------------------------------------------+
   | 6             | Diffusion term                                               |
   +---------------+--------------------------------------------------------------+
   | 7             | Thermal buoyancy term                                        |
   +---------------+--------------------------------------------------------------+
   | 8             | Chemical buoyancy term                                       |
   +---------------+--------------------------------------------------------------+
   | 9             | Pressure gradient term                                       |
   +---------------+--------------------------------------------------------------+
   | 10            | Sum of force terms: geostrophic balance                      |
   +---------------+--------------------------------------------------------------+
   | 11            | Sum of force terms: pressure, Coriolis and Lorentz           |
   +---------------+--------------------------------------------------------------+
   | 12            | Sum of force terms: pressure, buoyancy and Coriolis          |
   +---------------+--------------------------------------------------------------+
   | 13            | Sum of force terms: pressure, buoyancy, Coriolis and Lorentz |
   +---------------+--------------------------------------------------------------+
   | 14            | Sum of force terms: Lorentz/Coriolis                         |
   +---------------+--------------------------------------------------------------+
   | 15            | Sum of force terms: Pressure/Lorentz                         |
   +---------------+--------------------------------------------------------------+
   | 16            | Sum of force terms: Coriolis/Inertia/Archimedean             |
   +---------------+--------------------------------------------------------------+

This file can be read using :py:class:`MagicTs <magic.MagicTs>` with the following options:

   >>> # To stack all the dtVrms.TAG files of the current directory
   >>> ts = MagicTs(field='dtVrms', all=True)


.. _secdtBrmsFile:

``dtBrms.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS terms that enter the induction equation. This file is
written by the subroutine :f:subr:`dtBrms <rms/dtbrms()>`.

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
   | 1             | time                                                            |
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

.. _secphaseFile:

``phase.TAG``
-------------

This file contains several diagnostic related to phase field whenever this field is used
by MagIC. This file is calculated by the subroutine
:f:subr:`outPhase <outmisc_mod/outphase()>`.

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | time                                                            |
   +---------------+-----------------------------------------------------------------+
   | 2             | Average radius of the solidus                                   |
   +---------------+-----------------------------------------------------------------+
   | 3             | Average temperature at the solidus (should be close to tmelt)   |
   +---------------+-----------------------------------------------------------------+
   | 4             | Mean spherically-symmetric radius of the solidus                |
   +---------------+-----------------------------------------------------------------+
   | 5             | Mean spherically-symmetric temperature at the mean spherically  |
   |               | symmetric radius of the solidus                                 |
   +---------------+-----------------------------------------------------------------+
   | 6             | Minimum radius of the solidus                                   |
   +---------------+-----------------------------------------------------------------+
   | 7             | Maximum radius of the solidus                                   |
   +---------------+-----------------------------------------------------------------+
   | 8             | Volume of the solid phase                                       |
   +---------------+-----------------------------------------------------------------+
   | 9             | Kinetic energy of the solid phase (should be small)             |
   +---------------+-----------------------------------------------------------------+
   | 10            | Kinetic energy of the liquid phase                              |
   +---------------+-----------------------------------------------------------------+
   | 11            | Heat flux at the outer core boundary                            |
   +---------------+-----------------------------------------------------------------+
   | 12            | Heat flux at the inner core boundary                            |
   +---------------+-----------------------------------------------------------------+
   | 13            | Time variation of of temperature and phase field:               |
   |               | :math:`\frac{\partial}{\partial t}(T-St\Phi)`                   |
   +---------------+-----------------------------------------------------------------+
   | 14            | Maximum value of phase field (should not exceed one by much,    |
   |               | otherwise, Gibbs phenomenon is likely occurring)                |
   +---------------+-----------------------------------------------------------------+
   | 15            | Minimum value of phase field (should be close to zero,          |
   |               | otherwise, Gibbs phenomenon is likely occurring)                |
   +---------------+-----------------------------------------------------------------+

   >>> # To stack all the phase.TAG files of the current directory
   >>> ts = MagicTs(field='phase', all=True)

.. _secHemiFile:

``hemi.TAG``
-------------

This file contains diagnostics related to North/South hemisphericity in kinetic and
magnetic energies. This is based on Dietrich and Wicht (2013) work. The file is
calculated by the subroutine :f:subr:`outHemi <outmisc_mod/outhemi()>`.

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | time                                                            |
   +---------------+-----------------------------------------------------------------+
   | 2             | relative hemisphericity of :math:`|u_r|`                        |
   +---------------+-----------------------------------------------------------------+
   | 3             | relative hemisphericity of kinetic energy                       |
   +---------------+-----------------------------------------------------------------+
   | 4             | relative hemisphericity of :math:`|B_r|`                        |
   +---------------+-----------------------------------------------------------------+
   | 5             | relative hemisphericity of magnetic energy                      |
   +---------------+-----------------------------------------------------------------+
   | 6             | relative hemisphericity of :math:`|B_r|` at the CMB             |
   +---------------+-----------------------------------------------------------------+
   | 7             | total kinetic energy (to assess method accuracy)                |
   +---------------+-----------------------------------------------------------------+
   | 8             | total magnetic energy (to assess method accuracy)               |
   +---------------+-----------------------------------------------------------------+

   >>> # To stack all the hemi.TAG files of the current directory
   >>> ts = MagicTs(field='hemi', all=True)


``growth_sym.TAG`` and ``growth_asym.TAG``
------------------------------------------

Those files contain the time series of growth rate of different azimuthal wavenumbers ranging from :f:var:`m_min <m_min>` to :f:var:`m_max <_m_max>`. This file is produced when MagIC is used to compute the onset of convection, i.e. when :f:var:`mode=5 <mode>`. `growth_sym` corresponds to equatorially-symmetric mode, `growth_asym` to equatorially-asymmetric modes. Those files are produced by the routine :f:subr:`get_onset <outmisc_mod/get_onset()>`.

   +---------------+---------------------------------------------------------+
   | No. of column | Contents                                                |
   +===============+=========================================================+
   | 1             | time                                                    |
   +---------------+---------------------------------------------------------+
   | 2             | growth rate of the azimuthal wave number `m_min`        |
   +---------------+---------------------------------------------------------+
   | 3             | growth rate of the azimuthal wave number `m_min+1`      |
   +---------------+---------------------------------------------------------+
   | 4             | growth rate of the azimuthal wave number `m_min+2`      |
   +---------------+---------------------------------------------------------+
   | ...           | growth rate of the azimuthal wave number `m_max`        |
   +---------------+---------------------------------------------------------+

``drift_sym.TAG`` and ``drift_asym.TAG``
----------------------------------------

Those files contain the time series of drift frequency of different azimuthal wavenumbers ranging from :f:var:`m_min <m_min>` to :f:var:`m_max <_m_max>`. This file is produced when MagIC is used to compute the onset of convection, i.e. when :f:var:`mode=5 <mode>`. `drift_sym` corresponds to equatorially-symmetric modes, `drift_asym` to equatorially-asymmetric modes. Those files are produced by the routine :f:subr:`get_onset <outmisc_mod/get_onset()>`.

   +---------------+-------------------------------------------------------------+
   | No. of column | Contents                                                    |
   +===============+=============================================================+
   | 1             | time                                                        |
   +---------------+-------------------------------------------------------------+
   | 2             | drift frequency of the azimuthal wave number `m_min`        |
   +---------------+-------------------------------------------------------------+
   | 3             | drift frequency of the azimuthal wave number `m_min+1`      |
   +---------------+-------------------------------------------------------------+
   | 4             | drift frequency of the azimuthal wave number `m_min+2`      |
   +---------------+-------------------------------------------------------------+
   | ...           | drift frequency of the azimuthal wave number `m_max`        |
   +---------------+-------------------------------------------------------------+
