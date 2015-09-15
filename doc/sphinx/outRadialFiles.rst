Time and horizontally averaged files
====================================

.. _secEkinRFile:

``eKinR.TAG``
-------------

This file contains the time and horizontally averaged outer core kinetic energy along the radius. This file is calculated by the subroutine ``get_e_kin`` in ``kinetic_energy.f90``.

   +---------------+----------------------------------------------------------------+
   | No. of column | Contents                                                       |
   +===============+================================================================+
   | 1             | radial level                                                   |
   +---------------+----------------------------------------------------------------+
   | 2             | time and horizontally averaged poloidal energy                 |
   +---------------+----------------------------------------------------------------+
   | 3             | time and horizontally averaged axisymmetric poloidal energy    |
   +---------------+----------------------------------------------------------------+
   | 4             | time and horizontally averaged toroidal energy                 |
   +---------------+----------------------------------------------------------------+
   | 5             | time and horizontally averaged axisymmetric toroidal energy    |
   +---------------+----------------------------------------------------------------+
   | 6             | time and horizontally averaged poloidal energy,                |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 7             | time and horizontally averaged axisymmetric poloidal energy,   |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 8             | time and horizontally averaged toroidal energy,                |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 9             | time and horizontally averaged axisymmetric toroidal energy,   |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+


This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='eKinR')

.. _secEmagRfile:

``eMagR.TAG``
-------------

This file contains the time and horizontally averaged outer core magnetic energy along the radius. This file is calculated by the subroutine ``get_e_mag`` in ``magnetic_energy.f90``.

   +---------------+----------------------------------------------------------------+
   | No. of column | Contents                                                       |
   +===============+================================================================+
   | 1             | radial level                                                   |
   +---------------+----------------------------------------------------------------+
   | 2             | time and horizontally averaged poloidal energy                 |
   +---------------+----------------------------------------------------------------+
   | 3             | time and horizontally averaged axisymmetric poloidal energy    |
   +---------------+----------------------------------------------------------------+
   | 4             | time and horizontally averaged toroidal energy                 |
   +---------------+----------------------------------------------------------------+
   | 5             | time and horizontally averaged axisymmetric toroidal energy    |
   +---------------+----------------------------------------------------------------+
   | 6             | time and horizontally averaged poloidal energy,                |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 7             | time and horizontally averaged axisymmetric poloidal energy,   |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 8             | time and horizontally averaged toroidal energy,                |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 9             | time and horizontally averaged axisymmetric toroidal energy,   |
   |               | normalized by surface area at this radial level                |
   +---------------+----------------------------------------------------------------+
   | 10            | ratio between time-averaged dipole energy and                  |
   |               | time-averaged total energy                                     |
   +---------------+----------------------------------------------------------------+


This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='eMagR')

.. _secParRfile:

``parR.TAG``
------------

This file contains several time and horizontally averaged flow properties (magnetic Reynolds number, Rossby number, etc.). This file is calculated by the subroutine ``outPar`` in ``outPar.f90``.

   +---------------+----------------------------------------------------------------+
   | No. of column | Contents                                                       |
   +===============+================================================================+
   | 1             | radial level                                                   |
   +---------------+----------------------------------------------------------------+
   | 2             | Magnetic Reynolds number                                       |
   +---------------+----------------------------------------------------------------+
   | 3             | Local Rossby number (based on the mass-weighted velocity)      |
   +---------------+----------------------------------------------------------------+
   | 4             | Local Rossby number (based on the RMS velocity)                |
   +---------------+----------------------------------------------------------------+
   | 5             | Local flow length-scale (based on the mass-weighted velocity)  |
   +---------------+----------------------------------------------------------------+
   | 6             | Local flow length-scale based on the non-axisymmetric flow     |
   |               | components (based on the mass-weighted velocity)               |
   +---------------+----------------------------------------------------------------+
   | 7             | Local flow length-scale (based on the RMS velocity)            |
   +---------------+----------------------------------------------------------------+
   | 8             | Local flow length-scale based on the non-axisymmetric flow     |
   |               | components (based on the RMS velocity)                         |
   +---------------+----------------------------------------------------------------+

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='parR')


.. _secPowerRfile:

``powerR.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_power=.true. <varl_power>`

This file contains the time and horizontally averaged power input (Buoyancy power) and outputs (viscous and Ohmic heating). This file is calculated by the subroutine ``get_power`` in ``power.f90``.

   +---------------+--------------------------------------------------------------------+
   | No. of column | Contents                                                           |
   +===============+====================================================================+
   | 1             | radial level                                                       |
   +---------------+--------------------------------------------------------------------+
   | 2             | Buoyancy power: :math:`Ra\,g(r)\,\langle u_r T'\rangle_\phi`       |
   +---------------+--------------------------------------------------------------------+
   | 3             | Viscous dissipation: :math:`\langle(\nabla \times u)^2\rangle_\phi`|
   +---------------+--------------------------------------------------------------------+
   | 4             | Ohmic dissipation: :math:`\langle(\nabla \times B)^2\rangle_\phi`  |
   +---------------+--------------------------------------------------------------------+

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='powerR')
