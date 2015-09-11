Time and horizontally averaged files
====================================

.. _secEkinRFile:

``eKinR.TAG``
-------------

This file contains the time and horizontally averaged outer core kinetic energy along the radius. This file is calulated by the subroutine ``get_e_kin`` in ``kinetic_energy.f90``.

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
