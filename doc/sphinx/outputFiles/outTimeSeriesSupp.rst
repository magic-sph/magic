
Additional optional time-series outputs
=======================================

.. _secAMFile:

``AM.TAG``
-------------

.. note:: This file is **only** written when :ref:`l_AM=.true. <varl_AM>`

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the AM.TAG files of the current directory
   >>> ts = MagicTs(field='AM', all=True)


.. _secpowerFile:

``power.TAG``
-------------

.. note:: This file is **only** written when :ref:`l_power=.true. <varl_power>`

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the power.TAG files of the current directory
   >>> ts = MagicTs(field='power', all=True)


.. _secdtVrmsFile:

``dtVrms.TAG``
--------------

.. warning:: The RMS calculation is actually wrong in the current version. This 
             needs again to be ported from MagIC 3.44. This issue only affects 
             ``dtVrms.TAG``, though. A ticket has been opened on github regarding
	     this issue: https://github.com/magic-sph/magic/issues/1

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS force balance of the Navier Stokes equation. This file is
written by the subroutine ``dtVrms`` in the file ``outRMS.f90``.

   +---------------+--------------------------------------------------+
   | No. of column | Contents                                         |
   +===============+==================================================+
   | 1             | time                                             |
   +---------------+--------------------------------------------------+
   | 2             | Poloidal flow changes: inertia--advection        |
   +---------------+--------------------------------------------------+
   | 3             | Toroidal flow changes: inertia--advection        |
   +---------------+--------------------------------------------------+
   | 4             | Poloidal Coriolis force                          |
   +---------------+--------------------------------------------------+
   | 5             | Toroidal Coriolis force                          |
   +---------------+--------------------------------------------------+
   | 6             | Poloidal Lorentz force                           |
   +---------------+--------------------------------------------------+
   | 7             | Toroidal Lorentz force                           |
   +---------------+--------------------------------------------------+
   | 8             | Poloidal advection term                          |
   +---------------+--------------------------------------------------+
   | 9             | Toroidal advection term                          |
   +---------------+--------------------------------------------------+
   | 10            | Poloidal diffusion term                          |
   +---------------+--------------------------------------------------+
   | 11            | Toroidal diffusion term                          |
   +---------------+--------------------------------------------------+
   | 12            | Buoyancy term                                    |
   +---------------+--------------------------------------------------+
   | 13            | Pressure gradient term                           |
   +---------------+--------------------------------------------------+
   | 14            | Sum of force terms: geostrophic balance          |
   +---------------+--------------------------------------------------+
   | 15            | Sum of force terms: magnetostrophic balance      |
   +---------------+--------------------------------------------------+
   | 16            | Sum of force terms: Archemidian balance          |
   +---------------+--------------------------------------------------+


.. _secdtBrmsFile:

``dtBrms.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS terms that enter the induction equation. This file is
written by the subroutine ``dtBrms`` in the file ``out_RMS.f90``.

   +---------------+-------------------------------------------------------+
   | No. of column | Contents                                              |
   +===============+=======================================================+
   | 1             | time                                                  |
   +---------------+-------------------------------------------------------+
   | 2             | Changes in magnetic field (poloidal)                  |
   +---------------+-------------------------------------------------------+
   | 3             | Changes in magnetic field (toroidal)                  |
   +---------------+-------------------------------------------------------+
   | 4             | Poloidal strecthing term                              |
   +---------------+-------------------------------------------------------+
   | 5             | Toroidal strecthing term                              |
   +---------------+-------------------------------------------------------+
   | 6             | Poloidal field advection term                         |
   +---------------+-------------------------------------------------------+
   | 7             | Toroidal field advection term                         |
   +---------------+-------------------------------------------------------+
   | 8             | Poloidal diffusion term                               |
   +---------------+-------------------------------------------------------+
   | 9             | Toroidal diffusion term                               |
   +---------------+-------------------------------------------------------+
   | 10            | Omega effect / toroidal strecthing term               |
   +---------------+-------------------------------------------------------+
   | 11            | Omega effect                                          |
   +---------------+-------------------------------------------------------+
   | 12            | Poloidal field production (stretching+advection)      |
   +---------------+-------------------------------------------------------+
   | 13            | Toroidal field production (stretching+advection)      |
   +---------------+-------------------------------------------------------+


.. _secdtDrmsFile:

``dtDrms.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This files contains the RMS terms that enter the induction equation of the
dipole. This file is written by the subroutine ``dtBrms`` in the file
``out_RMS.f90``.

   +---------------+-------------------------------------------------------+
   | No. of column | Contents                                              |
   +===============+=======================================================+
   | 1             | time                                                  |
   +---------------+-------------------------------------------------------+
   | 2             | Dipole stretching                                     |
   +---------------+-------------------------------------------------------+
   | 3             | Dipole advection term                                 |
   +---------------+-------------------------------------------------------+
   | 4             | Dipole diffusion term                                 |
   +---------------+-------------------------------------------------------+

