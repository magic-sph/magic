
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


.. _secu_squareFile:

``u_square.TAG``
----------------

.. note:: This file is **only** written in anelastic models, i.e. either when
          :ref:`strat/=0 <varstrat>` or when :ref:`interior_model/="None" <varinterior_model>`

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the u_square.TAG files of the current directory
   >>> ts = MagicTs(field='u_square', all=True)

.. _secdriftFile:

``drift[V|B][D|Q].TAG``
-----------------------

.. note:: These files are **only** written when :ref:`l_drift=.true. <varl_drift>`

.. _secinerFile:

``iner[P|T].TAG``
-----------------------

.. note:: These files are **only** written when :ref:`l_iner=.true. <varl_iner>` and :ref:`minc = 1 <varMinc>`.

These files contain time series of spherical harmonic coefficients upto degree, :math:`l=6` at a radius :math:`r = (r_{cmb} - r_{icb})/2`. The ``inerP.TAG`` contains coefficients of the poloidal potential while the ``inerT.TAG`` contains coefficients of the toroidal potential. The oscillations of these coefficients can be analysed to look for inertial modes. As an example, the columns of the ``inerP.TAG`` look like follows:

  +--------------+-------------+
  | No. of column| Coefficient |
  +==============+=============+
  | 1            | w(1,1)      |
  +--------------+-------------+
  | 2            | w(2,1)      |
  +--------------+-------------+
  | 3            | w(2,2)      |
  +--------------+-------------+
  | 4            | w(3,1)      |
  +--------------+-------------+
  |             ...            |
  +--------------+-------------+
  | 20           | w(6,5)      |
  +--------------+-------------+
  | 21           | w(6,6)      |
  +--------------+-------------+

where ``w(l,m)`` is the spherical harmonic coefficient with degree :math:`l` and order :math:`m`, of the poloidal potential.

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

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the dtVrms.TAG files of the current directory
   >>> ts = MagicTs(field='dtVrms', all=True)


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

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the dtBrms.TAG files of the current directory
   >>> ts = MagicTs(field='dtBrms', all=True)


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


.. _secperpParFile:

``perpPar.TAG``
---------------

.. note:: This file is **only** written when :ref:`l_perpPar=.true. <varl_perpPar>`

This file can be read using :py:class:`magic.MagicTs` with the following options::
   >>> # To stack all the perpPar.TAG files of the current directory
   >>> ts = MagicTs(field='perpPar', all=True)


.. _secrBspecFiles:

``rB[r|p]Spec.TAG``
-------------------

.. note:: This file is **only** written when :ref:`l_rMagSpec=.true. <varl_rMagSpec>`

The calculations for done in the ``radial_spectra.f90`` files.
