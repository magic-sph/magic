.. _secOutNmlMisc:

RMS force balance
-----------------

  .. warning:: The RMS calculation is actually wrong in the current version. This needs again to be ported from MagIC 3.44. The RMS contributions to the induction equation are correct, though. A ticket has been opened on github regarding this issue: https://github.com/magic-sph/magic/issues/1

The code can compute the RMS of the force balance and the induction equation.

.. _varl_RMS:

* **l_RMS** (default ``l_RMS=.false.``) is a logical, which enables the calculation of RMS force balance, when set to ``.true.``.

* **l_RMStest** (default ``l_RMStest=.false.``) is a logical. This is a debug flag to check the consistency of the RMS calculation.


Additional possible diagnostics
-------------------------------

Geostrophy
++++++++++

.. _varl_par:

* **l_par** (default ``l_par=.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of the degree of geostrophy. These quantities are then stored in the columns 10-16 of the :ref:`misc.TAG <secMiscFile>` file.

Helicity
++++++++

.. _varl_hel:

* **l_hel** (default ``l_hel=.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of helicity (RMS, northern and southern hemisphere, etc.). The outputs are stored in the columns 6-9 of the :ref:`misc.TAG <secMiscFile>` file.

.. _varl_power:

Power budget
++++++++++++

* **l_power** (default ``l_power.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of input and output power (buoyancy, viscous and ohmic dissipations, torques). The time series are stored in ``power.TAG`` and the time-averaged radial profiles in :ref:`powerR.TAG <secPowerRfile>`.

.. _varl_fluxProfs:

Heat transport
++++++++++++++


* **l_fluxProfs** (default ``l_fluxProfs=.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of time-averaged radial heat flux profiles (conductive flux, convective flux, kinetic flux, viscous flux, Poynting flux and resistive flux). The time-averaged radial profiles are stored in the :ref:`fluxesR.TAG <secFluxesRfile>` file.

.. _varl_viscBcCalc:

Boundary layer analysis
+++++++++++++++++++++++

* **l_viscBcCalc** (default ``l_viscBcCalc=.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of time-averaged radial profiles that can be further use to determine the viscous and thermal boundary layer thicknesses: temperature, temperature variance, horizontal velocity, etc. The time-averaged radial profiles are stored in the :ref:`bLayersR.TAG <secBLayersRfile>` file.

.. _varl_perpPar:

Parallel/perpendicular decomposition
++++++++++++++++++++++++++++++++++++

* **l_perpPar** (default ``l_perpPar=.false.``) is a logical. When set to ``.true.``, this logical enables the decomposition of kinetic energy into components parallel and perpendicular to the rotation axis. The time series are stored in ``perpPar.TAG`` and the time-averaged radial profiles in :ref:`perpParR.TAG <secPerpParRfile>`.
