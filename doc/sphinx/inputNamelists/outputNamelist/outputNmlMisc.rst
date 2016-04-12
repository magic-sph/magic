.. _secOutNmlMisc:

RMS force balance
-----------------

The code can compute the RMS contributions of the different forces that
contribute to the Navier-Stokes equation and the the different terms that enter
the induction equation.

.. _varl_RMS:

* **l_RMS** (default :f:var:`l_RMS=.false. <l_rms>`) is a logical, which enables the calculation of RMS force balance, when set to ``.true.``. The outputs are stored in :ref:`dtVrms.TAG <secdtVrmsFile>` and :ref:`dtBrms.TAG <secdtBrmsFile>`.

* **rCut** (default :f:var:`rCut=0.0 <rcut>`) is a float. This is the thickness of the layer which is left out at both boundaries for the RMS calculation. ``rCut=0.075`` actually means that 7.5% below the CMB and above the ICB are disregarded in the force balance calculation.

* **rDea** (default  :f:var:`rDea=0.0 <rdea>`) is a float. This controls the dealiasing in RMS calculations. ``rDea=0.1`` means that the highest 10% of the Chebyshev modes are set to zero.


Additional possible diagnostics
-------------------------------

Geostrophy
++++++++++

.. _varl_par:

* **l_par** (default :f:var:`l_par=.false. <l_par>`) is a logical. When set to ``.true.``, this logical enables additional calculations (for instance the degree of geostrophy). The details of these calculations can be found in the subroutine ``getEgeos`` in the ``Egeos.f90`` file. These quantities are then stored in the columns 10-16 of the :ref:`geos.TAG <secGeosFile>` file.

* **l_corrMov** (default :f:var:`l_corrMov=.false. <l_corrmov>`) is a logical. When set to ``.true.``, this logical enables the calculation of a movie file that stores North/South correlation in the ``CVorz_mov.TAG`` file.

Helicity
++++++++

.. _varl_hel:

* **l_hel** (default :f:var:`l_hel=.false. <l_hel>`) is a logical. When set to ``.true.``, this logical enables the calculation of helicity (RMS, northern and southern hemisphere, etc.). The outputs are stored in the columns 6-9 of the :ref:`helicity.TAG <secHelicityFile>` file.

.. _varl_power:

Power budget
++++++++++++

* **l_power** (default :f:var:`l_power.false. <l_power>`) is a logical. When
  set to ``.true.``, this logical enables the calculation of input and output
  power (buoyancy, viscous and ohmic dissipations, torques). The time series
  are stored in :ref:`power.TAG <secPowerFile>` and :ref:`dtE.TAG <secdtEFile>` 
  and the time-averaged radial profiles in :ref:`powerR.TAG <secPowerRfile>`.

.. _varl_AM:

Angular momentum
++++++++++++++++

* **l_AM** (default :f:var:`l_AM=.false. <l_am>`) is a logical. When set to ``.true.``, this logical enables the calculation of angular momentum. The time series are stored in :ref:`AM.TAG <secAMFile>`.

.. _varl_drift:

Drift rates
+++++++++++

* **l_drift** (default :f:var:`l_drift=.false. <l_drift>`) is a logical. When set to ``.true.``, this logical enables the storage of some selected coefficients to allow the calculation of the drift rate. The time series are stored in :ref:`drift[V|B][DQ].TAG <secdriftFile>`.

.. _varl_iner:

Inertial modes
++++++++++++++

* **l_iner** (default :f:var:`l_iner=.false. <l_iner>`) is a logical. When set to ``.true.``, this logical enables the storage of some selected :math:`w(\ell, m)` at mid-shell (stored in :ref:`inerP.TAG <secinerFile>`) and :math:`z(\ell, m)` at mid-shell (stored in :ref:`inerT.TAG <secinerFile>`). Those files can be further used to identify inertial modes.

.. _varl_rMagSpec:

Radial spectra
++++++++++++++

* **l_rMagSpec** (default :f:var:`l_rMagSpec=.false <l_rmagspec>`) is a logical. When set to ``.true.``, the magnetic spectra for the first 6 spherical harmonic degree :math:`\ell` for all radii are stored at times of log ouputs. This produces the unformatted fortran files :ref:`rBrSpec.TAG <secrBspecFiles>` and :ref:`rBpSpec.TAG <secrBspecFiles>`.

* **l_DTrMagSpec** (default :f:var:`l_DTrMagSpec=.false <l_dtrmagspec>`) is a logical. When set to ``.true.``, the magnetic spectra of the magnetic field production terms for the first 6 spherical harmonic degree :math:`\ell` for all radii are stored at times of log ouputs. This produces the unformatted fortran files ``rBrProSpec.TAG``, ``rBrAdvSpec.TAG``, ``rBrDifSpec.TAG``, ``rBrDynSpec.TAG``, ``rBpProSpec.TAG``, ``rBpAdvSpec.TAG``, ``rBpDifSpec.TAG`` and ``rBpDynSpec.TAG``. All those files have exactly the same format as the :ref:`rBrSpec.TAG <secrBspecFiles>`.

.. _varl_fluxProfs:

Heat transport
++++++++++++++

* **l_fluxProfs** (default :f:var:`l_fluxProfs=.false. <l_fluxprofs>`) is a logical. When set to ``.true.``, this logical enables the calculation of time-averaged radial heat flux profiles (conductive flux, convective flux, kinetic flux, viscous flux, Poynting flux and resistive flux). The time-averaged radial profiles are stored in the :ref:`fluxesR.TAG <secFluxesRfile>` file.

.. _varl_viscBcCalc:

Boundary layer analysis
+++++++++++++++++++++++

* **l_viscBcCalc** (default :f:var:`l_viscBcCalc=.false. <l_viscbccalc>`) is a logical. When set to ``.true.``, this logical enables the calculation of time-averaged radial profiles that can be further use to determine the viscous and thermal boundary layer thicknesses: temperature, temperature variance, horizontal velocity, etc. The time-averaged radial profiles are stored in the :ref:`bLayersR.TAG <secBLayersRfile>` file.

.. _varl_perpPar:

Parallel/perpendicular decomposition
++++++++++++++++++++++++++++++++++++

* **l_perpPar** (default :f:var:`l_perpPar=.false. <l_perppar>`) is a logical. When set to ``.true.``, this logical enables the decomposition of kinetic energy into components parallel and perpendicular to the rotation axis. The time series are stored in :ref:`perpPar.TAG <secperpParFile>` and the time-averaged radial profiles in :ref:`perpParR.TAG <secPerpParRfile>`.

Potential vorticity
+++++++++++++++++++

* **l_PV** (default :f:var:`l_PV=.false. <l_pv>`) is a logical. When set to ``.true.``, this logical enables some potential vorticity diagnostics. At the end of the run, the results are stored in the the files ``PVZ.TAG`` and ``Vcy.TAG``.

Pressure
++++++++

* **l_PressGraph** (default :f:var:`l_PressGraph=.true. <l_pressgraph>`) is a logical. When set to ``.true.``, this logical enables the storage of pressure in the :ref:`graphic files <secGraphFile>`.

Time evolution of the m-spectra
+++++++++++++++++++++++++++++++

* **l_energy_modes** (default :f:var:`l_energy_modes=.false. <l_energy_modes>`) is a logical. When set to ``.true.``, this logical enables the storage of the time-evolution of the kinetic and magnetic energy spectra for a given range of spherical harmonic orders: :ref:`time spectra <secTimeSpectraFiles>`.

* **m_max_modes** (default :f:var:`m_max_modes=13 <m_max_modes>`) is an integer. This controls the maximum spherical harmonic order when :f:var:`l_energy_modes=.true. <l_energy_modes>`.
