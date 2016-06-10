Mantle and Inner Core Namelists
===============================

.. _secMantle:


Mantle Namelist
---------------

This namelist defines mantle properties 

* **conductance_ma** (default :f:var:`conductance_ma=0.0 <conductance_ma>`) is a real that defines the conductance (dimensionless) of the mantle.

.. _varnRotMa:

* **nRotMa** (default :f:var:`nRotMa=0 <nrotma>`) is an integer that defines the rotation of the mantle:

 +-----------+----------------------------------------------------------------------------------+
 | nRotMa=-1 | Mantle rotates with prescribed rate (see ``omega_ma1``  and ``omega_ma2`` below) |
 +-----------+----------------------------------------------------------------------------------+
 | nRotMa=0  | Fixed, non-rotating mantle                                                       |
 +-----------+----------------------------------------------------------------------------------+
 | nRotMa=1  | Mantle rotates according to torques                                              |
 +-----------+----------------------------------------------------------------------------------+

* **rho_ratio_ma** (default :f:var:`rho_ratio_ma=1 <rho_ratio_ma>`) is a real which gives the density of the mantle in terms of that of the outer core.

* **omega_ma1** (default :f:var:`omega_ma1=0.0 <omega_ma1>`) is a real which defines a mantle rotation rate (used when ``nRotMa=-1``).

* **omegaOsz_ma1** (default :f:var:`omegaOsz_ma1=0.0 <omegaosz_ma1>`) is a real which prescribes the oscillation frequency of the mantle rotation rate. In this case, ``omega_ma1`` is the amplitude of the oscillation.

* **tShift_ma1** (default :f:var:`tShift_ma1=0.0 <tshift_ma1>`) is a real which defines the time shift of the mantle rotation rate ``omega_ma1``.

* **omega_ma2** (default :f:var:`omega_ma2=0.0 <omega_ma2>`) is a real which defines a second mantle rotation rate.

* **omegaOsz_ma2** (default :f:var:`omegaOsz_ma2=0.0 <omegaosz_ma2>`) is a real which defines the oscillation frequency of the second mantle rotation rate ``omega_ma2``.

* **tShift_ma2** (default :f:var:`tShift_ma2=0.0 <tshift_ma2>`) is a real which defines the time shift for ``omega_ma2``.

The resultant prescribed mantle rotation rate is computed as:

.. code-block:: fortran

  omega_ma = omega_ma1*cos(omegaOsz_ma1*(time+tShift_ma1)) + &
             omega_ma2*cos(omegaOsz_ma2*(time+tShift_ma2))

The following defines the parameters when one wants to excite inertial modes in
the system artificially using a method similar to `Rieutord et. al 2012
<http://dx.doi.org/10.1103/PhysRevE.86.026304>`_ .

* **amp_RiMaSym** (default :f:var:`amp_RiMaSym=0.0 <amp_rimasym>`) is a real which defines the amplitude of forcing on the outer boundary for an equatorially symmetric mode

* **omega_RiMaSym** (default :f:var:`omega_RiMaSym=0.0 <omega_rimasym>`) is a real which defines the frequency of forcing on the outer boundary for an equatorially symmetric mode 

* **m_RiMaSym** (default :f:var:`m_RiMaSym=0.0 <m_rimasym>`) is an integer which defines the wavenumber of the equatorially symmetric mode one wants to excite

The following variables define the same for an equatorially anti-symmetric mode:

* **amp_RiMaAsym** (default :f:var:`amp_RiMaAsym=0.0 <amp_rimaasym>`)
* **omega_RiMaAsym** (default :f:var:`omega_RiMaAsym=0.0 <omega_rimaasym>`)
* **m_RiMaAsym** (default :f:var:`m_RiMaAsym=0.0 <m_rimaasym>`)

.. _secInnerCore:

Inner Core Namelist
-------------------

This namelist defines properties of the inner core

* **sigma_ratio** (default :f:var:`sigma_ratio=0.0 <sigma_ratio>`) is a real that defines the conductivity of the inner core with respect to the value of the outer core. ``sigma_ratio=0`` thus corresponds to a non-conducting inner core.

.. _varnRotIc:

* **nRotIc** (default :f:var:`nRotIc=0 <nrotic>`) is an integer that defines the rotation of the inner core. Behaves the same way as :f:var:`nRotMa <nrotma>` (above).

* **rho_ratio_ic** (default :f:var:`rho_ratio_ic=1.0 <rho_ratio_ic>`) is a real which defines the density of the inner core in terms of that of the outer core.

* **BIC** (default :f:var:`BIC=0.0 <bic>`) is a real which gives the imposed dipole field strength at the Inner Core Boundary. Having ``BIC > 0`` implies that the inner core acts as a dipole magnet - as implemented in the DTS experiment at Grenoble, France.

* **Variables prescribing rotation rate of inner core** The following variables are used to prescribe rotation rate of the inner core. They behave in the same way as the corresponding variables for the mantle. They are usd only when ``nRotIC=0``.

  - **omega_ic1** (default :f:var:`omega_ic1=0.0 <omega_ic>`)
  - **omegaOsz_ic1** (default :f:var:`omegaOsz_ic1=0.0 <omegaosz_ic1>`)
  - **tShift_ic1** (default :f:var:`tShift_ic1=0.0 <tshift_ic1>`)
  - **omega_ic2** (default :f:var:`omega_ic2=0.0 <omega_ic2>`)
  - **omegaOsz_ic2** (default :f:var:`omegaOsz_ic2=0.0 <omegaosz_ic2>`)
  - **tShift_ic2** (default :f:var:`tShift_ic2=0.0 <tshift_ic2>`)

As with the mantle, the resultant prescribed rotation rate for the inner core is computed as:

.. code-block:: fortran

  omega_ic = omega_ic1*cos(omegaOsz_ic1*(time+tShift_ic1)) + &
             omega_ic2*cos(omegaOsz_ic2*(time+tShift_ic2))

The following, as for the mantle namelist, is for artificially exciting inertial modes in the spherical shell, but for the inner boundary.

* **amp_RiIcSym** (default :f:var:`amp_RiIcSym=0.0 <amp_riicsym>`) is a real which defines the amplitude of forcing on the inner boundary for an equatorially symmetric mode

* **omega_RiIcSym** (default :f:var:`omega_RiIcSym=0.0 <omega_riicsym>`) is a real which defines the frequency of forcing on the inner boundary for an equatorially symmetric mode 

* **m_RiIcSym** (default :f:var:`m_RiIcSym=0.0 <m_riicsym>`) is an integer which defines the wavenumber of the equatorially symmetric mode one wants to excite

The following variables define the same for an equatorially anti-symmetric mode:

* **amp_RiIcAsym** (default :f:var:`amp_RiIcAsym=0.0 <amp_riicasym>`)
* **omega_RiIcAsym** (default :f:var:`omega_RiIcAsym=0.0 <omega_riicasym>`)
* **m_RiIcAsym** (default :f:var:`m_RiIcAsym=0.0 <m_riicasym>`)
