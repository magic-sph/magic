Time-averaged radial profiles
=============================

.. _secEkinRFile:

``eKinR.TAG``
-------------

This file contains the time and horizontally averaged outer core kinetic energy along the radius. This file is calculated by the subroutine :f:subr:`get_e_kin <kinetic_energy/get_e_kin()>`.

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


This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='eKinR')

.. _secEmagRfile:

``eMagR.TAG``
-------------

This file contains the time and horizontally averaged outer core magnetic energy along the radius. This file is calculated by the subroutine :f:subr:`get_e_mag <magnetic_energy/get_e_mag()>`.

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


This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='eMagR')

.. _secParRfile:

``parR.TAG``
------------

This file contains several time and horizontally averaged flow properties (magnetic Reynolds number, Rossby number, etc.). This file is calculated by the subroutine :f:var:`outPar <outpar_mod/outpar()>`.

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
   | 5             | Local flow length-scale                                        |
   +---------------+----------------------------------------------------------------+
   | 6             | Local flow length-scale based on the non-axisymmetric flow     |
   |               | components                                                     |
   +---------------+----------------------------------------------------------------+
   | 7             | Local flow length-scale based on the peak of the poloidal      |
   |               | kinetic energy                                                 |
   +---------------+----------------------------------------------------------------+
   | 8             | Standard deviation of magnetic Reynolds number                 |
   +---------------+----------------------------------------------------------------+
   | 9             | Standard deviation of local Rossby number (mass-weighted)      |
   +---------------+----------------------------------------------------------------+
   | 10            | Standard deviation of local Rossby number (RMS velocity)       |
   +---------------+----------------------------------------------------------------+
   | 11            | Standard deviation of convective lengthscale                   |
   +---------------+----------------------------------------------------------------+
   | 12            | Standard deviation of convective lengthscale (non-axi)         |
   +---------------+----------------------------------------------------------------+
   | 13            | Standard deviation of convective lengthscale (pol. peak)       |
   +---------------+----------------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='parR')


.. _secHeatRfile:

``heatR.TAG``
-------------

.. note:: This file is **only** written when an equation for the heat transport (temperature or entropy) is solved.

This file contains several time and horizontally averaged thermodynamic properties (temperature, pressure, entropy, etc.) and their variance. This file is calculated by the subroutine :f:var:`outHeat <outmisc_mod/outheat()>`.

   +---------------+------------------------------------------------------------+
   | No. of column | Contents                                                   |
   +===============+============================================================+
   | 1             | Radial level                                               |
   +---------------+------------------------------------------------------------+
   | 2             | Entropy (spherically-symetric contribution)                |
   +---------------+------------------------------------------------------------+
   | 3             | Temperature (spherically-symetric contribution)            |
   +---------------+------------------------------------------------------------+
   | 4             | Pressure (spherically-symetric contribution)               |
   +---------------+------------------------------------------------------------+
   | 5             | Density (spherically-symetric contribution)                |
   +---------------+------------------------------------------------------------+
   | 6             | Chemical composition (spherically-symetric contribution)   |
   +---------------+------------------------------------------------------------+
   | 7             | Standard deviation of entropy                              |
   +---------------+------------------------------------------------------------+
   | 8             | Standard deviation of temperature                          |
   +---------------+------------------------------------------------------------+
   | 9             | Standard deviation of pressure                             |
   +---------------+------------------------------------------------------------+
   | 10            | Standard deviation of density                              |
   +---------------+------------------------------------------------------------+
   | 11            | Standard deviation of chemical composition                 |
   +---------------+------------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='heatR')

.. _secPowerRfile:

``powerR.TAG``
--------------

.. note:: This file is **only** written when :ref:`l_power=.true. <varl_power>`

This file contains the time and horizontally averaged power input (Buoyancy power) and outputs (viscous and Ohmic heating). This file is calculated by the subroutine :f:subr:`get_power <power/get_power()>`.

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | radial level                                                    |
   +---------------+-----------------------------------------------------------------+
   | 2             | Buoyancy power: :math:`Ra\,g(r)\,\langle u_r T'\rangle_s`       |
   +---------------+-----------------------------------------------------------------+
   | 3             | Chemical power: :math:`Ra_\xi\,g(r)\,\langle u_r \xi'\rangle_s` |
   +---------------+-----------------------------------------------------------------+
   | 4             | Viscous dissipation: :math:`\langle(\sigma)^2\rangle_s`         |
   +---------------+-----------------------------------------------------------------+
   | 5             | Ohmic dissipation: :math:`\langle(\nabla \times B)^2\rangle_s`  |
   +---------------+-----------------------------------------------------------------+
   | 6             | Standard deviation of buoyancy power                            |
   +---------------+-----------------------------------------------------------------+
   | 7             | Standard deviation of chemical power                            |
   +---------------+-----------------------------------------------------------------+
   | 8             | Standard deviation of viscous dissipation                       |
   +---------------+-----------------------------------------------------------------+
   | 9             | Standard deviation of ohmic dissipation                         |
   +---------------+-----------------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='powerR')

.. _secFluxesRfile:

``fluxesR.TAG``
---------------

.. note:: This file is **only** written when :ref:`l_fluxProfs=.true. <varl_fluxProfs>`

This file contains the time and horizontally averaged heat flux carried out by several physical processes: conductive flux, convective flux, kinetic flux, viscous flux, Poynting flux and resistive flux. This file is calculated by the subroutine :f:subr:`outPar <outpar_mod/outpar()>`.

   .. tabularcolumns:: |l|p{12cm}|

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | radial level                                                    |
   +---------------+-----------------------------------------------------------------+
   | 2             | conductive flux:                                                |
   |               |    .. math:: {\cal F}_{cond} = -\frac{1}{Pr}\kappa\tilde{\rho}  |
   |               |              \tilde{T}\frac{\partial \langle s \rangle_s}       |
   |               |              {\partial r}                                       |
   +---------------+-----------------------------------------------------------------+
   | 3             | convective flux:                                                |
   |               |    .. math:: {\cal F}_{conv}= \tilde{\rho}\tilde{T} \langle     |
   |               |              s\,u_r \rangle_s+\frac{Pr\,Di}{E\,Ra}\langle       |
   |               |              p\,u_r \rangle_s                                   |
   +---------------+-----------------------------------------------------------------+
   | 4             | kinetic flux:                                                   |
   |               |    .. math:: {\cal F}_{kin}= \frac{1}{2}\frac{Pr\,Di}{Ra}       |
   |               |              \langle u_r (\tilde{\rho}u^2) \rangle_s            |
   +---------------+-----------------------------------------------------------------+
   | 5             | viscous flux:                                                   |
   |               |    .. math:: {\cal F}_{visc}= -\frac{Pr\,Di}{Ra}                |
   |               |              \langle \vec{u}\cdot S \rangle_s                   |
   +---------------+-----------------------------------------------------------------+
   | 6             | Poynting flux:                                                  |
   |               |    .. math:: {\cal F}_{poyn}= -\frac{Pr\,Di}{Ra\,E\,Pm}         |
   |               |              \langle (\vec{u}\times\vec{B})\times\vec{B}        |
   |               |              \rangle_s                                          |
   +---------------+-----------------------------------------------------------------+
   | 7             | resistive flux:                                                 |
   |               |    .. math:: {\cal F}_{poyn}= \frac{Pr\,Di}{Ra\,E\,Pm^2}        |
   |               |              \langle (\vec{\nabla}\times\vec{B})\times\vec{B}   |
   |               |              \rangle_s                                          |
   +---------------+-----------------------------------------------------------------+
   | 8             | Standard deviation of conductive flux                           |
   +---------------+-----------------------------------------------------------------+
   | 9             | Standard deviation of convective flux                           |
   +---------------+-----------------------------------------------------------------+
   | 10            | Standard deviation of kinetic flux                              |
   +---------------+-----------------------------------------------------------------+
   | 11            | Standard deviation of viscous flux                              |
   +---------------+-----------------------------------------------------------------+
   | 12            | Standard deviation of Poynting flux                             |
   +---------------+-----------------------------------------------------------------+
   | 13            | Standard deviation of resistive flux                            |
   +---------------+-----------------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='fluxesR')

.. _secBLayersRfile:

``bLayersR.TAG``
----------------

.. note:: This file is **only** written when :ref:`l_viscBcCalc=.true. <varl_viscBcCalc>`

This file contains several time and horizontally averaged profiles that can be further used to determine thermal and viscous boundary layers: entropy (or temperature), entropy variance, horizontal velocity, radial derivative of the horizontal velocity, thermal dissipation rate. This file is calculated by the subroutine :f:subr:`outPar <outpar_mod/outpar()>`.

   .. tabularcolumns:: |l|p{12cm}|

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | radial level                                                    |
   +---------------+-----------------------------------------------------------------+
   | 2             | entropy: :math:`\langle s \rangle_s`                            |
   +---------------+-----------------------------------------------------------------+
   | 3             | horizontal velocity:                                            |
   |               |    .. math:: u_h=\left\langle\sqrt{u_\theta^2+u_\phi^2}         |
   |               |              \right\rangle_s                                    |
   +---------------+-----------------------------------------------------------------+
   | 4             | radial derivative of the horizontal velocity:                   |
   |               |    .. math:: \partial u_h/\partial r                            |
   +---------------+-----------------------------------------------------------------+
   | 5             | thermal dissipation rate:                                       |
   |               |    .. math:: \epsilon_T=\langle (\nabla T)^2 \rangle_s          |
   +---------------+-----------------------------------------------------------------+
   | 6             | Standard deviation of entropy                                   |
   +---------------+-----------------------------------------------------------------+
   | 7             | Standard deviation of horizontal velocity :math:`u_h`           |
   +---------------+-----------------------------------------------------------------+
   | 8             | Standard deviation of the radial derivative of :math:`u_h`      |
   +---------------+-----------------------------------------------------------------+
   | 9             | Standard deviation of the thermal dissipation rate              |
   +---------------+-----------------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='bLayersR')

Additional analyses of the boundary layers can then be carried out using :py:class:`BLayers <magic.bLayers.BLayers>`:

   >>> bl = BLayers(iplot=True)

.. _secPerpParRfile:

``perpParR.TAG``
----------------

.. note:: This file is **only** written when :ref:`l_perpPar=.true. <varl_perpPar>`

This file contains several time and horizontally averaged profiles that decompose the kinetic energy into components parallel and perpendicular to the rotation axis. This file is calculated by the subroutine :f:subr:`outPerpPar <outpar_mod/outperppar()>`.

   .. tabularcolumns:: |l|p{12cm}|

   +---------------+-----------------------------------------------------------------+
   | No. of column | Contents                                                        |
   +===============+=================================================================+
   | 1             | radial level                                                    |
   +---------------+-----------------------------------------------------------------+
   | 2             | Total kinetic energy perpendicular to the rotation axis:        |
   |               |    .. math:: \frac{1}{2}\langle u_s^2+u_\phi^2 \rangle_s        |
   +---------------+-----------------------------------------------------------------+
   | 3             | Total kinetic energy parallel to the rotation axis:             |
   |               |    .. math:: \frac{1}{2}\langle u_z^2\rangle_s                  |
   +---------------+-----------------------------------------------------------------+
   | 4             | Axisymmetric kinetic energy perpendicular to the rotation axis  |
   +---------------+-----------------------------------------------------------------+
   | 5             | Axisymmetric kinetic energy parallel to the rotation axis       |
   +---------------+-----------------------------------------------------------------+
   | 6             | Standard deviation of energy perpendicular to the rotation axis |
   +---------------+-----------------------------------------------------------------+
   | 7             | Standard deviation of energy parallel to the rotation axis      |
   +---------------+-----------------------------------------------------------------+
   | 8             | Standard deviation of axisymmetric energy perpendicular to the  |
   |               | rotation axis                                                   |
   +---------------+-----------------------------------------------------------------+
   | 9             | Standard deviation of axisymmetric energy parallel to the       |
   |               | rotation axis                                                   |
   +---------------+-----------------------------------------------------------------+


This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='perpParR')
