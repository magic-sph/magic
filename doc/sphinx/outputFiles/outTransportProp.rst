Transport properties of the reference state
===========================================

These files define the radial transport properties of the reference state.
These arrays are calculated in the subroutines :f:subr:`radial
<radial_functions/radial()>` and :f:subr:`transportProperties
<radial_functions/transportproperties()>`.  The output files are written in the
subroutine :f:subr:`preCalc <precalculations/precalc()>`.


.. _secAnelFile:

``anel.TAG``
------------

.. note::
   This output is only calculated when an anelastic model is run, that is when
   :f:var:`l_anel=.true. <l_anel>` or :f:var:`l_anelastic_liquid=.true.
   <l_anelastic_liquid>`.


This file contains the radial profiles of the reference state (density,
temperature, gravity, etc.).

   +----------------+---------------------------------------------------------+
   | No. of column  |  Contents                                               |
   +================+=========================================================+
   | 1              | radial level: :math:`r`                                 |
   +----------------+---------------------------------------------------------+
   | 2              | temperature: :math:`\tilde{T}(r)`                       |
   +----------------+---------------------------------------------------------+
   | 3              | density: :math:`\tilde{\rho}(r)`                        |
   +----------------+---------------------------------------------------------+
   | 4              | radial derivative of the log of the density:            |
   |                | :math:`\beta={\rm d} \ln\tilde{\rho}/{\rm d} r`         |
   +----------------+---------------------------------------------------------+
   | 5              | radial derivative of :math:`\beta`:                     |
   |                | :math:`{\rm d} \beta/{\rm d} r`                         |
   +----------------+---------------------------------------------------------+
   | 6              | gravity: :math:`g(r)`                                   |
   +----------------+---------------------------------------------------------+
   | 7              | entropy gradient: :math:`{\rm d} s_0/{\rm d} r`         |
   +----------------+---------------------------------------------------------+
   | 8              | thermal diffusion operator:                             |
   |                | :math:`\nabla \cdot (K(r)\tilde{T}(r)\nabla s_0)`       |
   +----------------+---------------------------------------------------------+
   | 9              | inverse of the Gruneisen parameter :math`1/\Gamma`:     |
   |                | :math:`(\partial\ln\tilde{\rho}/\partial\ln\tilde{T})_S`|
   +----------------+---------------------------------------------------------+
   | 10             | radial derivative of the log of temperature:            |
   |                | :math:`\beta={\rm d} \ln\tilde{T}/{\rm d} r`            |
   +----------------+---------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='anel')
   >>> # print radius and density
   >>> print(rad.radius, rad.rho0)

.. _secVarCondFile:


``varCond.TAG``
---------------

.. note::
   This output is only calculated when the electrical conductivity varies with radius, i.e. when :ref:`nVarCond /= 0 <varnVarCond>`

This file contains the radial profiles of the electrical conductivity, the electrical
diffusivity and its radial derivative.

   +----------------+---------------------------------------------------------+
   | No. of column  |  Contents                                               |
   +================+=========================================================+
   | 1              | radial level: :math:`r`                                 |
   +----------------+---------------------------------------------------------+
   | 2              | electrical conductivity: :math:`\sigma(r)`              |
   +----------------+---------------------------------------------------------+
   | 3              | electrical diffusivity: :math:`\lambda(r)=1/\sigma(r)`  |
   +----------------+---------------------------------------------------------+
   | 4              | radial derivative of the electrical diffusivity:        |
   |                | :math:`{\rm d} \ln\lambda/{\rm d} r`                    |
   +----------------+---------------------------------------------------------+

This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='varCond')
   >>> print(rad.conduc) # Electrical conductivity

.. _secVarDiffFile:

``varDiff.TAG``
---------------

.. note::
   This output is only calculated when the thermal diffusivity varies with
   radius, i.e. when :ref:`nVarDiff /= 0 <varnVarDiff>`

This file contains the radial profiles of the thermal conductivity, the thermal
diffusivity and its radial derivative.

   +----------------+--------------------------------------------------------------+
   | No. of column  |  Contents                                                    |
   +================+==============================================================+
   | 1              | radial level: :math:`r`                                      |
   +----------------+--------------------------------------------------------------+
   | 2              | thermal conductivity: :math:`K(r)`                           |
   +----------------+--------------------------------------------------------------+
   | 3              | thermal diffusivity: :math:`\kappa(r)=K(r)/\tilde{\rho}(r)`  |
   +----------------+--------------------------------------------------------------+
   | 4              | radial derivative of the electrical diffusivity:             |
   |                | :math:`{\rm d} \ln\kappa/{\rm d} r`                          |
   +----------------+--------------------------------------------------------------+
   | 5              | Prandtl number: :math:`Pr(r)=\nu(r)/\kappa(r)`               |
   +----------------+--------------------------------------------------------------+


This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='varDiff')
   >>> print(rad.kappa) # Thermal diffusivity

.. _secVarViscFile:

``varVisc.TAG``
----------------

.. note::
   This output is only calculated when the kinematic viscosity varies with
   radius, i.e. when :ref:`nVarVisc /= 0 <varnVarVisc>`

This file contains the radial profiles of the dynamic viscosity, the kinematic
viscosity and its radial derivative.

   +----------------+--------------------------------------------------------------+
   | No. of column  |  Contents                                                    |
   +================+==============================================================+
   | 1              | radial level: :math:`r`                                      |
   +----------------+--------------------------------------------------------------+
   | 2              | dynamic viscosity: :math:`\mu(r)`                            |
   +----------------+--------------------------------------------------------------+
   | 3              | kinetmatic viscosity: :math:`\nu(r)=\mu(r)/\tilde{\rho}(r)`  |
   +----------------+--------------------------------------------------------------+
   | 4              | radial derivative of the kinematic viscosity:                |
   |                | :math:`{\rm d} \ln\nu/{\rm d} r`                             |
   +----------------+--------------------------------------------------------------+
   | 5              | Prandtl number: :math:`Pr(r)=\nu(r)/\kappa(r)`               |
   +----------------+--------------------------------------------------------------+
   | 6              | magnetic Prandtl number :math:`Pm(r)=\nu(r)/\lambda(r)`      |
   +----------------+--------------------------------------------------------------+


This file can be read using :py:class:`MagicRadial <magic.MagicRadial>` with the following options:

   >>> rad = MagicRadial(field='varVisc')
   >>> # print kinematic viscosity and Ekman
   >>> print(rad.kinVisc, rad.ekman)


.. _secMappingFile:

Nonlinear mapping of the Chebyshev grid
=======================================

``rNM.TAG``
-----------

.. note::
   This file is only written when :ref:`l_newmap=.true. <varl_newmap>`.

This file contains the profile of the radial mapping and its derivatives:


  +----------------+-------------------------------------------------------+
  | No. of column  | Contents                                              |
  +================+=======================================================+
  | 1              | Grid point index                                      |
  +----------------+-------------------------------------------------------+
  | 2              | Radius of a grid point                                |
  +----------------+-------------------------------------------------------+
  | 3              | First derivative of the mapping at a grid point       |
  +----------------+-------------------------------------------------------+
  | 4              | Second derivative of the mapping at a grid point      |
  +----------------+-------------------------------------------------------+
  | 5              | Third derivative of the mapping at a grid point       |
  +----------------+-------------------------------------------------------+

