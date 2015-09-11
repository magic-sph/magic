Transport properties for the reference state
============================================


.. _secAnelFile:

``anel.TAG``
------------

.. note::
   This output is only calculated when an anelastic model is run, that is when ``l_anel=.true.`` or ``l_anelastic_liquid=.true.``.

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='anel')

.. _secVarCondFile:

``varCond.TAG``
---------------

.. note::
   This output is only calculated when the electrical conductivity varies with radius, i.e. when ``nVarCond /= 0`` (see :ref:`nVarCond <varnVarCond>`)

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='varCond')

.. _secVarDiffFile:

``varDiff.TAG``
---------------

.. note::
   This output is only calculated when the thermal diffusivity varies with radius, i.e. when ``nVarDiff /= 0`` (see :ref:`nVarDiff <varnVarDiff>`)

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='varDiff')

.. _secVarViscFile:

``varVisc.TAG``
----------------

.. note::
   This output is only calculated when the kinematic viscosity varies with radius, i.e. when ``nVarVisc /= 0`` (see :ref:`nVarVisc <varnVarVisc>`)

This file can be read using :py:class:`magic.MagicRadial` with the following options::
   >>> ts = MagicRadial(field='varVisc')

