.. _secOutNmlCoeff:

Poloidal field potential at CMB
-------------------------------

The **cmb** outputs controls the output of poloidal field potential at the CMB :math:`b_{\ell m}(r=r_o)`: ``B_coeff_cmb.TAG``.

.. note:: This calculation is **only** enabled when ``l_cmb_field=.true.``.

Specific inputs
+++++++++++++++

* **l_cmb_field** (default ``l_cmb_field=.false.``) is a logical. It needs to be turned on to get ``cmb`` files computed.

* **l_dt_cmb_field** (default ``l_dt_cmb_field=.false.``) is a logical. When set to ``.true.``, it allows the calculation of the secular variation of the magnetic field at the CMB.

* **l_max_cmb** (default ``l_max_cmb=14``) is an integer. This is the maximum spherical harmonic degree :math:`\ell` stored in the ``cmb`` file, i.e. only :math:`\ell \leq \ell_{maxcmb}` are stored.

Standard inputs
+++++++++++++++

* **n_cmb_step** (default ``n_cmb_step=0``) is an integer. This is the number of timesteps between two ``cmb`` outputs.

* **n_cmbs** (default ``n_cmbs=0``) is an integer. This is the number of ``cmb`` outputs to be written.

* **t_cmb**  (default  ``t_cmb=-1.0 -1.0 ...``) is real array, which contains the times when ``cmb`` outputs are requested.

* **dt_cmb** (default ``dt_cmb=0.0``) is a real, which defines the time interval between ``cmb`` outputs.

* **t_cmb_start** (default ``t_cmb_start=0.0``) is a real, which defines the time to start writing ``cmb`` outputs.

* **t_cmb_stop** (default ``t_cmb_stop=0.0``) is a real, which defines the time to stop writing ``cmb`` outputs.


Potential at several depths
---------------------------

The **coeff_r** outputs controls the output of the potential at several depths: ``B_coeff_r*.TAG``, ``V_coeff_r*.TAG`` and ``T_coeff_r*.TAG`` are produced.

.. note:: This calculation is **only** enabled when ``l_r_field=.true.``.

Specific inputs
+++++++++++++++

* **l_r_field** (default ``l_r_field=.false.``) is a logical. It needs to be turned on to get ``r_field`` files computed.

* **l_r_fieldT** (default ``l_r_fieldT=.false.``) is a logical. When set to ``.true.``, the thermal field is also stored in a file named ``T_coeff_r*.TAG``.

* **l_max_r** (default ``l_max_r=l_max``) is an integer. This is the maximum spherical harmonic degree :math:`\ell` stored in the ``r_field`` file, i.e. only :math:`\ell \leq \ell_{maxcmb}` are stored.

Standard inputs
+++++++++++++++

* **n_r_field_step** (default ``n_r_field_step=0``) is an integer. This is the number of timesteps between two ``r_field`` outputs.

* **n_r_fields** (default ``n_r_fields=0``) is an integer. This is the number of ``r_field`` outputs to be written.

* **t_r_field**  (default  ``t_r_field=-1.0 -1.0 ...``) is real array, which contains the times when ``r_field`` outputs are requested.

* **dt_r_field** (default ``dt_r_field=0.0``) is a real, which defines the time interval between ``r_field`` outputs.

* **t_r_field_start** (default ``t_r_field_start=0.0``) is a real, which defines the time to start writing ``r_field`` outputs.

* **t_r_field_stop** (default ``t_r_field_stop=0.0``) is a real, which defines the time to stop writing ``r_field`` outputs.

