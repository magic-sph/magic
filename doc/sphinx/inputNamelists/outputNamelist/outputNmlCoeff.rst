.. _secOutNmlCoeff:

Poloidal magnetic field potential at CMB
----------------------------------------

The **cmb** outputs controls the output of poloidal field potential coefficients at the CMB :math:`b_{\ell m}(r=r_o)`: :ref:`B_coeff_cmb.TAG <secCmbFile>` up to a maximum spherical harmonic degree ``l_max_cmb``.

.. note:: This calculation is **only** enabled when ``l_cmb_field=.true.`` or when
          ``l_dt_cmb_field=.true.``.

Specific inputs
+++++++++++++++

.. _varl_cmb_field:

* **l_cmb_field** (default :f:var:`l_cmb_field=.false. <l_cmb_field>`) is a logical. It needs to be turned on to get ``cmb`` files computed.

* **l_dt_cmb_field** (default :f:var:`l_dt_cmb_field=.false. <l_dt_cmb_field>`) is a logical. When set to ``.true.``, it allows the calculation of the secular variation of the magnetic field at the CMB.

.. _varl_max_cmb:

* **l_max_cmb** (default :f:var:`l_max_cmb=14 <l_max_cmb>`) is an integer. This is the maximum spherical harmonic degree :math:`\ell` stored in :ref:`B_coeff_cmb.TAG <secCmbFile>`, i.e. only :math:`\ell \leq \ell_{maxcmb}` are stored. For example, the following input parameter means that the :ref:`B_coeff_cmb.TAG <secCmbFile>` file is stored up to a spherical harmonic degree of :math:`\ell`:

  .. code-block:: fortran
  
     l_cmb_field = .true.,
     l_max_cmb   = 20,



Standard inputs
+++++++++++++++

* **n_cmb_step** (default :f:var:`n_cmb_step=0 <n_cmb_step>`) is an integer. This is the number of timesteps between two ``cmb`` outputs.

* **n_cmbs** (default :f:var:`n_cmbs=0 <n_cmbs>`) is an integer. This is the number of ``cmb`` outputs to be written.

* **t_cmb**  (default :f:var:`t_cmb=-1.0 -1.0 ... <t_cmb>`) is real array, which contains the times when ``cmb`` outputs are requested.

* **dt_cmb** (default :f:var:`dt_cmb=0.0 <dt_cmb>`) is a real, which defines the time interval between ``cmb`` outputs.

* **t_cmb_start** (default :f:var:`t_cmb_start=0.0 <t_cmb_start>`) is a real, which defines the time to start writing ``cmb`` outputs.

* **t_cmb_stop** (default :f:var:`t_cmb_stop=0.0 <t_cmb_stop>`) is a real, which defines the time to stop writing ``cmb`` outputs.


Poloidal and toroidal potentials at several depths
--------------------------------------------------

The ``coeff_r#`` outputs controls the output of the poloidal and
toroidal potential coefficients at several depths up to a maximum spherical
harmonic degree ``l_max_r``. The files :ref:`B_coeff_r#.TAG <secBcoeffrFile>`
and :ref:`V_coeff_r#.TAG <secVcoeffrFile>` are written when
``l_r_field=.true.``. The file :ref:`T_coeff_r#.TAG <secTcoeffrFile>` is
written when ``l_r_fieldT=.true.``.

.. note:: This calculation is **only** enabled when ``l_r_field=.true.`` or when
          ``l_r_fieldT=.true.``.

Specific inputs
+++++++++++++++

.. _varl_r_field:

* **l_r_field** (default :f:var:`l_r_field=.false. <l_r_field>`) is a logical. It needs to be turned on to get ``r_field`` files computed.

.. _varl_r_fieldT:

* **l_r_fieldT** (default :f:var:`l_r_fieldT=.false. <l_r_fieldt>`) is a logical. When set to ``.true.``, the thermal field is also stored in a file named :ref:`T_coeff_r*.TAG <secTcoeffrFile>`.

.. _varl_max_r:

* **l_max_r** (default :f:var:`l_max_r=l_max <l_max_r>`) is an integer. This is the maximum spherical harmonic degree :math:`\ell` stored in the ``r_field`` file, i.e. only :math:`\ell \leq \ell_{maxcmb}` are stored.

There are two ways to specify the radial grid points where you want to store
the ``[B|V|T]_coeff_r#.TAG`` files. You can specify a stepping ``n_r_step``: in that
case **5** ``coeff_r#.TAG`` files will be stored at 5 different radial levels every
``n_r_step`` grid point:

  .. code-block:: fortran
  
     l_r_field = .true.,
     n_r_step  = 6,
     l_max_r   = 30,

  This will produces 5 files that contain the poloidal and toroidal potentials
  up to spherical harmonic degree :math:`\ell=30`: 

     - ``[B|V|T]_coeff_r1.TAG`` corresponds to the radial grid point with the index ``nR=6``.
     - ``[B|V|T]_coeff_r2.TAG`` to ``nR=12``.
     - ``[B|V|T]_coeff_r3.TAG`` to ``nR=18``.
     - ``[B|V|T]_coeff_r4.TAG`` to ``nR=24``.
     - ``[B|V|T]_coeff_r5.TAG`` to ``nR=30``.

.. _varn_r_step:

* **n_r_step** (default :f:var:`n_r_step=2 <n_r_step>`) is an integer. This specifies the stepping between two consecutive ``[B|V|T]_coeff_r#.TAG`` files.

Alternatively, the input array ``n_r_array`` can be used to specify the radial grid points you exactly want to store:

  .. code-block:: fortran
  
     l_r_field = .true.,
     n_r_array = 8, 24, 47,
     l_max_r   = 10,

  This will produces 3 files that contain the poloidal and toroidal potentials
  up to spherical harmonic degree :math:`\ell=10`: 

     - ``[B|V|T]_coeff_r1.TAG`` corresponds to the radial grid point with the index ``nR=8``.
     - ``[B|V|T]_coeff_r2.TAG`` to ``nR=24``.
     - ``[B|V|T]_coeff_r3.TAG`` to ``nR=47``.

.. _varn_r_array:

* **n_r_array** (default :f:var:`n_r_array=0 0 0 ... <n_r_array>`) a an integer array. You can specify the radial grid points (starting from ``n_r_cmb=1``) where you want to store the coefficients.

Standard inputs
+++++++++++++++

* **n_r_field_step** (default :f:var:`n_r_field_step=0 <n_r_field_step>`) is an integer. This is the number of timesteps between two ``r_field`` outputs.

* **n_r_fields** (default :f:var:`n_r_fields=0 <n_r_fields>`) is an integer. This is the number of ``r_field`` outputs to be written.

* **t_r_field**  (default :f:var:`t_r_field=-1.0 -1.0 ... <t_r_field>`) is real array, which contains the times when ``r_field`` outputs are requested.

* **dt_r_field** (default :f:var:`dt_r_field=0.0 <dt_r_field>`) is a real, which defines the time interval between ``r_field`` outputs.

* **t_r_field_start** (default :f:var:`t_r_field_start=0.0 <t_r_field_start>`) is a real, which defines the time to start writing ``r_field`` outputs.

* **t_r_field_stop** (default :f:var:`t_r_field_stop=0.0 <t_r_field_stop>`) is a real, which defines the time to stop writing ``r_field`` outputs.

