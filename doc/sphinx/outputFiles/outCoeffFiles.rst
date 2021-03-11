.. _secCoeffFiles:

Poloidal and toroidal potentials at given depths
================================================

These are fortran unformatted files which store time series of poloidal and
toroidal coefficients of different fields (magnetic field, velocity and
temperature) at specific depths.

In the following, :code:`time(j)`  is the time during the :math:`j^{th}` time
step, :code:`time(N)` being the last step. :code:`real` and :code:`imag` denote
real and imaginary parts, respectively, of spherical harmonic coefficients.
Also, the following notations will be used for the coefficients of potentials
(note that scalar fields like temperature do not have a poloidal/toroidal
decomposition):

    +----------------+-----------+------------+
    | Field          | Poloidal  | Toroidal   |
    +================+===========+============+
    | Magnetic       | :f:var:`b`| :f:var:`aj`|
    +----------------+-----------+------------+
    | Velocity       | :f:var:`w`| :f:var:`z` |
    +----------------+-----------+------------+
    | Temperature    |      :f:var:`s`        |
    +----------------+-----------+------------+


First and second derivatives are denoted with a differential notation. e.g:
:f:var:`dw` is the first derivative of :f:var:`w`, while :f:var:`ddb` is the second
derivative of :f:var:`b`.

.. _secCmbFile:

``B_coeff_cmb.TAG``
-------------------

.. note:: This file is **only** written when :ref:`l_cmb_field=.true. <varl_cmb_field>`

This file contains time series of spherical harmonic coefficients for the
poloidal potential of the magnetic field at the outer boundary (CMB) up to a
spherical harmonic degree given by :ref:`l_max_cmb <varl_max_cmb>`.
The detailed calculations are done in the subroutine :f:subr:`write_Bcmb
<out_coeff/write_bcmb()>`. The contents of the file look as follows:

 * **Header** The file header consists of the information: :ref:`l_max_cmb <varl_max_cmb>`, :ref:`minc <varMinc>` and the number of data points ``n_data``.
 * **Data** Each chunk of data after the header has the same pattern of ``time`` followed by a list of real and imaginary values of coefficients.

Thus, on a whole, the structure of the file looks like follows:

    .. code-block:: fortran

          !------------
          ! Line 1
          !------------

          l_max_cmb, minc, n_data

          !------------------------------
          ...

          !------------
          ! Line j + 1
          !------------

          time(j),
          real(b(l=1,m=0)), imag(b(l=1,m=0)),
          real(b(l=2,m=0)), imag(b(l=2,m=0)),
          ...
          real(b(l=l_max_cmb,m=l_max_cmb)), imag(b(l=l_max_cmb,m=l_max_cmb)),

          ...

          !-------------
          ! Line N + 1
          !-------------

          time(N),
          real(b(l=1,m=0)), imag(b(l=1,m=0)),
          real(b(l=2,m=0)), imag(b(l=2,m=0)),
          ...
          real(b(l=l_max_cmb,m=l_max_cmb)), imag(b(l=l_max_cmb,m=l_max_cmb))

This file can be read using :py:class:`MagicCoeffCmb <magic.coeff.MagicCoeffCmb>` with the following options:

   >>> # To stack the files B_cmb_coeff.testc to B_cmb_coeff.testf
   >>> cmb = MagicCoeffCmb(tag='test[c-f]')
   >>> # print Gauss coefficient for (\ell=10, m=3)
   >>> print(cmb.glm[:, cmb.idx[10, 3]])



.. _secCoeffrFiles:

Coefficients at desired radii
------------------------------

The following files **[B|V|T]_coeff_r#.TAG** save coefficients at specified
depths and are written by the subroutine :f:subr:`write_coeff_r
<out_coeff/write_coeff_r>`. See the section on :ref:`CMB and radial
coefficients <secOutNmlCoeff>` in the :ref:`ouput control namelist
<secOutputNml>` for details of specifying depth, using :ref:`n_r_step
<varn_r_step>` or :ref:`n_r_array <varn_r_array>` and desired maximum degree of
output :ref:`l_max_r <varl_max_r>`. A separate file for each desired radius
is written, numbered suitably as ``[B|V|T]_coeff_r1.TAG``,
``[B|V|T]_coeff_r2.TAG`` etc.


.. _secBcoeffrFile:

``B_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_field=.true. <varl_r_field>`.

This file contains output of time series of the spherical harmonic coefficients
of the poloidal and toroidal magnetic field potentials and the first and second
derivatives of the poloidal potential coefficients in the order :f:var:`b`,
:f:var:`db`, :f:var:`aj` and :f:var:`ddb`.  The output is for a specific
radius, :math:`r` up to degree :ref:`l_max_r <varl_max_cmb>`.

 * **Header** The file header consists of the information: :ref:`l_max_r
   <varl_max_r>`, :ref:`minc <varMinc>`,  the number of data points ``n_data``
   and the radius, :f:var:`r`.
 * **Data** Each chunk of data after the header contains the ``time`` at which
   the coefficients are stored, followed by the real and imaginary parts of:
   the poloidal coefficient ``b``, it's first derivative :f:var:`db`, the toroidal
   coefficient :f:var:`aj` and the second derivative of the poloidal coefficient
   :f:var:`ddb`.


The complete structure of the file looks like follows:

    .. code-block:: fortran

          !------------
          ! Line 1
          !------------

          l_max_r, minc, n_data, r

          !-------------------------------------------
          ...

          !------------
          ! Line j + 1
          !------------

          time(j),
          real(b(lm=1)), imag(b(lm=1)),
          real(b(lm=2)), imag(b(lm=2)),
          ...
          real(b(lm=lm_max)), imag(b(lm=lm_max)),
          real(db(lm=1)), imag(db(lm=1)),
          real(db(lm=2)), imag(db(lm=2)),
          ...
          real(db(lm=lm_max)), imag(db(lm=lm_max)),
          real(aj(lm=1)), imag(aj(lm=1)),
          real(aj(lm=2)), imag(aj(lm=2)),
          ...
          real(aj(lm=lm_max)), imag(aj(lm=lm_max)),
          real(ddb(lm=1)), imag(ddb(lm=1)),
          real(ddb(lm=1)), imag(ddb(lm=1)),
          ...
          real(ddb(lm=lm_max)), imag(ddb(lm=lm_max)),

          ...

          !------------
          ! Line N + 1
          !------------

          time(N),
          real(b(lm=1)), imag(b(lm=1)),
          real(b(lm=2)), imag(b(lm=2)),
          ...
          real(b(lm=lm_max)), imag(b(lm=lm_max)),
          real(db(lm=1)), imag(db(lm=1)),
          real(db(lm=2)), imag(db(lm=2)),
          ...
          real(db(lm=lm_max)), imag(db(lm=lm_max)),
          real(aj(lm=1)), imag(aj(lm=1)),
          real(aj(lm=2)), imag(aj(lm=2)),
          ...
          real(aj(lm=lm_max)), imag(aj(lm=lm_max)),
          real(ddb(lm=1)), imag(ddb(lm=1)),
          real(ddb(lm=1)), imag(ddb(lm=1)),
          ...
          real(ddb(lm=lm_max)), imag(ddb(lm=lm_max)),


This file can be read using :py:class:`MagicCoeffR <magic.coeff.MagicCoeffR>` with the following options:

   >>> # To stack the files B_coeff_r3.test* from the working directory
   >>> cr = MagicCoeffR(tag='test*', field='B', r=3)
   >>> # print the time and the poloidal potential for (\ell=3, m=3)
   >>> print(cr.time, cr.wlm[:, cr.idx[3, 3])



.. _secVcoeffrFile:

``V_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_field=.true. <varl_r_field>`

This file contains output of time series of the spherical harmonic coefficients
of the poloidal and toroidal velocity field potentials and the first
derivatives of the poloidal potential coefficients in the order :f:var:`w`,
:f:var:`dw`, and :f:var:`z`.  The output is for a specific radius, :math:`r` up
to degree :ref:`l_max_r <varl_max_cmb>`.

 * **Header** The file header consists of the information: :ref:`l_max_r
   <varl_max_r>`, :ref:`minc <varMinc>`,  the number of data points ``n_data``
   and the radius, :f:var:`r`.

 * **Data** Each chunk of data after the header contains the ``time`` at which
   the coefficients are stored, followed by the real and imaginary parts of:
   the poloidal coefficient :f:var:`w`, it's first derivative :f:var:`dw` and the
   toroidal coefficient :f:var:`z`.

The complete structure of the file looks like follows:

    .. code-block:: fortran

        !------------
        ! Line 1
        !------------

        l_max_r, minc, n_data, r

        !----------------------------------
        ...

        !------------
        ! Line j + 1
        !------------

        time(j),
        real(w(lm=1)), imag(w(lm=1)),
        real(w(lm=2)), imag(w(lm=2)),
        ...
        real(w(lm=lm_max)), imag(w(lm=lm_max)),
        real(dw(lm=1)), imag(dw(lm=1)),
        real(dw(lm=2)), imag(dw(lm=2)),
        ...
        real(dw(lm=lm_max)), imag(dw(lm=lm_max)),
        real(z(lm=1)), imag(z(lm=1)),
        real(z(lm=2)), imag(z(lm=2)),
        ...
        real(z(lm=lm_max)), imag(z(lm=lm_max)),

        ...

        !--------------
        ! Line N + 1
        !--------------

        time(N),
        real(w(lm=1)), imag(w(lm=1)),
        real(w(lm=2)), imag(w(lm=2)),
        ...
        real(w(lm=lm_max)), imag(w(lm=lm_max)),
        real(dw(lm=1)), imag(dw(lm=1)),
        real(dw(lm=2)), imag(dw(lm=2)),
        ...
        real(dw(lm=lm_max)), imag(dw(lm=lm_max)),
        real(z(lm=1)), imag(z(lm=1)),
        real(z(lm=2)), imag(z(lm=2)),
        ...
        real(z(lm=lm_max)), imag(z(lm=lm_max))

This file can be read using :py:class:`MagicCoeffR <magic.coeff.MagicCoeffR>` with the following options:

   >>> # To stack the files V_coeff_r3.test* from the working directory
   >>> cr = MagicCoeffR(tag='test*', field='V', r=3)
   >>> # print the poloidal and toroidal potentials for (\ell=6, m=0)
   >>> print(cr.wlm[:, cr.idx[6, 0]], cr.zlm[:, cr.idx[6, 0]])


.. _secTcoeffrFile:

``T_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_fieldT=.true. <varl_r_fieldT>`

This file contains output of time series of the spherical harmonic coefficients
of the temperature (or entropy) field. The output is for a specific radius,
:math:`r` up to degree :ref:`l_max_r <varl_max_cmb>`.

 * **Header** The file header consists of the information: :ref:`l_max_r
   <varl_max_r>`, :ref:`minc <varMinc>`,  the number of data points ``n_data``
   and the radius, :f:var:`r`.

 * **Data** Each chunk of data after the header contains the ``time`` at which
   the coefficients are stored, followed by the real and imaginary parts of the
   coefficient :f:var:`s`.

The complete structure of the file looks like follows:

    .. code-block:: fortran

        !------------
        ! Line 1
        !------------

        l_max_r, minc, n_data, r

        !---------------------------------

        ...

        !------------
        ! Line j + 1
        !------------

        time(j),
        real(s(lm=0)), imag(s(lm=0)),
        real(s(lm=1)), imag(s(lm=1)),
        real(s(lm=2)), imag(s(lm=2)),
        ...
        real(s(lm=lm_max)), imag(s(lm=lm_max)),

        !------------
        ! Line N + 1
        !------------

        time(N),
        real(s(lm=0)), imag(s(lm=0)),
        real(s(lm=1)), imag(s(lm=1)),
        real(s(lm=2)), imag(s(lm=2)),
        ...
        real(s(lm=lm_max)), imag(s(lm=lm_max)),

