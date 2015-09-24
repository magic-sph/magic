.. _secCoeffFiles:

Poloidal and toroidal potentials at given depths
================================================

These are fortran unformatted files which store time series of poloidal and toroidal coefficients of different fields (magnetic field, velocity and temeperature) at specific depths. Unformatted files are not directly human readable, and are used to store binary data and move it around without changing the internal representation. In fortran, the open, read and write operations for these files are performed as follows:

.. code-block:: fortran

  open(unit=4, file='test', form='unformatted')
  read(unit=4) readVar
  write(unit=n_out, iostat=ios) writeVar !Unformatted write

.. _secCmbFile:

``B_coeff_cmb.TAG``
-------------------

.. note:: This file is **only** written when :ref:`l_cmb_field=.true. <varl_cmb_field>` 

This file contains time series of spherical harmonic coefficients for the
poloidal part of the magnetic field at the outer boundary (CMB) upto a degree
given by :ref:`l_max_cmb <varl_max_cmb>`. The contents of the file look as
follows:

 * **Header** The file header consists of the information: :ref:`l_max_cmb <varl_max_cmb>`, :ref:`minc <varMinc>` and the number of data points - ``n_data``.
 * **Data** Each chunk of data after the header has the same pattern of ``time`` followed by a list of real and imaginary values of coefficients.

Thus, on a whole, the structure of the file looks like follows:

   .. code-block:: fortran
   
   	    !------------
   	    !-- Line 1
   	    !------------
            l_max_cmb, minc, n_data
            ...
   	    !------------
   	    !-- Line n
   	    !------------
            time(n), 
            real(w(l=1,m=0)),imag(w(l=1,m=0)),                  
            real(w(l=2,m=0)),imag(w(l=2,m=0)),                  
            ...
            real(w(l=l_max_cmb,m=l_max_cmb)),imag(w(l=l_max_cmb,m=l_max_cmb)),                  
            real(dw(l=1,m=0)),imag(dw(l=1,m=0)),                  
            real(dw(l=2,m=0)),imag(dw(l=2,m=0)),                  
            ...
            real(dw(l=l_max_cmb,m=l_max_cmb)),imag(dw(l=l_max_cmb,m=l_max_cmb)),                  
            real(z(l=1,m=0)),imag(z(l=1,m=0)),                  
            real(z(l=2,m=0)),imag(z(l=2,m=0)),                  
            ...
            real(z(l=l_max_cmb,m=l_max_cmb)),imag(z(l=l_max_cmb,m=l_max_cmb)),                  
   	    !------------
   	    !-- Line n+1
   	    !------------
            ...


The detailed calculations are done in the subroutine :f:subr:`write_Bcmb <out_coeff/write_bcmb()>`.


 +-----------------------------------------------------------------------------------------------------------+
 | Header: ``l_max_cmb, minc, n_data``                                                                       |
 +-----------------------------------------------------------------------------------------------------------+ 
 | :math:`t_1`, :math:`Re( b(\ell=0,m=0) ), Im( b(\ell=0,m=0) ),`                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=0) ), Im( b(\ell=1,m=0) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=1) ), Im( b(\ell=1,m=1) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=2,m=0) ), Im( b(\ell=2,m=0) ),`                                                    |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=\ell_{max,cmb},m=\ell_{max,cmb}) ), Im( b(\ell=\ell_{max,cmb},m=\ell_{max,cmb}) )` |
 +-----------------------------------------------------------------------------------------------------------+
 |                                           :math:`\cdots`                                                  |
 +-----------------------------------------------------------------------------------------------------------+
 | :math:`t_N`, :math:`Re( b(\ell=0,m=0) ), Im( b(\ell=0,m=0) ),`                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=0) ), Im( b(\ell=1,m=0) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=1) ), Im( b(\ell=1,m=1) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=2,m=0) ), Im( b(\ell=2,m=0) ),`                                                    |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=\ell_{max,cmb},m=\ell_{max,cmb}) ), Im( b(\ell=\ell_{max,cmb},m=\ell_{max,cmb}) )` |
 +-----------------------------------------------------------------------------------------------------------+ 

where :math:`t_j` is the time during the :math:`j^{th}` time step, :math:`t_N` being the last step. :math:`Re` and :math:`Im` denote real and imaginary parts, respectively, of :math:`b(\ell,m)` - the coefficient of the spherical harmonic with degree :math:`\ell` and order :math:`m`, for the poloidal potential of the magnetic field at the outer boundary (CMB).

.. _secBcoeffrFile:

``B_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_field=.true. <varl_r_field>`.

This file saves time and the poloidal and toroidal coefficients of the magnetic field - :math:`w,dw,z` at a specific radius, :math:`r` up to degree :ref:`l_max_r <varl_max_cmb>`.

 +-----------------------------------------------------------------------------------------------------------+
 | Header: ``l_max_r,minc,n_data,r``                                                                         |
 +-----------------------------------------------------------------------------------------------------------+ 
 | :math:`t_1`, :math:`Re( b(\ell=0,m=0) ), Im( b(\ell=0,m=0) ),`                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=0) ), Im( b(\ell=1,m=0) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=1) ), Im( b(\ell=1,m=1) ),`                                                    |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |      :math:`Re( db(\ell=0,m=0) ), Im( db(\ell=0,m=0) ),`                                                  |
 |                                                                                                           |
 |      :math:`Re( db(\ell=1,m=0) ), Im( db(\ell=1,m=0) ),`                                                  |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |  :math:`Re( ddb(\ell=\ell_{max,r},m=\ell_{max,r}) ), Im( ddb(\ell=\ell_{max,cmb},m=\ell_{max,r}) )`       |
 +-----------------------------------------------------------------------------------------------------------+
 |                                           :math:`\cdots`                                                  |
 +-----------------------------------------------------------------------------------------------------------+
 | :math:`t_N`, :math:`Re( b(\ell=0,m=0) ), Im( b(\ell=0,m=0) ),`                                            |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=0) ), Im( b(\ell=1,m=0) ),`                                                    |
 |                                                                                                           |
 |      :math:`Re( b(\ell=1,m=1) ), Im( b(\ell=1,m=1) ),`                                                    |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |      :math:`Re( db(\ell=0,m=0) ), Im( db(\ell=0,m=0) ),`                                                  |
 |                                                                                                           |
 |      :math:`Re( db(\ell=1,m=0) ), Im( db(\ell=1,m=0) ),`                                                  |
 |                                                                                                           |
 |                 :math:`\cdots`                                                                            |
 |                                                                                                           |
 |  :math:`Re( ddb(\ell=\ell_{max,r},m=\ell_{max,r}) ), Im( ddb(\ell=\ell_{max,cmb},m=\ell_{max,r}) )`       |
 +-----------------------------------------------------------------------------------------------------------+
  

.. _secVcoeffrFile:

``V_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_field=.true. <varl_r_field>`


.. _secTcoeffrFile:

``T_coeff_r#.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_r_fieldT=.true. <varl_r_fieldT>`
