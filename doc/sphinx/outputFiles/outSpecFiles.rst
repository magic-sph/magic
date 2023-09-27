.. _secSpecFiles:

Spectra
=======

.. _secKinSpecFile:

``kin_spec_#.TAG``
------------------

This file contains the kinetic energy spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+--------------------------------------------------------------+
   | No. of column | Contents                                                     |
   +===============+==============================================================+
   | 1             | degree / order                                               |
   +---------------+--------------------------------------------------------------+
   | 2             | Poloidal kinetic energy versus degree                        |
   +---------------+--------------------------------------------------------------+
   | 3             | Poloidal kinetic energy versus order                         |
   +---------------+--------------------------------------------------------------+
   | 4             | Toroidal kinetic energy versus degree                        |
   +---------------+--------------------------------------------------------------+
   | 5             | Toroidal kinetic energy versus order                         |
   +---------------+--------------------------------------------------------------+
   | 6             | Poloidal kinetic energy  at :math:`r=r_o-0.01` versus degree |
   +---------------+--------------------------------------------------------------+
   | 7             | Poloidal kinetic energy  at :math:`r=r_o-0.01` versus order  |
   +---------------+--------------------------------------------------------------+
   | 8             | Axisymmetric poloidal kinetic energy versus degree           |
   +---------------+--------------------------------------------------------------+
   | 9             | Axisymmetric toroidal kinetic energy versus degree           |
   +---------------+--------------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> sp = MagicSpectrum(field='ekin')

.. _secMagSpecFile:

``mag_spec_#.TAG``
------------------

This file contains the magnetic energy spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal magnetic energy in the outer core versus degree  |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal magnetic energy in the outer core versus order   |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal magnetic energy in the outer core versus degree  |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal magnetic energy in the outer core versus order   |
   +---------------+-----------------------------------------------------------+
   | 6             | Poloidal magnetic energy in the inner core versus degree  |
   +---------------+-----------------------------------------------------------+
   | 7             | Poloidal magnetic energy in the inner core versus order   |
   +---------------+-----------------------------------------------------------+
   | 8             | Toroidal magnetic energy in the inner core versus degree  |
   +---------------+-----------------------------------------------------------+
   | 9             | Toroidal magnetic energy in the inner core versus order   |
   +---------------+-----------------------------------------------------------+
   | 10            | Poloidal magnetic energy at the CMB versus degree         |
   +---------------+-----------------------------------------------------------+
   | 11            | Poloidal magnetic energy at the CMB versus order          |
   +---------------+-----------------------------------------------------------+
   | 12            | Poloidal magnetic energy at the CMB                       |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> sp = MagicSpectrum(field='emag')

.. _secu2SpecFile:

``u2_spec_#.TAG``
-----------------

.. note:: This file is **only** written in anelastic models, i.e. either when
          :ref:`strat/=0 <varstrat>` or when :ref:`interior_model/="None" <varinterior_model>`

This file contains the spectra of the square velocity. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal contribution per degree in the outer core        |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal contribution per order in the outer core         |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal contribution per degree in the outer core        |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal contribution per order in the outer core         |
   +---------------+-----------------------------------------------------------+
   | 6             | Axisymmetric poloidal contribution versus degree          |
   +---------------+-----------------------------------------------------------+
   | 7             | Axisymmetric toroidal contribution versus degree          |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``u2_spec_1.test``:
   >>> sp = MagicSpectrum(field='u2', ispec=1, tag='test')


.. _secTSpecFile:

``T_spec_#.TAG``
----------------

This file contains the temperature/entropy spectra, those are defined by taking the
square of temperature/entropy. It is written by the subroutine
:f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Square temperature/entropy versus degree                  |
   +---------------+-----------------------------------------------------------+
   | 3             | Square temperature/entropy versus order                   |
   +---------------+-----------------------------------------------------------+
   | 4             | Square temperature/entropy at the ICB versus degree       |
   +---------------+-----------------------------------------------------------+
   | 5             | Square temperature/entropy at the ICB versus order        |
   +---------------+-----------------------------------------------------------+
   | 6             | Square radial derivative of temperature/entropy at the ICB|
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 7             | Square radial derivative of temperature/entropy at the ICB|
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``T_spec_3.test_a``:
   >>> sp = MagicSpectrum(field='T', ispec=3, tag='test_a')

.. _secXiSpecFile:

``Xi_spec_#.TAG``
-----------------

This file contains the spectra of chemical composition, this is defined by taking the
square of chemical composition. It is written by the subroutine
:f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-------------------------------------------------------------+
   | No. of column | Contents                                                    |
   +===============+=============================================================+
   | 1             | degree / order                                              |
   +---------------+-------------------------------------------------------------+
   | 2             | Square chemical composition versus degree                   |
   +---------------+-------------------------------------------------------------+
   | 3             | Square chemical composition versus order                    |
   +---------------+-------------------------------------------------------------+
   | 4             | Square chemical composition at the ICB versus degree        |
   +---------------+-------------------------------------------------------------+
   | 5             | Square chemical composition at the ICB versus order         |
   +---------------+-------------------------------------------------------------+
   | 6             | Square radial derivative of chemical composition at the ICB |
   |               | versus degree                                               |
   +---------------+-------------------------------------------------------------+
   | 7             | Square radial derivative of chemical composition at the ICB |
   |               | versus order                                                |
   +---------------+-------------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``Xi_spec_3.test``:
   >>> sp = MagicSpectrum(field='Xi', ispec=3, tag='test')

.. _secPhaseSpecFile:

``Phase_spec_#.TAG``
--------------------

This file contains the phase field spectra, those are defined by taking the
square of phase field. It is written by the subroutine
:f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Square of the phase field as a function of degree         |
   +---------------+-----------------------------------------------------------+
   | 3             | Square of the phase field as a function of order          |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``Phase_spec_3.cheb``:
   >>> sp = MagicSpectrum(field='phase', ispec=3, tag='cheb')

.. _sec2DSpectra:

2D spectra ``2D_[kin|mag]_spec_#.TAG`` and ``2D_[kin|mag]_spec_ave.TAG``
---------------------------------------------------------------------------

.. note:: Those files are **only** written when :ref:`l_2D_spectra=.true. <varl_2D_spectra>`. The time-averaged files also require that :ref:`l_spec_avg=.true. <varl_spec_avg>`.

Those files contain 2-D spectra in the :math:`(r,\ell)` and in the
:math:`(r,m)` planes.  In other words, the poloidal and toroidal energies
versus degree :math:`\ell` or versus order :math:`m` are computed for all
radii. There are two kinds of those files that correspond to the
aforementioned spectra, namely **2D_kin_spec_#.TAG**, **2D_mag_spec_#.TAG**.
In case time-averages are requested, **2D_kin_spec_ave.TAG** and
**2D_mag_spec_ave.TAG** will also be stored. The calculations are done
in the subroutine
:f:subr:`spectrum <spectra/spectrum()>`. The structure of the output files
are same for these three outputs. They are stored as fortran unformatted files.

Unformatted files are not directly human readable, and are used to store binary
data and move it around without changing the internal representation. In
fortran, the open, read and write operations for these files are performed as follows:

.. code-block:: fortran

  open(unit=4, file='test', form='unformatted')
  read(unit=4) readVar
  write(unit=n_out, iostat=ios) writeVar ! Unformatted write

2D spectra files have the following structure:

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       version ! version number

       !-------------
       ! Line 2
       !-------------

       time, n_r_max, l_max, minc ! Time, resolution, max(\ell), azimuthal symmetry

       !-------------
       ! Line 3
       !-------------

       r(1), r(2), r(3), ..., r(n_r_max)                  ! Radius

       !-------------
       ! Line 4
       !-------------

       e_p_l(l=1,r=1), e_p_l(l=1,r=2), ..., e_p_l(l=1,r=n_r_max),      ! Poloidal energy
       ...                                                             ! versus degree
       e_p_l(l=l_max,r=1), e_p_l(l=l_max,r=2), ..., e_p_l(l=l_max,r=n_r_max),

       !-------------
       ! Line 5
       !-------------

       e_p_m(m=0,r=1), e_p_l(m=0,r=2), ..., e_p_l(m=1,r=n_r_max),      ! Poloidal energy
       ...                                                             ! versus order
       e_p_l(m=l_max,r=1), e_p_l(m=l_max,r=2), ..., e_p_l(m=l_max,r=n_r_max),

       !-------------
       ! Line 6
       !-------------

       e_t_l(l=1,r=1), e_t_l(l=1,r=2), ..., e_t_l(l=1,r=n_r_max),      ! Toroidal energy
       ...                                                             ! versus degree
       e_t_l(l=l_max,r=1), e_t_l(l=l_max,r=2), ..., e_t_l(l=l_max,r=n_r_max),

       !-------------
       ! Line 7
       !-------------

       e_t_m(m=0,r=1), e_t_l(m=0,r=2), ..., e_t_l(m=1,r=n_r_max),      ! Toroidal energy
       ...                                                             ! versus order
       e_t_l(m=l_max,r=1), e_t_l(m=l_max,r=2), ..., e_t_l(m=l_max,r=n_r_max),

       !-------------
       ! Line 8
       !-------------

       e_pa_l(l=1,r=1), e_pa_l(l=1,r=2), ..., e_pa_l(l=1,r=n_r_max),   ! Pol. axi. energy
       ...                                                             ! versus degree
       e_pa_l(l=l_max,r=1), e_pa_l(l=l_max,r=2), ..., e_pa_l(l=l_max,r=n_r_max),

       !-------------
       ! Line 9
       !-------------

       e_ta_l(l=1,r=1), e_ta_l(l=1,r=2), ..., e_ta_l(l=1,r=n_r_max),   ! Tor. axi. energy
       ...                                                             ! versus degree
       e_ta_l(l=l_max,r=1), e_ta_l(l=l_max,r=2), ..., e_ta_l(l=l_max,r=n_r_max),


Those files can be read using the python class :py:class:`MagicSpectrum2D <magic.MagicSpectrum2D>` with
the following options:

   >>> # Read the file 2D_mag_spec_3.ext
   >>> sp = MagicSpectrum2D(tag='ext', field='e_mag', ispec=3)
   >>> # Print e_pol_l and e_tor_m
   >>> print(sp.e_pol_l, sp.e_tor_m)


.. _secKinSpecAveFile:

``kin_spec_ave.TAG``
--------------------

.. note:: This file is **only** written when :ref:`l_spec_avg=.true. <varl_spec_avg>`


This file contains the time-average kinetic energy spectra as well as squared quantities
to allow a possible further reconstruction of the standard deviation.
This file is written by the subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+---------------------------------------------------------------------+
   | No. of column | Contents                                                            |
   +===============+=====================================================================+
   | 1             | degree / order                                                      |
   +---------------+---------------------------------------------------------------------+
   | 2             | Time-averaged poloidal kinetic energy versus degree                 |
   +---------------+---------------------------------------------------------------------+
   | 3             | Time-averaged poloidal kinetic energy versus order                  |
   +---------------+---------------------------------------------------------------------+
   | 4             | Time-averaged toroidal kinetic energy versus degree                 |
   +---------------+---------------------------------------------------------------------+
   | 5             | Time-averaged toroidal kinetic energy versus order                  |
   +---------------+---------------------------------------------------------------------+
   | 6             | Time-averaged axisymmetric poloidal kinetic energy versus degree    |
   +---------------+---------------------------------------------------------------------+
   | 7             | Time-averaged axisymmetric toroidal kinetic energy versus degree    |
   +---------------+---------------------------------------------------------------------+
   | 8             | Standard deviation of poloidal kinetic energy versus degree         |
   +---------------+---------------------------------------------------------------------+
   | 9             | Standard deviation of poloidal kinetic energy versus order          |
   +---------------+---------------------------------------------------------------------+
   | 10            | Standard deviation of toroidal kinetic energy versus degree         |
   +---------------+---------------------------------------------------------------------+
   | 11            | Standard deviation of toroidal kinetic energy versus order          |
   +---------------+---------------------------------------------------------------------+
   | 12            | Standard deviation of axisym. poloidal kinetic energy versus degree |
   +---------------+---------------------------------------------------------------------+
   | 13            | Standard deviation of axisym. toroidal kinetic energy versus degree |
   +---------------+---------------------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``kin_spec_ave.test``:
   >>> sp = MagicSpectrum(field='kin', ave=True, tag='test')

.. _secMagSpecAveFile:

``mag_spec_ave.TAG``
--------------------

.. note:: This file is **only** written when :ref:`l_spec_avg=.true. <varl_spec_avg>` and
          the run is magnetic

This file contains the time-average magnetic energy spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+------------------------------------------------------------------------+
   | No. of column | Contents                                                               |
   +===============+========================================================================+
   | 1             | degree / order                                                         |
   +---------------+------------------------------------------------------------------------+
   | 2             | Time-averaged poloidal magnetic energy in the outer core versus degree |
   +---------------+------------------------------------------------------------------------+
   | 3             | Time-averaged poloidal magnetic energy in the outer core versus order  |
   +---------------+------------------------------------------------------------------------+
   | 4             | Time-averaged toroidal magnetic energy in the outer core versus degree |
   +---------------+------------------------------------------------------------------------+
   | 5             | Time-averaged toroidal magnetic energy in the outer core versus order  |
   +---------------+------------------------------------------------------------------------+
   | 6             | Time-averaged poloidal magnetic energy at the CMB versus degree        |
   +---------------+------------------------------------------------------------------------+
   | 7             | Time-averaged poloidal magnetic energy at the CMB versus order         |
   +---------------+------------------------------------------------------------------------+
   | 8             | Standard deviation of the poloidal magnetic energy in the outer        |
   |               | core versus degree                                                     |
   +---------------+------------------------------------------------------------------------+
   | 9             | Standard deviation of the poloidal magnetic energy in the outer core   |
   |               | versus order                                                           |
   +---------------+------------------------------------------------------------------------+
   | 10            | Standard deviation of the toroidal magnetic energy in the outer core   |
   |               | versus degree                                                          |
   +---------------+------------------------------------------------------------------------+
   | 11            | Standard deviation of the toroidal magnetic energy in the outer core   |
   |               | versus order                                                           |
   +---------------+------------------------------------------------------------------------+
   | 12            | Standard deviation of the magnetic energy at the CMB                   |
   |               | versus degree                                                          |
   +---------------+------------------------------------------------------------------------+
   | 13            | Standard deviation of the magnetic energy at the CMB                   |
   |               | versus order                                                           |
   +---------------+------------------------------------------------------------------------+


This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``mag_spec_ave.test``:
   >>> sp = MagicSpectrum(field='mag', ave=True, tag='test')


.. _secTempSpecAveFile:

``T_spec_ave.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_spec_avg=.true. <varl_spec_avg>`

This file contains the time-averaged temperature/entropy spectra and their standard
deviation. It is written by the subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | Spherical harmonic degree/order                           |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged RMS temperature/entropy versus degree       |
   +---------------+-----------------------------------------------------------+
   | 3             | Time-averaged RMS temperature/entropy versus order        |
   +---------------+-----------------------------------------------------------+
   | 4             | Time-averaged RMS temperature/entropy at the ICB versus   |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 5             | Time-averaged RMS temperature/entropy at the ICB versus   |
   |               | order                                                     |
   +---------------+-----------------------------------------------------------+
   | 6             | Time-averaged temperature/entropy gradient at the ICB     |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 7             | Time-averaged temperature/entropy gradient at the ICB     |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+
   | 8             | Standard deviation of the temperature/entropy versus      |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 9             | Standard deviation of the temperature/entropy versus      |
   |               | order                                                     |
   +---------------+-----------------------------------------------------------+
   | 10            | Standard deviation of the temperature/entropy at the ICB  |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 11            | Standard deviation of the temperature/entropy at the ICB  |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+
   | 12            | Standard deviation of the temperature/entropy gradient    |
   |               | at the ICB  versus degree                                 |
   +---------------+-----------------------------------------------------------+
   | 13            | Standard deviation of the temperature/entropy gradient    |
   |               | at the ICB  versus order                                  |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``T_spec_ave.test``:
   >>> sp = MagicSpectrum(field='T', ave=True, tag='test')

.. _secXipecAveFile:

``Xi_spec_ave.TAG``
-------------------

.. note:: This file is **only** written when :ref:`l_spec_avg=.true. <varl_spec_avg>`

This file contains the time-averaged composition spectra and their standard
deviation. It is written by the subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | Spherical harmonic degree/order                           |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged squared composition versus degree           |
   +---------------+-----------------------------------------------------------+
   | 3             | Time-averaged squared composition versus order            |
   +---------------+-----------------------------------------------------------+
   | 4             | Time-averaged squared composition at the ICB versus       |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 5             | Time-averaged squared composition at the ICB versus       |
   |               | order                                                     |
   +---------------+-----------------------------------------------------------+
   | 6             | Time-averaged squared composition gradient at the ICB     |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 7             | Time-averaged squared composition gradient at the ICB     |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+
   | 8             | Standard deviation of the squared composition versus      |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 9             | Standard deviation of the squared composition versus      |
   |               | order                                                     |
   +---------------+-----------------------------------------------------------+
   | 10            | Standard deviation of the squared composition at the ICB  |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 11            | Standard deviation of the squared composition at the ICB  |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+
   | 12            | Standard deviation of the squared composition gradient    |
   |               | at the ICB  versus degree                                 |
   +---------------+-----------------------------------------------------------+
   | 13            | Standard deviation of the squared composition gradient    |
   |               | at the ICB  versus order                                  |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``Xi_spec_ave.cheb``:
   >>> sp = MagicSpectrum(field='Xi', ave=True, tag='cheb')

.. _secPhaseSpecAveFile:

``Phase_spec_ave.TAG``
-----------------------

.. note:: This file is **only** written when :ref:`l_spec_avg=.true. <varl_spec_avg>`

This file contains the time-averaged phase field spectra and their standard
deviation. It is written by the subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | Spherical harmonic degree/order                           |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged square phase field versus degree            |
   +---------------+-----------------------------------------------------------+
   | 3             | Time-averaged square phase field versus order             |
   +---------------+-----------------------------------------------------------+
   | 4             | Standard deviation of the squared phase field field       |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 5             | Standard deviation of the squared phase field field       |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``Phase_spec_ave.test``:
   >>> sp = MagicSpectrum(field='phase', ave=True, tag='test')


.. _secRMSSpectra:

``dtVrms_spec.TAG``
--------------------

.. note:: This file is **only** written when :ref:`l_RMS=.true. <varl_RMS>`

This file contains the time-averaged force balance spectra as well as their standard deviation.
The calculations are done in the subroutine :f:subr:`dtVrms <out_rms/dtvrms()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree + 1                                                |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged Inertia versus degree                       |
   +---------------+-----------------------------------------------------------+
   | 3             | Time-averaged Coriolis force versus degree                |
   +---------------+-----------------------------------------------------------+
   | 4             | Time-averaged Lorentz force versus degree                 |
   +---------------+-----------------------------------------------------------+
   | 5             | Time-averaged Advection term versus degree                |
   +---------------+-----------------------------------------------------------+
   | 6             | Time-averaged Viscous force versus degree                 |
   +---------------+-----------------------------------------------------------+
   | 7             | Time-averaged thermal Buoyancy versus degree              |
   +---------------+-----------------------------------------------------------+
   | 8             | Time-averaged chemical Buoyancy versus degree             |
   +---------------+-----------------------------------------------------------+
   | 9             | Time-averaged Pressure gradient versus degree             |
   +---------------+-----------------------------------------------------------+
   | 10             | Time-averaged Pressure/Coriolis balance versus degree    |
   +---------------+-----------------------------------------------------------+
   | 11            | Time-averaged Pressure/Coriolis/Lorentz balance versus    |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 12            | Time-averaged Pressure/Coriolis/Buoyancy balance versus   |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 13            | Time-averaged Pressure/Coriolis/Lorentz/Buoyancy balance  |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 14            | Time-averaged Coriolis/Lorentz balance versus degree      |
   +---------------+-----------------------------------------------------------+
   | 15            | Time-averaged Pressure/Lorentz balance versus degree      |
   +---------------+-----------------------------------------------------------+
   | 16            | Time-averaged Coriolis/Inertia/Buoyancy balance versus    |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 17            | Standard deviation of Inertia versus degree               |
   +---------------+-----------------------------------------------------------+
   | 18            | Standard deviation of Coriolis force versus degree        |
   +---------------+-----------------------------------------------------------+
   | 19            | Standard deviation of Lorentz force versus degree         |
   +---------------+-----------------------------------------------------------+
   | 20            | Standard deviation of Advection term versus degree        |
   +---------------+-----------------------------------------------------------+
   | 21            | Standard deviation of Viscous force versus degree         |
   +---------------+-----------------------------------------------------------+
   | 22            | Standard deviation of thermal Buoyancy versus degree      |
   +---------------+-----------------------------------------------------------+
   | 23            | Standard deviation of chemical Buoyancy versus degree     |
   +---------------+-----------------------------------------------------------+
   | 24            | Standard deviation of Pressure gradient versus degree     |
   +---------------+-----------------------------------------------------------+
   | 25            | Standard deviation of Pressure/Coriolis balance versus    |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 26            | Standard deviation of Pressure/Coriolis/Lorentz balance   |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 27            | Standard deviation of Pressure/Coriolis/Buoyancy balance  |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 28            | Standard deviation of Pressure/Coriolis/Lorentz/Buoyancy  |
   |               | balance versus degree                                     |
   +---------------+-----------------------------------------------------------+
   | 29            | Standard deviation of Coriolis/Lorentz balance versus     |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 30            | Standard deviation of Pressure/Lorentz balance versus     |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 31            | Standard deviation of Coriolis/Inertia/Buoyancy balance   |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+


This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``dtVrms_spec.test``:
   >>> sp = MagicSpectrum(field='dtVrms', tag='test')


.. _sec2DRMSSpectra:

2D force balance spectra ``2D_dtVrms_spec.TAG``
-----------------------------------------------

.. note:: Those files are **only** written when :ref:`l_RMS=.true. <varl_RMS>` and :ref:`l_2D_RMS=.true. <varl_2D_RMS>`.

Those files contain 2-D force balance spectra in the :math:`(r,\ell)` plane.
The calculations are done in the subroutine :f:subr:`dtVrms <out_rms/dtvrms()>`.
The output file is stored as a Fortran unformatted file.

The structure of the 2D force balance spectra files are as follows:

   .. code-block:: fortran

       !------------
       ! Line 1
       !------------

       version

       !-------------
       ! Line 2
       !-------------

       n_r_max, l_max ! radial resolution, max(\ell)

       !-------------
       ! Line 3
       !-------------

       r(1), r(2), r(3), ..., r(n_r_max)                  ! Radius

       !-------------
       ! Line 4
       !-------------

       Cor_l(l=1,r=1), Cor_l(l=1,r=2), ..., Cor_l(l=1,r=n_r_max),      ! Coriolis force
       ...                                                             ! versus degree
       Cor_l(l=l_max,r=1), Cor_l(l=l_max,r=2), ..., Cor_l(l=l_max,r=n_r_max),

       !-------------
       ! Line 5
       !-------------

       Adv_l ! Advection

       !-------------
       ! Line 6
       !-------------

       LF_l ! Lorentz force

       !-------------
       ! Line 7
       !-------------

       Buo_temp_l ! Thermal buoyancy

       !------------
       ! Line 8
       !------------

       Buo_xi_l ! Chemical buoyancy

       !-------------
       ! Line 9
       !-------------

       Pre_l ! Pressure

       !-------------
       ! Line 10
       !-------------

       Dif_l ! Viscosity

       !-------------
       ! Line 11
       !-------------

       Iner_l ! Inertia

       !-------------
       ! Line 12
       !-------------

       Geo_l ! Sum of force terms: geostrophic balance

       !-------------
       ! Line 13
       !-------------

       Mag_l ! Sum of force terms: pressure, Coriolis and Lorentz

       !-------------
       ! Line 14
       !-------------

       Arc_l ! Sum of force terms: pressure, buoyancy and Coriolis

       !-------------
       ! Line 15
       !-------------

       ArcMag_l ! Sum of force terms: pressure, buoyancy, Coriolis and Lorentz

       !-------------
       ! Line 16
       !-------------

       CIA_l ! Sum of force terms Coriolis/Inertia/Archimedean

       !-------------
       ! Line 17
       !-------------

       CLF_l ! Sum of force terms Coriolis/Lorentz

       !-------------
       ! Line 18
       !-------------

       PLF_l ! Sum of force terms Pression/Lorentz


Those files can be read using the python class :py:class:`MagicSpectrum2D <magic.MagicSpectrum2D>` with the following options:

   >>> # Read the file 2D_dtVrms_spec.ext
   >>> sp = MagicSpectrum2D(tag='ext', field='dtVrms')
   >>> # Print Cor_l
   >>> print(sp.Cor_l)

.. _secTimeSpectraFiles:

2D spectra `am_[kin|mag]_[pol|tor].TAG`
---------------------------------------

Those files contain the time evolution of the poloidal and toroidal kinetic and
magnetic spectra for a given range of spherical harmonic orders :math:`m`.
There are four kinds of those files that correspond to the aforementioned
spectra, namely **am_kin_pol.TAG**, **am_kin_tor.TAG**, **am_mag_pol.TAG** and
**am_mag_tor.TAG**. The calculations are done in the subroutine
:f:subr:`get_amplitude <spectra/get_amplitude()>`. The structure of the output
files is the same for the four outputs (fortran unformatted
files):

   .. code-block:: fortran

       !-------------
       ! Line 1
       !-------------

       time(t=0), e_p_m(m=0,t=0), e_p_m(m=1,t=0), ..., e_p_m(m=m_max_modes,t=0)

       ...

       !-------------
       ! Line N
       !-------------

       time(t=N), e_p_m(m=0,t=N), e_p_m(m=1,t=N), ..., e_p_m(m=m_max_modes,t=N)

       ...


Those files can be read using the python class :py:class:`MagicTs
<magic.MagicTs>` with the following options:

   >>> # Read the file am_mag_pol.ext
   >>> ts = MagicTs(field='am_mag_pol', tag='ext')
   >>> # Print the time
   >>> print(ts.time)
   >>> # Print the energy content in m=11 for all times
   >>> print(ts.coeffs[:, 11])
