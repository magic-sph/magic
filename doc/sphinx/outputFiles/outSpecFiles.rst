Spectra
=======

.. _secKinSpecFile:

kin_spec_#.TAG
--------------

This files contains the magnetic spectra. This file is written by the subroutine ``spectrum`` in the file ``s_spectrum.f90``.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal magnetic energy per degree in the outer core     |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal magnetic energy per order in the outer core      |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal magnetic energy per degree in the outer core     |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal magnetic energy per order in the outer core      |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`magic.MagicSpectrum` with the following options::
   >>> sp = MagicSpectrum(field='ekin')

.. _secMagSpecFile:

mag_spec_#.TAG
--------------

This files contains the magnetic spectra. This file is written by the subroutine ``spectrum`` in the file ``s_spectrum.f90``.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal magnetic energy per degree in the outer core     |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal magnetic energy per order in the outer core      |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal magnetic energy per degree in the outer core     |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal magnetic energy per order in the outer core      |
   +---------------+-----------------------------------------------------------+
   | 6             | Poloidal magnetic energy per degree in the inner core     |
   +---------------+-----------------------------------------------------------+
   | 7             | Poloidal magnetic energy per order in the inner core      |
   +---------------+-----------------------------------------------------------+
   | 8             | Toroidal magnetic energy per degree in the inner core     |
   +---------------+-----------------------------------------------------------+
   | 9             | Toroidal magnetic energy per order in the inner core      |
   +---------------+-----------------------------------------------------------+
   | 10            | Poloidal magnetic energy per degree at the CMB            |
   +---------------+-----------------------------------------------------------+
   | 11            | Poloidal magnetic energy per order at the CMB             |
   +---------------+-----------------------------------------------------------+
   | 12            | Poloidal magnetic energy at the CMB                       |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`magic.MagicSpectrum` with the following options::
   >>> sp = MagicSpectrum(field='emag')
                                             
.. _secTSpecFile:

T_spec_#.TAG
------------
