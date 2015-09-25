.. _secSpecFiles:

Spectra
=======

.. _secKinSpecFile:

``kin_spec_#.TAG``
------------------

This file contains the magnetic spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal kinetic energy per degree in the outer core      |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal kinetic energy per order in the outer core       |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal kinetic energy per degree in the outer core      |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal kinetic energy per order in the outer core       |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> sp = MagicSpectrum(field='ekin')

.. _secMagSpecFile:

``mag_spec_#.TAG``
------------------

This file contains the magnetic spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

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

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``u2_spec_1.test``:
   >>> sp = MagicSpectrum(field='u2', ispec=1, tag='test')

                                             
.. _secTSpecFile:

``T_spec_#.TAG``
----------------

This file contains the temperature/entropy spectra. It is written by the subroutine
:f:subr:`spectrum_temp <spectra/spectrum_temp()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | RMS temperature/entropy versus degree                     |
   +---------------+-----------------------------------------------------------+
   | 3             | RMS temperature/entropy versus order                      |
   +---------------+-----------------------------------------------------------+
   | 4             | RMS temperature/entropy at the ICB versus degree          |
   +---------------+-----------------------------------------------------------+
   | 5             | RMS temperature/entropy at the ICB versus order           |
   +---------------+-----------------------------------------------------------+
   | 6             | RMS radial derivative of temperature/entropy at the ICB   |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 7             | RMS radial derivative of temperature/entropy at the ICB   |
   |               | versus order                                              |
   +---------------+-----------------------------------------------------------+

