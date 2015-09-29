.. _secSpecFiles:

Spectra
=======

.. _secKinSpecFile:

``kin_spec_#.TAG``
------------------

This file contains the kinetic energy spectra. This file is written by the
subroutine :f:subr:`spectrum <spectra/spectrum()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Poloidal kinetic energy versus degree                     |
   +---------------+-----------------------------------------------------------+
   | 3             | Poloidal kinetic energy versus order                      |
   +---------------+-----------------------------------------------------------+
   | 4             | Toroidal kinetic energy versus degree                     |
   +---------------+-----------------------------------------------------------+
   | 5             | Toroidal kinetic energy versus order                      |
   +---------------+-----------------------------------------------------------+

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


.. _secKinSpecAveFile:

``kin_spec_ave.TAG``
--------------------

.. note:: This file is **only** written when :ref:`l_average=.true. <varl_average>`


This file contains the time-average kinetic energy spectra as well as squared quantities
to allow a possible further reconstruction of the standard deviation. 
This file is written by the subroutine :f:subr:`spectrum_average <spectra/spectrum_average()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | degree / order                                            |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged poloidal kinetic energy versus degree       |
   +---------------+-----------------------------------------------------------+
   | 3             | Time-averaged poloidal kinetic energy versus order        |
   +---------------+-----------------------------------------------------------+
   | 4             | Time-averaged toroidal kinetic energy versus degree       |
   +---------------+-----------------------------------------------------------+
   | 5             | Time-averaged toroidal kinetic energy versus order        |
   +---------------+-----------------------------------------------------------+
   | 6             | Time-averaged poloidal kinetic energy square versus degree|
   +---------------+-----------------------------------------------------------+
   | 7             | Time-averaged poloidal kinetic energy square versus order |
   +---------------+-----------------------------------------------------------+
   | 8             | Time-averaged toroidal kinetic energy square versus degree|
   +---------------+-----------------------------------------------------------+
   | 9             | Time-averaged toroidal kinetic energy square versus order |
   +---------------+-----------------------------------------------------------+

This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``kin_spec_ave.test``:
   >>> sp = MagicSpectrum(field='kin', ave=True, tag='test')

.. _secMagSpecAveFile:

``mag_spec_ave.TAG``
--------------------

.. note:: This file is **only** written when :ref:`l_average=.true. <varl_average>` and
          the run is magnetic

This file contains the time-average magnetic energy spectra. This file is written by the
subroutine :f:subr:`spectrum_average <spectra/spectrum_average()>`.

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
   | 8             | Time-averaged poloidal magnetic energy in the outer core + its standard|
   |               | deviation versus degree                                                |
   +---------------+------------------------------------------------------------------------+
   | 9             | Time-averaged poloidal magnetic energy in the outer core - its standard|
   |               | deviation versus degree                                                |
   +---------------+------------------------------------------------------------------------+
   | 10            | Time-averaged poloidal magnetic energy in the outer core + its standard|
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+
   | 11            | Time-averaged poloidal magnetic energy in the outer core - its standard|
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+
   | 12            | Time-averaged toroidal magnetic energy in the outer core + its standard|
   |               | deviation versus degree                                                |
   +---------------+------------------------------------------------------------------------+
   | 13            | Time-averaged toroidal magnetic energy in the outer core - its standard|
   |               | deviation versus degree                                                |
   +---------------+------------------------------------------------------------------------+
   | 14            | Time-averaged toroidal magnetic energy in the outer core + its standard|
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+
   | 15            | Time-averaged toroidal magnetic energy in the outer core - its standard|
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+
   | 16            | Time-averaged poloidal magnetic energy at the CMB + its standard       |
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+
   | 17            | Time-averaged poloidal magnetic energy at the CMB - its standard       |
   |               | deviation versus order                                                 |
   +---------------+------------------------------------------------------------------------+


This file can be read using :py:class:`MagicSpectrum <magic.MagicSpectrum>` with the following options:

   >>> # To read the file ``mag_spec_ave.test``:
   >>> sp = MagicSpectrum(field='mag', ave=True, tag='test')


.. _secTempSpecAveFile:

``T_spec_ave.TAG``
------------------

.. note:: This file is **only** written when :ref:`l_average=.true. <varl_average>`

This file contains the time-averaged temperature/entropy spectra and their standard
deviation. It is written by the subroutine :f:subr:`spectrum_temp_average <spectra/spectrum_temp_average()>`.

   +---------------+-----------------------------------------------------------+
   | No. of column | Contents                                                  |
   +===============+===========================================================+
   | 1             | Spherical harmonic degree                                 |
   +---------------+-----------------------------------------------------------+
   | 2             | Time-averaged RMS temperature/entropy versus degree       |
   +---------------+-----------------------------------------------------------+
   | 3             | Standard deviation of the temperature/entropy versus      |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 4             | Time-averaged RMS temperature/entropy at the ICB versus   |
   |               | degree                                                    |
   +---------------+-----------------------------------------------------------+
   | 5             | Standard deviation of the temperature/entropy at the ICB  |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 6             | Time-averaged temperature/entropy gradient at the ICB     |
   |               | versus degree                                             |
   +---------------+-----------------------------------------------------------+
   | 7             | Standard deviation of the temperature/entropy gradient    |
   |               | at the ICB  versus degree                                 |
   +---------------+-----------------------------------------------------------+


.. _sec2DSpectra:

2D spectra
----------
