Output files
############

While some information of a run is written into ``STDOUT`` to monitor its
progress, most outputs are printed into dedicated files identified by the
:ref:`TAG <varTag>` extension. Most of the information found in ``STDOUT`` is also written
to the log-file called :ref:`log.TAG <secLogFile>`. In addition, this file contains all input
parameters, truncation, information on other output files, and some results
like the time averaged energies (for :ref:`l_average=.true. <varl_average>`). The python
functions that can be used to process the results can extract the informations
concerning the run from :ref:`log.TAG <secLogFile>`.  Other output files are organised in
columns or lines (datasets). Their meaning is explained below. The number of
graphic files (``G_#.TAG``), restart files (``rst_#.tag``), and spectrum files
(``(kin|mag)_spec_#.TAG``) are determined by the input parameters in the
namelist :ref:`&output <secOutputNml>`. A new file is produced for each output time. If time
averaging has been chosen the time-averaged graphic, potential and spectra
files will have the number 0 and a prefix of the form ``_ave``. All other files
are numbered consecutively starting with 1. Times at which these files have
been written can be found in :ref:`log.tag <secLogFile>`. Other files (except the log-file)
contain time series. The frequency of storage for the kinetic and magnetic
energy files, the :ref:`rot.TAG <secRotFile>` file, the :ref:`dipole.TAG
<secDipoleFile>` file, the :ref:`par.TAG <secParFile>` file and the
:ref:`misc.TAG <secMiscFile>` file are determined by the log-times. Output
times for cmb files (``B_coeff_cmb.TAG``), coeff files
(``(T|V|B)_coeff_r.TAG``), potential files, torsional oscillations files and
movie frames are chosen independently (see above). Those files are stored as
unformatted files described below.


.. toctree::

   logFile.rst
   outTimeSeries.rst
   outRadialFiles.rst
   outTransportProp.rst
   outSpecFiles.rst
