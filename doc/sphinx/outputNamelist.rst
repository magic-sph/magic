.. _secOutputNml:

Output control namelist
=======================

This namelist contains all the parameters that can be adjusted to control the outputs and diagnostics calculated by the code.


Determining times for output
----------------------------

There are four different ways to control at which time step a specific output should be written. Outputs are generally distributed over the total calculation intervall unless an output time interval is defined by a start time ``t_start`` and a stop time ``t_stop``. If no ``t_start`` is provided, the start time of the calculation is used. If no ``t_stop`` is provided or ``t_stop>t_start`` the total calculation interval is assumed

   1. **Prescribed number of outputs**. The outputs are distributed evenly over the total calculation interval so that the number of timesteps between two outputs is always the same, with the possible exception of the first interval. Last output is written for the last time step, and to compensate the interval before the first output may be longer. However, if ``t_stop`` is provided, the outputs are distributed evenly over the interval [``t_stop``,``t_start``] with equal times intervals between them. 

   .. note:: These input variables are usually named with a pattern that follows ``n_outputName``, for instance, ``n_graphs``, ``n_rsts``, ``n_specs``, ``n_logs``, etc.
             
             In case you want to make use ot a specific time interval, the input variables follow a pattern of the form ``t_outputName_start``, ``t_outputName_stop``. For instance, ``t_graph_start``, ``t_graph_stop``, ``t_log_start``, ``t_log_stop``, ``t_spec_start``, ``t_spec_stop``, etc.

   ..

   2. **User-defined interval between two outputs, given in number of time steps**. Again the last output is performed at the end of the run and a compensation may take place at the beginning.

   .. note:: These input variables are usually named with a pattern that follows ``n_outputName_step``, for instance, ``n_graph_step``, ``n_rst_step``, ``n_spec_step``, ``n_log_step``, ``n_movie_step``, etc.

   ..

   3. **Defined time interval between two outputs**.

   .. note:: These input variables are usually named with a pattern that follows ``dt_outputName``, for instance, :ref:`dt_graph <dtGraphVar>`, ``dt_rst``, ``dt_spec``, ``dt_log``, ``dt_movie``, etc.

   ..

   4. **User-defined times for output**. By default 5000 different times can be defined for each output type. This can be increased by increasing n_time_hits in the file ``output_data.f90``. While the first three possibilities can only be used alternatively, the fourth one can be employed in addition to one of the two others.

   .. note:: These input variables are usually named with a pattern that follows ``t_outputName``, for instance, ``t_graph``, ``t_rst``, ``t_spec``, ``t_log``, ``t_movie``, etc.

   ..


An important parameter in this context is ``l_true_time``. If this is set to ``.true.``, the time steps of the program are modified to meet a desired output time. This forces a recalculation of the inversion matricies and therefore requires some additional computing time. When ``l_true_time=.false.``, the values at the timestep closest to the desired output time are chosen. Since the timesteps are generally small, this option suffices for most applications.

* **l_true_time** (default ``l_true_time=.false.``) is a logical. It causes the code to change time steps to exactly meet the requested output times.

* **l_save_out** (default ``l_save_out=.false.``) is a logical. When set to ``.true.``, the diagnostic files will be safely opened and closed before and after any outputs. When set to ``.false.``, the diagnostic files will be opened before the first iteration timestep and close at the end of the run. This may cost some computing time, but guarantees that only minimal information is lost in case of a crash.

* **lVerbose** (default ``lVerbose=.false.``) is a logical. When set to ``.true.``, the code displays a lot of debugging informations.

  .. warning:: Never set ``lVerbose`` to ``.true.`` for a production run!

* **runid** (default, ``runid="MAGIC default run"``) is a character string. This can be used to briefly describe your run. This information is then for instance stored in the header of the graphic files.


Standard time-series outputs
----------------------------

The **log** outputs controls the output of all the default time series of the
file: kinetic and magnetic energies (:ref:`e_kin.TAG <secEkinFile>`,
:ref:`e_mag_oc.TAG <secEmagocFile>` and :ref:`e_mag_ic.TAG <secEmagicFile>`
files), dipole information (:ref:`dipole.TAG <secDipoleFile>` file), rotation
(:ref:`rot.TAG <secRotFile>`) parameters (:ref:`par.TAG <secParFile>`) and
various additional diagnostics (:ref:`misc.TAG <secMiscFile>`):

* **n_log_step** (default ``n_log_step=50``) is an integer. This is the number of timesteps between two log outputs.

  .. warning:: Be careful: when using too small ``n_log_step``, the disk access will dramatically increases, thus decreasing the code performance.

* **n_logs** (default ``n_logs=0``) is an integer. This is the number of log-information sets to be written.

* **t_log**  (default  ``t_log=-1.0 -1.0 ...``) is real array, which contains the times when log outputs are requested.

* **dt_log** (default ``dt_log=0.0``) is a real, which defines the time interval between log outputs.

* **t_log_start** (default ``t_log_start=0.0``) is a real, which defines the time to start writing log outputs.

* **t_log_stop** (default ``t_log_stop=0.0``) is a real, which defines the time to stop writing log outputs.


Restart files
-------------

The **rst** outputs controls the output of restart files (``rst_t_#.TAG``) (i.e. check points in time from which the code could be restarted):

* **n_rst_step** (default ``n_rst_step=0``) is an integer. This is the number of timesteps between two restart files.

* **n_rsts** (default ``n_rsts=1``) is an integer. This is the number of restart files to be written.

* **t_rst**  (default  ``t_rst=-1.0 -1.0 ...``) is real array, which contains the times when restart files are requested.

* **dt_rst** (default ``dt_rst=0.0``) is a real, which defines the time interval between restart files.

* **t_rst_start** (default ``t_rst_start=0.0``) is a real, which defines the time to start writing restart files.

* **t_rst_stop** (default ``t_rst_stop=0.0``) is a real, which defines the time to stop writing restart files.

* **n_stores** (default ``n_stores=0``) is an integer. This is another way of requesting a certain number of restart files. However, instead of creating each time a new restart file, if ``n_stores > n_rsts``  the restart file is overwritten, which can possibly help saving some disk space.

.. warning:: The ``rst`` files can become quite big and writting them too frequently will slow down the code. Except for very special use, the default set up should be sufficient.


Graphic files
-------------

The **graph** outputs controls the output of graphic files (``G_#.TAG`` and ``G_t_#.TAG``) which contain a snapshot the entropy, the velocity field and the magnetic fields:

* **n_graph_step** (default ``n_graph_step=0``) is an integer. This is the number of timesteps between two graphic files.

* **n_graphs** (default ``n_graphss=1``) is an integer. This is the number of graphic files to be written.

* **t_graph**  (default  ``t_graph=-1.0 -1.0 ...``) is real array, which contains the times when graphic files are requested.

.. _dtGraphVar:

* **dt_graph** (default ``dt_graph=0.0``) is a real, which defines the time interval between graphic files.

* **t_graph_start** (default ``t_graph_start=0.0``) is a real, which defines the time to start writing graphic files.

* **t_graph_stop** (default ``t_graph_stop=0.0``) is a real, which defines the time to stop writing graphic files.



Spectra
-------

The **spec** outputs controls the output of spectra: kinetic energy spectra (``kin_spec_#.TAG``), magnetic energy spectra (``mag_spec_#.TAG``) and thermal spectra (``T_spec_#.TAG``):

* **n_spec_step** (default ``n_spec_step=0``) is an integer. This is the number of timesteps between two spectra.

* **n_specs** (default ``n_specs=0``) is an integer. This is the number of spectra to be written.

* **t_spec**  (default  ``t_spec=-1.0 -1.0 ...``) is real array, which contains the times when spectra are requested.

* **dt_spec** (default ``dt_spec=0.0``) is a real, which defines the time interval between spectra.

* **t_spec_start** (default ``t_spec_start=0.0``) is a real, which defines the time to start writing spectra.

* **t_spec_stop** (default ``t_spec_stop=0.0``) is a real, which defines the time to stop writing spectra.


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

The **cmb** outputs controls the output of the potential at several depths: ``B_coeff_r*.TAG``, ``V_coeff_r*.TAG`` and ``T_coeff_r*.TAG`` are produced.

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



Movie files
-----------


The **movie** outputs controls the output of movie files (``*_mov.TAG``). 

.. note:: This calculation is **only** enabled when ``l_movie=.true.``.

Specific inputs
+++++++++++++++

* **l_movie** (default ``l_movie=.false.``) is a logical. It needs to be turned on to get movie computed.

* **movie** (default ``movie=' ', ' ', ...``) is a character string array. It contains the name of the movies one wants to compute.

Standard inputs
+++++++++++++++

* **n_movie_step** (default ``n_movie_step=0``) is an integer. This is the number of timesteps between two movie outputs.

* **n_movies** (default ``n_moviess=1``) is an integer. This is the number of movie outputs to be written.

* **t_movie**  (default  ``t_movie=-1.0 -1.0 ...``) is real array, which contains the times when movie outputs are requested.

* **dt_movie** (default ``dt_movie=0.0``) is a real, which defines the time interval between movie outputs.

* **t_movie_start** (default ``t_movie_start=0.0``) is a real, which defines the time to start writing movie outputs.

* **t_movie_stop** (default ``t_movie_stop=0.0``) is a real, which defines the time to stop writing movie outputs.



Field Averages
--------------

The code can perform on-the-fly time-averaging of entropy, velocity field and magnetic field. Respective graphic output and spectra are written into the corresponding files (with ``G_ave.TAG``, ``kin_spec_avec.TAG``). The time-averaged energies are written into the log file.

.. _varl_average:

* **l_average** (default ``l_average=.false.``) is a logical, which enables the time-averaging of fields when set to ``.true.``.

  .. warning:: Time-averaging has a large memory imprint as it requires the storage of 3-D arrays. Be careful, when using large truncations.

RMS force balance
-----------------

  .. warning:: The RMS calculation is actually wrong in the current version. This needs again to be ported from MagIC 3.44. The RMS contributions to the induction equation are correct, though. A ticket has been opened on github regarding this issue: https://github.com/magic-sph/magic/issues/1

The code can compute the RMS of the force balance and the induction equation.

* **l_RMS** (default ``l_RMS=.false.``) is a logical, which enables the calculation of RMS force balance, when set to ``.true.``.

* **l_RMStest** (default ``l_RMStest=.false.``) is a logical. This is a debug flag to check the consistency of the RMS calculation.

Torsional oscillations
----------------------


Miscellaneous
-------------

Helicity
++++++++

* **l_hel** (default ``l_hel=.false.``) is a logical. When set to ``.true.``, this logical enables the calculation of helicity (RMS, northern and southern hemisphere, etc.). The outputs are stored in the :ref:`misc.TAG <secMiscFile>` file.

.. _varl_power:

Power budget
++++++++++++

* **l_power** (default ``l_power.false.``) is a logical. When set to ``.true.``, this logical enables the calculation if input and output power (buoyancy, viscous and ohmic dissipations, torques). The time series are stored in ``power.TAG`` and the time-averaged radial profiles in :ref:`powerR.TAG <secPowerRfile>`.
