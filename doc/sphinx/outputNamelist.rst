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

   .. note:: These input variables are usually named with a pattern that follows ``dt_outputName``, for instance, ``dt_graph``, ``dt_rst``, ``dt_spec``, ``dt_log``, ``dt_movie``, etc.

   ..

   4. **User-defined times for output**. By default 5000 different times can be defined for each output type. This can be increased by increasing n_time_hits in the file ``output_data.f90``. While the first three possibilities can only be used alternatively, the fourth one can be employed in addition to one of the two others.

   .. note:: These input variables are usually named with a pattern that follows ``t_outputName``, for instance, ``t_graph``, ``t_rst``, ``t_spec``, ``t_log``, ``t_movie``, etc.

   ..


An important parameter in this context is ``l_true_time``. If this is set to ``.true.``, the time steps of the program are modified to meet a desired output time. This forces a recalculation of the inversion matricies and therefore requires some additional computing time. When ``l_true_time=.false.``, the values at the timestep closest to the desired output time are chosen. Since the timesteps are generally small, this option suffices for most applications.

* **l_true_time** (default ``l_true_time=.false.``) is a logical. It causes the code to change time steps to exactly meet the requested output times.

* **l_save_out** (default ``l_save_out=.false.``) is a logical. When set to ``.true.``, the diagnostic files will be safely opened and closed before and after any outputs. When set to ``.false.``, the diagnostic files will be opened before the first iteration timestep and close at the end of the run.

* **lVerbose** (default ``lVerbose=.false.``) is a logical. When set to ``.true.``, the code displays a lot of debugging informations.

  .. warning:: Never set ``lVerbose`` to ``.true.`` for a production run!


Standard time-series outputs
----------------------------

The **log** outputs controls the output of all the default time series of the file: kinetic and magnetic energies (``e_kin.TAG``, ``e_mag_oc.TAG`` and ``e_mag_ic.TAG`` files), dipole information (``dipole.TAG`` file), parameters (``par.TAG``) and various additional diagnostics (``misc.TAG``):

* **n_log_step** (default ``n_log_step=50``) is an integer. This is the number of timesteps between two log outputs.

  .. warning:: Be careful: when using too small ``n_log_step``, the disk access will dramatically increases, thus decreasing the code performance.

* **n_logs** (default ``n_logs=0``) is an integer. This is the number of log-information sets to be written.

* **t_log**  (default  ``t_log=0.0 0.0``) is real array, which contains the times when log outputs are requested.

* **dt_log** (default ``dt_log=0.0``) is a real, which defines the time interval between log outputs.

* **t_log_start** (default ``t_log_start=0.0``) is a real, which defines the time to start writing log outputs.

* **t_log_stop** (default ``t_log_stop=0.0``) is a real, which defines the time to stop writing log outputs.

Restart files
-------------

The **rst** outputs controls the output of restart files (``rst_t_#.TAG``):

* **n_rst_step** (default ``n_rst_step=0``) is an integer. This is the number of timesteps between two restart files.

* **n_rsts** (default ``n_rsts=1``) is an integer. This is the number of restart files to be written.

* **t_rst**  (default  ``t_rst=0.0 0.0``) is real array, which contains the times when restart files are requested.

* **dt_rst** (default ``dt_rst=0.0``) is a real, which defines the time interval between restart files.

* **t_rst_start** (default ``t_rst_start=0.0``) is a real, which defines the time to start writing restart files.

* **t_spec_stop** (default ``t_rst_stop=0.0``) is a real, which defines the time to stop writing restart files.

.. warning:: The ``rst`` files can become quite big and writting them too frequently will slow down the code. Except for very special use, the default set up should be sufficient.


Spectra
-------

The **spec** outputs controls the output of spectra: kinetic energy spectra (``kin_spec_#.TAG``), magnetic energy spectra (``mag_spec_#.TAG``) and thermal spectra (``T_spec_#.TAG``):

* **n_spec_step** (default ``n_spec_step=0``) is an integer. This is the number of timesteps between two spectra.

* **n_specs** (default ``n_specs=0``) is an integer. This is the number of spectra to be written.

* **t_spec**  (default  ``t_spec=0.0 0.0``) is real array, which contains the times when spectra are requested.

* **dt_spec** (default ``dt_spec=0.0``) is a real, which defines the time interval between spectra.

* **t_spec_start** (default ``t_spec_start=0.0``) is a real, which defines the time to start writing spectra.

* **t_spec_stop** (default ``t_spec_stop=0.0``) is a real, which defines the time to stop writing spectra.

