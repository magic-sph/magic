.. _secOutputNml:

Output control namelist
=======================

This namelist contains all the parameters that can be adjusted to control the outputs and diagnostics calculated by the code.

There are four different ways to control at which time step a specific output
should be written. Outputs are generally distributed over the total calculation
intervall unless an output time interval is defined by a start time ``t_start``
and a stop time ``t_stop``. If no ``t_start`` is provided, the start time of
the calculation is used. If no ``t_stop`` is provided or ``t_stop>t_start`` the
total calculation interval is assumed

   1. **Prescribed number of outputs**. The outputs are distributed evenly over the total calculation interval so that the number of timesteps between two outputs is always the same, with the possible exception of the first interval. Last output is written for the last time step, and to compensate the interval before the first output may be longer. However, if ``t_stop`` is provided, the outputs are distributed evenly over the interval [``t_stop``, ``t_start``] with equal times intervals between them. 

   .. note:: These input variables are usually named with a pattern that follows 
             ``n_outputName``, for instance, :ref:`n_graphs <varn_graphs>`, 
	     :ref:`n_rsts <varn_rsts>`, :ref:`n_specs <varn_specs>`, 
	     :ref:`n_logs <varn_logs>`, etc.
             
             In case you want to make use ot a specific time interval, the input
	     variables follow a pattern of the form ``t_outputName_start``, 
	     ``t_outputName_stop``. For instance, 
	     :ref:`t_graph_start <vart_graph_start>`, 
	     :ref:`t_graph_stop <vart_graph_stop>`, :ref:`t_log_start <vart_log_start>`,
	     :ref:`t_log_stop <vart_log_stop>`, :ref:`t_spec_start <vart_spec_start>`,
	     :ref:`t_spec_stop <vart_spec_stop>`, etc.

   ..

   2. **User-defined interval between two outputs, given in number of time steps**. Again the last output is performed at the end of the run and a compensation may take place at the beginning.

   .. note:: These input variables are usually named with a pattern that follows
             ``n_outputName_step``, for instance, 
	     :ref:`n_graph_step <varn_graph_step>`, :ref:`n_rst_step <varn_rst_step>`,
	     :ref:`n_spec_step <varn_spec_step>`, :ref:`n_log_step <varn_log_step>`,
	     :ref:`n_movie_step <varn_movie_step>`, etc.

   ..

   3. **Defined time interval between two outputs**.

   .. note:: These input variables are usually named with a pattern that follows 
             ``dt_outputName``, for instance, :ref:`dt_graph <vardt_Graph>`, 
	     :ref:`dt_rst <vardt_rst>`, :ref:`dt_spec <vardt_spec>`, 
	     :ref:`dt_log <vardt_log>`, :ref:`dt_movie <vardt_movie>`, etc.

   ..

   4. **User-defined times for output**. By default 5000 different times can be defined for each output type. This can be increased by increasing n_time_hits in the file ``output_data.f90``. While the first three possibilities can only be used alternatively, the fourth one can be employed in addition to one of the two others.

   .. note:: These input variables are usually named with a pattern that follows 
             ``t_outputName``, for instance, :ref:`t_graph <vart_graph>`, 
	     :ref:`t_rst <vart_rst>`, :ref:`t_spec <vart_spec>`,
	     :ref:`t_log <vart_log>`, :ref:`t_movie <vart_movie>`, etc.

   ..


An important parameter in this context is :ref:`l_true_time <varl_true_time>`.
If this is set to ``.true.``, the time steps of the program are modified to
meet a desired output time. This forces a recalculation of the inversion
matricies and therefore requires some additional computing time. When
``l_true_time=.false.``, the values at the timestep closest to the desired
output time are chosen. Since the timesteps are generally small, this option
suffices for most applications.


.. _varl_true_time:

* **l_true_time** (default :f:var:`l_true_time=.false. <l_true_time>`) is a logical. It causes the code to change time steps to exactly meet the requested output times.

The different possible outputs control parameters are then extensively described in the following pages:

.. topic:: Possible outputs

      1. :ref:`Control standard/common outputs <secOutNmlStd>`

      2. :ref:`CMB and radial coefficients <secOutNmlCoeff>` 

      3. :ref:`Storage of potentials in spectral space <secOutNmlPot>`

      4. :ref:`Torsional oscillations diagnostics <secOutNmlTO>`

      5. :ref:`Additional possible diagnostics <secOutNmlMisc>`


.. toctree::
   :hidden:

   outputNamelist/outputNmlStd.rst
   outputNamelist/outputNmlCoeff.rst
   outputNamelist/outputNmlPot.rst
   outputNamelist/outputNmlTO.rst
   outputNamelist/outputNmlMisc.rst


Generic options
---------------

* **l_save_out** (default :f:var:`l_save_out=.false. <l_save_out>`) is a logical. When set to ``.true.``, the diagnostic files will be safely opened and closed before and after any outputs. When set to ``.false.``, the diagnostic files will be opened before the first iteration timestep and close at the end of the run. This may cost some computing time, but guarantees that only minimal information is lost in case of a crash.

* **lVerbose** (default :f:var:`lVerbose=.false. <lverbose>`) is a logical. When set to ``.true.``, the code displays a lot of debugging informations.

  .. warning:: Never set :f:var:`lVerbose` to ``.true.`` for a production run!

* **runid** (default, :f:var:`runid="MAGIC default run" <runid>`) is a character string. This can be used to briefly describe your run. This information is then for instance stored in the header of the graphic files.
