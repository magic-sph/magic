.. _secOutNmlStd:

Standard time-series outputs
----------------------------

The **log** outputs controls the output of all the default time series of the
file: kinetic and magnetic energies (:ref:`e_kin.TAG <secEkinFile>`,
:ref:`e_mag_oc.TAG <secEmagocFile>` and :ref:`e_mag_ic.TAG <secEmagicFile>`
files), dipole information (:ref:`dipole.TAG <secDipoleFile>` file), rotation
(:ref:`rot.TAG <secRotFile>`) parameters (:ref:`par.TAG <secParFile>`) and
various additional diagnostics (:ref:`misc.TAG <secMiscFile>`):

.. _varn_log_step:

* **n_log_step** (default ``n_log_step=50``) is an integer. This is the number of timesteps between two log outputs.

  .. warning:: Be careful: when using too small ``n_log_step``, the disk access will dramatically increases, thus decreasing the code performance.

.. _varn_logs:

* **n_logs** (default ``n_logs=0``) is an integer. This is the number of log-information sets to be written.

.. _vart_log:

* **t_log**  (default  ``t_log=-1.0 -1.0 ...``) is real array, which contains the times when log outputs are requested.

.. _vardt_log:

* **dt_log** (default ``dt_log=0.0``) is a real, which defines the time interval between log outputs.

.. _vart_log_start:

* **t_log_start** (default ``t_log_start=0.0``) is a real, which defines the time to start writing log outputs.

.. _vart_log_stop:

* **t_log_stop** (default ``t_log_stop=0.0``) is a real, which defines the time to stop writing log outputs.


Restart files
-------------

The **rst** outputs controls the output of restart files (:ref:`rst_t_#.TAG <secRestartFile>`) (i.e. check points in time from which the code could be restarted):

.. _varn_rst_step:

* **n_rst_step** (default ``n_rst_step=0``) is an integer. This is the number of timesteps between two restart files.

.. _varn_rsts:

* **n_rsts** (default ``n_rsts=1``) is an integer. This is the number of restart files to be written.

.. _vart_rst:

* **t_rst**  (default  ``t_rst=-1.0 -1.0 ...``) is real array, which contains the times when restart files are requested.

.. _vardt_rst:

* **dt_rst** (default ``dt_rst=0.0``) is a real, which defines the time interval between restart files.

* **t_rst_start** (default ``t_rst_start=0.0``) is a real, which defines the time to start writing restart files.

* **t_rst_stop** (default ``t_rst_stop=0.0``) is a real, which defines the time to stop writing restart files.

* **n_stores** (default ``n_stores=0``) is an integer. This is another way of requesting a certain number of restart files. However, instead of creating each time a new restart file, if ``n_stores > n_rsts``  the restart file is overwritten, which can possibly help saving some disk space.

.. warning:: The ``rst`` files can become quite big and writting them too frequently will slow down the code. Except for very special use, the default set up should be sufficient.


Graphic files
-------------

The **graph** outputs controls the output of graphic files (:ref:`G_#.TAG <secGraphFile>`) which contain a snapshot the entropy, the velocity field and the magnetic fields:

.. _varn_graph_step:

* **n_graph_step** (default ``n_graph_step=0``) is an integer. This is the number of timesteps between two graphic files.

.. _varn_graphs:

* **n_graphs** (default ``n_graphss=1``) is an integer. This is the number of graphic files to be written.

.. _vart_graph:

* **t_graph**  (default  ``t_graph=-1.0 -1.0 ...``) is real array, which contains the times when graphic files are requested.

.. _vardt_graph:

* **dt_graph** (default ``dt_graph=0.0``) is a real, which defines the time interval between graphic files.

.. _vart_graph_start:

* **t_graph_start** (default ``t_graph_start=0.0``) is a real, which defines the time to start writing graphic files.

.. _vart_graph_stop:

* **t_graph_stop** (default ``t_graph_stop=0.0``) is a real, which defines the time to stop writing graphic files.



Spectra
-------

The **spec** outputs controls the output of spectra: kinetic energy spectra (:ref:`kin_spec_#.TAG <secKinSpecFile`), magnetic energy spectra (:ref:`mag_spec_#.TAG <secMagSpecFile>`) and thermal spectra (:ref:`T_spec_#.TAG <secTSpecFile>`):

.. _varn_spec_step:

* **n_spec_step** (default ``n_spec_step=0``) is an integer. This is the number of timesteps between two spectra.

.. _varn_specs:

* **n_specs** (default ``n_specs=0``) is an integer. This is the number of spectra to be written.

.. _vart_spec:

* **t_spec**  (default  ``t_spec=-1.0 -1.0 ...``) is real array, which contains the times when spectra are requested.

.. _vardt_spec:

* **dt_spec** (default ``dt_spec=0.0``) is a real, which defines the time interval between spectra.

.. _vart_spec_start:

* **t_spec_start** (default ``t_spec_start=0.0``) is a real, which defines the time to start writing spectra.

.. _vart_spec_stop:

* **t_spec_stop** (default ``t_spec_stop=0.0``) is a real, which defines the time to stop writing spectra.


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

.. _varn_movie_step:

* **n_movie_step** (default ``n_movie_step=0``) is an integer. This is the number of timesteps between two movie outputs.

* **n_movies** (default ``n_moviess=1``) is an integer. This is the number of movie outputs to be written.

.. _vart_movie:

* **t_movie**  (default  ``t_movie=-1.0 -1.0 ...``) is real array, which contains the times when movie outputs are requested.

.. _vardt_movie:

* **dt_movie** (default ``dt_movie=0.0``) is a real, which defines the time interval between movie outputs.

* **t_movie_start** (default ``t_movie_start=0.0``) is a real, which defines the time to start writing movie outputs.

* **t_movie_stop** (default ``t_movie_stop=0.0``) is a real, which defines the time to stop writing movie outputs.


Field Averages
--------------

The code can perform on-the-fly time-averaging of entropy, velocity field and magnetic field. Respective graphic output and spectra are written into the corresponding files (with :ref:`G_ave.TAG <secGraphFile>`, :ref:`kin_spec_avec.TAG <secKinSpecFile>`). The time-averaged energies are written into the :ref:`log.TAG <secLogFile>` file.

.. _varl_average:

* **l_average** (default ``l_average=.false.``) is a logical, which enables the time-averaging of fields when set to ``.true.``.

  .. warning:: Time-averaging has a large memory imprint as it requires the storage of 3-D arrays. Be careful, when using large truncations.

