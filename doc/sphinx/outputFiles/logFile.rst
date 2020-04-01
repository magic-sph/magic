.. _secLogFile:

Log file: ``log.TAG``
=====================

This is a text file contains information about the run, including many of the
things which are printed to ``STDOUT``. It has the following information in order
of appearance:

   * **Code version**: the version of the code
   
   * **Parallelization**: information about number of MPI ranks being used,
     blocking information of OpenMP chunks and processor load balancing
   
   * **Namelists**: displays values of all namelist variables. The ones input by
     the user should have the input values while the rest of them are set to their
     default values.
   
   * **Mode** The mode of the run - self-consistent/kinematic dynamo, convection,
     couette flow etc. See the :ref:`control namelist <secControlNml>` for more
     information about :ref:`mode <varmode>`.
   
   * **Grid parameters**: information about the grid sizes and truncation being
     used. More information about this in the :ref:`grid namelist <secGridNml>`.
     If a new grid, different from that in the restart file is used, then a
     comparison is shown between old and new grid parameters and the user is
     informed that the data is being mapped from the old to the new grid.
   
   * **Progress**: information about the progress of the run for every 10% of the
     run and the mean wall time for time step.
   
   * **Writing of graphic, movie, restart and spectra files**: displays the time
     step and tells the user whenever a :ref:`G_#.TAG <secGraphFile>`,
     :ref:`checkpoint_#.TAG <secRestartFile>` or :ref:`spectra <secSpecFiles>` file or a
     :ref:`movie frame <secMovieFile>` is written disk.
   
   * **Energies**: gives kinetic and magnetic energies (total, poloidal, toroidal,
     total density) at the end of the run.
   
   * **Time averages**:  this part gives time averaged kinetic and magnetic
     energies (total, poloidal, toroidal, total density) and time averaged
     parameters (Rm, Elsass, Rol etc.). If :ref:`l_spec_avg=.true. <varl_spec_avg>`,
     this section also provides information about average spectra being written.
     If :ref:`l_average=.true. <varl_average>`, it is additionally mentioned that 
     time averaged graphic files are written.
   
   * **Wall times**: this is the last part of the log file and it provides
     information about the mean wall time for running different parts of the code.
     These values can be used to judge the speed and scaling capabilities of
     your computer.


Most of these informations can be parsed and stored into a python class using
:py:class:`MagicSetup <magic.MagicSetup>`:

>>> # read log.N0m2
>>> stp = MagicSetup(nml='log.N0m2')
>>> print(stp.ek, stp.prmag) # print Ekman and magnetic Prandtl numbers
>>> print(stp.l_max) # print l_max
