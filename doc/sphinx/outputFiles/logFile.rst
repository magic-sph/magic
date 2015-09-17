.. _secLogFile:

Log file: ``log.TAG``
=====================

This is a text file contains information about the run, including many of the things which are printed to STDOUT. It has the following information in order of appearance:

* **Code version** The version of the code

* **Parallelization** Information about number of MPI ranks being used, blocking information of OpenMP chunks and processor load balancing

* **Namelist** Displays values of all namelist variables. The ones input by the user should have the input values while the rest of them are set to their default values.

* **Mode** The mode of the run - self-consistent/kinematic dynamo, convection, couette flow etc. See the :ref:`control namelist <secControlNml>` for more information about ``mode``.

* **Grid parameters** Information about the grid sizes and truncation being used. More information about this in the :ref:`grid namelist <secGridNml>`. If a new grid, different from that in the restart file is used, then a comparison is shown between old and new grid parameters and the user is informed that the data is being mapped from the old to the new grid.

* **Progress** Information about the progress of the run for every 10% of the run and the mean wall time for time step.

* **Writing of graphic, movie, restart and spectra files** Displays the time step and tells the user whenever a graphic,  restart or spectrum file or a movie frame is written to the disk.

* **Energies** Gives kinetic and magnetic energies (total, poloidal, toroidal, total density) at the end of the run.

* **Time averages**  This part gives time averaged kinetic and magnetic energies (total, poloidal, toroidal, total density) and time averaged parameters (Rm, Elsass, Rol etc.). If ``l_average=.true.``, this section also provides information about average spectra and graphic files being written.

* **Wall times** This is the last part of the log file and it provides information about the mean wall time for running different parts of the code.
