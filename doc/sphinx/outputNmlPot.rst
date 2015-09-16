.. _secOutNmlPot:

Flow poloidal and toroidal potentials in spectral space
-------------------------------------------------------

Specific inputs
+++++++++++++++

* **l_storeVpot** (default ``l_storeVpot=.false.``) is a logical. It needs to be turned on to store the flow poloidal and toroidal potentials. It then writes the following file: ``Vpot_#.TAG``.

Standard inputs
+++++++++++++++

* **n_Vpot_step** (default ``n_Vpot_step=0``) is an integer. This is the number of timesteps between two ``Vpot`` outputs.

* **n_Vpots** (default ``n_Vpots=1``) is an integer. This is the number of ``Vpot`` outputs to be written.

* **t_Vpot**  (default  ``t_Vpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Vpot`` outputs are requested.

* **dt_Vpot** (default ``dt_Vpot=0.0``) is a real, which defines the time interval between ``Vpot`` outputs.

* **t_Vpot_start** (default ``t_Vpot_start=0.0``) is a real, which defines the time to start writing ``Vpot`` outputs.

* **t_Vpot_stop** (default ``t_Vpot_stop=0.0``) is a real, which defines the time to stop writing ``Vpot`` outputs.


Magnetic field poloidal and toroidal potentials in spectral space
-----------------------------------------------------------------

Specific inputs
+++++++++++++++

* **l_storeBpot** (default ``l_storeBpot=.false.``) is a logical. It needs to be turned on to store the magnetic field poloidal and toroidal potentials. It then writes the following file: ``Bpot_#.TAG``.

Standard inputs
+++++++++++++++

* **n_Bpot_step** (default ``n_Bpot_step=0``) is an integer. This is the number of timesteps between two ``Bpot`` outputs.

* **n_Bpots** (default ``n_Bpots=1``) is an integer. This is the number of ``Bpot`` outputs to be written.

* **t_Bpot**  (default  ``t_Bpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Bpot`` outputs are requested.

* **dt_Bpot** (default ``dt_Bpot=0.0``) is a real, which defines the time interval between ``Bpot`` outputs.

* **t_Bpot_start** (default ``t_Bpot_start=0.0``) is a real, which defines the time to start writing ``Bpot`` outputs.

* **t_Bpot_stop** (default ``t_Bpot_stop=0.0``) is a real, which defines the time to stop writing ``Bpot`` outputs.

Entropy/temperature in spectral space
-------------------------------------

Specific inputs
+++++++++++++++

* **l_storeTpot** (default ``l_storeTpot=.false.``) is a logical. It needs to be turned on to store the entropy. It then writes the following file: ``Tpot_#.TAG``.

Standard inputs
+++++++++++++++

* **n_Tpot_step** (default ``n_Tpot_step=0``) is an integer. This is the number of timesteps between two ``Tpot`` outputs.

* **n_Tpots** (default ``n_Tpots=1``) is an integer. This is the number of ``Tpot`` outputs to be written.

* **t_Tpot**  (default  ``t_Tpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Tpot`` outputs are requested.

* **dt_Tpot** (default ``dt_Tpot=0.0``) is a real, which defines the time interval between ``Tpot`` outputs.

* **t_Tpot_start** (default ``t_Tpot_start=0.0``) is a real, which defines the time to start writing ``Tpot`` outputs.

* **t_Tpot_stop** (default ``t_Tpot_stop=0.0``) is a real, which defines the time to stop writing ``Tpot`` outputs.

