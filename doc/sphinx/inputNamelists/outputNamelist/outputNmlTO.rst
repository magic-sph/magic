.. _secOutNmlTO:

Torsional oscillations (``TO``)
-------------------------------

Specific inputs
+++++++++++++++

* **l_TO** (default ``l_TO=.false.``) is a logical. It needs to be turned on to compute the torsional oscillations diagnostics (``TO``) computed.


Standard inputs
+++++++++++++++

* **n_TO_step** (default ``n_TO_step=0``) is an integer. This is the number of timesteps between two ``TO`` outputs.

* **n_TOs** (default ``n_TOs=1``) is an integer. This is the number of ``TO`` outputs to be written.

* **t_TO**  (default  ``t_TO=-1.0 -1.0 ...``) is real array, which contains the times when ``TO`` outputs are requested.

* **dt_TO** (default ``dt_TO=0.0``) is a real, which defines the time interval between ``TO`` outputs.

* **t_TO_start** (default ``t_TO_start=0.0``) is a real, which defines the time to start writing ``TO`` outputs.

* **t_TO_stop** (default ``t_TO_stop=0.0``) is a real, which defines the time to stop writing ``TO`` outputs.


