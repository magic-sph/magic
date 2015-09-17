.. _secOutNmlTO:

Torsional oscillations (``TO``)
-------------------------------

Specific inputs
+++++++++++++++

* **l_TO** (default ``l_TO=.false.``) is a logical. It needs to be turned on to compute the torsional oscillations diagnostics (``TO``) computed.

* **l_TOmovie** (default ``l_TOmovie=.false``) is a logical. It needs to be turned on to store the :ref:`TO_movie.TAG <secTO_movieFile>` files.

* **sDens** (default ``sDens=1.0``) is a float. It gives the relative point density of the cylindrical grid (in the radial direction).

* **zDens** (default ``zDens=1.0``) is a float. It gives the relative point density of the cylindrical grid (in the vertical direction).


Standard inputs
+++++++++++++++

* **n_TO_step** (default ``n_TO_step=0``) is an integer. This is the number of timesteps between two ``TO`` outputs.

* **n_TOs** (default ``n_TOs=1``) is an integer. This is the number of ``TO`` outputs to be written.

* **t_TO**  (default  ``t_TO=-1.0 -1.0 ...``) is real array, which contains the times when ``TO`` outputs are requested.

* **dt_TO** (default ``dt_TO=0.0``) is a real, which defines the time interval between ``TO`` outputs.

* **t_TO_start** (default ``t_TO_start=0.0``) is a real, which defines the time to start writing ``TO`` outputs.

* **t_TO_stop** (default ``t_TO_stop=0.0``) is a real, which defines the time to stop writing ``TO`` outputs.


* **n_TOZ_step** (default ``n_TOZ_step=0``) is an integer. This is the number of timesteps between two ``TO`` outputs.

* **n_TOZs** (default ``n_TOZs=1``) is an integer. This is the number of ``TO`` outputs to be written.

* **t_TOZ**  (default  ``t_TOZ=-1.0 -1.0 ...``) is real array, which contains the times when ``TO`` outputs are requested.

* **dt_TOZ** (default ``dt_TOZ=0.0``) is a real, which defines the time interval between ``TO`` outputs.

* **t_TOZ_start** (default ``t_TOZ_start=0.0``) is a real, which defines the time to start writing ``TO`` outputs.

* **t_TOZ_stop** (default ``t_TOZ_stop=0.0``) is a real, which defines the time to stop writing ``TO`` outputs.


* **n_TOmovie_step** (default ``n_TOmovie_step=0``) is an integer. This is the number of timesteps between two ``TO_mov`` outputs.

* **n_TOmovies** (default ``n_TOmovies=1``) is an integer. This is the number of ``TO_mov`` outputs to be written.

* **t_TOmovie**  (default  ``t_TOmovie=-1.0 -1.0 ...``) is real array, which contains the times when ``TO_mov`` outputs are requested.

* **dt_TOmovie** (default ``dt_TOmovie=0.0``) is a real, which defines the time interval between ``TO_mov`` outputs.

* **t_TOmovie_start** (default ``t_TOmovie_start=0.0``) is a real, which defines the time to start writing ``TO_mov`` outputs.

* **t_TOmovie_stop** (default ``t_TOmovie_stop=0.0``) is a real, which defines the time to stop writing ``TO_mov`` outputs.

