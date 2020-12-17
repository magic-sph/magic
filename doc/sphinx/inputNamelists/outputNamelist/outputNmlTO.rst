.. _secOutNmlTO:

Torsional oscillations (``TO``)
-------------------------------

Specific inputs
+++++++++++++++

.. _varl_TO:

* **l_TO** (default :f:var:`l_TO=.false. <l_to>`) is a logical. It needs to be turned on to compute the torsional oscillations diagnostics (``TO``) computed.

.. _varl_TOmovie:

* **l_TOmovie** (default :f:var:`l_TOmovie=.false <l_tomovie>`) is a logical. It needs to be turned on to store the :ref:`TO_movie.TAG <secTO_movieFile>` files.

.. _varsDens:

* **sDens** (default :f:var:`sDens=1.0 <sdens>`) is a float. It gives the relative point density of the cylindrical grid (in the radial direction).

.. _varzDens:

* **zDens** (default :f:var:`zDens=1.0 <zdens>`) is a float. It gives the relative point density of the cylindrical grid (in the vertical direction).


Standard inputs
+++++++++++++++

* **n_TO_step** (default :f:var:`n_TO_step=0 <n_to_step>`) is an integer. This is the number of timesteps between two ``TO`` outputs.

* **n_TOs** (default :f:var:`n_TOs=1 <n_tos>`) is an integer. This is the number of ``TO`` outputs to be written.

* **t_TO**  (default  :f:var:`t_TO=-1.0 -1.0 ...<t_to>`) is real array, which contains the times when ``TO`` outputs are requested.

* **dt_TO** (default :f:var:`dt_TO=0.0 <dt_to>`) is a real, which defines the time interval between ``TO`` outputs.

* **t_TO_start** (default :f:var:`t_TO_start=0.0 <t_to_start>`) is a real, which defines the time to start writing ``TO`` outputs.

* **t_TO_stop** (default :f:var:`t_TO_stop=0.0 <t_to_stop>`) is a real, which defines the time to stop writing ``TO`` outputs.

* **n_TOmovie_step** (default :f:var:`n_TOmovie_step=0 <n_tomovie_step>`) is an integer. This is the number of timesteps between two ``TO_mov`` outputs.

* **n_TOmovie_frames** (default :f:var:`n_TOmovies=1 <n_tomovie_frames>`) is an integer. This is the number of ``TO_mov`` outputs to be written.

* **t_TOmovie**  (default  :f:var:`t_TOmovie=-1.0 -1.0 ... <t_tomovie>`) is real array, which contains the times when ``TO_mov`` outputs are requested.

* **dt_TOmovie** (default ``dt_TOmovie=0.0``) is a real, which defines the time interval between ``TO_mov`` outputs.

* **t_TOmovie_start** (default :f:var:`t_TOmovie_start=0.0 <t_tomovie_start>`) is a real, which defines the time to start writing ``TO_mov`` outputs.

* **t_TOmovie_stop** (default :f:var:`t_TOmovie_stop=0.0 <t_tomovie_stop>`) is a real, which defines the time to stop writing ``TO_mov`` outputs.

