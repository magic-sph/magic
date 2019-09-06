.. _secOutNmlPot:

Poloidal and toroidal potentials in spectral and radial space
----------------------------------------------------------------

The **[V|B|T]_lmr** outputs controls the output of potential files
(:ref:`V_lmr_#.TAG <secPotFiles>`, :ref:`B_lmr_#.TAG <secPotFiles>` and 
:ref:`T_lmr_#.TAG <secPotFiles>`). These are files contain
the poloidal and toroidal flow and magnetic field potentials (and entropy/temperature)
written in spectral and radial spaces (for instance :code:`w(lm_max, n_r_max)`).
These files can be quite handy since they can be possibly used to reconstruct any
quantity in the spectral space or in the physical space you may be interested in.


Standard inputs
+++++++++++++++

.. _varn_pot_step:

* **n_pot_step** (default :f:var:`n_pot_step=0 <n_pot_step>`) is an integer. This is the number of timesteps between two  ``[V|B|T|Xi]_lmr`` outputs.

.. _varn_pots:

* **n_pots** (default :f:var:`n_pots=1 <n_pots>`) is an integer. This is the number of ``[V|B|T|Xi]_lmr`` outputs to be written.

.. _vart_pot:

* **t_pot**  (default  :f:var:`t_pot=-1.0 -1.0 ... <t_pot>`) is real array, which contains the times when  ``[V|B|T|Xi]_lmr`` outputs are requested.

* **dt_pot** (default :f:var:`dt_pot=0.0 <dt_pot>`) is a real, which defines the time interval between two ``[V|B|T|Xi]_lmr`` outputs.

* **t_pot_start** (default :f:var:`t_pot_start=0.0 <t_pot_start>`) is a real, which defines the time to start writing ``[V|B|T|Xi]_lmr`` outputs.

* **t_pot_stop** (default :f:var:`t_pot_stop=0.0 <t_pot_stop>`) is a real, which defines the time to stop writing ``[V|B|T|Xi]_lmr`` outputs.
