.. _secOutNmlPot:

Poloidal and toroidal potentials in spectral and Chebyshev space
----------------------------------------------------------------

The **[V|B|T]pot** outputs controls the output of potential files
(``Vpot_#.TAG``, ``Bpot_#.TAG`` and ``Tpot_#.TAG``). These are files that contain
the poloidal and toroidal flow and magnetic field potentials (and entropy/temperature)
written in spectral and Chebyshev spaces (for instance ``w(lm_max, n_cheb_max)``).
These files can be pretty useful since they can be possibly used to reconstruct any
quantity in the spectral space or in the physical space you may be interested in.

Specific inputs
+++++++++++++++

They are two ways to store those files. The first option is to use
``l_storePot=.true.`` and the corresponding time control parameters (
:ref:`n_pot_step <varn_pot_step>`, :ref:`t_pot <vart_pot>`, :ref:`n_pots
<varn_pots>`, etc.). In that case the three files ``Vpot_#.TAG``,
``Bpot_#.TAG`` and ``Tpot_#.TAG`` will be stored. The following example
will create new ``Vpot_#.TAG``, ``Bpot_#.TAG`` and ``Tpot_#.TAG`` files
every 1000 time steps:

  .. code-block :: fortran

     l_storePot = .true.,
     n_pot_step = 1000, 

* **l_storePot** (default ``l_storePot=.false.``) is a logical. It needs to be turned on to store all the potentials in three different files: ``Vpot_#.TAG``, ``Bpot_#.TAG`` and
``Tpot_#.TAG``.

The second option is control separately the writing of the three files using
the three logicals ``l_storeVpot``, ``l_storeBpot`` and ``l_storeTpot`` and their
corresponding time control parameters. The following example wrill create a new
``Vpot_#.TAG`` file every 1000 time steps and a new ``Bpot_#.TAG`` file every
3000 time steps (no ``Tpot_#.TAG`` files are stored in that case):

  .. code-block:: fortran

     l_storeVpot = .true.,
     n_Vpot_step = 1000, 
     l_storeBpot = .true.,
     n_Bpot_step = 3000, 
     l_storeTpot = .false.,

* **l_storeVpot** (default ``l_storeVpot=.false.``) is a logical. It needs to be turned on to store the flow poloidal and toroidal potentials. It then writes the ``Vpot_#.TAG`` file.

* **l_storeBpot** (default ``l_storeBpot=.false.``) is a logical. It needs to be turned on to store the magnetic field poloidal and toroidal potentials. It then writes the  ``Bpot_#.TAG`` file.

* **l_storeTpot** (default ``l_storeTpot=.false.``) is a logical. It needs to be turned on to store the entropy. It then writes the ``Tpot_#.TAG`` file.

Standard inputs
+++++++++++++++

.. _varn_pot_step:

* **n_pot_step** (default ``n_pot_step=0``) is an integer. This is the number of timesteps between two  ``[V|B|P]pot`` outputs.

.. _varn_pots:

* **n_pots** (default ``n_pots=1``) is an integer. This is the number of ``[V|B|P]pot`` outputs to be written.

.. _vart_pot:

* **t_pot**  (default  ``t_pot=-1.0 -1.0 ...``) is real array, which contains the times when  ``[V|B|P]pot`` outputs are requested.

* **dt_pot** (default ``dt_pot=0.0``) is a real, which defines the time interval between two ``[V|B|P]pot`` outputs.

* **t_pot_start** (default ``t_pot_start=0.0``) is a real, which defines the time to start writing ``[V|B|P]pot`` outputs.

* **t_pot_stop** (default ``t_pot_stop=0.0``) is a real, which defines the time to stop writing ``[V|B|P]pot`` outputs.

* **n_Vpot_step** (default ``n_Vpot_step=0``) is an integer. This is the number of timesteps between two ``Vpot`` outputs.

* **n_Vpots** (default ``n_Vpots=1``) is an integer. This is the number of ``Vpot`` outputs to be written.

* **t_Vpot**  (default  ``t_Vpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Vpot`` outputs are requested.

* **dt_Vpot** (default ``dt_Vpot=0.0``) is a real, which defines the time interval between ``Vpot`` outputs.

* **t_Vpot_start** (default ``t_Vpot_start=0.0``) is a real, which defines the time to start writing ``Vpot`` outputs.

* **t_Vpot_stop** (default ``t_Vpot_stop=0.0``) is a real, which defines the time to stop writing ``Vpot`` outputs.

* **n_Bpot_step** (default ``n_Bpot_step=0``) is an integer. This is the number of timesteps between two ``Bpot`` outputs.

* **n_Bpots** (default ``n_Bpots=1``) is an integer. This is the number of ``Bpot`` outputs to be written.

* **t_Bpot**  (default  ``t_Bpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Bpot`` outputs are requested.

* **dt_Bpot** (default ``dt_Bpot=0.0``) is a real, which defines the time interval between ``Bpot`` outputs.

* **t_Bpot_start** (default ``t_Bpot_start=0.0``) is a real, which defines the time to start writing ``Bpot`` outputs.

* **t_Bpot_stop** (default ``t_Bpot_stop=0.0``) is a real, which defines the time to stop writing ``Bpot`` outputs.

* **n_Tpot_step** (default ``n_Tpot_step=0``) is an integer. This is the number of timesteps between two ``Tpot`` outputs.

* **n_Tpots** (default ``n_Tpots=1``) is an integer. This is the number of ``Tpot`` outputs to be written.

* **t_Tpot**  (default  ``t_Tpot=-1.0 -1.0 ...``) is real array, which contains the times when ``Tpot`` outputs are requested.

* **dt_Tpot** (default ``dt_Tpot=0.0``) is a real, which defines the time interval between ``Tpot`` outputs.

* **t_Tpot_start** (default ``t_Tpot_start=0.0``) is a real, which defines the time to start writing ``Tpot`` outputs.

* **t_Tpot_stop** (default ``t_Tpot_stop=0.0``) is a real, which defines the time to stop writing ``Tpot`` outputs.

