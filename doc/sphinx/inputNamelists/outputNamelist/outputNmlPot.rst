.. _secOutNmlPot:

Poloidal and toroidal potentials in spectral and radial space
----------------------------------------------------------------

The **[V|B|T]_lmr** outputs controls the output of potential files
(:ref:`V_lmr_#.TAG <secVpotFile>`, :ref:`B_lmr_#.TAG <secBpotFile>` and 
:ref:`T_lmr_#.TAG <secTpotFile>`). These are files contain
the poloidal and toroidal flow and magnetic field potentials (and entropy/temperature)
written in spectral and radial spaces (for instance :code:`w(lm_max, n_r_max)`).
These files can be pretty useful since they can be possibly used to reconstruct any
quantity in the spectral space or in the physical space you may be interested in.

Specific inputs
+++++++++++++++

They are two ways to store those files. The first option is to use
:f:var:`l_storePot=.true. <l_storepot>` and the corresponding time control
parameters ( :ref:`n_pot_step <varn_pot_step>`, :ref:`t_pot <vart_pot>`,
:ref:`n_pots <varn_pots>`, etc.). In that case the three files :ref:`V_lmr_#.TAG
<secVpotFile>`, :ref:`B_lmr_#.TAG <secBpotFile>` and :ref:`Tpot_#.TAG
<secTpotFile>` will be stored. The following example will create new
:ref:`V_lmr_#.TAG <secVpotFile>`, :ref:`Bpot_#.TAG <secBpotFile>` and
:ref:`T_lmr_#.TAG <secTpotFile>` files every 1000 time steps:

  .. code-block :: fortran

     l_storePot = .true.,
     n_pot_step = 1000, 

.. _varl_storePot:

* **l_storePot** (default :f:var:`l_storePot=.false. <l_storepot>`) is a
  logical. It needs to be turned on to store all the potentials in three
  different files: :ref:`V_lmr_#.TAG <secVpotFile>`, :ref:`B_lmr_#.TAG <secBpotFile>` 
  and :ref:`T_lmr_#.TAG <secTpotFile>`.

The second option is control separately the writing of the three files using
the three logicals :f:var:`l_storeVpot <l_storevpot>`, :f:var:`l_storeBpot
<l_storebpot>` and :f:var:`l_storeTpot <l_storetpot>` and their corresponding
time control parameters. The following example wrill create a new
:ref:`V_lmr_#.TAG <secVpotFile>` file every 1000 time steps and a new
:ref:`B_lmr_#.TAG <secBpotFile>` file every 3000 time steps (no :ref:`T_lmr_#.TAG
<secTpotFile>` files are stored in that case):

  .. code-block:: fortran

     l_storeVpot = .true.,
     n_Vpot_step = 1000, 
     l_storeBpot = .true.,
     n_Bpot_step = 3000, 
     l_storeTpot = .false.,

.. _varl_storeVpot:


* **l_storeVpot** (default :f:var:`l_storeVpot=.false. <l_storevpot>`) is a
  logical. It needs to be turned on to store the flow poloidal and toroidal
  potentials. It then writes the :ref:`V_lmr_#.TAG <secVpotFile>` file.

.. _varl_storeBpot:

* **l_storeBpot** (default :f:var:`l_storeBpot=.false. <l_storebpot>`) is a
  logical. It needs to be turned on to store the magnetic field poloidal and
  toroidal potentials. It then writes the :ref:`B_lmr_#.TAG <secBpotFile>` file.

.. _varl_storeTpot:

* **l_storeTpot** (default :f:var:`l_storeTpot=.false. <l_storetpot>`) is a
  logical. It needs to be turned on to store the entropy. It then writes the
  :ref:`T_lmr_#.TAG <secTpotFile>` file.

Standard inputs
+++++++++++++++

.. _varn_pot_step:

* **n_pot_step** (default :f:var:`n_pot_step=0 <n_pot_step>`) is an integer. This is the number of timesteps between two  ``[V|B|P]_lmr`` outputs.

.. _varn_pots:

* **n_pots** (default :f:var:`n_pots=1 <n_pots>`) is an integer. This is the number of ``[V|B|P]_lmr`` outputs to be written.

.. _vart_pot:

* **t_pot**  (default  :f:var:`t_pot=-1.0 -1.0 ... <t_pot>`) is real array, which contains the times when  ``[V|B|P]_lmr`` outputs are requested.

* **dt_pot** (default :f:var:`dt_pot=0.0 <dt_pot>`) is a real, which defines the time interval between two ``[V|B|P]_lmr`` outputs.

* **t_pot_start** (default :f:var:`t_pot_start=0.0 <t_pot_start>`) is a real, which defines the time to start writing ``[V|B|P]_lmr`` outputs.

* **t_pot_stop** (default :f:var:`t_pot_stop=0.0 <t_pot_stop>`) is a real, which defines the time to stop writing ``[V|B|P]_lmr`` outputs.

* **n_Vpot_step** (default :f:var:`n_Vpot_step=0 <n_vpot_step>`) is an integer. This is the number of timesteps between two ``V_lmr`` outputs.

* **n_Vpots** (default :f:var:`n_Vpots=1 <n_vpots>`) is an integer. This is the number of ``Vlmr`` outputs to be written.

* **t_Vpot**  (default  :f:var:`t_Vpot=-1.0 -1.0 ...<t_vpot>`) is real array, which contains the times when ``V_lmr`` outputs are requested.

* **dt_Vpot** (default :f:var:`dt_Vpot=0.0 <dt_vpot>`) is a real, which defines the time interval between ``V_lmr`` outputs.

* **t_Vpot_start** (default :f:var:`t_Vpot_start=0.0 <t_vpot_start>`) is a real, which defines the time to start writing ``V_lmr`` outputs.

* **t_Vpot_stop** (default :f:var:`t_Vpot_stop=0.0 <t_vpot_stop>`) is a real, which defines the time to stop writing ``V_lmr`` outputs.

* **n_Bpot_step** (default :f:var:`n_Bpot_step=0 <n_bpot_step>`) is an integer. This is the number of timesteps between two ``B_lmr`` outputs.

* **n_Bpots** (default :f:var:`n_Bpots=1 <n_bpots>`) is an integer. This is the number of ``B_lmr`` outputs to be written.

* **t_Bpot**  (default  :f:var:`t_Bpot=-1.0 -1.0 ... <t_bpot>`) is real array, which contains the times when ``B_lmr`` outputs are requested.

* **dt_Bpot** (default :f:var:`dt_Bpot=0.0 <dt_bpot>`) is a real, which defines the time interval between ``B_lmr`` outputs.

* **t_Bpot_start** (default :f:var:`t_Bpot_start=0.0 <t_bpot_start>`) is a real, which defines the time to start writing ``B_lmr`` outputs.

* **t_Bpot_stop** (default :f:var:`t_Bpot_stop=0.0 <t_bpot_stop>`) is a real, which defines the time to stop writing ``B_lmr`` outputs.

* **n_Tpot_step** (default :f:var:`n_Tpot_step=0 <n_tpot_step>`) is an integer. This is the number of timesteps between two ``T_lmr`` outputs.

* **n_Tpots** (default :f:var:`n_Tpots=1 <n_tpots>`) is an integer. This is the number of ``T_lmr`` outputs to be written.

* **t_Tpot**  (default  :f:var:`t_Tpot=-1.0 -1.0 ... <t_tpot>`) is real array, which contains the times when ``T_lmr`` outputs are requested.

* **dt_Tpot** (default :f:var:`dt_Tpot=0.0 <dt_tpot>`) is a real, which defines the time interval between ``T_lmr`` outputs.

* **t_Tpot_start** (default :f:var:`t_Tpot_start=0.0 <t_tpot_start>`) is a real, which defines the time to start writing ``T_lmr`` outputs.

* **t_Tpot_stop** (default :f:var:`t_Tpot_stop=0.0 <t_tpot_stop>`) is a real, which defines the time to stop writing ``T_lmr`` outputs.

