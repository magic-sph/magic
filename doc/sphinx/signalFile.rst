Interactive communication with the code using ``signal.TAG``
############################################################

It is possible to interactively communicate with the MagIC code **during a run**,
using a file which is systematically created when the simulation starts, called
**signal.TAG**. By default, this file contains only the word ``NOT`` and does
nothing to the simulation. Replacing ``NOT`` by one of the following allowed
keywords will have some influence on the outputs or possibly force the code to
terminate its execution:

   * **END**: Changing the word ``NOT`` to ``END`` will cause the code to
     finish after the current time step and write all the outputs as if it was
     programmed to finish at that time from the start. This will thus normally
     produce the :ref:`checkpoint_end.TAG <secRestartFile>` file that will possibly
     allow you to continue this run later at your convenience.


   * **GRA**: Changing the word ``NOT`` to ``GRA`` will cause the code to produce
     a graphic ouptut file :ref:`G_#.TAG <secGraphFile>`. The keyword will be
     automatically restored to ``NOT`` once the graphic file has been produced.


   * **RST**: Changing the word ``NOT`` to ``RST`` will cause the code to produce
     a restart file :ref:`checkpoint_t#.TAG <secRestartFile>`. The keyword will then be
     restored to ``NOT`` once the restart file has been written.


   * **SPE**: Changing the word ``NOT`` to ``SPE`` will cause the code to produce
     spectra :ref:`kin_spec_#.TAG <secKinSpecFile>` (and possibly 
     :ref:`mag_spec_#.TAG <secMagSpecFile>` and `T_spec_#.TAG <secTSpecFile>` depending
     if the run is magnetic or not, or if it solves a temperature/entropy equation).
     Once the spectra files have been written, the keyword will be automatically replaced
     by ``NOT``.

   * **POT**: Changing the word ``NOT`` to ``POT`` will cause the code to produce
     the potential files :ref:`V_lmr_#.TAG <secVpotFile>` (and possibly 
     :ref:`B_lmr_#.TAG <secBpotFile>` and `T_lmr_#.TAG <secTpotFile>` depending
     if the run is magnetic or not, or if it solves a temperature/entropy equation).
     Once the potential files have been written, the keyword will be automatically replaced
     by ``NOT``.


   .. note:: Those keywords are **case-insensitive**.


Instead of editing the file with your favorite editor to specify the requested
keyword, we recommand using instead the shell command ``echo`` to avoid some
possible crash during the code execution when writing into the ``signal.TAG``
file. For instance, if you want a :ref:`graphic output file <secGraphFile>`, just use
the following command (adapted to your current :ref:`TAG <varTAG>`):

   .. code-block:: bash

      $ echo GRA > signal.TAG
